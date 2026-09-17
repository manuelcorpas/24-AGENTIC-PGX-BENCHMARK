"""
Blinded CPIC pharmacogenomics study - 110 cases.

Pipeline
--------
  genotype string  ->  diplotype  ->  CPIC phenotype  ->  CPIC recommendation

Design
------
* Star-allele calling is done by CONSTRAINT SOLVING over the CPIC allele-definition
  table (the same logic PharmCAT uses): enumerate all ordered pairs of alleles whose
  definitions reproduce the observed genotype at every typed locus, discard any allele
  that has a defining variant at an UNTYPED locus (it cannot be confirmed), then rank
  surviving candidates by parsimony (fewest non-reference alleles, then fewest
  defining variants).
* Non-SNV inputs (copy number, TA-repeat, HLA serotype, mtDNA, named G6PD haplotypes)
  are handled by explicit, individually documented rules.
* Phenotype comes from the CPIC `diplotype` table (diplotype -> generesult + lookupkey).
* Recommendation comes from the CPIC `recommendation` table, matched on `lookupkey`.
  When a guideline is multi-gene and only one gene was typed, the untyped partner gene
  is set to CPIC's own "No Result" key rather than being assumed normal.
* Where CPIC publishes no machine-readable recommendation (warfarin dose table,
  ondansetron, irinotecan, PEG-IFN/IFNL3) the provided guideline corpus is used and the
  source is recorded per case.

Sources are recorded per gene in SOURCES (see METHODS.md).
"""
from __future__ import annotations

import itertools
import json
import os
import re
from collections import defaultdict

import pandas as pd
import requests

CPIC = "https://api.cpicpgx.org/v1/"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cpic_cache.json")

GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "CYP2B6", "CYP3A5", "CYP4F2", "DPYD", "G6PD",
         "HLA-A", "HLA-B", "IFNL3", "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT",
         "UGT1A1", "VKORC1", "CFTR"]

DRUGNAMES = ["warfarin", "codeine", "clopidogrel", "voriconazole", "citalopram", "omeprazole",
             "tamoxifen", "paroxetine", "amitriptyline", "efavirenz", "fluorouracil",
             "azathioprine", "mercaptopurine", "thioguanine", "carbamazepine", "ondansetron",
             "phenytoin", "celecoxib", "tacrolimus", "capecitabine", "simvastatin",
             "atorvastatin", "irinotecan", "atazanavir", "rasburicase", "primaquine",
             "abacavir", "oxcarbazepine", "allopurinol", "peginterferon alfa-2a",
             "peginterferon alfa-2b", "ivacaftor", "gentamicin", "sevoflurane"]

# case drug label -> CPIC drug used for the recommendation lookup
DRUG_PROXY = {"thiopurines": "mercaptopurine", "aminoglycosides": "gentamicin",
              "volatile-anaesthetics": "sevoflurane", "PEG-IFN-α": "peginterferon alfa-2a"}


# --------------------------------------------------------------------------- fetch
def _get(table, **params):
    params.setdefault("limit", 100000)
    r = requests.get(CPIC + table, params=params, timeout=180)
    r.raise_for_status()
    return r.json()


def fetch_cpic(force=False):
    """Download the CPIC reference tables we need (cached to disk)."""
    if os.path.exists(CACHE) and not force:
        return json.load(open(CACHE))
    ginlist = "in.(" + ",".join(f'"{g}"' for g in GENES) + ")"
    data = {
        "retrieved": pd.Timestamp.utcnow().isoformat(),
        "allele_definition": _get(
            "allele_definition",
            select=("id,genesymbol,name,structuralvariation,matchesreferencesequence,"
                    "allele_location_value(variantallele,sequence_location(dbsnpid,name))"),
            genesymbol=ginlist),
        "allele": _get("allele", select="genesymbol,name,clinicalfunctionalstatus,activityvalue",
                       genesymbol=ginlist),
        "drug": _get("drug", select="drugid,name"),
        "pair": _get("pair", select="genesymbol,drugid,cpiclevel,pgxtesting"),
    }
    dmap = {d["name"].lower(): d["drugid"] for d in data["drug"]}
    ids = [dmap[n] for n in DRUGNAMES if n in dmap]
    data["recommendation"] = _get(
        "recommendation",
        select=("id,drugid,implications,drugrecommendation,classification,phenotypes,"
                "activityscore,allelestatus,lookupkey,population,comments"),
        drugid="in.(" + ",".join(f'"{i}"' for i in ids) + ")")
    data["diplotype"] = {}
    for g in GENES:
        if g in ("RYR1", "G6PD", "CFTR"):        # huge tables; fetched filtered below
            continue
        data["diplotype"][g] = _get("diplotype", genesymbol=f"eq.{g}",
                                    select="genesymbol,diplotype,generesult,lookupkey,totalactivityscore")
    json.dump(data, open(CACHE, "w"))
    return data


# ------------------------------------------------------------------- allele model
class AlleleModel:
    """rsID-indexed CPIC allele definitions for one gene."""

    REF = {"CYP2D6": "*1", "CYP2C19": "*38", "CYP2C9": "*1", "CYP2B6": "*1", "CYP3A5": "*1",
           "TPMT": "*1", "NUDT15": "*1", "SLCO1B1": "*1", "UGT1A1": "*1", "CYP4F2": "*1",
           "DPYD": "Reference", "G6PD": "B (reference)",
           "VKORC1": "rs9923231 reference (C)", "IFNL3": "rs12979860 reference (C)",
           "MT-RNR1": "Reference", "RYR1": "Reference",
           "CFTR": "ivacaftor non-responsive CFTR sequence"}

    def __init__(self, gene, allele_definitions, functions):
        self.gene = gene
        self.defs = {}
        self.sv = set()
        for ad in allele_definitions:
            if ad["genesymbol"] != gene:
                continue
            v = {}
            for lv in ad["allele_location_value"]:
                sl = lv["sequence_location"]
                key = sl["dbsnpid"] or sl["name"]
                v[key] = lv["variantallele"]
            self.defs[ad["name"]] = v
            if ad["structuralvariation"]:
                self.sv.add(ad["name"])
        self.ref = self.REF[gene]
        self.refvals = self.defs[self.ref]
        self.func = {n: functions.get((gene, n), ("", "")) for n in self.defs}

    # CPIC allele definitions use IUPAC ambiguity codes where sub-alleles differ
    IUPAC = {"R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
             "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

    def value(self, allele, locus):
        return self.defs[allele].get(locus, self.refvals.get(locus))

    def expand(self, allele, locus):
        """Set of bases an allele definition permits at a locus."""
        v = self.value(allele, locus)
        if v is None:
            return set()
        if len(v) == 1 and v in self.IUPAC:
            return set(self.IUPAC[v])
        return {v}

    def defining_loci(self, allele):
        """Loci where this allele DEFINITELY differs from the reference allele.

        A position carrying an IUPAC ambiguity code that still permits the reference
        base is not a defining variant - the allele's sub-alleles disagree there.
        """
        out = set()
        for l in self.defs[allele]:
            if l not in self.refvals:
                continue
            if not (self.expand(self.ref, l) & self.expand(allele, l)):
                out.add(l)
        return out

    def solve(self, genotypes, ploidy=2):
        """Call a diplotype (or a single allele when ploidy=1, i.e. hemizygous X / mtDNA).

        genotypes: {rsid: 'A/B'} for ploidy 2, {rsid: 'A'} for ploidy 1.

        Every non-structural allele in the gene is a candidate. An allele is compatible
        if its definition matches the observation at every TYPED locus; defining variants
        at untyped loci are unobserved and therefore unconfirmed. Candidates are ranked by
          1. number of UNCONFIRMED defining variants (assumption burden -> prefer 0),
          2. number of non-reference alleles (prefer the reference allele),
          3. total number of defining variants,
        which reproduces the tag-SNP reading a clinical laboratory applies to a targeted
        panel while never silently inventing a variant at a locus that was assayed.
        """
        typed = list(genotypes)
        obs = {l: sorted(genotypes[l].split("/")) for l in typed}
        # candidates: non-structural alleles that actually carry a variant definition
        # (CPIC lists hybrid/combination haplotypes such as CYP2D6 '*13+*1' with no
        #  positional definition at all; they would match any observation vacuously)
        cands = sorted(a for a in self.defs
                       if a not in self.sv and (self.defs[a] or a == self.ref))
        unconf = {a: len(self.defining_loci(a) - set(typed)) for a in cands}
        hits = []
        def ok1(a, l):
            return obs[l][0] in self.expand(a, l)

        def ok2(a, b, l):
            x, y = obs[l]
            ea, eb = self.expand(a, l), self.expand(b, l)
            return (x in ea and y in eb) or (y in ea and x in eb)

        if ploidy == 1:
            for a in cands:
                if all(ok1(a, l) for l in typed):
                    hits.append((unconf[a], a != self.ref, len(self.defining_loci(a)), a, None))
        else:
            for a, b in itertools.combinations_with_replacement(cands, 2):
                if all(ok2(a, b, l) for l in typed):
                    hits.append((unconf[a] + unconf[b],
                                 (a != self.ref) + (b != self.ref),
                                 len(self.defining_loci(a)) + len(self.defining_loci(b)), a, b))
        hits.sort()
        out = []
        for u, nv, nd, a, b in hits:
            name = a if b is None else (f"{a}/{a}" if a == b else f"{a}/{b}")
            out.append({"call": name, "unconfirmed": u, "n_variant_alleles": nv,
                        "n_defining": nd})
        return out


# --------------------------------------------------------------------- genotype IO
RS_RE = re.compile(r"(rs\d+)\s+([ACGTRYacgt])\s*/\s*([ACGTRYacgt])")
RS_HEMI_RE = re.compile(r"(rs\d+)\s+([ACGT])(?![/\w])")


def parse_rs_genotypes(text):
    """Extract {rsid: 'A/B'} diploid calls from a free-text genotype description."""
    return {m.group(1): "/".join(sorted([m.group(2).upper(), m.group(3).upper()]))
            for m in RS_RE.finditer(text)}


def parse_rs_hemizygous(text):
    return {m.group(1): m.group(2).upper() for m in RS_HEMI_RE.finditer(text)}


# ------------------------------------------------------------------- phenotype IO
def phenotype_for(diplotable, gene, diplotype):
    t = diplotable.get(gene, {})
    for key in (diplotype, "/".join(reversed(diplotype.split("/")))):
        if key in t:
            return t[key]["generesult"], t[key]["lookupkey"]
    return None, None


def recommendation_for(recs, drug, lookupkey):
    """Exact match of the CPIC lookupkey dict."""
    out = [r for r in recs if r["drug"] == drug and r["lookupkey"] == lookupkey]
    return out


if __name__ == "__main__":
    raise SystemExit("import this module from run_study.py")
