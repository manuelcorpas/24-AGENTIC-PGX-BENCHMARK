#!/usr/bin/env python3
"""
Blinded pharmacogenomics study - CPIC interpretation pipeline
=============================================================

Input : cases_blinded.json          (110 cases: id, gene, drug, genotype)
        cpic_guideline_corpus.json  (CPIC guideline text, 21 genes, snapshot 2026-05-09)

Output: claude_science_answers.csv          case_id,diplotype,phenotype,recommendation,abstained
        table_diplotype_to_phenotype.csv    diplotype -> phenotype lookup actually used
        table_phenotype_to_recommendation.csv  phenotype -> CPIC recommendation actually used
        table_variant_to_allele.csv         rsID/base -> star allele definitions used
        case_audit.csv                      full per-case audit trail (strand, rule, source, lethal class)
        sources_and_versions.csv            source + version for every gene

Design
------
Three deterministic stages, each grounded in a named, versioned source:

  (1) genotype string -> diplotype
      Variant bases are matched against CPIC allele-definition tables
      (sequence_location + allele_location_value). CPIC stores variant bases on the
      FORWARD GENOMIC strand (its `chromosomelocation` HGVS), while several of the
      supplied genotype strings are written in GENE orientation. Orientation is
      therefore resolved per variant: the literal bases are tried first, and if they
      are not valid alleles at that position the reverse complement is tried. The
      orientation used is recorded for every case.
      Star alleles are assigned from a variant->allele table derived from CPIC by the
      minimal-definition rule (the allele with the fewest defining positions that
      carries the observed base), with documented overrides where a limited SNV panel
      is genuinely ambiguous.

  (2) diplotype -> phenotype
      Exact lookup in the CPIC `diplotype` table (both allele orders tried), which
      carries `generesult` (phenotype) and `lookupkey` (activity score or phenotype
      string) for every CPIC-defined diplotype.

  (3) phenotype -> recommendation
      Exact lookup in the CPIC `recommendation` table on drugid + lookupkey. For
      multi-gene guidelines (thiopurines TPMT+NUDT15, TCAs CYP2D6+CYP2C19, phenytoin
      CYP2C9+HLA-B, MH RYR1+CACNA1S) the gene that was not genotyped is set to
      CPIC's own "No Result" key rather than being assumed normal.
      Three CPIC guidelines carry no machine-readable recommendation rows in the API
      (ondansetron/tropisetron, warfarin, IFNL3/peginterferon); for those the
      recommendation is taken verbatim from the supplied guideline corpus and is
      labelled as corpus-derived.

No published benchmark, paper or repository reporting answers to these cases was
consulted or searched for.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "cpic_api_cache")
API_ROOT = "https://api.cpicpgx.org/v1"
os.makedirs(CACHE, exist_ok=True)


# --------------------------------------------------------------------------------------
# CPIC data API client (on-disk cache so the run is reproducible and re-runnable)
# --------------------------------------------------------------------------------------
def api(path: str, **params) -> list:
    key = hashlib.sha256((path + "|" + json.dumps(params, sort_keys=True)).encode()).hexdigest()[:24]
    fp = os.path.join(CACHE, f"{path.replace('/', '_')}__{key}.json")
    if os.path.exists(fp):
        with open(fp) as fh:
            return json.load(fh)
    q = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"{API_ROOT}/{path}" + (f"?{q}" if q else "")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                data = json.loads(r.read())
            with open(fp, "w") as fh:
                json.dump(data, fh)
            return data
        except Exception as exc:                      # transient API/network failure
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"CPIC API failed for {url}: {last}")


def complement(base: str) -> str:
    return "".join({"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, b) for b in base[::-1])


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------
with open(os.path.join(BASE, "cases_blinded.json")) as fh:
    CASES = json.load(fh)
with open(os.path.join(BASE, "cpic_guideline_corpus.json")) as fh:
    CORPUS = json.load(fh)

CORPUS_GENES = CORPUS["genes"]
CORPUS_TAG = f"corpus {CORPUS['schema_version']} snapshot {CORPUS['generated_at_utc']}"

# gene symbol as used in the case file -> CPIC gene symbol
GENE_SYMBOL = {
    "HLA-B*57:01": "HLA-B", "HLA-B*15:02": "HLA-B", "HLA-B*58:01": "HLA-B", "HLA-A*31:01": "HLA-A",
}
# drug label as used in the case file -> CPIC drug name used for the recommendation lookup
DRUG_LOOKUP = {
    "thiopurines": "azathioprine",          # TPMT/NUDT15 thiopurine guideline, AZA row
    "PEG-IFN-α": "peginterferon alfa-2a",
    "aminoglycosides": "gentamicin",        # MT-RNR1 guideline: one recommendation for all AGs
    "volatile-anaesthetics": "desflurane",  # MH guideline: one recommendation for all volatiles
}
# CPIC recommendation `population` to use when a guideline is stratified, and why
POPULATION = {
    "clopidogrel": ("CVI ACS PCI ", "ACS/PCI - the indication with the life-threatening failure mode"),
    "voriconazole": ("adults", "adult dosing"),
    "phenytoin": ("PHT naive", "new prescription"),
    "carbamazepine": ("CBZ naive", "new prescription"),
    "oxcarbazepine": ("OXC naive", "new prescription"),
}


# --------------------------------------------------------------------------------------
# Stage 0: pull the CPIC reference tables needed for these 21 genes / 36 pairs
# --------------------------------------------------------------------------------------
CPIC_GENES = ["CFTR", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "DPYD", "G6PD",
              "HLA-A", "HLA-B", "IFNL3", "MT-RNR1", "NUDT15", "RYR1", "CACNA1S", "SLCO1B1", "TPMT",
              "UGT1A1", "VKORC1"]
gene_in = "in.(" + ",".join(f'"{x}"' for x in CPIC_GENES) + ")"

DRUGS = api("drug", limit=5000)
DRUG_ID = {}
for d in DRUGS:
    DRUG_ID.setdefault(d["name"].lower(), d["drugid"])
GUIDELINES = {x["id"]: x for x in api("guideline", limit=500)}
PAIRS = api("pair", limit=2000)
SEQLOC = api("sequence_location", genesymbol=gene_in, limit=20000)
ALLELE_DEF = {d["id"]: d for d in api("allele_definition", genesymbol=gene_in, limit=5000)}
ALLELES = api("allele", genesymbol=gene_in, limit=5000)
RECS = api("recommendation", limit=5000)

LOC_BY_RS = {}
for s in SEQLOC:
    if s.get("dbsnpid"):
        LOC_BY_RS.setdefault(s["dbsnpid"], []).append(s)


def variant_alleles_at(locid: int) -> list:
    """All CPIC allele definitions carrying a value at one sequence location."""
    return api("allele_location_value", locationid=f"eq.{locid}", limit=5000)


def n_defining(defid: int) -> int:
    return len(api("allele_location_value", alleledefinitionid=f"eq.{defid}", limit=5000))


# ---- variant -> star allele table, derived from CPIC by the minimal-definition rule ----
# Overrides record where a single-SNV panel is genuinely ambiguous and which allele the
# CPIC guideline table is keyed on (the phenotype is unchanged by the choice in every
# case below, and the ambiguity is reported in the methods note).
ALLELE_OVERRIDE = {
    ("CYP2B6", "rs3745274", "T"): ("*6", "rs3745274 (516G>T) alone is shared by *6 and *9; "
                                         "the CPIC efavirenz table is keyed on *6. Both are "
                                         "decreased function, so the phenotype is identical. "
                                         "785A>G (rs2279343) was not reported."),
    ("CYP2D6", "rs1065852", "A"): ("*10", "rs1065852 (100C>T) alone: *10 is the minimal CPIC "
                                          "allele carrying it (4180G>C not reported)."),
}


def build_variant_allele_table(genes_rs: dict) -> dict:
    """{(gene, rsid, forward_base): (allele_name, note)} from the CPIC allele definitions."""
    table = {}
    for gene, rsids in genes_rs.items():
        for rs in rsids:
            for loc in LOC_BY_RS.get(rs, []):
                if loc["genesymbol"] != gene:
                    continue
                rows = variant_alleles_at(loc["id"])
                by_base = {}
                for r in rows:
                    base = r["variantallele"]
                    if base not in "ACGT":            # IUPAC/indel codes: not used here
                        continue
                    defid = r["alleledefinitionid"]
                    name = ALLELE_DEF[defid]["name"]
                    by_base.setdefault(base, []).append((n_defining(defid), name))
                for base, cands in by_base.items():
                    cands.sort()
                    ov = ALLELE_OVERRIDE.get((gene, rs, base))
                    if ov:
                        table[(gene, rs, base)] = ov
                    else:
                        table[(gene, rs, base)] = (
                            cands[0][1],
                            f"minimal CPIC allele definition carrying {rs}={base} "
                            f"({cands[0][0]} defining position(s))")
    return table


RS_BY_GENE = {
    "CYP2D6": ["rs3892097", "rs1065852"],
    "CYP2C19": ["rs4244285", "rs12248560"],
    "CYP2C9": ["rs1057910"],
    "CYP2B6": ["rs3745274", "rs2279343"],
    "CYP3A5": ["rs776746"],
    "CYP4F2": ["rs2108622"],
    "DPYD": ["rs3918290", "rs56038477"],
    "TPMT": ["rs1800460", "rs1142345"],
    "NUDT15": ["rs116855232"],
    "SLCO1B1": ["rs4149056"],
    "G6PD": ["rs1050828", "rs1050829", "rs5030868"],
    "VKORC1": ["rs9923231"],
    "IFNL3": ["rs12979860"],
    "CFTR": ["rs75527207"],
}
VAR2ALLELE = build_variant_allele_table(RS_BY_GENE)

# reference (non-variant) base per (gene, rsid), on the forward genomic strand
REF_BASE = {}
for gene, rsids in RS_BY_GENE.items():
    for rs in rsids:
        for loc in LOC_BY_RS.get(rs, []):
            if loc["genesymbol"] != gene:
                continue
            m = re.search(r"g\.\d+([ACGT])>([ACGT])", loc["chromosomelocation"] or "")
            if m:
                REF_BASE[(gene, rs)] = m.group(1)

# reference allele name per gene (what a non-variant chromosome is called by CPIC)
REF_ALLELE = {"CYP2D6": "*1", "CYP2C19": "*1", "CYP2C9": "*1", "CYP2B6": "*1", "CYP3A5": "*1",
              "CYP4F2": "*1", "TPMT": "*1", "NUDT15": "*1", "SLCO1B1": "*1", "UGT1A1": "*1",
              "DPYD": "Reference", "RYR1": "Reference", "G6PD": "B (reference)"}


# --------------------------------------------------------------------------------------
# Stage 1: genotype string -> diplotype
# --------------------------------------------------------------------------------------
def parse_calls(genotype: str) -> dict:
    """{rsid: [observed bases]} from strings like 'rs3892097 C/T', 'rs1050828 T', 'rs8175347 TA6/TA7'."""
    calls = {}
    for rs, call in re.findall(r"(rs\d+)\s+((?:TA\d+|[ACGT])(?:\s*/\s*(?:TA\d+|[ACGT]))?)", genotype):
        calls.setdefault(rs, []).append([b.strip() for b in call.split("/")])
    return calls


def resolve_orientation(gene: str, rs: str, bases: list) -> tuple:
    """Return (forward_genomic_bases, orientation_label). CPIC stores forward-genomic bases."""
    valid = {b for (gn, r, b) in VAR2ALLELE if gn == gene and r == rs}
    ref = REF_BASE.get((gene, rs))
    if ref:
        valid = valid | {ref}
    if not valid:
        return bases, "unknown (no CPIC allele values at this position)"
    if all(b in valid for b in bases):
        return bases, "forward-genomic (as given)"
    flipped = [complement(b) for b in bases]
    if all(b in valid for b in flipped):
        return flipped, "gene-orientation input -> complemented to forward-genomic"
    return bases, "UNRESOLVED"


class Call:
    """One resolved diplotype with its provenance."""

    def __init__(self, diplotype, rule, orientation="n/a", caveat="", abstain=False, phenotype_key=None):
        self.diplotype = diplotype
        self.rule = rule
        self.orientation = orientation
        self.caveat = caveat
        self.abstain = abstain
        self.phenotype_key = phenotype_key or diplotype   # diplotype string used for CPIC lookup


# CPIC writes the reference allele LAST for the genes whose alleles are named by HGVS
REF_LAST = {"DPYD", "RYR1", "G6PD"}


def star_pair(gene: str, alleles: list) -> str:
    """Order a diplotype the way CPIC writes it in its diplotype table."""
    ref = REF_ALLELE.get(gene)
    a, b = alleles
    if gene in REF_LAST:
        if ref and a == ref and b != ref:
            a, b = b, a
    elif ref and b == ref and a != ref:
        a, b = b, a
    return f"{a}/{b}"


def call_snv_gene(gene: str, genotype: str) -> Call:
    """Generic SNV panel -> diplotype for the star-allele genes."""
    calls = parse_calls(genotype)
    ref_allele = REF_ALLELE.get(gene, "*1")

    if not calls:                                          # "All <gene> SNPs homozygous reference"
        if re.search(r"homozygous reference|risk variants reference|no deficiency variants|no variant",
                     genotype, re.I):
            return Call(star_pair(gene, [ref_allele, ref_allele]),
                        "no variant detected on the tested panel -> reference/reference",
                        caveat="A negative SNV panel cannot exclude untested alleles "
                               "(or, for CYP2D6, gene deletions/duplications); CPIC reports "
                               "such a result as the reference diplotype.")
        return Call("", "no parseable variant call", abstain=True,
                    caveat="genotype string carries no interpretable variant call")

    # ---- resolve every call to forward-genomic bases and count variant copies ----------
    copies, orients, notes = {}, [], []
    ploidy = 2
    for rs, callsets in calls.items():
        bases = callsets[-1]
        fwd, orient = resolve_orientation(gene, rs, bases)
        orients.append(f"{rs}: {orient}")
        if orient == "UNRESOLVED":
            return Call("", f"could not resolve {rs} bases {bases} against CPIC allele values",
                        orientation="; ".join(orients), abstain=True,
                        caveat=f"{rs} call {'/'.join(bases)} matches neither the forward-genomic "
                               f"nor the complemented CPIC allele values at that position")
        if len(fwd) == 1:
            ploidy = 1                                        # hemizygous call
        ref = REF_BASE.get((gene, rs))
        var = [b for b in fwd if b != ref]
        if var:
            copies[rs] = (len(var), var[0])
    orient_s = "; ".join(orients)

    # ---- no variant on the tested panel ------------------------------------------------
    if not copies:
        dip = ref_allele if ploidy == 1 else star_pair(gene, [ref_allele, ref_allele])
        return Call(dip, "no variant allele detected at any tested position -> reference",
                    orientation=orient_s,
                    caveat="A negative SNV panel cannot exclude untested alleles (or, for CYP2D6, "
                           "gene deletions/duplications); CPIC reports such a result as the "
                           "reference diplotype.")

    # ---- one variant position -----------------------------------------------------------
    if len(copies) == 1:
        rs, (n, base) = next(iter(copies.items()))
        name, note = VAR2ALLELE[(gene, rs, base)]
        notes.append(note)
        chroms = [name] * n + [ref_allele] * (ploidy - n)
        dip = chroms[0] if ploidy == 1 else star_pair(gene, chroms)
        return Call(dip, "; ".join(notes), orientation=orient_s)

    # ---- two or more variant positions: resolve phase against CPIC haplotypes ------------
    cis = cis_allele(gene, [(rs, b) for rs, (_, b) in copies.items()])
    counts = {rs: n for rs, (n, _) in copies.items()}
    if cis and len(set(counts.values())) == 1:
        n = next(iter(counts.values()))
        notes.append(f"all observed variants are admitted by the CPIC definition of {cis}; "
                     f"called in cis (phase was not provided)")
        chroms = [cis] * n + [ref_allele] * (ploidy - n)
        dip = chroms[0] if ploidy == 1 else star_pair(gene, chroms)
        return Call(dip, "; ".join(notes), orientation=orient_s)

    chroms = []
    for rs, (n, base) in copies.items():
        name, note = VAR2ALLELE[(gene, rs, base)]
        notes.append(note)
        chroms.extend([name] * n)
    if len(chroms) == ploidy:
        notes.append("variant positions are not explained by a single CPIC haplotype; "
                     "assigned to separate chromosomes (in trans)")
        dip = chroms[0] if ploidy == 1 else star_pair(gene, chroms)
        return Call(dip, "; ".join(notes), orientation=orient_s)
    return Call("", f"variant pattern cannot be resolved into {ploidy} allele(s): {chroms}",
                abstain=True, orientation=orient_s,
                caveat="The observed combination of variant copy numbers is not explained by any "
                       "single CPIC haplotype and does not partition into the expected number of "
                       "chromosomes; phase or additional positions would be required.")


IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def allele_size(gene: str, name: str) -> int:
    return min([n_defining(d) for d, v in ALLELE_DEF.items()
                if v["name"] == name and v["genesymbol"] == gene] or [10 ** 6])


def cis_allele(gene: str, variants: list):
    """
    Name of the CPIC allele whose definition admits the observed variant base at EVERY
    position given, i.e. the haplotype that would explain all of them in cis.
    `variants` is [(rsid, forward_genomic_base), ...].
    IUPAC ambiguity codes in a CPIC definition (e.g. CYP2D6*4 carries R at rs1065852)
    are treated as admitting either base, because the sub-alleles differ there.
    """
    sets = []
    for rs, base in variants:
        locs = [l for l in LOC_BY_RS.get(rs, []) if l["genesymbol"] == gene]
        if not locs:
            return None
        names = {ALLELE_DEF[r["alleledefinitionid"]]["name"]
                 for r in variant_alleles_at(locs[0]["id"])
                 if base in IUPAC.get(r["variantallele"], "")}
        sets.append(names)
    shared = set.intersection(*sets) if sets else set()
    shared = {nm for nm in shared if nm not in ("Reference", REF_ALLELE.get(gene))}
    if not shared:
        return None
    return sorted(shared, key=lambda nm: (allele_size(gene, nm), nm))[0]


# ---- gene-specific handlers for the non-SNV-panel genes -------------------------------
def call_cyp2d6(genotype: str) -> Call:
    if re.search(r"\*1xN duplication|duplication detected", genotype, re.I):
        return Call("*1/*1xN",
                    "copy-number assay reports a duplicated normal-function allele with no "
                    "variant SNV -> *1/*1xN",
                    caveat="N is not resolved by the assay; every xN>=2 gives activity score >=3.0, "
                           "so the CPIC phenotype (ultrarapid metabolizer) is the same. Phenotype "
                           "looked up as *1/*1x2 (activity score 3.0).",
                    phenotype_key="*1/*1x2")
    return call_snv_gene("CYP2D6", genotype)


def call_ugt1a1(genotype: str) -> Call:
    m = re.findall(r"TA(\d)", genotype)
    if len(m) == 2:
        alleles = ["*28" if n == "7" else "*1" for n in m]
        return Call(star_pair("UGT1A1", alleles),
                    "UGT1A1 promoter (TA)n repeat: TA7 = *28 (rs8175347/rs3064744), TA6 = *1")
    return Call("", "unparsed UGT1A1 repeat genotype", abstain=True)


def call_hla(gene_label: str, genotype: str) -> Call:
    allele = gene_label.split("*", 1)[1]
    if re.search(r"negative|no risk allele", genotype, re.I):
        return Call(f"*{allele} negative", "HLA typing negative")
    if re.search(r"positive", genotype, re.I):
        note = ""
        if "rs2395029" in genotype:
            note = ("Called from the HCP5 rs2395029 tag SNP, not direct HLA-B typing. The tag is "
                    "in strong LD with HLA-B*57:01 but is a proxy; confirmatory HLA-B typing is "
                    "the reference method.")
        return Call(f"*{allele} positive", "HLA typing positive", caveat=note)
    return Call("", "unparsed HLA typing result", abstain=True)


def call_g6pd(genotype: str) -> Call:
    """G6PD is X-linked: hemizygous males get one allele, females two."""
    calls = parse_calls(genotype)
    fwd = {}
    orients = []
    single_base_only = all(len(cs) == 1 for sets in calls.values() for cs in sets)
    for rs, callsets in calls.items():
        allb = [b for cs in callsets for b in cs]
        f, orient = resolve_orientation("G6PD", rs, allb)
        orients.append(f"{rs}: {orient}")
        fwd[rs] = set(f)
    orient_s = "; ".join(orients)
    med = "Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham"
    a_minus = "A- 202A_376G"
    ref = "B (reference)"
    var202 = "T" in fwd.get("rs1050828", set())
    var376 = "C" in fwd.get("rs1050829", set())
    var563 = "A" in fwd.get("rs5030868", set())
    female_het = bool(re.search(r"heterozygous female", genotype, re.I)) or \
        ("C" in fwd.get("rs1050828", set()) and "T" in fwd.get("rs1050828", set()))
    # A single base per position is a single-X (hemizygous) call; asserting a homozygous
    # female diplotype would claim more than the input states.
    hemizygous = (bool(re.search(r"hemizygous", genotype, re.I)) or single_base_only) \
        and not female_het

    if var563:
        var = med
        rule = "rs5030868 c.563C>T = G6PD Mediterranean (WHO class II, deficient)"
    elif var202 and var376:
        var = a_minus
        rule = "rs1050828 c.202G>A + rs1050829 c.376A>G in cis = G6PD A- 202A/376G (deficient)"
    elif not (var202 or var376 or var563):
        return Call(ref, "no deficiency variant detected on the tested panel", orientation=orient_s,
                    caveat="A negative panel cannot exclude untested G6PD variants.")
    else:
        return Call("", "partial G6PD variant pattern", abstain=True, orientation=orient_s)

    if female_het:
        return Call(f"{var}/{ref}", rule + "; heterozygous female", orientation=orient_s,
                    caveat="Heterozygous females have variable enzyme activity because of "
                           "X-inactivation; CPIC assigns phenotype 'Variable' and requires "
                           "measured enzyme activity.")
    if hemizygous:
        return Call(var, rule + "; hemizygous male (single X allele)", orientation=orient_s,
                    caveat="Case text allows either a hemizygous male or a homozygous female; "
                           "both give the same CPIC phenotype (Deficient).")
    return Call(f"{var}/{var}", rule + "; homozygous", orientation=orient_s)


def call_cftr(genotype: str) -> Call:
    if "rs75527207" in genotype or "G551D" in genotype:
        return Call("G551D positive (heterozygous or homozygous)",
                    "rs75527207 A allele = G551D (gating class III); case states G/A or A/A",
                    caveat="Het (G551D/other) and hom (G551D/G551D) cannot be distinguished from "
                           "the given string, and CPIC gives the same result for both "
                           "('Homozygous or Heterozygous G551D-CFTR, rs75527207 genotype AA or "
                           "AG'); phenotype looked up as G551D/G551D.",
                    phenotype_key="G551D/G551D")
    if "rs113993960" in genotype or "F508del" in genotype:
        return Call("F508del/F508del",
                    "rs113993960 del/del = F508del homozygous",
                    caveat="CPIC classifies F508del/F508del under 'ivacaftor non-responsive CFTR "
                           "sequence'; looked up under that CPIC diplotype label.",
                    phenotype_key="ivacaftor non-responsive CFTR sequence/"
                                  "ivacaftor non-responsive CFTR sequence")
    return Call("", "unparsed CFTR genotype", abstain=True)


def call_mtrnr1(genotype: str) -> Call:
    if re.search(r"m\.1555A>G", genotype):
        return Call("m.1555A>G", "homoplasmic m.1555A>G (mitochondrial; no diploid diplotype)")
    if re.search(r"m\.1555A\b", genotype):
        return Call("m.1555A (reference)", "no MT-RNR1 risk variant detected",
                    caveat="CPIC assigns 'normal risk' when no MT-RNR1 variant is found; the "
                           "reference has no diplotype row in the CPIC table, so the phenotype "
                           "was taken from the CPIC recommendation key for normal risk.",
                    phenotype_key="__MTRNR1_NORMAL__")
    return Call("", "unparsed MT-RNR1 genotype", abstain=True)


def call_ryr1(genotype: str) -> Call:
    if re.search(r"no pathogenic variant", genotype, re.I):
        return Call("Reference/Reference", "no RYR1 pathogenic variant detected",
                    caveat="A negative RYR1 result does not exclude malignant hyperthermia "
                           "susceptibility; CPIC calls this 'Uncertain Susceptibility'.")
    if re.search(r"c\.1840C>T|Arg614Cys|RYR1 mutation", genotype, re.I):
        return Call("c.1840C>T/Reference",
                    "RYR1 p.Arg614Cys (c.1840C>T), a CPIC MH-associated pathogenic variant, "
                    "heterozygous (autosomal dominant)")
    return Call("", "unparsed RYR1 genotype", abstain=True)


def call_vkorc1(genotype: str) -> Call:
    calls = parse_calls(genotype)
    bases = calls.get("rs9923231", [[]])[-1]
    fwd, orient = resolve_orientation("VKORC1", "rs9923231", bases)
    gene_or = [complement(b) for b in fwd]                 # report in the -1639 G>A convention
    dip = "/".join(f"-1639{b}" for b in gene_or)
    return Call(dip, "rs9923231 -1639G>A promoter genotype", orientation=f"rs9923231: {orient}",
                caveat="CPIC assigns no metaboliser phenotype to VKORC1; it is a warfarin "
                       "dose-sensitivity genotype used inside the CPIC dosing table.",
                phenotype_key="__VKORC1__:" + "".join(gene_or))


def call_cyp4f2(genotype: str) -> Call:
    c = call_snv_gene("CYP4F2", genotype)
    c.caveat = ("CPIC assigns no metaboliser phenotype to CYP4F2 and publishes no CYP4F2-specific "
                "dose action; it appears only as a minor term in some warfarin dosing algorithms.")
    c.phenotype_key = "__CYP4F2__:" + c.diplotype
    return c


def call_ifnl3(genotype: str) -> Call:
    calls = parse_calls(genotype)
    bases = calls.get("rs12979860", [[]])[-1]
    fwd, orient = resolve_orientation("IFNL3", "rs12979860", bases)
    dip = "rs12979860 " + "/".join(fwd)
    return Call(dip, "rs12979860 genotype (IFNL3/IL28B); IFNL3 has no star-allele diplotype "
                     "system and CPIC keys this guideline directly on the rs12979860 genotype",
                orientation=f"rs12979860: {orient}",
                phenotype_key="__IFNL3__:" + "".join(sorted(fwd)))


def make_diplotype(case: dict) -> Call:
    gene, gt = case["gene"], case["genotype"]
    if gene.startswith("HLA-"):
        return call_hla(gene, gt)
    return {
        "CYP2D6": lambda: call_cyp2d6(gt),
        "UGT1A1": lambda: call_ugt1a1(gt),
        "G6PD": lambda: call_g6pd(gt),
        "CFTR": lambda: call_cftr(gt),
        "MT-RNR1": lambda: call_mtrnr1(gt),
        "RYR1": lambda: call_ryr1(gt),
        "VKORC1": lambda: call_vkorc1(gt),
        "CYP4F2": lambda: call_cyp4f2(gt),
        "IFNL3": lambda: call_ifnl3(gt),
    }.get(gene, lambda: call_snv_gene(gene, gt))()


# --------------------------------------------------------------------------------------
# Stage 2: diplotype -> phenotype (CPIC diplotype table)
# --------------------------------------------------------------------------------------
# phenotypes for the genes CPIC does not phenotype-classify, taken from the guideline corpus
CORPUS_PHENOTYPE = {
    "__VKORC1__:AA": ("VKORC1 -1639 A/A (low warfarin dose requirement; no CPIC metaboliser "
                      "phenotype assigned)", None),
    "__VKORC1__:GG": ("VKORC1 -1639 G/G (reference; normal warfarin dose requirement; no CPIC "
                      "metaboliser phenotype assigned)", None),
    "__CYP4F2__:*3/*3": ("CYP4F2 *3/*3 (no CPIC phenotype assigned for CYP4F2)", None),
    "__CYP4F2__:*1/*1": ("CYP4F2 *1/*1 (no CPIC phenotype assigned for CYP4F2)", None),
    "__IFNL3__:CC": ("Favorable response genotype", None),
    "__IFNL3__:TT": ("Unfavorable response genotype", None),
    "__MTRNR1_NORMAL__": ("normal risk of aminoglycoside-induced hearing loss",
                          "normal risk of aminoglycoside-induced hearing loss"),
}


def lookup_phenotype(gene: str, key: str) -> tuple:
    """(phenotype, lookupkey_value, source) for a diplotype string."""
    if key in CORPUS_PHENOTYPE:
        ph, lk = CORPUS_PHENOTYPE[key]
        return ph, lk, CORPUS_TAG
    cpic_gene = GENE_SYMBOL.get(gene, gene)
    orders = [key]
    if "/" in key:
        a, b = key.split("/", 1)
        orders.append(f"{b}/{a}")
    for dp in orders:
        rows = api("diplotype", genesymbol=f"eq.{cpic_gene}",
                   diplotype=f"eq.{urllib.parse.quote(dp, safe='')}", limit=5)
        if rows:
            r = rows[0]
            lk = (r["lookupkey"] or {}).get(cpic_gene)
            return r["generesult"], lk, f"CPIC API diplotype table ({cpic_gene})"
    return None, None, None


# --------------------------------------------------------------------------------------
# Stage 3: phenotype -> CPIC recommendation
# --------------------------------------------------------------------------------------
RECS_BY_DRUG = {}
for r in RECS:
    RECS_BY_DRUG.setdefault(r["drugid"], []).append(r)

# CPIC's own key for a gene in a multi-gene guideline that was not genotyped
NO_RESULT = "No Result"

# Recommendations for the three CPIC guidelines that carry no rows in the data API.
# Text is condensed from the guideline tables in cpic_guideline_corpus.json.
CORPUS_RECS = {
    # ondansetron / tropisetron and CYP2D6  (PharmGKB PA166161954)
    ("ondansetron", "Ultrarapid Metabolizer"):
        ("Select an alternative antiemetic not predominantly metabolised by CYP2D6 "
         "(e.g. granisetron)", "Moderate"),
    ("ondansetron", "Normal Metabolizer"):
        ("Initiate therapy with the recommended starting dose", "Strong"),
    ("ondansetron", "Intermediate Metabolizer"):
        ("Insufficient evidence of clinical impact; initiate therapy with the recommended "
         "starting dose", "No recommendation"),
    ("ondansetron", "Poor Metabolizer"):
        ("Insufficient evidence of clinical impact; initiate therapy with the recommended "
         "starting dose", "No recommendation"),
    # IFNL3 and peginterferon-alpha-based regimens  (PharmGKB PA166110235) - guideline RETIRED
    ("peginterferon alfa-2a", "Favorable response genotype"):
        ("Approximately 70% chance of SVR after 48 weeks of PEG-IFN-alpha + ribavirin "
         "(~90% with a protease-inhibitor combination); weighs in favour of using a "
         "PEG-IFN-alpha-containing regimen", "Strong"),
    ("peginterferon alfa-2a", "Unfavorable response genotype"):
        ("Approximately 30% chance of SVR after 48 weeks of PEG-IFN-alpha + ribavirin "
         "(~60% with a protease-inhibitor combination); consider the implications before "
         "initiating a PEG-IFN-alpha-containing regimen", "Strong"),
}

# UGT1A1 / irinotecan is a CPIC *level A* pair for which CPIC has published NO guideline.
# The supplied corpus carries DPWG and RNPGx guidance instead; it is quoted here and
# explicitly attributed, so that no non-CPIC text is ever presented as a CPIC recommendation.
NON_CPIC_RECS = {
    ("irinotecan", "Poor Metabolizer"):
        ("no CPIC recommendation exists (CPIC level A pair, no CPIC guideline published). "
         "DPWG: start at 70% of the normal dose and increase as tolerated, guided by neutrophil "
         "count; RNPGx: for doses >250 mg/m2 reduce the initial dose by 30%", "DPWG/RNPGx"),
    ("irinotecan", "Intermediate Metabolizer"):
        ("no CPIC recommendation exists (CPIC level A pair, no CPIC guideline published). "
         "DPWG: no action needed for this gene-drug interaction", "DPWG"),
    ("irinotecan", "Normal Metabolizer"):
        ("no CPIC recommendation exists (CPIC level A pair, no CPIC guideline published). "
         "DPWG/RNPGx: no dose adjustment", "DPWG/RNPGx"),
}

# warfarin: CPIC publishes a dose table (mg/day), not a phenotype -> action row.
# Keys are the CPIC table cells from the corpus (VKORC1 -1639 genotype x CYP2C9 diplotype).
WARFARIN_TABLE = {                       # (VKORC1, CYP2C9) -> mg/day
    ("GG", "*1/*1"): "5-7", ("GG", "*1/*2"): "5-7", ("GG", "*1/*3"): "3-4",
    ("GG", "*2/*2"): "3-4", ("GG", "*2/*3"): "3-4", ("GG", "*3/*3"): "0.5-2",
    ("GA", "*1/*1"): "5-7", ("GA", "*1/*2"): "3-4", ("GA", "*1/*3"): "3-4",
    ("GA", "*2/*2"): "3-4", ("GA", "*2/*3"): "0.5-2", ("GA", "*3/*3"): "0.5-2",
    ("AA", "*1/*1"): "3-4", ("AA", "*1/*2"): "3-4", ("AA", "*1/*3"): "0.5-2",
    ("AA", "*2/*2"): "0.5-2", ("AA", "*2/*3"): "0.5-2", ("AA", "*3/*3"): "0.5-2",
}


def warfarin_recommendation(gene: str, diplotype: str) -> tuple:
    """CPIC warfarin dosing with only one of CYP2C9/VKORC1/CYP4F2 genotyped."""
    if gene == "CYP2C9":
        doses = {v for (vk, c9), v in WARFARIN_TABLE.items() if c9 == diplotype}
        if len(doses) == 1:
            d = doses.pop()
            return (f"use a validated pharmacogenetic dosing algorithm; CPIC dose range "
                    f"{d} mg/day for CYP2C9 {diplotype} at any VKORC1 genotype", "Strong")
        conds = ", ".join(f"VKORC1 {vk} {WARFARIN_TABLE[(vk, diplotype)]} mg/day"
                          for vk in ("GG", "GA", "AA"))
        return (f"use a validated pharmacogenetic dosing algorithm including CYP2C9 and VKORC1; "
                f"VKORC1 was not genotyped ({conds})", "Strong")
    if gene == "VKORC1":
        vk = diplotype.replace("-1639", "").replace("/", "")
        conds = ", ".join(f"CYP2C9 {c9} {WARFARIN_TABLE[(vk, c9)]} mg/day"
                          for c9 in ("*1/*1", "*1/*2", "*1/*3", "*2/*2", "*2/*3", "*3/*3"))
        return (f"use a validated pharmacogenetic dosing algorithm including CYP2C9 and VKORC1; "
                f"CYP2C9 was not genotyped ({conds})", "Strong")
    return ("no CYP4F2-specific dose action is recommended by CPIC; dose with a validated "
            "pharmacogenetic algorithm based on CYP2C9 and VKORC1 (CYP4F2*3 carriers may need a "
            "slightly higher dose)", "n/a")


def condense(text: str, limit: int = 300) -> str:
    """Trim CPIC recommendation text to its leading sentence(s), verbatim."""
    t = re.sub(r"\s+", " ", text or "").strip()
    if len(t) <= limit:
        return t
    cut = t[:limit]
    for sep in (". ", "; "):
        if sep in cut:
            return cut[:cut.rfind(sep) + 1].strip()
    return cut.rstrip() + "..."


def lookup_recommendation(case: dict, phenotype: str, lookupkey) -> tuple:
    """(recommendation_action, classification, source, full_text)."""
    gene = GENE_SYMBOL.get(case["gene"], case["gene"])
    drug = case["drug"]
    cpic_drug = DRUG_LOOKUP.get(drug, drug)

    if gene in ("CYP2C9", "VKORC1", "CYP4F2") and drug == "warfarin":
        act, cls = warfarin_recommendation(gene, case["_diplotype"])
        return act, cls, f"{CORPUS_TAG} (CPIC warfarin dose table; no API rows)", act

    ckey = (cpic_drug, phenotype)
    if ckey in CORPUS_RECS:
        act, cls = CORPUS_RECS[ckey]
        return act, cls, f"{CORPUS_TAG} (guideline table; no API rows)", act
    if ckey in NON_CPIC_RECS:
        act, org = NON_CPIC_RECS[ckey]
        return act, "n/a (not CPIC)", f"{CORPUS_TAG} ({org} guideline; CPIC has none)", act

    did = DRUG_ID.get(cpic_drug.lower())
    rows = RECS_BY_DRUG.get(did, [])
    if not rows or lookupkey is None:
        return None, None, None, None

    want_pop, _ = POPULATION.get(cpic_drug, (None, None))
    cands = []
    for r in rows:
        lk = r["lookupkey"] or {}
        if lk.get(gene) != lookupkey:
            continue
        others = {k: v for k, v in lk.items() if k != gene}
        if others and any(v != NO_RESULT for v in others.values()):
            continue                      # keep only rows where untested genes are "No Result"
        if want_pop and r["population"].strip() != want_pop.strip():
            continue
        cands.append(r)
    if not cands:                         # single-gene guideline, or no No-Result row exists
        for r in rows:
            lk = r["lookupkey"] or {}
            if lk.get(gene) == lookupkey and len(lk) == 1:
                if want_pop and r["population"].strip() != want_pop.strip():
                    continue
                cands.append(r)
    if not cands:
        return None, None, None, None
    r = cands[0]
    gname = GUIDELINES.get(r["guidelineid"], {}).get("name", "CPIC")
    gver = GUIDELINES.get(r["guidelineid"], {}).get("version", "?")
    src = f"CPIC API recommendation table: {gname} (v{gver}), population '{r['population'].strip()}'"
    return condense(r["drugrecommendation"]), r["classification"], src, r["drugrecommendation"]


# --------------------------------------------------------------------------------------
# Lethal class: standard dosing in this genotype risks death or a life-threatening event
# --------------------------------------------------------------------------------------
def lethal_class(case: dict, phenotype: str) -> tuple:
    gene, drug, ph = case["gene"], case["drug"], (phenotype or "")
    rules = [
        (gene == "DPYD" and drug in ("fluorouracil", "capecitabine") and "Poor Metaboliz" in ph,
         "lethal", "Near-complete DPD deficiency: standard fluoropyrimidine dosing causes "
                   "early-onset, potentially fatal mucositis, myelosuppression and neurotoxicity."),
        (gene == "DPYD" and drug in ("fluorouracil", "capecitabine") and "Intermediate" in ph,
         "lethal", "Partial DPD deficiency: full-dose fluoropyrimidine carries a documented risk "
                   "of severe, sometimes fatal toxicity; CPIC mandates a 50% dose reduction."),
        (gene == "TPMT" and "Poor Metaboliz" in ph,
         "lethal", "TPMT poor metaboliser: conventional thiopurine doses cause profound, "
                   "potentially fatal myelosuppression."),
        (gene == "NUDT15" and "Poor Metaboliz" in ph,
         "lethal", "NUDT15 poor metaboliser: conventional thiopurine doses cause severe "
                   "leukopenia/myelosuppression."),
        (gene == "CYP2D6" and drug == "codeine" and "Ultrarapid" in ph,
         "lethal", "CYP2D6 ultrarapid metaboliser: excess morphine formation from codeine has "
                   "caused fatal respiratory depression, including in children and breastfed infants."),
        (gene.startswith("HLA-B*57:01") and "positive" in ph,
         "lethal", "HLA-B*57:01 positive: abacavir hypersensitivity syndrome, fatal on rechallenge."),
        (gene.startswith("HLA-B*15:02") and "positive" in ph,
         "lethal", "HLA-B*15:02 positive: high risk of carbamazepine/oxcarbazepine "
                   "Stevens-Johnson syndrome / toxic epidermal necrolysis, which is often fatal."),
        (gene.startswith("HLA-A*31:01") and "positive" in ph,
         "lethal", "HLA-A*31:01 positive: risk of carbamazepine DRESS and SJS/TEN."),
        (gene.startswith("HLA-B*58:01") and "positive" in ph,
         "lethal", "HLA-B*58:01 positive: allopurinol SJS/TEN/DRESS, with mortality up to ~25%."),
        (gene == "G6PD" and drug == "rasburicase" and ph in ("Deficient", "Deficient with CNSHA"),
         "lethal", "G6PD deficiency: rasburicase causes severe acute haemolysis and "
                   "methaemoglobinaemia, which can be fatal; it is contraindicated."),
        (gene == "G6PD" and drug == "primaquine" and ph in ("Deficient", "Deficient with CNSHA"),
         "lethal", "G6PD deficiency: primaquine causes acute haemolytic anaemia that can be "
                   "life-threatening."),
        (gene == "RYR1" and "Malignant Hyperthermia" in ph,
         "lethal", "RYR1 malignant hyperthermia susceptibility: volatile anaesthetics or "
                   "succinylcholine can trigger a fulminant, rapidly fatal MH crisis."),
        (gene == "CYP2C19" and drug == "clopidogrel" and ("Poor Metaboliz" in ph or "Intermediate" in ph),
         "lethal", "Reduced clopidogrel activation after ACS/PCI: markedly increased risk of stent "
                   "thrombosis, myocardial infarction and cardiovascular death."),
        (gene == "CYP2C9" and drug == "warfarin" and "Poor Metaboliz" in ph,
         "lethal", "CYP2C9 poor metaboliser: standard warfarin loading causes over-anticoagulation "
                   "and potentially fatal haemorrhage."),
        (gene == "UGT1A1" and drug == "irinotecan" and "Poor Metaboliz" in ph,
         "lethal", "UGT1A1*28/*28: full-dose irinotecan causes severe neutropenia and diarrhoea "
                   "that can be fatal."),
        (gene == "SLCO1B1" and drug == "simvastatin" and ph in ("Poor Function", "Decreased Function"),
         "lethal", "Impaired OATP1B1 transport: high simvastatin exposure risks myopathy and "
                   "rhabdomyolysis, which can be fatal."),
        (gene == "MT-RNR1" and "increased risk" in ph,
         "severe-irreversible", "m.1555A>G: aminoglycosides cause permanent bilateral hearing "
                                "loss - disabling and irreversible rather than lethal."),
        (gene == "CYP2C19" and drug == "voriconazole" and ph in ("Ultrarapid Metabolizer", "Rapid Metabolizer"),
         "severe-efficacy", "Subtherapeutic voriconazole in invasive fungal infection: treatment "
                            "failure carries high mortality."),
        (gene == "CYP2D6" and drug == "tamoxifen" and ("Poor Metaboliz" in ph or "Intermediate" in ph),
         "severe-efficacy", "Low endoxifen: increased breast-cancer recurrence risk (an efficacy "
                            "failure, not acute toxicity)."),
        (gene == "CYP2C9" and drug == "phenytoin" and "Poor Metaboliz" in ph,
         "severe-toxicity", "Phenytoin accumulation: ataxia, encephalopathy and cardiovascular "
                            "toxicity at standard loading doses."),
    ]
    for hit, cls, why in rules:
        if hit:
            return cls, why
    return "standard", ""


# --------------------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------------------
def cpic_level(case: dict) -> tuple:
    gene = GENE_SYMBOL.get(case["gene"], case["gene"])
    cpic_drug = DRUG_LOOKUP.get(case["drug"], case["drug"])
    did = DRUG_ID.get(cpic_drug.lower())
    for p in PAIRS:
        if p["genesymbol"] == gene and p["drugid"] == did:
            gname = GUIDELINES.get(p["guidelineid"], {}).get("name", "no CPIC guideline")
            gver = GUIDELINES.get(p["guidelineid"], {}).get("version", "n/a")
            return p["cpiclevel"], f"{gname} (v{gver})"
    return "unknown", "no CPIC pair found"


def main():
    rows, audit = [], []
    for case in CASES:
        call = make_diplotype(case)
        case["_diplotype"] = call.diplotype
        level, guideline = cpic_level(case)

        phenotype = lookupkey = ph_src = None
        if not call.abstain and call.diplotype:
            phenotype, lookupkey, ph_src = lookup_phenotype(case["gene"], call.phenotype_key)

        action = cls = rec_src = full = None
        if phenotype:
            case["_phenotype"] = phenotype
            action, cls, rec_src, full = lookup_recommendation(case, phenotype, lookupkey)

        lclass, lwhy = lethal_class(case, phenotype)

        # ---- abstention -------------------------------------------------------------
        abstain = False
        reasons = []
        if call.abstain or not call.diplotype:
            abstain = True
            reasons.append(f"genotype not mappable: {call.rule}")
        elif phenotype is None:
            abstain = True
            reasons.append("no CPIC phenotype assignment found for this diplotype")

        if abstain:
            rec_text = "ABSTAIN: " + "; ".join(reasons)
        else:
            drug = case["drug"]
            if action is None:
                if level == "A" and guideline.startswith("no CPIC guideline"):
                    rec_text = (f"{drug}: no CPIC recommendation exists - CPIC lists "
                                f"{case['gene']}/{drug} as a level A pair but has not published a "
                                f"guideline for it, so no CPIC action can be quoted")
                else:
                    rec_text = (f"{drug}: no CPIC recommendation row matched this phenotype "
                                f"({phenotype})")
            else:
                rec_text = f"{drug}: {action}"
            if level == "Retired":
                rec_text += " [CPIC guideline RETIRED - no longer a level A pair]"

        rows.append({
            "case_id": case["id"],
            "diplotype": call.diplotype or "",
            "phenotype": phenotype or "",
            "recommendation": rec_text,
            "abstained": "true" if abstain else "false",
        })
        audit.append({
            "case_id": case["id"], "gene": case["gene"], "drug": case["drug"],
            "genotype_given": case["genotype"],
            "strand_orientation": call.orientation,
            "diplotype": call.diplotype, "diplotype_rule": call.rule,
            "diplotype_lookup_key": call.phenotype_key,
            "phenotype": phenotype or "", "cpic_lookupkey": lookupkey or "",
            "phenotype_source": ph_src or "",
            "recommendation_action": action or "", "classification": cls or "",
            "recommendation_source": rec_src or "",
            "recommendation_full_text": re.sub(r"\s+", " ", full or ""),
            "cpic_level": level, "cpic_guideline": guideline,
            "lethal_class": lclass, "lethal_rationale": lwhy,
            "abstained": "true" if abstain else "false",
            "caveat": call.caveat,
        })

    def write(fn, fieldnames, data):
        with open(os.path.join(BASE, fn), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(data)
        print(f"wrote {fn}: {len(data)} rows")

    write("claude_science_answers.csv",
          ["case_id", "diplotype", "phenotype", "recommendation", "abstained"], rows)
    write("case_audit.csv", list(audit[0].keys()), audit)

    # ---- table 2a: diplotype -> phenotype --------------------------------------------
    d2p, seen = [], set()
    for a in audit:
        k = (a["gene"], a["diplotype"])
        if not a["diplotype"] or k in seen:
            continue
        seen.add(k)
        d2p.append({"gene": a["gene"], "diplotype": a["diplotype"],
                    "cpic_diplotype_key_used": a["diplotype_lookup_key"],
                    "phenotype": a["phenotype"], "cpic_lookupkey": a["cpic_lookupkey"],
                    "source": a["phenotype_source"]})
    write("table_diplotype_to_phenotype.csv", list(d2p[0].keys()), sorted(
        d2p, key=lambda r: (r["gene"], r["diplotype"])))

    # ---- table 2b: phenotype -> recommendation ---------------------------------------
    p2r, seen = [], set()
    for a in audit:
        k = (a["gene"], a["drug"], a["phenotype"])
        if not a["phenotype"] or k in seen:
            continue
        seen.add(k)
        p2r.append({"gene": a["gene"], "drug": a["drug"], "phenotype": a["phenotype"],
                    "cpic_lookupkey": a["cpic_lookupkey"],
                    "recommendation": a["recommendation_action"],
                    "classification": a["classification"],
                    "cpic_level": a["cpic_level"], "source": a["recommendation_source"]})
    write("table_phenotype_to_recommendation.csv", list(p2r[0].keys()), sorted(
        p2r, key=lambda r: (r["gene"], r["drug"], r["phenotype"])))

    # ---- table 2c: variant -> allele -------------------------------------------------
    v2a = []
    for (gene, rs, base), (name, note) in sorted(VAR2ALLELE.items()):
        loc = [l for l in LOC_BY_RS.get(rs, []) if l["genesymbol"] == gene][0]
        v2a.append({"gene": gene, "rsid": rs, "cpic_variant_name": loc["name"],
                    "chromosome_hgvs": loc["chromosomelocation"], "gene_hgvs": loc["genelocation"],
                    "forward_genomic_base": base, "reference_base": REF_BASE.get((gene, rs), ""),
                    "allele": name, "basis": note,
                    "source": "CPIC API allele_definition + allele_location_value + sequence_location"})
    write("table_variant_to_allele.csv", list(v2a[0].keys()), v2a)

    # ---- sources and versions --------------------------------------------------------
    src = []
    for gene in sorted({GENE_SYMBOL.get(c["gene"], c["gene"]) for c in CASES}):
        alle = [a for a in ALLELES if a["genesymbol"] == gene]
        corpus_gene = [k for k in CORPUS_GENES if GENE_SYMBOL.get(k, k) == gene]
        cg = CORPUS_GENES[corpus_gene[0]] if corpus_gene else None
        pairs_g = [p for p in PAIRS if p["genesymbol"] == gene]
        src.append({
            "gene": gene,
            "cpic_api_allele_table_version": max([a["version"] for a in alle], default="n/a"),
            "cpic_api_alleles_n": len(alle),
            "cpic_guidelines": "; ".join(sorted({
                f"{GUIDELINES[p['guidelineid']]['name']} (v{GUIDELINES[p['guidelineid']]['version']})"
                for p in pairs_g if p["guidelineid"] in GUIDELINES})) or "none",
            "cpic_levels_for_case_drugs": "; ".join(sorted({
                f"{drugmap_name(p['drugid'])}={p['cpiclevel']}" for p in pairs_g
                if drugmap_name(p["drugid"]) in case_drugs_for(gene)})),
            "corpus_source": "; ".join(s["name"] for s in cg["sources"]) if cg else "not in corpus",
            "corpus_sha256": cg["sha256"] if cg else "",
            "corpus_retrieved": "; ".join(sorted({s["retrieved_date"] for s in cg["sources"]})) if cg else "",
        })
    write("sources_and_versions.csv", list(src[0].keys()), src)

    n_abs = sum(1 for r in rows if r["abstained"] == "true")
    print(f"\n{len(rows)} cases; {n_abs} abstentions")
    for a in audit:
        if a["lethal_class"] == "lethal":
            print(f"  LETHAL-CLASS {a['case_id']} {a['gene']}/{a['drug']}: {a['phenotype']}")
    return rows, audit


_DRUG_NAME = {}


def drugmap_name(did):
    if not _DRUG_NAME:
        for d in DRUGS:
            _DRUG_NAME[d["drugid"]] = d["name"]
    return _DRUG_NAME.get(did, did)


def case_drugs_for(gene):
    out = set()
    for c in CASES:
        if GENE_SYMBOL.get(c["gene"], c["gene"]) == gene:
            out.add(DRUG_LOOKUP.get(c["drug"], c["drug"]))
    return out


if __name__ == "__main__":
    pass
    main()
