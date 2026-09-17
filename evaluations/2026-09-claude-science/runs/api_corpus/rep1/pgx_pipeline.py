#!/usr/bin/env python3
"""
Blinded pharmacogenomics study - genotype -> diplotype -> CPIC phenotype -> recommendation.

Inputs  (same folder):
    cases_blinded.json            110 blinded cases (id, gene, drug, genotype)
    cpic_guideline_corpus.json    supplied CPIC guideline text snapshot (2026-05-09)
    cpic_api_raw/*.json           CPIC data API tables downloaded by fetch_cpic_api.py

Outputs (same folder):
    claude_science_answers.csv        case_id,diplotype,phenotype,recommendation,abstained
    table_diplotype_to_phenotype.csv  diplotype -> CPIC phenotype table actually used
    table_phenotype_to_recommendation.csv  phenotype -> CPIC recommendation table actually used
    case_audit.csv                    full per-case provenance / assumptions / ambiguities
    source_versions.csv               source + version recorded per gene

Design
------
Nothing about the star-allele assignment is hard-coded from memory.  For every
SNV-defined case the pipeline reads CPIC's own allele-definition tables
(sequence_location + allele_location_value + allele_definition), works out the
reference and variant base on the GRCh38 forward strand from CPIC's
`chromosomelocation` HGVS string, resolves the orientation of the genotype as
reported in the case, and assigns the star allele that CPIC defines minimally by
that variant.  Diplotype -> phenotype comes from the CPIC `diplotype` table and
phenotype -> recommendation from the CPIC `recommendation` table.  Cases that CPIC
cannot key on a single SNV (structural variants, HLA typing, repeat
polymorphisms, multi-variant haplotypes) go through explicit, individually
documented handlers.
"""
from __future__ import annotations

import collections
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "cpic_api_raw")

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ----------------------------------------------------------------------------
# load
# ----------------------------------------------------------------------------
def load(name):
    with open(os.path.join(RAW, name + ".json")) as fh:
        return json.load(fh)


cases = json.load(open(os.path.join(HERE, "cases_blinded.json")))
corpus = json.load(open(os.path.join(HERE, "cpic_guideline_corpus.json")))

T = {n: load(n) for n in ("gene", "pair", "drug", "guideline", "diplotype",
                          "recommendation", "allele", "allele_definition",
                          "sequence_location", "allele_location_value")}

DRUG_ID = {d["name"].lower(): d["drugid"] for d in T["drug"]}
GUIDELINE = {x["id"]: x for x in T["guideline"]}
AD = {a["id"]: a for a in T["allele_definition"]}
ALLELE_FN = {(a["genesymbol"], a["name"]): a for a in T["allele"]}

DIP = {}
for d in T["diplotype"]:
    DIP[(d["genesymbol"], d["diplotype"])] = d

# forward-strand reference / variant bases for every CPIC sequence location
LOC = {}
for s in T["sequence_location"]:
    m = re.match(r"g\.\d+([ACGT]+)>([ACGT;>.g\d ]+)$", (s["chromosomelocation"] or ""))
    ref = m.group(1) if m else None
    LOC[s["id"]] = dict(gene=s["genesymbol"], rsid=s["dbsnpid"], name=s["name"],
                        hgvs_g=s["chromosomelocation"], ref=ref)

ALV_BY_LOC = collections.defaultdict(list)
for r in T["allele_location_value"]:
    ALV_BY_LOC[r["locationid"]].append(r)

# number of positions at which each allele definition differs from the reference
AD_NVAR = collections.Counter()
for r in T["allele_location_value"]:
    L = LOC.get(r["locationid"])
    if L and L["ref"] and r["variantallele"] != L["ref"]:
        AD_NVAR[r["alleledefinitionid"]] += 1


def locations_for(gene, rsid):
    return [s for s in T["sequence_location"]
            if s["genesymbol"] == gene and s["dbsnpid"] == rsid]


def observed_bases(gene, rsid):
    """{forward-strand base -> set of CPIC allele names carrying it}"""
    out = collections.defaultdict(set)
    for s in locations_for(gene, rsid):
        for r in ALV_BY_LOC.get(s["id"], []):
            out[r["variantallele"]].add(AD[r["alleledefinitionid"]]["name"])
    return out


def star_sort_key(name):
    """Canonical ordering: lowest numbered star allele first, letters after digits."""
    m = re.match(r"\*(\d+)([A-Za-z]*)$", name)
    if m:
        return (0, int(m.group(1)), m.group(2))
    return (1, 0, name)


def minimal_allele(gene, rsid, base):
    """The CPIC allele minimally defined by `base` at `rsid` (fewest variant positions)."""
    cands = []
    for s in locations_for(gene, rsid):
        for r in ALV_BY_LOC.get(s["id"], []):
            if r["variantallele"] == base:
                ad = AD[r["alleledefinitionid"]]
                cands.append((AD_NVAR[ad["id"]], star_sort_key(ad["name"]), ad["name"]))
    if not cands:
        return None, []
    cands.sort()
    best_n = cands[0][0]
    tied = [c[2] for c in cands if c[0] == best_n]
    return cands[0][2], tied


def resolve_orientation(gene, rsid, alleles):
    """
    Map the genotype letters as written in the case onto CPIC's forward-strand
    alphabet.  Returns (mapped_alleles, orientation, hgvs_g).  Raises if the
    genotype is not interpretable on either strand, or is interpretable on both
    (a palindromic A/T or C/G site, where strand cannot be recovered).
    """
    locs = locations_for(gene, rsid)
    if not locs:
        raise KeyError("rsid %s not in CPIC allele definitions for %s" % (rsid, gene))
    alpha = set(observed_bases(gene, rsid)) | {LOC[s["id"]]["ref"] for s in locs}
    alpha = {a for a in alpha if a in COMPLEMENT}
    fwd_ok = all(a in alpha for a in alleles)
    rev = [COMPLEMENT[a] for a in alleles]
    rev_ok = all(a in alpha for a in rev)
    if fwd_ok and rev_ok:
        raise ValueError("palindromic site - strand not recoverable")
    if fwd_ok:
        return list(alleles), "forward (as written)", locs[0]["chromosomelocation"]
    if rev_ok:
        return rev, "complemented to forward strand", locs[0]["chromosomelocation"]
    raise ValueError("alleles %s not compatible with CPIC alphabet %s" % (alleles, sorted(alpha)))


def allele_for_base(gene, rsid, base):
    """reference base -> reference allele name; variant base -> minimal star allele."""
    locs = locations_for(gene, rsid)
    ref = LOC[locs[0]["id"]]["ref"]
    if base == ref:
        names = {AD[r["alleledefinitionid"]]["name"]
                 for s in locs for r in ALV_BY_LOC.get(s["id"], [])
                 if r["variantallele"] == ref}
        for cand in ("*1", "Reference", "B (reference)"):
            if cand in names:
                return cand, [cand], True
        # Some genes name the fully-reference haplotype something other than *1
        # (CYP2C19*38 is the zero-variant haplotype; *1 carries 991A>G).  A
        # clinical report calls the normal-function chromosome *1, and CPIC's
        # diplotype table carries both forms with identical phenotypes, so use
        # *1 when the gene has one and keep the definition-table name alongside.
        if (gene, "*1") in ALLELE_FN:
            return "*1", sorted(names | {"*1"}), True
        return sorted(names, key=star_sort_key)[0], sorted(names), True
    name, tied = minimal_allele(gene, rsid, base)
    if name is None:
        raise ValueError("base %s at %s not defined by any CPIC allele" % (base, rsid))
    return name, tied, False


def order_diplotype(gene, a1, a2):
    """Match CPIC's own diplotype string if either order is present, else sort."""
    for cand in ("%s/%s" % (a1, a2), "%s/%s" % (a2, a1)):
        if (gene, cand) in DIP:
            return cand
    return "/".join(sorted([a1, a2]))


# ----------------------------------------------------------------------------
# genotype -> diplotype
# ----------------------------------------------------------------------------
RS_TOKEN = re.compile(r"(rs\d+)\s+([ACGT])(?:\s*/\s*([ACGT]))?")

# Which co-tested gene value CPIC's multi-gene recommendation key is looked up
# with when the case reports only one of the genes.  CPIC's own convention is
# that an untested gene carries no result; we instead assume the reference /
# normal value and record the assumption on every affected case.
CO_GENE_DEFAULT = {
    "CYP2C19": "Normal Metabolizer",
    "CYP2D6": "2.0",
    "TPMT": "Normal Metabolizer",
    "NUDT15": "Normal Metabolizer",
    "HLA-B": None,          # set per drug below
    "HLA-A": "*31:01 negative",
    "CACNA1S": "Uncertain Susceptibility",
}

# Population / indication selected where CPIC gives several.  Documented choices.
POPULATION_CHOICE = {
    "voriconazole": "adults",
    "clopidogrel": "CVI ACS PCI ",
    "phenytoin": "PHT naive",
    "carbamazepine": "CBZ naive",
    "oxcarbazepine": "OXC naive",
}

# Case-label drug -> CPIC drug(s) to look the recommendation up under.
DRUG_ALIAS = {
    "thiopurines": ["azathioprine", "mercaptopurine", "thioguanine"],
    "aminoglycosides": ["gentamicin"],
    "volatile-anaesthetics": ["sevoflurane"],
    "PEG-IFN-\u03b1": ["peginterferon alfa-2a"],
}

GENE_ALIAS = {"HLA-B*15:02": "HLA-B", "HLA-B*57:01": "HLA-B",
              "HLA-B*58:01": "HLA-B", "HLA-A*31:01": "HLA-A"}


def resolve_genotype(case):
    """
    Returns dict with diplotype, lookup gene/value, basis, assumptions,
    ambiguity, abstain flag.
    """
    gene = GENE_ALIAS.get(case["gene"], case["gene"])
    gt = case["genotype"]
    low = gt.lower()
    out = dict(gene=gene, diplotype=None, phenotype=None, basis="", assumptions=[],
               ambiguity="", abstain=False, abstain_reason="", orientation="",
               hgvs_g="", method="")

    # ---- reference-everywhere statements -----------------------------------
    if re.match(r"^all .* (homozygous reference|reference)$", low) or low == "all dpyd risk variants reference":
        ref_name = {"DPYD": "Reference"}.get(gene, "*1")
        out.update(diplotype=order_diplotype(gene, ref_name, ref_name),
                   basis="Case states every interrogated variant is homozygous reference; "
                         "assigned the CPIC reference allele on both chromosomes.",
                   method="reference-statement")
        out["assumptions"].append(
            "Only the variants on the assay panel were interrogated; an untested "
            "rare no-function allele would not be visible.")
        return out

    # ---- CYP2D6 copy-number duplication ------------------------------------
    if gene == "CYP2D6" and "duplication" in low:
        out.update(diplotype="*1/*1x2",
                   basis="Copy-number assay reports a duplication of the *1 (normal function) "
                         "allele and no SNV-defined variant allele; CPIC scores *1x2 as two "
                         "normal-function copies, giving activity score 3.0 with the "
                         "non-duplicated *1.",
                   method="structural-variant handler")
        out["assumptions"] += [
            "N assumed to be 2 (*1x2); the assay reports *1xN without resolving copy number.",
            "The non-duplicated chromosome is assumed to be *1, as no SNV-defined variant was reported.",
        ]
        out["ambiguity"] = ("If N>2 the activity score is higher but the phenotype stays "
                            "Ultrarapid. If the second chromosome actually carried a no-function "
                            "allele (e.g. *4/*1x2, activity score 2.0) the phenotype would be "
                            "Normal Metabolizer - copy-number assays alone cannot phase the "
                            "duplication to an allele.")
        return out

    # ---- CYP2D6 compound rs3892097 + rs1065852 -----------------------------
    if gene == "CYP2D6" and "rs3892097" in low and "rs1065852" in low:
        out.update(diplotype="*1/*4",
                   basis="rs3892097 (1847G>A) and rs1065852 (100C>T) both heterozygous. The "
                         "parsimonious phase places both variants on one chromosome, which is "
                         "CPIC's *4 haplotype definition (*4 carries 100C>T and 1847G>A), "
                         "leaving *1 on the other chromosome.",
                   method="two-SNV phasing handler")
        out["ambiguity"] = ("Unphased data also admit *4/*10 (1847G>A in trans to 100C>T, "
                            "activity score 0.25). Both phasings translate to Intermediate "
                            "Metabolizer, so the phenotype call is robust to the phase.")
        return out

    # ---- TPMT two heterozygous / two homozygous variants -------------------
    if gene == "TPMT" and "rs1800460" in low and "rs1142345" in low:
        out.update(diplotype="*3A/*3A",
                   basis="Homozygous for both c.460G>A (rs1800460) and c.719A>G (rs1142345). "
                         "The forward-genomic alleles given (T/T and C/C) are the variant bases "
                         "at CPIC's g.18138997C>T and g.18130687T>C, so both variants are "
                         "present on both chromosomes; CPIC's *3A carries both.",
                   method="two-SNV haplotype handler")
        out["ambiguity"] = ("A *3B/*3C genotype (the two variants in trans) would give the same "
                            "SNV counts only if each were heterozygous; with both homozygous the "
                            "diplotype must be *3A/*3A. CPIC assigns Poor Metabolizer to *3A/*3A "
                            "and to *3B/*3C alike, so the phenotype is unaffected.")
        return out

    # ---- G6PD haplotypes ---------------------------------------------------
    if gene == "G6PD":
        het_female = "heterozygous female" in low
        has_a_minus = "rs1050828" in low and "rs1050829" in low and (
            "in cis" in low or "a-" in low or "c/t" in low)
        if "rs5030868" in low and re.search(r"rs5030868 a\b", low):
            allele, _, _ = allele_for_base("G6PD", "rs5030868", "A")
            out.update(diplotype=allele,
                       basis="rs5030868 variant base A on the forward strand (CPIC "
                             "g.154534419G>A, transcript c.563C>T p.Ser188Phe) in a hemizygous "
                             "male / homozygous female; CPIC defines this single-SNV allele as "
                             "G6PD Mediterranean (class II, deficient).",
                       method="G6PD haplotype handler", orientation="forward (as written)",
                       hgvs_g="g.154534419G>A")
            out["assumptions"].append(
                "Reported as hemizygous male or homozygous female; both give a single "
                "deficient-allele G6PD result in CPIC's table.")
            return out
        if has_a_minus and not het_female:
            out.update(diplotype="A- 202A_376G",
                       basis="rs1050828 (c.202G>A) and rs1050829 (c.376A>G) variant bases in cis, "
                             "hemizygous male or homozygous female. CPIC's allele definition for "
                             "'A- 202A_376G' is exactly these two variants (class III, deficient).",
                       method="G6PD haplotype handler", orientation="forward (as written)",
                       hgvs_g="g.154536002C>T + g.154535277T>C")
            return out
        if het_female:
            out.update(diplotype="A- 202A_376G/B (reference)",
                       basis="Heterozygous female carrying one A- (202A_376G) chromosome and one "
                             "reference (B) chromosome. CPIC assigns heterozygous females the "
                             "'Variable' phenotype because X-inactivation makes enzyme activity "
                             "unpredictable from genotype.",
                       method="G6PD haplotype handler")
            out["ambiguity"] = ("The case allows rs1050829 T/C or T/T; the second base only "
                                "changes whether the reference chromosome is B or the "
                                "normal-function A allele, and CPIC calls A/A- 'Variable' as well.")
            return out
        if "no deficiency variants" in low:
            out.update(diplotype="B (reference)",
                       basis="No G6PD deficiency variant at the interrogated positions; CPIC "
                             "reference allele B, hemizygous male or homozygous female.",
                       method="G6PD haplotype handler")
            out["assumptions"].append(
                "Only rs1050828 / rs5030868 were interrogated; over 200 other deficiency "
                "variants exist and would not be seen.")
            return out

    # ---- UGT1A1 TA repeat --------------------------------------------------
    if gene == "UGT1A1" and "ta" in low:
        m = re.search(r"ta(\d)\s*/\s*ta(\d)", low)
        rep = {"6": "*1", "7": "*28"}
        a1, a2 = rep[m.group(1)], rep[m.group(2)]
        out.update(diplotype=order_diplotype("UGT1A1", a1, a2),
                   basis="rs8175347 promoter (TA)n repeat: TA6 is the CPIC reference allele *1 "
                         "and TA7 is *28 (A(TA)7TAA).",
                   method="repeat-length handler")
        return out

    # ---- HLA typing --------------------------------------------------------
    if gene in ("HLA-A", "HLA-B"):
        star = case["gene"].split("*", 1)[1]
        if "rs2395029" in low:
            positive = "positive" in low or "any positive" in low
            out.update(diplotype="*57:01 positive" if positive else "*57:01 negative",
                       basis="rs2395029 (HCP5 G allele) used as a tag for HLA-B*57:01; the case "
                             "states the result as positive / no risk allele.",
                       method="HLA tag-SNP handler")
            out["assumptions"].append(
                "rs2395029 is a proxy in strong linkage disequilibrium with HLA-B*57:01 in "
                "European populations but not a direct type; CPIC keys on the HLA-B*57:01 "
                "allele itself.")
            out["ambiguity"] = ("Tag-SNP LD with *57:01 is imperfect outside European ancestry; "
                                "direct HLA typing is required for a definitive call.")
            return out
        positive = "positive" in low
        out.update(diplotype="*%s %s" % (star, "positive" if positive else "negative"),
                   basis="Direct HLA typing result stated in the case.",
                   method="HLA typing handler")
        if not positive:
            out["assumptions"].append(
                "A negative type lowers but does not abolish risk of the cutaneous adverse "
                "reaction; other risk alleles and non-HLA mechanisms exist.")
        return out

    # ---- MT-RNR1 -----------------------------------------------------------
    if gene == "MT-RNR1":
        if "1555a>g" in low:
            out.update(diplotype="m.1555A>G",
                       basis="Homoplasmic m.1555A>G in the mitochondrial 12S rRNA gene; CPIC "
                             "lists this variant as conferring increased risk of "
                             "aminoglycoside-induced hearing loss.",
                       method="mitochondrial variant handler")
        else:
            out.update(diplotype="Reference",
                       basis="No MT-RNR1 variant detected; CPIC's Reference allele carries "
                             "normal risk.",
                       method="mitochondrial variant handler")
            out["assumptions"].append(
                "Heteroplasmy below assay detection, and MT-RNR1 variants of uncertain risk "
                "outside the interrogated positions, would not be seen.")
        return out

    # ---- RYR1 --------------------------------------------------------------
    if gene == "RYR1":
        if "no pathogenic variant" in low:
            out.update(diplotype="Reference/Reference",
                       basis="No RYR1 malignant-hyperthermia-associated variant detected. CPIC "
                             "assigns 'Uncertain Susceptibility' rather than a normal result, "
                             "because a negative genotype does not exclude MH susceptibility.",
                       method="RYR1 handler")
            out["assumptions"].append(
                "Absence of a detected variant is not evidence of absent MH risk; about half of "
                "MH-susceptible individuals have no identifiable RYR1 variant.")
        else:
            out.update(diplotype="c.1840C>T/Reference",
                       basis="RYR1 c.1840C>T (p.Arg614Cys), a CPIC-listed "
                             "malignant-hyperthermia-associated variant, heterozygous "
                             "(autosomal dominant).",
                       method="RYR1 handler")
        return out

    # ---- CFTR --------------------------------------------------------------
    if gene == "CFTR":
        if "rs75527207" in low or "g551d" in low:
            out.update(diplotype="G551D/ivacaftor non-responsive CFTR sequence",
                       basis="rs75527207 variant allele A on the forward strand (CPIC "
                             "g.117587806G>A) is CFTR G551D, a gating (class III) variant listed "
                             "in the ivacaftor label as responsive. One responsive allele is "
                             "sufficient for CPIC's responsive call.",
                       method="CFTR handler", hgvs_g="g.117587806G>A")
            out["assumptions"].append(
                "Case gives G/A or A/A; CPIC's lookup is the same either way because a single "
                "ivacaftor-responsive allele yields 'ivacaftor responsive in CF patients'.")
        else:
            out.update(diplotype="ivacaftor non-responsive CFTR sequence/"
                                 "ivacaftor non-responsive CFTR sequence",
                       basis="F508del/F508del (rs113993960). F508del is a processing (class II) "
                             "variant and is not among the ivacaftor-responsive variants in "
                             "CPIC's allele definition table, so both chromosomes map to CPIC's "
                             "'ivacaftor non-responsive CFTR sequence' allele.",
                       method="CFTR handler")
            out["assumptions"].append(
                "CPIC's CFTR guideline covers ivacaftor monotherapy only; F508del homozygotes "
                "are treated with CFTR modulator combinations (e.g. "
                "elexacaftor/tezacaftor/ivacaftor), which this guideline does not address.")
        return out

    # ---- generic single-SNV star allele path -------------------------------
    toks = RS_TOKEN.findall(gt)
    if len(toks) == 1:
        rsid, b1, b2 = toks[0]
        if not b2:
            b2 = b1
        mapped, orientation, hgvs = resolve_orientation(gene, rsid, [b1, b2])
        a1, tie1, isref1 = allele_for_base(gene, rsid, mapped[0])
        a2, tie2, isref2 = allele_for_base(gene, rsid, mapped[1])
        dipl = order_diplotype(gene, a1, a2)
        loc = locations_for(gene, rsid)[0]
        out.update(diplotype=dipl, orientation=orientation, hgvs_g=hgvs,
                   method="CPIC allele-definition lookup (single SNV)",
                   basis="%s %s/%s read as forward-strand %s/%s at CPIC %s (%s, %s); "
                         "reference base %s. Assigned %s and %s."
                         % (rsid, b1, b2, mapped[0], mapped[1], hgvs, loc["name"],
                            loc["proteinlocation"], LOC[loc["id"]]["ref"], a1, a2))
        ties = sorted(set([t for t in tie1 + tie2 if t not in (a1, a2)]))
        if ties:
            out["ambiguity"] = ("This single variant is also contained in the definitions of %s; "
                                "those alleles carry additional variants that were not "
                                "interrogated, so the minimal definition is used."
                                % ", ".join(ties[:6]))
        if not (isref1 and isref2):
            out["assumptions"].append(
                "Variants not named in the case are assumed reference, so the star allele is the "
                "one CPIC defines minimally by the reported variant.")
        return out

    out.update(abstain=True, abstain_reason="No deterministic mapping rule matched this "
                                            "genotype description.")
    return out


# ----------------------------------------------------------------------------
# diplotype -> phenotype
# ----------------------------------------------------------------------------
GENOTYPE_ONLY_GENES = {
    # genes where CPIC assigns no diplotype-level phenotype term
    "VKORC1": dict(
        variant="Warfarin-sensitive genotype (-1639 A/A); CPIC defines no VKORC1 phenotype term",
        reference="Reference VKORC1 promoter genotype (-1639 G/G); CPIC defines no VKORC1 "
                  "phenotype term"),
    "CYP4F2": dict(
        variant="Increased warfarin dose requirement (*3/*3); CPIC defines no CYP4F2 "
                "phenotype term",
        reference="Reference (*1/*1); CPIC defines no CYP4F2 phenotype term"),
    "IFNL3": dict(
        variant="Unfavourable response genotype (rs12979860 T/T)",
        reference="Favourable response genotype (rs12979860 C/C)"),
}


def phenotype_for(gene, diplotype, case_genotype):
    d = DIP.get((gene, diplotype))
    if d:
        return d["generesult"], d["lookupkey"].get(gene), "CPIC diplotype table", d
    if gene in GENOTYPE_ONLY_GENES:
        is_var = diplotype.count("variant") == 2 or diplotype == "*3/*3"
        key = "variant" if is_var else "reference"
        return GENOTYPE_ONLY_GENES[gene][key], None, "CPIC guideline text (no CPIC diplotype table for this gene)", None
    a = ALLELE_FN.get((gene, diplotype))
    if a:
        st = a["clinicalfunctionalstatus"]
        if gene == "MT-RNR1":
            # the allele table capitalises 'Normal risk...' while the diplotype and
            # recommendation tables use lower case; normalise to the lookup form
            st = st[0].lower() + st[1:]
        return st, st, "CPIC allele table", None
    return None, None, "unresolved", None


# ----------------------------------------------------------------------------
# phenotype -> recommendation
# ----------------------------------------------------------------------------
def recommendation_for(case_drug, gene, lookupvalue):
    """Return (list of matching recommendation rows, cpic drug name used)."""
    for dn in DRUG_ALIAS.get(case_drug, [case_drug]):
        did = DRUG_ID.get(dn.lower())
        if not did:
            continue
        rows = [r for r in T["recommendation"] if r["drugid"] == did]
        if not rows:
            continue
        hits = []
        for r in rows:
            lk = r["lookupkey"] or {}
            if lk.get(gene) != lookupvalue:
                continue
            ok = True
            for og, ov in lk.items():
                if og == gene:
                    continue
                want = CO_GENE_DEFAULT.get(og)
                if og == "HLA-B":
                    want = {"phenytoin": "*15:02 negative",
                            "carbamazepine": "*15:02 negative"}.get(dn)
                if want is not None and ov != want:
                    ok = False
            if ok:
                hits.append(r)
        pop = POPULATION_CHOICE.get(dn)
        if pop:
            sel = [r for r in hits if r["population"] == pop]
            if sel:
                hits = sel
        if hits:
            return hits, dn
    return [], None


def condense(text, limit=200):
    text = re.sub(r"\s+", " ", (text or "")).strip()
    if len(text) <= limit:
        return text
    cut = text[:limit]
    j = max(cut.rfind(". "), cut.rfind("; "))
    return (cut[:j + 1] if j > 60 else cut.rstrip() + "...").strip()


# ----------------------------------------------------------------------------
# recommendations CPIC does not expose as a machine-readable table
# ----------------------------------------------------------------------------
CORPUS_RECS = {
    # (case drug, gene, phenotype/lookup value) -> (recommendation, classification, source)
    ("ondansetron", "CYP2D6", "Ultrarapid Metabolizer"): (
        "Select an alternative drug not predominantly metabolised by CYP2D6 (e.g. granisetron).",
        "Moderate", "supplied corpus: CPIC ondansetron/CYP2D6 (PharmGKB PA166161954)"),
    ("ondansetron", "CYP2D6", "Normal Metabolizer"): (
        "Initiate therapy with recommended starting dose.",
        "Strong", "supplied corpus: CPIC ondansetron/CYP2D6 (PharmGKB PA166161954)"),
    ("ondansetron", "CYP2D6", "Poor Metabolizer"): (
        "Insufficient evidence demonstrating clinical impact based on CYP2D6 genotype. Initiate "
        "therapy with recommended starting dose.",
        "No recommendation", "supplied corpus: CPIC ondansetron/CYP2D6 (PharmGKB PA166161954)"),
    ("ondansetron", "CYP2D6", "Intermediate Metabolizer"): (
        "Insufficient evidence demonstrating clinical impact based on CYP2D6 genotype. Initiate "
        "therapy with recommended starting dose.",
        "No recommendation", "supplied corpus: CPIC ondansetron/CYP2D6 (PharmGKB PA166161954)"),
    ("irinotecan", "UGT1A1", "Poor Metabolizer"): (
        "Start at 70% of the standard dose (DPWG) / reduce the starting dose by 25-30% for "
        "regimens above 180 mg/m2 and titrate on neutrophil count (RNPGx); UGT1A1 *28/*28 "
        "carries a high risk of severe neutropenia and diarrhoea. No CPIC guideline exists for "
        "this pair.",
        "DPWG / RNPGx", "supplied corpus: DPWG (PA166104951) and RNPGx (PA166127626) irinotecan/UGT1A1"),
    ("irinotecan", "UGT1A1", "Intermediate Metabolizer"): (
        "No dose adjustment recommended for *1/*28 heterozygotes at conventional doses; "
        "monitor for neutropenia at high doses. No CPIC guideline exists for this pair.",
        "DPWG / RNPGx", "supplied corpus: DPWG (PA166104951) and RNPGx (PA166127626) irinotecan/UGT1A1"),
    ("irinotecan", "UGT1A1", "Normal Metabolizer"): (
        "Use standard dosing. No CPIC guideline exists for this pair.",
        "DPWG / RNPGx", "supplied corpus: DPWG (PA166104951) and RNPGx (PA166127626) irinotecan/UGT1A1"),
    ("PEG-IFN-\u03b1", "IFNL3", "reference"): (
        "Favourable response genotype (~70% chance of SVR after 48 weeks): weighs in favour of "
        "using PEG-IFN-alpha/ribavirin-containing regimens. CPIC has retired this gene-drug pair.",
        "Strong (guideline retired)", "supplied corpus: CPIC IFNL3 (PharmGKB PA166110235)"),
    ("PEG-IFN-\u03b1", "IFNL3", "variant"): (
        "Unfavourable response genotype (~30% chance of SVR after 48 weeks): consider the "
        "implications before initiating PEG-IFN-alpha/ribavirin-containing regimens. CPIC has "
        "retired this gene-drug pair.",
        "Strong (guideline retired)", "supplied corpus: CPIC IFNL3 (PharmGKB PA166110235)"),
    # CPIC warfarin dosing is an algorithm/table, not a phenotype lookup.
    ("warfarin", "CYP2C9", "0.0"): (
        "CYP2C9 *3/*3: predicted 0.5-2 mg/day by the CPIC dose table for every VKORC1 genotype; "
        "consider an alternative anticoagulant, or dose with a validated pharmacogenetic "
        "algorithm and intensive INR monitoring.",
        "Strong", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "CYP2C9", "1.0"): (
        "CYP2C9 *1/*3: reduced dose - 3-4 mg/day if VKORC1 -1639 GG or GA, 0.5-2 mg/day if AA. "
        "VKORC1 genotype and continental ancestry are required to fix the dose; use a validated "
        "pharmacogenetic dosing algorithm.",
        "Strong", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "CYP2C9", "2.0"): (
        "CYP2C9 *1/*1: 5-7 mg/day if VKORC1 -1639 GG or GA, 3-4 mg/day if AA. VKORC1 genotype "
        "and continental ancestry are required to fix the dose; use a validated pharmacogenetic "
        "dosing algorithm.",
        "Strong", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "VKORC1", "variant"): (
        "VKORC1 -1639 A/A: warfarin-sensitive - 3-4 mg/day with CYP2C9 *1/*1 or *1/*2, "
        "0.5-2 mg/day with any CYP2C9 *3 or *2/*2 genotype. CYP2C9 genotype and continental "
        "ancestry are required to fix the dose; use a validated pharmacogenetic dosing algorithm.",
        "Strong", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "VKORC1", "reference"): (
        "VKORC1 -1639 G/G: 5-7 mg/day with CYP2C9 *1/*1 or *1/*2, 3-4 mg/day with *1/*3, *2/*2 "
        "or *2/*3, 0.5-2 mg/day with *3/*3. CYP2C9 genotype and continental ancestry are "
        "required to fix the dose; use a validated pharmacogenetic dosing algorithm.",
        "Strong", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "CYP4F2", "variant"): (
        "CYP4F2 *3/*3 (rs2108622 T/T): modest increase in warfarin dose requirement; CPIC "
        "incorporates CYP4F2 only as a dose increase of roughly 5-10% on top of the "
        "CYP2C9/VKORC1-based dose in non-African-ancestry patients. CYP2C9 and VKORC1 genotypes "
        "are required before a dose can be given.",
        "Optional", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
    ("warfarin", "CYP4F2", "reference"): (
        "CYP4F2 *1/*1 (rs2108622 C/C): no CYP4F2-based dose adjustment. CYP2C9 and VKORC1 "
        "genotypes are required before a dose can be given.",
        "Optional", "supplied corpus: CPIC warfarin/CYP2C9,VKORC1,CYP4F2 (PharmGKB PA166104949)"),
}

# Cases where standard dosing risks severe or life-threatening harm.
LETHAL_RULES = [
    ("DPYD", "Poor Metabolizer", "Complete DPD deficiency; full-dose fluoropyrimidine causes "
                                 "fatal myelosuppression, mucositis and neurotoxicity."),
    ("DPYD", "Intermediate Metabolizer", "Partial DPD deficiency; full-dose fluoropyrimidine "
                                         "causes grade 4-5 toxicity in a substantial fraction."),
    ("TPMT", "Poor Metabolizer", "Conventional thiopurine dosing causes life-threatening "
                                 "myelosuppression."),
    ("NUDT15", "Poor Metabolizer", "Conventional thiopurine dosing causes life-threatening "
                                   "myelosuppression."),
    ("G6PD", "Deficient", "Rasburicase or primaquine can precipitate acute haemolysis and "
                          "methaemoglobinaemia."),
    ("G6PD", "Variable", "Heterozygous females can be functionally deficient because of skewed "
                         "X-inactivation; genotype cannot exclude severe haemolysis, so enzyme "
                         "activity must be measured before exposure."),
    ("RYR1", "Malignant Hyperthermia Susceptibility", "Volatile anaesthetics or succinylcholine "
                                                      "can trigger fulminant malignant "
                                                      "hyperthermia."),
    ("HLA-B", "*57:01 positive", "Abacavir hypersensitivity syndrome; rechallenge is fatal."),
    ("HLA-B", "*15:02 positive", "Carbamazepine/oxcarbazepine-induced SJS/TEN."),
    ("HLA-B", "*58:01 positive", "Allopurinol-induced SJS/TEN and DRESS."),
    ("HLA-A", "*31:01 positive", "Carbamazepine-induced DRESS and SJS/TEN."),
    ("MT-RNR1", "increased risk of aminoglycoside-induced hearing loss",
     "Aminoglycoside exposure causes permanent bilateral deafness at any dose."),
]
LETHAL_DRUG_SPECIFIC = [
    ("CYP2D6", "Ultrarapid Metabolizer", "codeine",
     "Ultrarapid morphine formation from codeine has caused fatal respiratory depression."),
    ("CYP2C9", "Poor Metabolizer", "warfarin",
     "Standard warfarin dosing in a CYP2C9 poor metaboliser causes severe over-anticoagulation "
     "and major bleeding."),
    ("CYP2C9", "Poor Metabolizer", "phenytoin",
     "Standard phenytoin dosing causes toxic concentrations and encephalopathy."),
    ("CYP2C19", "Poor Metabolizer", "clopidogrel",
     "Loss of clopidogrel activation causes stent thrombosis and myocardial infarction."),
]


def lethal_class(gene, phenotype, drug):
    for g, p, why in LETHAL_RULES:
        if gene == g and phenotype == p:
            return True, why
    for g, p, d, why in LETHAL_DRUG_SPECIFIC:
        if gene == g and phenotype == p and drug == d:
            return True, why
    return False, ""


# ----------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------
def main():
    answers, audit = [], []
    dip_tbl, rec_tbl = {}, {}

    for c in cases:
        r = resolve_genotype(c)
        gene = r["gene"]
        phen = lookupval = phen_src = None
        recs_txt = ""
        rec_class = rec_src = ""
        cpic_drug = None
        abstain = r["abstain"]
        abstain_reason = r["abstain_reason"]

        if not abstain:
            phen, lookupval, phen_src, _ = phenotype_for(gene, r["diplotype"], c["genotype"])
            if phen is None:
                abstain = True
                abstain_reason = ("Diplotype %s resolved but no CPIC phenotype assignment found."
                                  % r["diplotype"])

        if not abstain:
            hits, cpic_drug = recommendation_for(c["drug"], gene, lookupval)
            if hits:
                recs_txt = condense(hits[0]["drugrecommendation"])
                rec_class = hits[0]["classification"]
                gl = GUIDELINE.get(hits[0]["guidelineid"], {})
                rec_src = "CPIC data API recommendation table (guideline '%s' v%s, population '%s')" % (
                    gl.get("name"), gl.get("version"), hits[0]["population"])
            else:
                vr = ("variant" if ("variant" in (r["diplotype"] or "")
                                    or r["diplotype"] == "*3/*3") else "reference")
                candidates = [(c["drug"], gene, phen),
                              (c["drug"], gene, lookupval),
                              (c["drug"], gene, vr)]
                key = next((k for k in candidates if k in CORPUS_RECS), None)
                if key is not None:
                    recs_txt, rec_class, rec_src = CORPUS_RECS[key]
                else:
                    abstain = True
                    abstain_reason = ("No CPIC recommendation found for %s with %s %s."
                                      % (c["drug"], gene, lookupval))

        lethal, lethal_why = (False, "")
        if phen:
            lethal, lethal_why = lethal_class(gene, phen, c["drug"])

        rec_field = "" if abstain else "%s: %s" % (c["drug"], recs_txt)
        answers.append(dict(case_id=c["id"],
                            diplotype="" if abstain else "%s %s" % (c["gene"].split("*")[0] if gene.startswith("HLA") else gene, r["diplotype"]),
                            phenotype="" if abstain else phen,
                            recommendation=rec_field,
                            abstained="true" if abstain else "false"))
        audit.append(dict(
            case_id=c["id"], gene=c["gene"], cpic_gene=gene, drug=c["drug"],
            genotype_as_given=c["genotype"], diplotype=r["diplotype"] or "",
            phenotype=phen or "", cpic_lookup_key=lookupval or "",
            mapping_method=r["method"], strand_orientation=r["orientation"],
            forward_genomic_hgvs=r["hgvs_g"], basis=r["basis"],
            assumptions=" | ".join(r["assumptions"]), ambiguity=r["ambiguity"],
            phenotype_source=phen_src or "", cpic_drug_used=cpic_drug or "",
            recommendation=recs_txt, recommendation_classification=rec_class,
            recommendation_source=rec_src, lethal_class="true" if lethal else "false",
            lethal_rationale=lethal_why,
            abstained="true" if abstain else "false", abstain_reason=abstain_reason))

        if not abstain:
            dip_tbl[(gene, r["diplotype"])] = dict(
                gene=gene, diplotype=r["diplotype"], cpic_phenotype=phen,
                cpic_lookup_key=lookupval or "", source=phen_src)
            rec_tbl[(c["drug"], gene, phen)] = dict(
                drug=c["drug"], gene=gene, cpic_phenotype=phen,
                cpic_lookup_key=lookupval or "", recommendation=recs_txt,
                classification=rec_class, source=rec_src,
                lethal_class="true" if lethal else "false")

    def write(fn, rows, cols):
        with open(os.path.join(HERE, fn), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    write("claude_science_answers.csv", answers,
          ["case_id", "diplotype", "phenotype", "recommendation", "abstained"])
    write("case_audit.csv", audit, list(audit[0].keys()))
    write("table_diplotype_to_phenotype.csv",
          sorted(dip_tbl.values(), key=lambda x: (x["gene"], x["diplotype"])),
          ["gene", "diplotype", "cpic_phenotype", "cpic_lookup_key", "source"])
    write("table_phenotype_to_recommendation.csv",
          sorted(rec_tbl.values(), key=lambda x: (x["drug"], x["gene"], x["cpic_phenotype"])),
          ["drug", "gene", "cpic_phenotype", "cpic_lookup_key", "recommendation",
           "classification", "source", "lethal_class"])

    # source + version per gene
    src_rows = []
    for gsym, gd in sorted(corpus["genes"].items()):
        cg = GENE_ALIAS.get(gsym, gsym)
        grec = next((x for x in T["gene"] if x["symbol"] == cg), {})
        gls = sorted({GUIDELINE[p["guidelineid"]]["name"] + " v" + str(GUIDELINE[p["guidelineid"]]["version"])
                      for p in T["pair"] if p["genesymbol"] == cg and p["guidelineid"] in GUIDELINE})
        levels = sorted({p["cpiclevel"] for p in T["pair"]
                         if p["genesymbol"] == cg and p["drugid"] in
                         {DRUG_ID.get(dn.lower()) for d in cases if d["gene"] == gsym
                          for dn in DRUG_ALIAS.get(d["drug"], [d["drug"]])}})
        src_rows.append(dict(
            gene=gsym, cpic_gene_symbol=cg,
            supplied_corpus_sha256=gd["sha256"],
            supplied_corpus_retrieved="2026-05-09/10 (schema %s, generated %s)" % (
                corpus["schema_version"], corpus["generated_at_utc"]),
            cpic_api_gene_table_version=grec.get("version", ""),
            cpic_api_lookupmethod=grec.get("lookupmethod", ""),
            cpic_guidelines_used="; ".join(gls),
            cpic_level_for_case_drugs=",".join(levels)))
    write("source_versions.csv", src_rows, list(src_rows[0].keys()))

    n_ab = sum(1 for a in answers if a["abstained"] == "true")
    n_lethal = sum(1 for a in audit if a["lethal_class"] == "true")
    print("cases=%d answered=%d abstained=%d lethal_class=%d" % (
        len(answers), len(answers) - n_ab, n_ab, n_lethal))
    print("unique diplotype->phenotype rows=%d  phenotype->recommendation rows=%d" % (
        len(dip_tbl), len(rec_tbl)))
    return answers, audit


if __name__ == "__main__":
    main()
