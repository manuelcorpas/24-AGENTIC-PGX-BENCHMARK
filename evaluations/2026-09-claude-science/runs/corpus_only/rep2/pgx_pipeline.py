"""
Blinded pharmacogenomics study - genotype -> diplotype -> phenotype -> CPIC recommendation.

Method: rule-based pipeline, no network access.
  Stage 1  parse the free-text genotype description into per-locus variant dosages
           (or into a non-SNV descriptor: CNV, HLA typing call, mtDNA variant, TA repeat).
  Stage 2  assemble a diplotype from allele definitions (rsID -> star allele), honouring
           known cis relationships (CYP2D6 *4 carries rs1065852; TPMT *3A = *3B+*3C;
           G6PD A- = c.202A + c.376G in cis).
  Stage 3  translate diplotype -> phenotype using CPIC activity scores / function pairs.
  Stage 4  look up phenotype + drug -> recommendation in a table transcribed from the
           supplied CPIC guideline corpus (cpic_guideline_corpus.json, snapshot 2026-05-09).

Sources: cpic_guideline_corpus.json only (21 genes). Allele *definitions* (which rsID
defines which star allele, and on which strand the case coordinates are reported) are not
in the corpus and come from standard pharmacogene nomenclature; every one is recorded
explicitly in ALLELE_DEFS below so it can be audited.
"""
import json
import os
import re
import csv
import hashlib
from collections import OrderedDict

# --------------------------------------------------------------------------------------
# Stage 2 inputs: allele definitions.
# base_ref / base_var are given on the FORWARD GENOMIC strand, matching the convention the
# case descriptions state explicitly ("forward-genomic"). Genes on the minus strand
# (CYP2D6, TPMT, DPYD, G6PD, CYP4F2, VKORC1) therefore have forward bases that are the
# complement of the classic cDNA notation; the cDNA change is given for traceability.
# --------------------------------------------------------------------------------------
ALLELE_DEFS = {
    # gene: { rsid: dict }
    "CYP2D6": {
        "rs3892097": dict(allele="*4", cdna="c.1846G>A", base_ref="C", base_var="T",
                          function="no function", activity=0.0, strand="minus"),
        "rs1065852": dict(allele="*10-defining / also carried in cis on *4",
                          cdna="c.100C>T (p.P34S)", base_ref="G", base_var="A",
                          function="decreased function", activity=0.25, strand="minus"),
    },
    "CYP2C19": {
        "rs4244285": dict(allele="*2", cdna="c.681G>A", base_ref="G", base_var="A",
                          function="no function", activity=0.0, strand="plus"),
        "rs12248560": dict(allele="*17", cdna="c.-806C>T", base_ref="C", base_var="T",
                           function="increased function", activity=1.5, strand="plus"),
    },
    "CYP2C9": {
        "rs1057910": dict(allele="*3", cdna="c.1075A>C (p.I359L)", base_ref="A", base_var="C",
                          function="no function", activity=0.0, strand="plus"),
        "rs1799853": dict(allele="*2", cdna="c.430C>T (p.R144C)", base_ref="C", base_var="T",
                          function="decreased function", activity=0.5, strand="plus"),
    },
    "CYP3A5": {
        # The case text annotates the direction itself: "A/A (homozygous *1)".
        "rs776746": dict(allele="*3", cdna="g.6986A>G", base_ref="A", base_var="G",
                         function="no function", activity=0.0, strand="as-reported"),
    },
    "CYP2B6": {
        "rs3745274": dict(allele="*6 (516G>T component)", cdna="c.516G>T (p.Q172H)",
                          base_ref="G", base_var="T",
                          function="decreased function", activity=0.5, strand="plus"),
        "rs2279343": dict(allele="*4", cdna="c.785A>G (p.K262R)", base_ref="A", base_var="G",
                          function="increased function", activity=1.5, strand="plus"),
    },
    "DPYD": {
        "rs3918290": dict(allele="c.1905+1G>A (*2A)", cdna="c.1905+1G>A", base_ref="C",
                          base_var="T", function="no function", activity=0.0, strand="minus"),
        "rs56038477": dict(allele="c.1236G>A / HapB3 (c.1129-5923C>G)", cdna="c.1236G>A",
                           base_ref="C", base_var="T",
                           function="decreased function", activity=0.5, strand="minus"),
    },
    "TPMT": {
        "rs1800460": dict(allele="*3B component", cdna="c.460G>A (p.A154T)", base_ref="C",
                          base_var="T", function="no function", activity=0.0, strand="minus"),
        "rs1142345": dict(allele="*3C component", cdna="c.719A>G (p.Y240C)", base_ref="T",
                          base_var="C", function="no function", activity=0.0, strand="minus"),
    },
    "NUDT15": {
        "rs116855232": dict(allele="*3", cdna="c.415C>T (p.R139C)", base_ref="C", base_var="T",
                            function="no function", activity=0.0, strand="plus"),
    },
    "SLCO1B1": {
        "rs4149056": dict(allele="*15 (521T>C; *5 on the *1a background)",
                          cdna="c.521T>C (p.V174A)", base_ref="T", base_var="C",
                          function="no function", activity=0.0, strand="plus"),
    },
    "UGT1A1": {
        # TA-repeat length polymorphism, reported directly as TA6/TA7 rather than as bases.
        "rs8175347": dict(allele="*28 (TA7)", cdna="TA repeat in TATA box", base_ref="TA6",
                          base_var="TA7", function="decreased function", activity=0.5,
                          strand="plus"),
    },
    "G6PD": {
        "rs5030868": dict(allele="Mediterranean (class II)", cdna="c.563C>T (p.S188F)",
                          base_ref="G", base_var="A", function="deficient (class II)",
                          activity=None, strand="minus"),
        "rs1050828": dict(allele="A- component (c.202G>A)", cdna="c.202G>A (p.V68M)",
                          base_ref="C", base_var="T", function="deficient (class III)",
                          activity=None, strand="minus"),
        "rs1050829": dict(allele="A component (c.376A>G)", cdna="c.376A>G (p.N126D)",
                          base_ref="T", base_var="C", function="non-deficient (class IV)",
                          activity=None, strand="minus"),
    },
    "CYP4F2": {
        "rs2108622": dict(allele="*3", cdna="c.1297G>A (p.V433M)", base_ref="C", base_var="T",
                          function="decreased function", activity=None, strand="minus"),
    },
    "VKORC1": {
        # The corpus dose table itself uses -1639 G/A nomenclature, and the cases match it.
        "rs9923231": dict(allele="-1639A (haplotype A)", cdna="c.-1639G>A", base_ref="G",
                          base_var="A", function="reduced VKORC1 expression", activity=None,
                          strand="as-reported (corpus convention)"),
    },
    "IFNL3": {
        "rs12979860": dict(allele="T (unfavourable)", cdna="rs12979860 C>T", base_ref="C",
                           base_var="T", function="unfavourable response", activity=None,
                           strand="plus"),
    },
    "CFTR": {
        "rs75527207": dict(allele="G551D", cdna="c.1652G>A (p.G551D)", base_ref="G",
                           base_var="A", function="gating defect, ivacaftor-responsive",
                           activity=None, strand="plus"),
        "rs113993960": dict(allele="F508del", cdna="c.1521_1523delCTT (p.F508del)",
                            base_ref="CTT", base_var="del",
                            function="processing defect, ivacaftor-non-responsive",
                            activity=None, strand="plus"),
    },
    "HLA-B*57:01": {
        "rs2395029": dict(allele="HCP5 G, tags HLA-B*57:01", cdna="HCP5 g.-", base_ref="T",
                          base_var="G", function="risk-allele tag", activity=None,
                          strand="plus"),
    },
}

# --------------------------------------------------------------------------------------
# Stage 1: parse the free-text genotype description.
# --------------------------------------------------------------------------------------
RS_TOKEN = re.compile(r"(rs\d+)\s+([A-Za-z0-9]+(?:/[A-Za-z0-9]+)?)")


def parse_genotype(gene, text):
    """Return (kind, payload).

    kind is one of:
      'snv'        payload = {rsid: n_variant_copies}  (0/1/2, or 1 for hemizygous X)
      'cnv'        payload = {'allele': '*1xN'}
      'hla_pos' / 'hla_neg'
      'mito'       payload = {'variant': 'm.1555A>G'} or {}
      'ryr1_pos' / 'ryr1_neg'
      'reference'  all assayed loci reference
    """
    t = text.strip()
    low = t.lower()

    # --- non-SNV descriptors, checked first ---
    if "xn duplication" in low or "copy-number assay" in low:
        return "cnv", {"allele": "*1xN"}
    if low.startswith("hla-") and "positive" in low:
        return "hla_pos", {}
    if low.startswith("hla-") and "negative" in low:
        return "hla_neg", {}
    if low.startswith("m.1555a>g") or "m.1555a>g" in low:
        return "mito", {"variant": "m.1555A>G"}
    if re.match(r"^m\.1555a\b", low):
        return "mito", {}
    if "ryr1" in low and ("no pathogenic" in low or "negative" in low):
        return "ryr1_neg", {}
    if "ryr1 mutation" in low or ("ryr1" in low and "c.1840c>t" in low):
        return "ryr1_pos", {"variant": "RYR1 c.1840C>T p.(Arg614Cys)"}
    if re.search(r"all .*(homozygous reference|risk variants reference)", low) or \
       low in ("all dpyd risk variants reference",):
        return "reference", {}
    if "f508del homozygous" in low:
        return "snv", {"dosages": {"rs113993960": 2}, "meta": {}}

    # --- rsID tokens ---
    defs = ALLELE_DEFS.get(gene, {})
    dosages = OrderedDict()
    hemizygous = "hemizygous" in low
    het_female = "heterozygous female" in low
    for rsid, alleles in RS_TOKEN.findall(t):
        if rsid not in defs:
            continue
        d = defs[rsid]
        if "/" in alleles:
            calls = alleles.split("/")
        else:
            calls = [alleles]            # hemizygous male, or a phased single haplotype call
        n_var = sum(1 for c in calls if c.upper() == str(d["base_var"]).upper())
        # "G/T or G/G (any positive)" style: presence of the risk base is what matters
        if "any positive" in low and n_var > 0:
            n_var = 1
        prev = dosages.get(rsid)
        # a locus may be named twice ("rs1050828 C/T + ... or T/T"); keep the first (primary) call
        if prev is None:
            dosages[rsid] = n_var
    if not dosages:
        return "unparsed", {"text": t}
    meta = {"hemizygous": hemizygous, "het_female": het_female}
    return "snv", {"dosages": dosages, "meta": meta}


# --------------------------------------------------------------------------------------
# Stage 2: diplotype assembly. One resolver per gene; each returns
# (diplotype, activity_score_or_None, note).
# --------------------------------------------------------------------------------------
def _pair(a, b):
    """Order a diplotype string conventionally (*1 first, then ascending)."""
    def key(x):
        m = re.match(r"\*(\d+)", x)
        return (0 if x == "*1" else 1, int(m.group(1)) if m else 99, x)
    a, b = sorted([a, b], key=key)
    return f"{a}/{b}"


def resolve_CYP2D6(kind, payload):
    if kind == "cnv":
        return "*1/*1xN", 3.0, ("Duplication of a normal-function allele detected by CNV assay. "
                                "Assumes the non-duplicated allele is normal function and N=2; "
                                "the CNV assay does not phase the duplication to a star allele.")
    if kind == "reference":
        return "*1/*1", 2.0, ""
    d = payload["dosages"]
    n4 = d.get("rs3892097", 0)
    n100 = d.get("rs1065852", 0)
    if n4 == 2:
        return "*4/*4", 0.0, ""
    if n4 == 1 and n100 == 1:
        # rs1065852 (c.100C>T) is carried in cis on the *4 haplotype, so het/het is *1/*4 and
        # NOT two decreased-function alleles. *4/*10 would require rs1065852 to be homozygous.
        return "*1/*4", 1.0, ("rs1065852 c.100C>T is part of the *4 haplotype; heterozygosity at "
                              "both loci is explained by a single *4 allele in cis. *4/*10 is "
                              "excluded because it would make rs1065852 homozygous.")
    if n4 == 1:
        return "*1/*4", 1.0, ""
    if n100 == 2:
        return "*10/*10", 0.5, ""
    if n100 == 1:
        return "*1/*10", 1.25, ""
    return "*1/*1", 2.0, ""


def resolve_CYP2C19(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    d = payload["dosages"]
    n2, n17 = d.get("rs4244285", 0), d.get("rs12248560", 0)
    if n2 == 2:
        return "*2/*2", None, ""
    if n2 == 1 and n17 == 1:
        return "*2/*17", None, "Phase unknown; CPIC assigns *2/*17 to intermediate metabolizer."
    if n2 == 1:
        return "*1/*2", None, ""
    if n17 == 2:
        return "*17/*17", None, ""
    if n17 == 1:
        return "*1/*17", None, ""
    return "*1/*1", None, ""


def resolve_CYP2C9(kind, payload):
    if kind == "reference":
        return "*1/*1", 2.0, ""
    d = payload["dosages"]
    n3, n2 = d.get("rs1057910", 0), d.get("rs1799853", 0)
    alleles = ["*3"] * n3 + ["*2"] * n2
    alleles += ["*1"] * (2 - len(alleles))
    act = {"*1": 1.0, "*2": 0.5, "*3": 0.0}
    return _pair(alleles[0], alleles[1]), sum(act[a] for a in alleles), ""


def resolve_CYP3A5(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    n3 = payload["dosages"].get("rs776746", 0)
    return {0: "*1/*1", 1: "*1/*3", 2: "*3/*3"}[n3], None, ""


def resolve_CYP2B6(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    d = payload["dosages"]
    n6, n4 = d.get("rs3745274", 0), d.get("rs2279343", 0)
    if n6 == 2:
        return "*6/*6", None, ("c.516G>T homozygous. Reported as *6/*6 by the usual panel "
                               "convention; strictly, c.516G>T without c.785A>G defines *9, so "
                               "*9/*9 is the alternative call. Both alleles are decreased "
                               "function, so the phenotype is unchanged.")
    if n6 == 1:
        return "*1/*6", None, ("c.516G>T heterozygous. Reported as *1/*6 by panel convention; "
                               "*1/*9 is the alternative call if c.785A>G is reference. Both are "
                               "decreased function, so the phenotype is unchanged.")
    if n4 == 1:
        return "*1/*4", None, ("c.785A>G heterozygous with c.516G>T reference defines *4 "
                               "(increased function).")
    if n4 == 2:
        return "*4/*4", None, ""
    return "*1/*1", None, ""


def resolve_DPYD(kind, payload):
    if kind == "reference":
        return "c.[=];[=]", 2.0, ""
    d = payload["dosages"]
    n2a, nb3 = d.get("rs3918290", 0), d.get("rs56038477", 0)
    if n2a == 2:
        return "c.[1905+1G>A];[1905+1G>A]", 0.0, ""
    if n2a == 1:
        return "c.[1905+1G>A];[=]", 1.0, ""
    if nb3 == 2:
        return "c.[1129-5923C>G];[1129-5923C>G]", 1.0, ""
    if nb3 == 1:
        return "c.[1129-5923C>G];[=] (HapB3, tagged by c.1236G>A)", 1.5, \
               "rs56038477 (c.1236G>A) is in complete linkage with the HapB3 intronic variant."
    return "c.[=];[=]", 2.0, ""


def resolve_TPMT(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    d = payload["dosages"]
    n460, n719 = d.get("rs1800460", 0), d.get("rs1142345", 0)
    if n460 == 2 and n719 == 2:
        return "*3A/*3A", None, "*3A carries both c.460G>A and c.719A>G in cis."
    if n460 == 1 and n719 == 1:
        return "*1/*3A", None, "Both changes in cis define one *3A allele."
    if n460 >= 1 and n719 == 0:
        return "*1/*3B" if n460 == 1 else "*3B/*3B", None, "c.460G>A alone defines *3B."
    if n719 >= 1 and n460 == 0:
        return "*1/*3C" if n719 == 1 else "*3C/*3C", None, "c.719A>G alone defines *3C."
    return "*1/*1", None, ""


def resolve_NUDT15(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    n = payload["dosages"].get("rs116855232", 0)
    return {0: "*1/*1", 1: "*1/*3", 2: "*3/*3"}[n], None, ""


def resolve_SLCO1B1(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    n = payload["dosages"].get("rs4149056", 0)
    note = ("rs4149056 alone cannot distinguish *5 from *15 (which additionally carries "
            "c.388A>G); reported as *15 by convention.")
    return {0: "*1/*1", 1: "*1/*15", 2: "*15/*15"}[n], None, "" if n == 0 else note


def resolve_UGT1A1(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    n = payload["dosages"].get("rs8175347", 0)
    return {0: "*1/*1", 1: "*1/*28", 2: "*28/*28"}[n], None, ""


def resolve_G6PD(kind, payload):
    d = payload["dosages"]
    meta = payload.get("meta", {})
    n_med = d.get("rs5030868", 0)
    n_202 = d.get("rs1050828", 0)
    n_376 = d.get("rs1050829", 0)
    het_f = meta.get("het_female", False)

    def ploidy(allele):
        # G6PD is X-linked; unless the case states a heterozygous female, the same variant
        # call is hemizygous in males and homozygous in females, with the same phenotype.
        return f"{allele} (hemizygous male) or {allele}/{allele} (homozygous female)"

    if n_med:
        return ploidy("Mediterranean"), None, \
            "c.563C>T = G6PD Mediterranean, WHO class II (deficient)."
    if n_202 and n_376:
        if het_f:
            return "B/A-", None, ("Heterozygous female: one non-deficient (class IV) and one "
                                  "deficient A- allele; X-inactivation makes activity variable.")
        return ploidy("A-"), None, \
            ("c.202G>A and c.376A>G in cis = G6PD A- (A-202A_376G), WHO class III (deficient).")
    if n_202 and not n_376:
        if het_f:
            return "B/A-", None, "Heterozygous female carrying c.202G>A."
        return ploidy("c.202G>A"), None, "c.202G>A without c.376A>G."
    return ploidy("B"), None, \
        "No deficiency variant detected at the assayed loci (class IV, B allele)."


def resolve_CYP4F2(kind, payload):
    if kind == "reference":
        return "*1/*1", None, ""
    n = payload["dosages"].get("rs2108622", 0)
    return {0: "*1/*1", 1: "*1/*3", 2: "*3/*3"}[n], None, ""


def resolve_VKORC1(kind, payload):
    n = payload["dosages"].get("rs9923231", 0)
    return {0: "-1639G/G", 1: "-1639G/A", 2: "-1639A/A"}[n], None, ""


def resolve_IFNL3(kind, payload):
    n = payload["dosages"].get("rs12979860", 0)
    return {0: "rs12979860 C/C", 1: "rs12979860 C/T", 2: "rs12979860 T/T"}[n], None, ""


def resolve_CFTR(kind, payload):
    d = payload["dosages"]
    if d.get("rs113993960"):
        return "F508del/F508del", None, ""
    n = d.get("rs75527207", 0)
    if n == 2:
        return "G551D/G551D", None, ""
    if n == 1:
        return "G551D/other (heterozygous) or G551D/G551D", None, \
            ("The case reports rs75527207 as G/A or A/A, so G551D is present but the "
             "heterozygous and homozygous states are not distinguished, and the second CFTR "
             "allele is not reported. Ivacaftor is recommended for homozygous or heterozygous "
             "G551D, so the recommendation is unaffected.")
    return "no ivacaftor-responsive variant detected", None, ""


def resolve_HLA_B5701(kind, payload):
    if kind == "snv":
        n = payload["dosages"].get("rs2395029", 0)
        if n:
            return "HLA-B*57:01 positive (*57:01/*X)", None, \
                   "rs2395029 (HCP5) G is a tag SNP for HLA-B*57:01, not direct HLA typing."
        return "HLA-B*57:01 negative (*X/*X)", None, \
               "rs2395029 T/T; tag SNP negative, not direct HLA typing."
    return ("HLA-B*57:01 positive (*57:01/*X)" if kind == "hla_pos"
            else "HLA-B*57:01 negative (*X/*X)"), None, ""


def resolve_HLA_generic(gene):
    def f(kind, payload):
        allele = gene.split("*", 1)[1]
        if kind == "hla_pos":
            return f"{gene} positive (*{allele}/*X)", None, ""
        return f"{gene} negative (*X/*X)", None, ""
    return f


def resolve_MTRNR1(kind, payload):
    if kind == "mito" and payload.get("variant"):
        return "m.1555A>G (homoplasmic)", None, ""
    return "m.1555A (no increased-risk variant at the assayed position)", None, \
           ("Only m.1555 was reported; m.1095T>C and m.1494C>T were not assessed and would "
            "also confer increased risk.")


def resolve_RYR1(kind, payload):
    if kind == "ryr1_pos":
        return "RYR1 c.1840C>T p.(Arg614Cys), heterozygous", None, \
               "Listed as an EMHG malignant-hyperthermia causative variant in the CPIC corpus."
    return "no EMHG causative RYR1/CACNA1S variant detected", None, ""


RESOLVERS = {
    "CYP2D6": resolve_CYP2D6, "CYP2C19": resolve_CYP2C19, "CYP2C9": resolve_CYP2C9,
    "CYP3A5": resolve_CYP3A5, "CYP2B6": resolve_CYP2B6, "DPYD": resolve_DPYD,
    "TPMT": resolve_TPMT, "NUDT15": resolve_NUDT15, "SLCO1B1": resolve_SLCO1B1,
    "UGT1A1": resolve_UGT1A1, "G6PD": resolve_G6PD, "CYP4F2": resolve_CYP4F2,
    "VKORC1": resolve_VKORC1, "IFNL3": resolve_IFNL3, "CFTR": resolve_CFTR,
    "HLA-B*57:01": resolve_HLA_B5701, "HLA-B*15:02": resolve_HLA_generic("HLA-B*15:02"),
    "HLA-A*31:01": resolve_HLA_generic("HLA-A*31:01"),
    "HLA-B*58:01": resolve_HLA_generic("HLA-B*58:01"),
    "MT-RNR1": resolve_MTRNR1, "RYR1": resolve_RYR1,
}

# --------------------------------------------------------------------------------------
# Stage 3: diplotype -> phenotype.
# CYP2D6 uses the activity-score bands printed in the corpus codeine/paroxetine tables
# (the CPIC standardised translation: UM >2.25, NM 1.25-2.25, IM 0<AS<1.25, PM 0).
# That translation is applied to ALL CYP2D6 drugs, including ondansetron and amitriptyline
# whose corpus tables still print the superseded 1.0-2.0 "normal metabolizer" band; both
# of those annotations carry the note "CPIC has updated the CYP2D6 genotype-phenotype
# translation in all supplemental files since the guideline was published".
# --------------------------------------------------------------------------------------
def fmt_as(a):
    return f"{a:g}"


def phenotype(gene, diplotype, activity):
    """Return (phenotype_label, pheno_key)."""
    if gene == "CYP2D6":
        a = activity
        if a == 0:
            return f"Poor metabolizer (activity score {fmt_as(a)})", "PM"
        if a < 1.25:
            return f"Intermediate metabolizer (activity score {fmt_as(a)})", "IM"
        if a <= 2.25:
            return f"Normal metabolizer (activity score {fmt_as(a)})", "NM"
        return f"Ultrarapid metabolizer (activity score {fmt_as(a)})", "UM"

    if gene == "CYP2C19":
        k = {"*2/*2": "PM", "*1/*2": "IM", "*2/*17": "IM", "*1/*1": "NM",
             "*1/*17": "RM", "*17/*17": "UM"}[diplotype]
        return {"PM": "Poor metabolizer", "IM": "Intermediate metabolizer",
                "NM": "Normal metabolizer", "RM": "Rapid metabolizer",
                "UM": "Ultrarapid metabolizer"}[k], k

    if gene == "CYP2C9":
        a = activity
        if a == 2.0:
            return "Normal metabolizer (activity score 2)", "NM"
        if a == 1.5:
            return "Intermediate metabolizer (activity score 1.5)", "IM_AS1.5"
        if a == 1.0:
            return "Intermediate metabolizer (activity score 1)", "IM_AS1"
        return f"Poor metabolizer (activity score {fmt_as(a)})", "PM"

    if gene == "CYP3A5":
        k = {"*1/*1": "EM", "*1/*3": "IM", "*3/*3": "PM"}[diplotype]
        return {"EM": "Extensive metabolizer (CYP3A5 expresser)",
                "IM": "Intermediate metabolizer (CYP3A5 expresser)",
                "PM": "Poor metabolizer (CYP3A5 non-expresser)"}[k], k

    if gene == "CYP2B6":
        k = {"*1/*1": "NM", "*1/*6": "IM", "*6/*6": "PM", "*1/*4": "RM", "*4/*4": "UM"}[diplotype]
        return {"NM": "Normal metabolizer", "IM": "Intermediate metabolizer",
                "PM": "Poor metabolizer", "RM": "Rapid metabolizer",
                "UM": "Ultrarapid metabolizer"}[k], k

    if gene == "DPYD":
        a = activity
        if a == 2.0:
            return "Normal metabolizer (activity score 2)", "NM"
        if a in (1.0, 1.5):
            return f"Intermediate metabolizer (activity score {fmt_as(a)})", "IM"
        if a == 0.5:
            return "Poor metabolizer (activity score 0.5)", "PM_AS0.5"
        return "Poor metabolizer (activity score 0)", "PM_AS0"

    if gene in ("TPMT", "NUDT15"):
        if diplotype == "*1/*1":
            return "Normal metabolizer", "NM"
        a, b = diplotype.split("/")
        nvar = sum(1 for x in (a, b) if x != "*1")
        if nvar == 1:
            return "Intermediate metabolizer", "IM"
        return "Poor metabolizer", "PM"

    if gene == "SLCO1B1":
        k = {"*1/*1": "NORMAL", "*1/*15": "DECREASED", "*15/*15": "POOR"}[diplotype]
        return {"NORMAL": "Normal function", "DECREASED": "Decreased function",
                "POOR": "Poor function"}[k], k

    if gene == "UGT1A1":
        k = {"*1/*1": "NM", "*1/*28": "IM", "*28/*28": "PM"}[diplotype]
        return {"NM": "Normal (extensive) metabolizer", "IM": "Intermediate metabolizer",
                "PM": "Poor metabolizer"}[k], k

    if gene == "G6PD":
        if diplotype.startswith("B ") or diplotype in ("B", "B/B"):
            return "Normal", "NORMAL"
        if diplotype == "B/A-":
            return "Variable", "VARIABLE"
        return "Deficient (WHO class II-III)", "DEFICIENT"

    if gene == "CYP4F2":
        if diplotype == "*3/*3":
            return "Decreased CYP4F2 function (*3/*3)", "DECREASED"
        if diplotype == "*1/*3":
            return "Decreased CYP4F2 function (*1/*3)", "DECREASED"
        return "Normal CYP4F2 function", "NORMAL"

    if gene == "VKORC1":
        if diplotype == "-1639A/A":
            return "Increased warfarin sensitivity (low dose requirement)", "AA"
        if diplotype == "-1639G/A":
            return "Intermediate warfarin sensitivity", "GA"
        return "Normal warfarin sensitivity (reference)", "GG"

    if gene == "IFNL3":
        if diplotype == "rs12979860 C/C":
            return "Favorable response genotype", "FAV"
        return "Unfavorable response genotype", "UNFAV"

    if gene == "CFTR":
        if diplotype.startswith("G551D"):
            return "Ivacaftor-responsive CFTR gating variant (G551D)", "RESPONSIVE"
        if diplotype == "F508del/F508del":
            return "Unlikely to respond to ivacaftor (F508del homozygous)", "NONRESPONSIVE"
        return "No ivacaftor-responsive variant detected", "NONE"

    if gene.startswith("HLA-"):
        pos = "positive" in diplotype
        if gene == "HLA-B*57:01":
            return ("High risk of abacavir hypersensitivity" if pos
                    else "Very low risk of abacavir hypersensitivity"), ("POS" if pos else "NEG")
        if gene == "HLA-B*58:01":
            return ("Significantly increased risk of allopurinol SCAR" if pos
                    else "Low or reduced risk of allopurinol SCAR"), ("POS" if pos else "NEG")
        return (f"Carrier of {gene}: increased risk of severe cutaneous adverse reaction" if pos
                else f"Non-carrier of {gene}: normal or reduced risk"), ("POS" if pos else "NEG")

    if gene == "MT-RNR1":
        if diplotype.startswith("m.1555A>G"):
            return "Increased risk of aminoglycoside-induced hearing loss", "POS"
        return "Normal risk of aminoglycoside-induced hearing loss", "NEG"

    if gene == "RYR1":
        if diplotype.startswith("RYR1 c."):
            return "Malignant hyperthermia susceptible", "MHS"
        return "Uncertain susceptibility (no causative variant detected)", "UNCERTAIN"

    raise KeyError(gene)


# --------------------------------------------------------------------------------------
# Stage 4: phenotype + drug -> recommendation, transcribed from the corpus.
# ABSTAIN marks a gene-drug pair for which the supplied corpus contains no recommendation.
# --------------------------------------------------------------------------------------
ABSTAIN = "ABSTAIN"

RECS = {
    ("CYP2D6", "codeine"): {
        "UM": "codeine: avoid; risk of life-threatening morphine toxicity - if an opioid is needed use a non-tramadol opioid (e.g. morphine) or a non-opioid analgesic",
        "NM": "codeine: use label-recommended age- or weight-specific dosing",
        "IM": "codeine: use label-recommended age- or weight-specific dosing; if no response and an opioid is warranted, consider a non-tramadol opioid",
        "PM": "codeine: avoid because of diminished analgesia; if an opioid is warranted use a non-tramadol opioid",
    },
    ("CYP2D6", "tamoxifen"): {
        "UM": "tamoxifen: initiate standard-of-care dosing (20 mg/day); avoid moderate and strong CYP2D6 inhibitors",
        "NM": "tamoxifen: initiate standard-of-care dosing (20 mg/day); avoid moderate and strong CYP2D6 inhibitors",
        "IM": "tamoxifen: consider an aromatase inhibitor (plus ovarian function suppression if premenopausal); if an aromatase inhibitor is contraindicated consider tamoxifen 40 mg/day; avoid CYP2D6 inhibitors",
        "PM": "tamoxifen: avoid - recommend alternative hormonal therapy (aromatase inhibitor, plus ovarian function suppression if premenopausal); if an aromatase inhibitor is contraindicated, tamoxifen 40 mg/day may be considered",
    },
    ("CYP2D6", "ondansetron"): {
        "UM": "ondansetron: select an alternative antiemetic not predominantly metabolised by CYP2D6 (e.g. granisetron)",
        "NM": "ondansetron: initiate therapy with the recommended starting dose",
        "IM": "ondansetron: insufficient evidence of clinical impact - initiate therapy with the recommended starting dose",
        "PM": "ondansetron: insufficient evidence of clinical impact - initiate therapy with the recommended starting dose",
    },
    ("CYP2D6", "paroxetine"): {
        "UM": "paroxetine: select an alternative drug not predominantly metabolised by CYP2D6",
        "NM": "paroxetine: initiate therapy with the recommended starting dose",
        "IM": "paroxetine: consider a lower starting dose and a slower titration schedule than for normal metabolizers",
        "PM": "paroxetine: consider a 50% reduction in starting dose, a slower titration schedule and a 50% lower maintenance dose",
    },
    ("CYP2D6", "amitriptyline"): {
        "UM": "amitriptyline: avoid because of potential lack of efficacy; consider a drug not metabolised by CYP2D6 (CYP2C19 genotype not supplied)",
        "NM": "amitriptyline: initiate therapy with the recommended starting dose (CYP2C19 genotype not supplied)",
        "IM": "amitriptyline: consider a 25% reduction of the recommended starting dose with therapeutic drug monitoring (CYP2C19 genotype not supplied)",
        "PM": "amitriptyline: avoid because of potential for side effects; if a TCA is warranted consider a 50% dose reduction with therapeutic drug monitoring (CYP2C19 genotype not supplied)",
    },
    ("CYP2C19", "clopidogrel"): {
        "UM": "clopidogrel: if considering clopidogrel, use at standard dose (75 mg/day)",
        "RM": "clopidogrel: if considering clopidogrel, use at standard dose (75 mg/day)",
        "NM": "clopidogrel: use at standard dose (75 mg/day)",
        "IM": "clopidogrel: avoid standard-dose clopidogrel if possible - use prasugrel or ticagrelor at standard dose if no contraindication",
        "PM": "clopidogrel: avoid if possible - use prasugrel or ticagrelor at standard dose if no contraindication",
    },
    ("CYP2C19", "voriconazole"): {
        "UM": "voriconazole: choose an alternative agent not dependent on CYP2C19 metabolism (isavuconazole, liposomal amphotericin B, posaconazole)",
        "RM": "voriconazole: choose an alternative agent not dependent on CYP2C19 metabolism (adult recommendation; in children initiate standard dosing and titrate by therapeutic drug monitoring)",
        "NM": "voriconazole: initiate therapy with recommended standard-of-care dosing",
        "IM": "voriconazole: initiate therapy with recommended standard-of-care dosing",
        "PM": "voriconazole: choose an alternative agent not dependent on CYP2C19 metabolism; if voriconazole is judged most appropriate, use a lower than standard dose with careful therapeutic drug monitoring",
    },
    ("CYP2C19", "citalopram"): {
        "UM": "citalopram: consider a clinically appropriate antidepressant not predominantly metabolised by CYP2C19; if citalopram is used and efficacy is inadequate at standard maintenance dosing, consider titrating to a higher maintenance dose",
        "RM": "citalopram: initiate with the recommended starting dose; if response is inadequate consider a higher maintenance dose or an alternative antidepressant not predominantly metabolised by CYP2C19",
        "NM": "citalopram: initiate therapy with the recommended starting dose",
        "IM": "citalopram: initiate with the recommended starting dose; consider a slower titration schedule and a lower maintenance dose than for normal metabolizers",
        "PM": "citalopram: consider an antidepressant not predominantly metabolised by CYP2C19; if citalopram is used, lower starting dose, slower titration and 50% reduction of the standard maintenance dose (FDA maximum 20 mg/day in CYP2C19 poor metabolizers because of QT prolongation)",
    },
    ("CYP2C19", "omeprazole"): {
        "UM": "omeprazole: increase the starting daily dose by 100% (may be divided) and monitor for efficacy",
        "RM": "omeprazole: initiate the standard starting daily dose; consider a 50-100% increase for H. pylori infection or erosive oesophagitis; monitor for efficacy",
        "NM": "omeprazole: initiate the standard starting daily dose; consider a 50-100% increase for H. pylori infection or erosive oesophagitis; monitor for efficacy",
        "IM": "omeprazole: initiate the standard starting daily dose; for chronic therapy (>12 weeks) with efficacy achieved, consider a 50% dose reduction",
        "PM": "omeprazole: initiate the standard starting daily dose; for chronic therapy (>12 weeks) with efficacy achieved, consider a 50% dose reduction",
    },
    ("CYP2C9", "warfarin"): {
        "NM": "warfarin: use a genotype-guided dosing algorithm; CPIC dose table for *1/*1 gives 5-7 mg/day with VKORC1 -1639 GG or GA and 3-4 mg/day with AA (VKORC1 genotype not supplied)",
        "IM_AS1.5": "warfarin: use a genotype-guided dosing algorithm; CPIC dose table for *1/*2 gives 5-7 mg/day with VKORC1 -1639 GG and 3-4 mg/day with GA or AA (VKORC1 genotype not supplied)",
        "IM_AS1": "warfarin: reduced dose required - CPIC dose table for *1/*3 gives 3-4 mg/day with VKORC1 -1639 GG or GA and 0.5-2 mg/day with AA (VKORC1 genotype not supplied); use a genotype-guided algorithm with close INR monitoring",
        "PM": "warfarin: markedly reduced dose required - CPIC dose table gives 0.5-2 mg/day for *3/*3 at every VKORC1 -1639 genotype; consider an alternative anticoagulant, or use a genotype-guided algorithm with close INR monitoring",
    },
    ("CYP2C9", "phenytoin"): {
        "NM": "phenytoin: no adjustment from typical dosing; adjust by therapeutic drug monitoring (HLA-B*15:02 not supplied - CPIC advises avoiding phenytoin in *15:02 carriers)",
        "IM_AS1.5": "phenytoin: no adjustment from typical dosing; adjust by therapeutic drug monitoring (HLA-B*15:02 not supplied)",
        "IM_AS1": "phenytoin: use a typical loading dose, then approximately 25% less than the typical maintenance dose, adjusted by therapeutic drug monitoring (HLA-B*15:02 not supplied)",
        "PM": "phenytoin: use a typical loading dose, then approximately 50% less than the typical maintenance dose, adjusted by therapeutic drug monitoring (HLA-B*15:02 not supplied - CPIC advises avoiding phenytoin in *15:02 carriers)",
    },
    ("CYP2C9", "celecoxib"): {
        "NM": "celecoxib: initiate therapy with the recommended starting dose at the lowest effective dosage for the shortest duration",
        "IM_AS1.5": "celecoxib: initiate therapy with the recommended starting dose; monitor for adverse events",
        "IM_AS1": "celecoxib: initiate therapy with the lowest recommended starting dose and titrate upward with caution, monitoring blood pressure and kidney function",
        "PM": "celecoxib: initiate at 25-50% of the lowest recommended starting dose, titrating no sooner than steady state (>=8 days), or use an alternative NSAID not metabolised by CYP2C9 (e.g. naproxen, sulindac, aspirin, ketorolac)",
    },
    ("CYP3A5", "tacrolimus"): {
        "EM": "tacrolimus: increase the starting dose 1.5-2 times the recommended starting dose (not exceeding 0.3 mg/kg/day) and use therapeutic drug monitoring",
        "IM": "tacrolimus: increase the starting dose 1.5-2 times the recommended starting dose (not exceeding 0.3 mg/kg/day) and use therapeutic drug monitoring",
        "PM": "tacrolimus: initiate therapy with the standard recommended dose and use therapeutic drug monitoring",
    },
    ("CYP2B6", "efavirenz"): {
        "UM": "efavirenz: initiate with standard dosing (600 mg/day)",
        "RM": "efavirenz: initiate with standard dosing (600 mg/day)",
        "NM": "efavirenz: initiate with standard dosing (600 mg/day)",
        "IM": "efavirenz: consider initiating with a decreased dose of 400 mg/day",
        "PM": "efavirenz: consider initiating with a decreased dose of 400 or 200 mg/day",
    },
    ("DPYD", "fluorouracil"): {
        "NM": "fluorouracil: no change from label-recommended dosage and administration",
        "IM": "fluorouracil: reduce the starting dose by 50%, then titrate by toxicity or therapeutic drug monitoring",
        "PM_AS0.5": "fluorouracil: avoid fluorouracil and fluoropyrimidine-based regimens; if no suitable alternative exists, use a strongly reduced dose with early therapeutic drug monitoring",
        "PM_AS0": "fluorouracil: avoid fluorouracil and fluoropyrimidine prodrug-based regimens (complete DPD deficiency - risk of fatal toxicity)",
    },
    ("DPYD", "capecitabine"): {
        "NM": "capecitabine: no change from label-recommended dosage and administration",
        "IM": "capecitabine: reduce the starting dose by 50%, then titrate by toxicity or therapeutic drug monitoring",
        "PM_AS0.5": "capecitabine: avoid capecitabine and fluoropyrimidine-based regimens; if no suitable alternative exists, use a strongly reduced dose with early therapeutic drug monitoring",
        "PM_AS0": "capecitabine: avoid capecitabine and fluoropyrimidine prodrug-based regimens (complete DPD deficiency - risk of fatal toxicity)",
    },
    ("TPMT", "azathioprine"): {
        "NM": "azathioprine: start at the normal starting dose and adjust by disease-specific guidelines",
        "IM": "azathioprine: consider a lower than standard starting dose, with adjustment based on tolerance and myelosuppression monitoring",
        "PM": "azathioprine: consider an alternative immunosuppressant; if a thiopurine is used it must be at a drastically reduced dose with close haematological monitoring",
    },
    ("NUDT15", "thiopurines"): {
        "NM": "thiopurines (azathioprine/mercaptopurine/thioguanine): start at the normal starting dose",
        "IM": "thiopurines (azathioprine/mercaptopurine/thioguanine): consider a lower than standard starting dose, with adjustment based on tolerance and myelosuppression monitoring",
        "PM": "thiopurines (azathioprine/mercaptopurine/thioguanine): consider an alternative immunosuppressant; if a thiopurine is used it must be at a drastically reduced dose with close haematological monitoring",
    },
    ("SLCO1B1", "simvastatin"): {
        "NORMAL": "simvastatin: prescribe the desired starting dose and adjust by disease-specific guidelines",
        "DECREASED": "simvastatin: prescribe an alternative statin depending on the desired potency; if simvastatin is warranted limit the dose to <20 mg/day",
        "POOR": "simvastatin: prescribe an alternative statin depending on the desired potency (high myopathy risk)",
    },
    ("SLCO1B1", "atorvastatin"): {
        "NORMAL": "atorvastatin: prescribe the desired starting dose and adjust by disease-specific guidelines",
        "DECREASED": "atorvastatin: prescribe <=40 mg as a starting dose; be aware of increased myopathy risk, especially at 40 mg",
        "POOR": "atorvastatin: prescribe <=20 mg as a starting dose; if >20 mg is needed consider rosuvastatin or combination therapy",
    },
    ("UGT1A1", "atazanavir"): {
        "NM": "atazanavir: no need to avoid atazanavir on the basis of UGT1A1 genotype",
        "IM": "atazanavir: no need to avoid atazanavir; inform the patient that jaundice is unlikely with this genotype",
        "PM": "atazanavir: consider an alternative agent, particularly where jaundice would be of concern to the patient (high likelihood of bilirubin-related discontinuation)",
    },
    ("UGT1A1", "irinotecan"): {
        "NM": "irinotecan: no action needed for this gene-drug interaction (DPWG/RNPGx; the corpus contains no CPIC guideline for irinotecan)",
        "IM": "irinotecan: no action needed for this gene-drug interaction (DPWG/RNPGx; the corpus contains no CPIC guideline for irinotecan)",
        "PM": "irinotecan: start at 70% of the normal dose and increase as tolerated, guided by neutrophil count; an intensified dose (>=240 mg/m2) is contraindicated (DPWG/RNPGx; the corpus contains no CPIC guideline for irinotecan)",
    },
    ("G6PD", "rasburicase"): {
        "NORMAL": "rasburicase: no reason to withhold rasburicase on the basis of G6PD status",
        "DEFICIENT": "rasburicase: contraindicated (risk of acute haemolytic anaemia and methaemoglobinaemia); alternatives include allopurinol",
        "VARIABLE": "rasburicase: measure G6PD enzyme activity before use to confirm normal status; alternatives include allopurinol",
    },
    ("G6PD", "primaquine"): {
        "NORMAL": "primaquine: no reason to avoid primaquine on the basis of G6PD status",
        "DEFICIENT": "primaquine: avoid standard 14-day anti-relapse dosing; only use per expert consensus regimens (0.75 mg/kg weekly for 8 weeks with monitoring for haemolysis, or a single 0.25 mg/kg gametocytocidal dose)",
        "VARIABLE": "primaquine: measure G6PD enzyme activity and treat according to the activity-based phenotype",
    },
    ("HLA-B*57:01", "abacavir"): {
        "POS": "abacavir: not recommended (high risk of abacavir hypersensitivity reaction)",
        "NEG": "abacavir: use per standard dosing guidelines",
    },
    ("HLA-B*15:02", "carbamazepine"): {
        "POS": "carbamazepine: if carbamazepine-naive, do not use carbamazepine (high risk of SJS/TEN); avoid other aromatic anticonvulsants with caution",
        "NEG": "carbamazepine: use per standard dosing guidelines",
    },
    ("HLA-B*15:02", "oxcarbazepine"): {
        "POS": "oxcarbazepine: if oxcarbazepine-naive, do not use oxcarbazepine (high risk of SJS/TEN)",
        "NEG": "oxcarbazepine: use per standard dosing guidelines",
    },
    ("HLA-A*31:01", "carbamazepine"): {
        "POS": "carbamazepine: if carbamazepine-naive and alternatives are available, do not use carbamazepine (increased risk of SJS/TEN, DRESS and maculopapular exanthema); if no alternative exists, use with increased clinical monitoring and stop at the first cutaneous reaction",
        "NEG": "carbamazepine: use per standard dosing guidelines",
    },
    ("HLA-B*58:01", "allopurinol"): {
        "POS": "allopurinol: contraindicated (significantly increased risk of severe cutaneous adverse reaction)",
        "NEG": "allopurinol: use per standard dosing guidelines",
    },
    ("MT-RNR1", "aminoglycosides"): {
        "POS": "aminoglycosides: avoid unless the high risk of permanent hearing loss is outweighed by infection severity and the lack of a safe alternative; if unavoidable, use the lowest dose and shortest duration with frequent audiometry",
        "NEG": "aminoglycosides: use at standard doses for the shortest feasible course with therapeutic drug monitoring and regular hearing assessment",
    },
    ("RYR1", "volatile-anaesthetics"): {
        "MHS": "volatile anaesthetics and succinylcholine: relatively contraindicated - use a non-triggering anaesthetic technique (total intravenous anaesthesia), have dantrolene available",
        "UNCERTAIN": "volatile anaesthetics and succinylcholine: a negative result does not exclude malignant hyperthermia susceptibility; let clinical findings, family history and further testing guide use",
    },
    ("IFNL3", "PEG-IFN-\u03b1"): {
        "FAV": "PEG-IFN-alpha (with ribavirin): favourable response genotype - approximately 70% chance of SVR after 48 weeks, or approximately 90% with a protease-inhibitor combination; weighs in favour of an interferon-containing regimen",
        "UNFAV": "PEG-IFN-alpha (with ribavirin): unfavourable response genotype - approximately 30% chance of SVR after 48 weeks, or approximately 60% with a protease-inhibitor combination; consider the implications before initiating an interferon-containing regimen",
    },
    ("CFTR", "ivacaftor"): {
        "RESPONSIVE": "ivacaftor: use according to the product label (G551D is a responsive gating variant)",
        "NONRESPONSIVE": "ivacaftor: not recommended - F508del homozygotes are unlikely to respond to ivacaftor monotherapy",
        "NONE": "ivacaftor: no responsive CFTR variant detected; ivacaftor is recommended only for patients with responsive variants",
    },
    ("VKORC1", "warfarin"): {
        "AA": "warfarin: increased sensitivity, lower dose required - CPIC dose table gives 3-4 mg/day for -1639 AA with CYP2C9 *1/*1 and 0.5-2 mg/day with *1/*3, *2/*2, *2/*3 or *3/*3 (CYP2C9 genotype not supplied); use a genotype-guided algorithm",
        "GA": "warfarin: CPIC dose table gives 5-7 mg/day for -1639 GA with CYP2C9 *1/*1 (CYP2C9 genotype not supplied); use a genotype-guided algorithm",
        "GG": "warfarin: CPIC dose table gives 5-7 mg/day for -1639 GG with CYP2C9 *1/*1 or *1/*2 (CYP2C9 genotype not supplied); use a genotype-guided algorithm",
    },
    ("CYP4F2", "warfarin"): {
        "DECREASED": ABSTAIN,
        "NORMAL": ABSTAIN,
    },
}

ABSTAIN_REASON = {
    ("CYP4F2", "warfarin"): ("The supplied CPIC warfarin excerpt indexes its dose table on "
                             "VKORC1 -1639 and CYP2C9 only; it contains no CYP4F2 row and no "
                             "CYP4F2-specific dosing statement, so no recommendation can be "
                             "derived from the corpus for a CYP4F2-only genotype."),
}

# --------------------------------------------------------------------------------------
# Harm class. "lethal" = standard dosing in this genotype carries a documented risk of a
# fatal or immediately life-threatening event. "severe" = serious or irreversible harm
# that is not usually fatal. Everything else is "routine".
# --------------------------------------------------------------------------------------
HARM = {
    ("CYP2D6", "codeine", "UM"): ("lethal", "Ultrarapid morphine formation; fatal respiratory depression reported in children and in breastfed infants of UM mothers."),
    ("CYP2C19", "clopidogrel", "PM"): ("lethal", "Loss of antiplatelet effect; stent thrombosis and fatal myocardial infarction."),
    ("CYP2C19", "clopidogrel", "IM"): ("lethal", "Reduced active metabolite; increased risk of stent thrombosis and fatal cardiovascular events."),
    ("CYP2C19", "voriconazole", "UM"): ("lethal", "Subtherapeutic exposure; failure of therapy for invasive fungal infection."),
    ("CYP2C19", "voriconazole", "RM"): ("lethal", "Subtherapeutic exposure; failure of therapy for invasive fungal infection."),
    ("CYP2C9", "warfarin", "PM"): ("lethal", "Standard 5 mg induction in *3/*3 causes marked over-anticoagulation and major haemorrhage."),
    ("CYP2C9", "warfarin", "IM_AS1"): ("lethal", "Standard induction over-anticoagulates *1/*3; major bleeding risk."),
    ("VKORC1", "warfarin", "AA"): ("lethal", "Low dose requirement; standard induction causes over-anticoagulation and major haemorrhage."),
    ("DPYD", "fluorouracil", "PM_AS0"): ("lethal", "Complete DPD deficiency; standard 5-FU is reliably fatal without dose avoidance."),
    ("DPYD", "fluorouracil", "PM_AS0.5"): ("lethal", "Near-complete DPD deficiency; risk of fatal toxicity."),
    ("DPYD", "fluorouracil", "IM"): ("lethal", "Partial DPD deficiency; severe or fatal myelosuppression and mucositis at full dose."),
    ("DPYD", "capecitabine", "PM_AS0"): ("lethal", "Complete DPD deficiency; standard capecitabine is reliably fatal without dose avoidance."),
    ("DPYD", "capecitabine", "PM_AS0.5"): ("lethal", "Near-complete DPD deficiency; risk of fatal toxicity."),
    ("DPYD", "capecitabine", "IM"): ("lethal", "Partial DPD deficiency; severe or fatal toxicity at full dose."),
    ("TPMT", "azathioprine", "PM"): ("lethal", "Life-threatening myelosuppression at conventional thiopurine doses."),
    ("TPMT", "azathioprine", "IM"): ("lethal", "Substantial risk of severe myelosuppression at conventional doses."),
    ("NUDT15", "thiopurines", "PM"): ("lethal", "Severe, life-threatening myelosuppression; tolerated dose is roughly 10% of normal."),
    ("NUDT15", "thiopurines", "IM"): ("lethal", "Increased risk of severe myelosuppression at conventional doses."),
    ("HLA-B*15:02", "carbamazepine", "POS"): ("lethal", "SJS/TEN with mortality up to ~30% for TEN."),
    ("HLA-B*15:02", "oxcarbazepine", "POS"): ("lethal", "SJS/TEN with substantial mortality."),
    ("HLA-A*31:01", "carbamazepine", "POS"): ("lethal", "SJS/TEN and DRESS, both potentially fatal."),
    ("HLA-B*58:01", "allopurinol", "POS"): ("lethal", "Allopurinol-induced SCAR (DRESS/SJS/TEN) with reported mortality up to ~25%."),
    ("HLA-B*57:01", "abacavir", "POS"): ("lethal", "Abacavir hypersensitivity syndrome; fatal on re-challenge."),
    ("G6PD", "rasburicase", "DEFICIENT"): ("lethal", "Acute intravascular haemolysis and methaemoglobinaemia; rasburicase is contraindicated."),
    ("G6PD", "primaquine", "DEFICIENT"): ("lethal", "Acute haemolytic anaemia, which can be severe enough to be fatal at standard 14-day dosing."),
    ("RYR1", "volatile-anaesthetics", "MHS"): ("lethal", "Malignant hyperthermia; untreated mortality approaches 70%."),
    ("UGT1A1", "irinotecan", "PM"): ("lethal", "Life-threatening neutropenia and diarrhoea at full dose (DPWG lists death among the clinical effects)."),
    ("SLCO1B1", "simvastatin", "POOR"): ("lethal", "Markedly increased simvastatin acid exposure; rhabdomyolysis can be fatal."),
    ("G6PD", "primaquine", "VARIABLE"): ("severe", "X-inactivation mosaicism means a heterozygous female can be functionally deficient; standard 14-day primaquine could cause significant haemolysis if enzyme activity is not measured first."),
    ("G6PD", "rasburicase", "VARIABLE"): ("severe", "A heterozygous female can be functionally deficient; rasburicase without a prior enzyme activity measurement risks acute haemolysis."),
    ("MT-RNR1", "aminoglycosides", "POS"): ("severe", "Permanent bilateral sensorineural hearing loss after even a single dose; irreversible but not fatal."),
    ("SLCO1B1", "simvastatin", "DECREASED"): ("severe", "Increased myopathy risk."),
    ("SLCO1B1", "atorvastatin", "POOR"): ("severe", "Increased atorvastatin exposure and myopathy risk."),
    ("CYP2C9", "phenytoin", "PM"): ("severe", "Phenytoin accumulation causing ataxia, encephalopathy and nystagmus; rarely fatal."),
    ("CYP2C19", "citalopram", "PM"): ("severe", "Higher exposure with QT prolongation risk; the FDA caps the dose at 20 mg/day."),
    ("CYP2D6", "amitriptyline", "PM"): ("severe", "TCA accumulation with cardiac conduction and anticholinergic toxicity."),
    ("CYP2D6", "tamoxifen", "PM"): ("severe", "Reduced endoxifen; higher breast-cancer recurrence risk (harm from under-treatment)."),
    ("CYP2C19", "voriconazole", "PM"): ("severe", "Supratherapeutic troughs with hepatotoxicity and neurotoxicity."),
    ("CFTR", "ivacaftor", "NONRESPONSIVE"): ("severe", "No benefit; harm is the opportunity cost of ineffective therapy."),
}


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------
def run(indir, outdir):
    with open(os.path.join(indir, "cases_blinded.json")) as fh:
        cases = json.load(fh)
    with open(os.path.join(indir, "cpic_guideline_corpus.json")) as fh:
        corpus = json.load(fh)

    rows, audit = [], []
    for c in cases:
        gene, drug, gtxt = c["gene"], c["drug"], c["genotype"]
        kind, payload = parse_genotype(gene, gtxt)
        if kind == "unparsed":
            rows.append(dict(case_id=c["id"], diplotype="not assigned",
                             phenotype="not assigned",
                             recommendation=f"{drug}: abstained - genotype description could "
                                            f"not be parsed", abstained="true"))
            audit.append(dict(case_id=c["id"], gene=gene, drug=drug, genotype=gtxt,
                              parsed_kind=kind, diplotype="", activity_score="",
                              phenotype="", pheno_key="", harm_class="", harm_reason="",
                              abstained="true",
                              note="genotype string not matched by any parser rule"))
            continue

        diplotype, activity, note = RESOLVERS[gene](kind, payload)
        pheno, key = phenotype(gene, diplotype, activity)
        rec = RECS[(gene, drug)][key]
        abstained = rec is ABSTAIN
        if abstained:
            reason = ABSTAIN_REASON[(gene, drug)]
            rec_text = f"{drug}: abstained - {reason}"
        else:
            rec_text = rec
        harm, harm_reason = HARM.get((gene, drug, key), ("routine", ""))

        rows.append(dict(case_id=c["id"], diplotype=diplotype, phenotype=pheno,
                         recommendation=rec_text,
                         abstained="true" if abstained else "false"))
        audit.append(dict(case_id=c["id"], gene=gene, drug=drug, genotype=gtxt,
                          parsed_kind=kind, diplotype=diplotype,
                          activity_score="" if activity is None else fmt_as(activity),
                          phenotype=pheno, pheno_key=key, harm_class=harm,
                          harm_reason=harm_reason,
                          abstained="true" if abstained else "false", note=note))

    os.makedirs(outdir, exist_ok=True)
    _write(os.path.join(outdir, "claude_science_answers.csv"),
           ["case_id", "diplotype", "phenotype", "recommendation", "abstained"], rows)
    _write(os.path.join(outdir, "case_audit_trail.csv"),
           ["case_id", "gene", "drug", "genotype", "parsed_kind", "diplotype",
            "activity_score", "phenotype", "pheno_key", "harm_class", "harm_reason",
            "abstained", "note"], audit)

    # --- table 1: allele definitions actually used ---
    adefs = []
    for gene, d in ALLELE_DEFS.items():
        for rsid, v in d.items():
            adefs.append(dict(gene=gene, rsid=rsid, allele=v["allele"], cdna_change=v["cdna"],
                              forward_ref_allele=v["base_ref"],
                              forward_variant_allele=v["base_var"],
                              gene_strand=v["strand"], allele_function=v["function"],
                              cpic_activity_value="" if v["activity"] is None
                              else fmt_as(v["activity"])))
    _write(os.path.join(outdir, "allele_definitions.csv"),
           ["gene", "rsid", "allele", "cdna_change", "forward_ref_allele",
            "forward_variant_allele", "gene_strand", "allele_function",
            "cpic_activity_value"], adefs)

    # --- table 2: diplotype -> phenotype, enumerated over the diplotypes reached ---
    seen, d2p = set(), []
    for a in audit:
        k = (a["gene"], a["diplotype"])
        if k in seen or not a["diplotype"]:
            continue
        seen.add(k)
        d2p.append(dict(gene=a["gene"], diplotype=a["diplotype"],
                        activity_score=a["activity_score"], phenotype=a["phenotype"],
                        phenotype_key=a["pheno_key"]))
    d2p.sort(key=lambda r: (r["gene"], r["diplotype"]))
    _write(os.path.join(outdir, "diplotype_to_phenotype.csv"),
           ["gene", "diplotype", "activity_score", "phenotype", "phenotype_key"], d2p)

    # --- table 3: phenotype -> recommendation, the full lookup table ---
    p2r = []
    for (gene, drug), m in RECS.items():
        for key, text in m.items():
            p2r.append(dict(gene=gene, drug=drug, phenotype_key=key,
                            recommendation=("ABSTAIN: " + ABSTAIN_REASON[(gene, drug)])
                            if text is ABSTAIN else text,
                            harm_class=HARM.get((gene, drug, key), ("routine", ""))[0],
                            corpus_source=_source_names(corpus, gene, drug)))
    p2r.sort(key=lambda r: (r["gene"], r["drug"], r["phenotype_key"]))
    _write(os.path.join(outdir, "phenotype_to_recommendation.csv"),
           ["gene", "drug", "phenotype_key", "recommendation", "harm_class",
            "corpus_source"], p2r)

    # --- provenance ---
    prov = dict(corpus_schema_version=corpus["schema_version"],
                corpus_generated_at_utc=corpus["generated_at_utc"],
                n_genes=len(corpus["genes"]), n_cases=len(cases),
                cases_sha256=_sha(os.path.join(indir, "cases_blinded.json")),
                corpus_sha256=_sha(os.path.join(indir, "cpic_guideline_corpus.json")),
                n_abstained=sum(1 for r in rows if r["abstained"] == "true"),
                n_lethal_class=sum(1 for a in audit if a["harm_class"] == "lethal"),
                n_severe_class=sum(1 for a in audit if a["harm_class"] == "severe"))
    with open(os.path.join(outdir, "run_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    return rows, audit, prov


def _source_names(corpus, gene, drug):
    gd = corpus["genes"].get(gene, {})
    out = []
    for s in gd.get("sources", []):
        if any(drug.lower() in x.lower() for x in s.get("drug_scope", [])):
            out.append(f"{s['source_org']}:{s['pa_id']}")
    return ";".join(out) if out else "none for this drug in corpus"


def _write(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(65536), b""):
            h.update(blk)
    return h.hexdigest()


if __name__ == "__main__":
    import sys
    d = sys.argv[1]
    o = sys.argv[2] if len(sys.argv) > 2 else d
    rr, aa, pp = run(d, o)
    print(json.dumps(pp, indent=2))
