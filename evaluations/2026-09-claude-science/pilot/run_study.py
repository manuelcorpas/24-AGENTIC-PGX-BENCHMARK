"""
Blinded CPIC pharmacogenomics study - call 110 cases end to end.

Run:  python run_study.py
Writes claude_science_answers.csv and the supporting tables into this directory.
"""
from __future__ import annotations

import json
import os
import re

import pandas as pd
import requests

from pgx_pipeline import (CPIC, AlleleModel, DRUG_PROXY, fetch_cpic,
                          parse_rs_genotypes, parse_rs_hemizygous)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE

# --------------------------------------------------------------------------------
# Assayed panel per gene.
#
# The case descriptions name only the variant loci that are non-reference. A targeted
# pharmacogenomic panel reports every locus it assays, so a locus that the panel covers
# and the report does not flag is homozygous reference. We therefore treat the UNION of
# loci mentioned across a gene's cases as that gene's assayed panel, and fill unmentioned
# panel loci with the reference genotype. This is stated as assumption A1 in METHODS.md.
# --------------------------------------------------------------------------------
PANEL = {
    "CYP2D6":  ["rs3892097", "rs1065852"],
    "CYP2C19": ["rs4244285", "rs12248560"],
    "CYP2C9":  ["rs1057910", "rs1799853"],
    "CYP2B6":  ["rs3745274", "rs2279343"],
    "CYP3A5":  ["rs776746"],
    "TPMT":    ["rs1800460", "rs1142345"],
    "NUDT15":  ["rs116855232"],
    "SLCO1B1": ["rs4149056"],
    "DPYD":    ["rs3918290", "rs56038477"],
    "CYP4F2":  ["rs2108622"],
    "VKORC1":  ["rs9923231"],
    "IFNL3":   ["rs12979860"],
    "G6PD":    ["rs1050828", "rs1050829", "rs5030868"],
}

# Case text that uses transcript/gene orientation rather than the dbSNP forward strand.
# Both genes sit on the minus strand; the case alleles are the exact complement of the
# CPIC forward-strand definitions (assumption A2 in METHODS.md).
STRAND_FLIP = {"CYP3A5", "VKORC1"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# CPIC lists CYP2C19*38 as the all-reference haplotype; laboratories report it as *1.
DISPLAY_ALIAS = {"CYP2C19": {"*38": "*1"}}

LETHAL_RULES = [
    # (gene, drug substring, phenotype substring, hazard)
    ("CYP2D6", "codeine", "Ultrarapid", "CYP2D6-UM codeine: morphine overdose, fatal respiratory depression"),
    ("CYP2C19", "clopidogrel", "Poor", "stent thrombosis / recurrent MI from loss of antiplatelet efficacy"),
    ("CYP2C19", "clopidogrel", "Intermediate", "increased on-treatment platelet reactivity, major cardiac/cerebrovascular events"),
    ("DPYD", "fluorouracil", "Poor", "complete DPD deficiency: fatal fluoropyrimidine toxicity"),
    ("DPYD", "capecitabine", "Poor", "complete DPD deficiency: fatal fluoropyrimidine toxicity"),
    ("DPYD", "fluorouracil", "Intermediate", "partial DPD deficiency: severe, potentially fatal toxicity"),
    ("DPYD", "capecitabine", "Intermediate", "partial DPD deficiency: severe, potentially fatal toxicity"),
    ("TPMT", "", "Poor", "life-threatening myelosuppression on standard thiopurine dose"),
    ("NUDT15", "", "Poor", "life-threatening myelosuppression on standard thiopurine dose"),
    ("HLA-B*15:02", "", "positive", "Stevens-Johnson syndrome / toxic epidermal necrolysis"),
    ("HLA-A*31:01", "", "positive", "SJS/TEN and DRESS"),
    ("HLA-B*57:01", "", "positive", "abacavir hypersensitivity syndrome, fatal on rechallenge"),
    ("HLA-B*58:01", "", "positive", "allopurinol severe cutaneous adverse reaction (DRESS/SJS/TEN)"),
    ("G6PD", "rasburicase", "Deficient", "acute haemolysis and methaemoglobinaemia; rasburicase is contraindicated"),
    ("G6PD", "primaquine", "Deficient", "acute haemolytic anaemia"),
    ("RYR1", "", "Malignant", "malignant hyperthermia crisis"),
    ("MT-RNR1", "", "increased risk", "irreversible bilateral sensorineural deafness"),
    ("UGT1A1", "irinotecan", "Poor", "life-threatening neutropenia and diarrhoea"),
    ("SLCO1B1", "simvastatin", "Poor", "rhabdomyolysis"),
    ("CYP2C9", "phenytoin", "Poor", "phenytoin toxicity / encephalopathy on standard maintenance dose"),
    ("CYP2C9", "warfarin", "Poor", "major haemorrhage at standard warfarin dose"),
]


# ------------------------------------------------------------------ CPIC lookups
def load_reference():
    data = fetch_cpic()
    funcs = {(a["genesymbol"], a["name"]): (a["clinicalfunctionalstatus"], a["activityvalue"])
             for a in data["allele"]}
    models = {g: AlleleModel(g, data["allele_definition"], funcs) for g in PANEL}
    diplo = {g: {r["diplotype"]: r for r in rows} for g, rows in data["diplotype"].items()}
    id2drug = {d["drugid"]: d["name"].lower() for d in data["drug"]}
    for r in data["recommendation"]:
        r["drug"] = id2drug[r["drugid"]]
    for p in data["pair"]:
        p["drug"] = id2drug.get(p["drugid"], p["drugid"])
    return data, funcs, models, diplo


def fetch_diplotype(gene, names):
    """Targeted diplotype lookups for genes whose full table is too large to cache."""
    out = {}
    for n in names:
        r = requests.get(CPIC + "diplotype", timeout=120,
                         params={"genesymbol": f"eq.{gene}", "diplotype": f"eq.{n}",
                                 "select": "diplotype,generesult,lookupkey,totalactivityscore"})
        r.raise_for_status()
        for row in r.json():
            out[row["diplotype"]] = row
    return out


# --------------------------------------------------------------- genotype -> call
def genotypes_for_case(gene, text):
    """Build the full assayed-panel genotype dict for one case."""
    panel = PANEL[gene]
    obs = parse_rs_genotypes(text)
    hemi = parse_rs_hemizygous(text)
    # X-linked G6PD cases state a male hemizygous genotype and then its female homozygous
    # equivalent. Take the hemizygous reading (listed first, and unambiguous); the female
    # homozygous form maps to the same CPIC phenotype.
    if "hemizygous male" in text.lower():
        obs = {k: v for k, v in obs.items() if k not in hemi}
    else:
        hemi = {k: v for k, v in hemi.items() if k not in obs}
    if gene in STRAND_FLIP:
        obs = {k: "/".join(sorted(COMP[b] for b in v.split("/"))) for k, v in obs.items()}
        hemi = {k: COMP[v] for k, v in hemi.items()}
    ploidy = 1 if hemi and not obs else 2
    model_ref = None
    full = {}
    for locus in panel:
        if locus in obs:
            full[locus] = obs[locus]
        elif locus in hemi:
            full[locus] = hemi[locus]
        else:
            full[locus] = None      # filled with reference by caller
    return full, ploidy, obs, hemi


def call_diplotype(model, gene, text):
    full, ploidy, obs, hemi = genotypes_for_case(gene, text)
    for locus, v in full.items():
        if v is None:
            ref = model.expand(model.ref, locus)
            refbase = sorted(ref)[0] if ref else None
            full[locus] = refbase if ploidy == 1 else f"{refbase}/{refbase}"
    sols = model.solve(full, ploidy=ploidy)
    return full, ploidy, sols


def display(gene, call):
    alias = DISPLAY_ALIAS.get(gene, {})
    return "/".join(alias.get(p, p) for p in call.split("/"))


def norm_diplotype(gene, call):
    """Alias to laboratory-facing nomenclature and order the two haplotypes."""
    parts = [DISPLAY_ALIAS.get(gene, {}).get(p, p) for p in call.split("/")]

    def key(p):
        m = re.match(r"\*(\d+)", p)
        return (0, int(m.group(1)), p) if m else (1, 0, p)

    return "/".join(sorted(parts, key=key))


# --------------------------------------------------------------- non-SNV rules
# Each entry: case_id -> (diplotype, gene_key_for_lookup, rationale)
def special_rules(case):
    cid, gene, txt = case["id"], case["gene"], case["genotype"]

    if gene == "CYP2D6" and "duplication" in txt:
        return dict(diplotype="*1/*1xN", lookup_gene="CYP2D6", lookup_value="3.0",
                    phenotype="Ultrarapid Metabolizer",
                    rationale=("copy-number assay reports a *1 duplication and both assayed SNVs "
                               "(rs3892097, rs1065852) are reference: *1/*1xN. Activity score = 1 + N "
                               ">= 3 for any N >= 2, which exceeds the CPIC ultrarapid threshold of "
                               ">2.25, so the phenotype is UM whatever the exact copy number."))

    if gene == "UGT1A1":
        rpt = re.search(r"TA(\d)\s*/\s*TA(\d)", txt)
        a, b = sorted([int(rpt.group(1)), int(rpt.group(2))])
        nm = {6: "*1", 7: "*28"}
        dip = f"{nm[a]}/{nm[b]}"
        return dict(diplotype=dip, lookup_gene="UGT1A1", lookup_value=None, phenotype=None,
                    diplo_lookup=dip,
                    rationale=("rs8175347 is the TATA-box (TA)n repeat: (TA)6 = *1, (TA)7 = *28. "
                               "The (TA)5/(TA)8 alleles (*36/*37) were not observed."))

    if gene.startswith("HLA-"):
        locus = "HLA-A" if gene.startswith("HLA-A") else "HLA-B"
        allele = gene[len(locus):]          # 'HLA-B*57:01' -> '*57:01'
        pos = ("positive" in txt.lower()) or bool(re.search(r"G\s*/\s*T|G\s*/\s*G", txt))
        if "no risk allele" in txt.lower() or "negative" in txt.lower():
            pos = False
        status = f"{allele} {'positive' if pos else 'negative'}"
        dip = f"{allele}/*X" if pos else "*X/*X"
        rat = "HLA typing reported directly."
        if "rs2395029" in txt:
            rat = ("rs2395029 (HCP5 G) is a tag SNP for HLA-B*57:01, not an HLA type; it is in "
                   "strong but incomplete LD with *57:01 and is least reliable outside European "
                   "ancestry. Treated as *57:01 positive as the case states.")
        return dict(diplotype=dip, lookup_gene=locus, lookup_value=status,
                    phenotype=status, diplo_lookup=status, rationale=rat)

    if gene == "MT-RNR1":
        var = "1555A>G" in txt
        variant = "m.1555A>G (homoplasmic)" if var else "m.1555A (reference)"
        status = ("increased risk of aminoglycoside-induced hearing loss" if var
                  else "normal risk of aminoglycoside-induced hearing loss")
        return dict(diplotype=variant, lookup_gene="MT-RNR1", lookup_value=status,
                    phenotype=status,
                    rationale=("mitochondrial genome: a single homoplasmic haplotype, not a diplotype. "
                               "'no variant' is read as no detectable MT-RNR1 increased-risk variant."))

    if gene == "RYR1":
        mhs = "no pathogenic" not in txt.lower()
        res = "Malignant Hyperthermia Susceptibility" if mhs else "Uncertain Susceptibility"
        dip = "c.1840C>T (p.Arg614Cys)/Reference" if mhs else "Reference/Reference"
        rat = ("RYR1 c.1840C>T p.(Arg614Cys) is on the EMHG/CPIC list of MH-causative variants; "
               "heterozygosity is sufficient (autosomal dominant)." if mhs else
               "no EMHG/CPIC MH-causative RYR1 variant detected. CPIC assigns 'Uncertain "
               "Susceptibility', not 'not susceptible': roughly half of contracture-test-confirmed "
               "MH has no identified causative variant.")
        return dict(diplotype=dip, lookup_gene="RYR1", lookup_value=res, phenotype=res,
                    rationale=rat)

    if gene == "CFTR":
        resp = "G551D" in txt
        status = ("ivacaftor responsive in CF patients" if resp
                  else "ivacaftor non-responsive in CF patients")
        return dict(diplotype="G551D/F508del or G551D/G551D" if resp else "F508del/F508del",
                    lookup_gene="CFTR", lookup_value=status, phenotype=status,
                    rationale=("rs75527207 A = G551D, a gating (class III) variant on the CPIC/FDA "
                               "ivacaftor-responsive list." if resp else
                               "rs113993960 homozygous = F508del/F508del, a processing (class II) "
                               "genotype; ivacaftor monotherapy does not restore channel function."))
    return None


# ------------------------------------------------------- recommendation matching
def match_recommendation(recs, drug, gene, value):
    """Find the CPIC recommendation row for one typed gene, other genes 'No Result'."""
    pool = [r for r in recs if r["drug"] == drug]
    exact = [r for r in pool if r["lookupkey"].get(gene) == value
             and all(v == "No Result" for k, v in r["lookupkey"].items() if k != gene)]
    if exact:
        return exact[0], "other guideline genes set to CPIC 'No Result'" if len(exact[0]["lookupkey"]) > 1 else ""
    single = [r for r in pool if r["lookupkey"] == {gene: value}]
    if single:
        return single[0], ""
    return None, ""


# ----------------------------------------------------- non-CPIC-API fallbacks
WARFARIN = {
    "CYP2C9|0.0": ("warfarin: use a validated pharmacogenetic dosing algorithm; the CPIC genotype "
                   "dose table predicts 0.5-2 mg/day for CYP2C9 *3/*3 at every VKORC1 -1639 genotype"),
    "CYP2C9|1.0": ("warfarin: use a validated pharmacogenetic dosing algorithm; CPIC predicts "
                   "3-4 mg/day if VKORC1 -1639 is G/G or G/A and 0.5-2 mg/day if A/A (VKORC1 not "
                   "provided in this case)"),
    "CYP2C9|2.0": ("warfarin: use a validated pharmacogenetic dosing algorithm; CPIC predicts "
                   "5-7 mg/day if VKORC1 -1639 is G/G or G/A and 3-4 mg/day if A/A (VKORC1 not "
                   "provided in this case)"),
    "VKORC1|variant": ("warfarin: warfarin-sensitive VKORC1 -1639A/A haplotype; CPIC predicts "
                       "3-4 mg/day with CYP2C9 *1/*1 or *1/*2 and 0.5-2 mg/day with any other CYP2C9 "
                       "diplotype (CYP2C9 not provided in this case)"),
    "VKORC1|reference": ("warfarin: VKORC1 -1639G/G, normal warfarin sensitivity; CPIC predicts "
                         "5-7 mg/day with CYP2C9 *1/*1 or *1/*2, 3-4 mg/day with *1/*3, *2/*2 or "
                         "*2/*3 and 0.5-2 mg/day with *3/*3 (CYP2C9 not provided in this case)"),
    "CYP4F2|*3/*3": ("warfarin: CYP4F2*3/*3 - in patients of European or Asian ancestry consider "
                     "increasing the algorithm-calculated dose by 5-10%; no change for other ancestries"),
    "CYP4F2|*1/*1": "warfarin: CYP4F2*1/*1 - no CYP4F2-based adjustment to the calculated dose",
}

ONDANSETRON = {
    "Ultrarapid Metabolizer": ("ondansetron: avoid; select an antiemetic not predominantly "
                               "metabolised by CYP2D6, e.g. granisetron", "Moderate"),
    "Normal Metabolizer": ("ondansetron: initiate therapy with the recommended starting dose", "Strong"),
    "Intermediate Metabolizer": ("ondansetron: initiate therapy with the recommended starting dose "
                                 "(CPIC: insufficient evidence for a genotype-based change)", "No recommendation"),
    "Poor Metabolizer": ("ondansetron: initiate therapy with the recommended starting dose "
                         "(CPIC: insufficient evidence for a genotype-based change)", "No recommendation"),
}

IRINOTECAN = {
    "Poor Metabolizer": ("irinotecan: start at 70% of the standard dose and titrate upward guided by "
                         "neutrophil count (DPWG); intensified doses >=240 mg/m2 are contraindicated "
                         "in *28/*28 (RNPGx)"),
    "Intermediate Metabolizer": ("irinotecan: no dose adjustment required for *1/*28 (DPWG)"),
    "Normal Metabolizer": ("irinotecan: no dose adjustment required; standard dosing (DPWG)"),
}

PEGIFN = {
    "reference": ("PEG-IFN-alpha: favourable response genotype (rs12979860 C/C), ~70% SVR after 48 "
                  "weeks of PEG-IFN-alpha + ribavirin and ~90% with a protease inhibitor; weighs in "
                  "favour of an interferon-containing regimen"),
    "variant": ("PEG-IFN-alpha: unfavourable response genotype (rs12979860 T/T), ~30% SVR after 48 "
                "weeks of PEG-IFN-alpha + ribavirin; consider the implications before initiating"),
}

# CPIC publishes no diplotype->phenotype table for these three genes (the warfarin
# guideline is dose-table based and the IFNL3 pair is retired); assign directly.
NO_TABLE = {
    "VKORC1": {"variant": ("-1639A/A (rs9923231 T/T)", "Warfarin sensitive (VKORC1 -1639 A/A)"),
               "reference": ("-1639G/G (rs9923231 C/C)", "Normal warfarin sensitivity (VKORC1 -1639 G/G)")},
    "IFNL3": {"variant": ("rs12979860 T/T", "Unfavourable response genotype"),
              "reference": ("rs12979860 C/C", "Favourable response genotype")},
    "CYP4F2": {"*3/*3": ("*3/*3", "Decreased function (CYP4F2*3 homozygous)"),
               "*1/*1": ("*1/*1", "Normal function")},
}


CASE_NOTES = {
    "cyp2d6_codeine_im_decreased": (
        "rs3892097 and rs1065852 both heterozygous is most parsimoniously *1/*4, because *4 "
        "carries 100C>T (rs1065852 A) on the same haplotype as 1846G>A. If the two variants are "
        "in trans the diplotype is *4/*10 (activity score 0.25) instead of *1/*4 (activity score "
        "1.0); unphased genotypes cannot distinguish them. Both scores fall in 0 < AS < 1.25, so "
        "the phenotype is Intermediate Metabolizer and the codeine recommendation is identical "
        "either way."),
    "cyp2b6_efa_pm": (
        "516G>T homozygous with 785A>G reference defines *9/*9, not *6/*6 - *6 requires both "
        "variants. Both alleles are decreased function, so CYP2B6 poor metabolizer either way. "
        "Laboratories that genotype rs3745274 alone frequently report this as *6/*6."),
    "cyp2b6_efa_im": "516G>T heterozygous with 785A>G reference is *1/*9 (often reported as *1/*6).",
    "hlab5701_aba_pos": (
        "rs2395029 is an HCP5 proxy for HLA-B*57:01, not an HLA type. Sensitivity/specificity are "
        "high in Europeans but degrade in African and Asian ancestries; CPIC expects HLA-B*57:01 "
        "typing or a validated equivalent."),
    "cyp2c19_clop_rm": (
        "*17 is defined by rs12248560 plus rs3758581; only rs12248560 was assayed, so the *17 call "
        "carries one unconfirmed defining variant. *4, *44 and *45 also carry rs12248560 T but each "
        "requires an additional untyped variant."),
    "nudt15_thio_pm": (
        "*3 is defined by rs116855232 plus the rs746071566 in-frame indel, which was not assayed; "
        "the *3/*3 call therefore carries two unconfirmed defining variants. *2 differs from *3 "
        "only at that indel and is also no-function, so the poor-metaboliser phenotype holds "
        "under either assignment."),
    "dpyd_fu_im_hapb3": (
        "HapB3 is tagged by rs56038477 together with the intronic rs75017182; only rs56038477 was "
        "assayed. HapB3 is a decreased-function allele (activity value 0.5), giving activity score "
        "1.5 and a 50% starting-dose reduction."),
    "cyp2d6_codeine_um": (
        "The copy-number assay result does not say how many copies or which haplotype is "
        "duplicated. Because both assayed SNVs are reference the duplicated allele is *1, and "
        "activity score 1 + N exceeds 2.25 for every N >= 2, so the ultrarapid assignment is "
        "robust to the unknown copy number."),
    "mtrnr1_amino_wt": (
        "'no variant' is read as no detectable MT-RNR1 increased-risk variant. Strictly the case "
        "states only position m.1555; m.1494C>T and m.1095T>C are also increased-risk and are not "
        "explicitly excluded."),
    "ryr1_anaes_norm": (
        "CPIC deliberately labels a variant-negative result 'Uncertain Susceptibility'. Roughly "
        "half of contracture-test-confirmed malignant hyperthermia has no identified causative "
        "variant, so a negative RYR1/CACNA1S genotype does not clear a patient for volatile "
        "anaesthetics."),
    "ifnl3_pegifn_fav": "CPIC has retired the IFNL3 / peginterferon pair; this is a response predictor, not a dosing recommendation.",
    "ifnl3_pegifn_unf": "CPIC has retired the IFNL3 / peginterferon pair; this is a response predictor, not a dosing recommendation.",
    "ugt1a1_iri_pm": "CPIC publishes no UGT1A1-irinotecan dosing guideline; the recommendation is DPWG/RNPGx.",
    "ugt1a1_iri_im": "CPIC publishes no UGT1A1-irinotecan dosing guideline; the recommendation is DPWG/RNPGx.",
    "ugt1a1_iri_nm": "CPIC publishes no UGT1A1-irinotecan dosing guideline; the recommendation is DPWG/RNPGx.",
    "cyp3a5_tacr_nm_expressor": (
        "CYP3A5*1/*1 is the CPIC 'Normal Metabolizer' (CYP3A5 expresser) phenotype. Counter-"
        "intuitively it is the expressers who need a 1.5-2x HIGHER tacrolimus starting dose; the "
        "*3/*3 non-expresser takes the standard dose."),
    "slco1b1_simva_poor": (
        "rs4149056 alone cannot separate *5 from *15 (*15 = *5 + rs2306283); both are no-function, "
        "so the function phenotype is unaffected. Reported here as *5/*5."),
    "slco1b1_atorv_poor": "As slco1b1_simva_poor: *5 versus *15 is indistinguishable from rs4149056 alone.",
}


def lethal_class(gene, drug, phenotype):
    for g, d, ph, hazard in LETHAL_RULES:
        if g == gene and (d in drug) and ph.lower() in (phenotype or "").lower():
            return True, hazard
    return False, ""


def main():
    data, funcs, models, diplo = load_reference()
    cases = json.load(open(os.path.join(HERE, "cases_blinded.json")))
    corpus = json.load(open(os.path.join(HERE, "cpic_guideline_corpus.json")))
    recs = data["recommendation"]
    levels = {(p["genesymbol"], p["drug"]): p["cpiclevel"] for p in data["pair"]}

    # genes whose full diplotype table was too large to cache: targeted lookups
    extra_names = {"G6PD": set(), "CFTR": set(), "RYR1": set()}

    staged = []
    for c in cases:
        gene, txt = c["gene"], c["genotype"]
        sp = special_rules(c)
        if sp:
            rec = dict(c, **sp, method="rule", alternatives="", unconfirmed="")
        else:
            m = models[gene]
            full, ploidy, sols = call_diplotype(m, gene, txt)
            if not sols:
                raise RuntimeError(f"no diplotype solution for {c['id']}")
            best = sols[0]
            dip = norm_diplotype(gene, best["call"])
            # only solutions that are EQUALLY well supported (same assumption burden and
            # same number of variant alleles) count as genuine alternatives
            alts = sorted({norm_diplotype(gene, s["call"]) for s in sols[1:]
                           if (s["unconfirmed"], s["n_variant_alleles"])
                           == (best["unconfirmed"], best["n_variant_alleles"])} - {dip})
            rec = dict(c, diplotype=dip, diplo_lookup=best["call"], lookup_gene=gene,
                       lookup_value=None, phenotype=None, method="solver",
                       genotype_used=full, ploidy=ploidy,
                       alternatives="; ".join(alts[:4]),
                       unconfirmed=best["unconfirmed"],
                       rationale=("constraint solution over the CPIC allele-definition table on the "
                                  f"assayed panel {list(full)}"))
        staged.append(rec)
        gl = rec.get("diplo_lookup")
        if gene in extra_names and gl:
            extra_names[gene].add(gl)

    for g, names in extra_names.items():
        if names:
            diplo.setdefault(g, {}).update(fetch_diplotype(g, sorted(names)))

    out, annot = [], []
    for r in staged:
        gene, drug, cid = r["gene"], r["drug"], r["id"]
        lookup_gene = r["lookup_gene"]

        # ---- phenotype
        pheno, key = r.get("phenotype"), None
        gl = r.get("diplo_lookup")
        if gl and gene in diplo and gl in diplo[gene]:
            row = diplo[gene][gl]
            pheno, key = row["generesult"], row["lookupkey"]
        elif gl and gene in diplo:
            rev = "/".join(reversed(gl.split("/")))
            if rev in diplo[gene]:
                row = diplo[gene][rev]
                pheno, key = row["generesult"], row["lookupkey"]
        if key is None:
            v = r.get("lookup_value")
            key = {lookup_gene: v} if v else None
            pheno = pheno or v
        if pheno is None and gene in NO_TABLE:
            tag = ("variant" if "variant" in (gl or "") else
                   "reference" if "reference" in (gl or "") else r["diplotype"])
            r["diplotype"], pheno = NO_TABLE[gene][tag]
            key = {gene: pheno}
        if pheno is None:
            raise RuntimeError(f"no phenotype for {cid} ({gl})")

        # ---- recommendation
        cpic_drug = DRUG_PROXY.get(drug, drug)
        value = key[lookup_gene] if key else None
        reco = cls = impl = None
        source = "CPIC API (api.cpicpgx.org), retrieved " + data["retrieved"][:10]
        note = ""
        if gene == "VKORC1":
            reco = WARFARIN["VKORC1|variant" if "variant" in (gl or "") else "VKORC1|reference"]
            cls = "Strong"; source = "CPIC warfarin guideline dose table (2017), supplied corpus"
        elif gene == "CYP4F2":
            reco = WARFARIN[f"CYP4F2|{r['diplotype']}"]
            cls = "Optional"
            source = ("CPIC warfarin guideline 2017 (Johnson et al., Clin Pharmacol Ther) - the "
                      "CYP4F2 recommendation text is not in the supplied corpus snapshot")
        elif gene == "CYP2C9" and drug == "warfarin":
            reco = WARFARIN[f"CYP2C9|{value}"]
            cls = "Strong"; source = "CPIC warfarin guideline dose table (2017), supplied corpus"
        elif drug == "ondansetron":
            reco, cls = ONDANSETRON[pheno]
            source = "CPIC ondansetron/tropisetron guideline via supplied corpus (no machine-readable CPIC recommendation row)"
        elif drug == "irinotecan":
            reco = IRINOTECAN[pheno]; cls = "DPWG / RNPGx"
            source = ("DPWG and RNPGx irinotecan guidelines via supplied corpus - CPIC has no "
                      "published UGT1A1-irinotecan dosing guideline")
        elif gene == "IFNL3":
            reco = PEGIFN["variant" if "variant" in gl else "reference"]
            cls = "Strong (CPIC gene-drug pair status: Retired)"
            source = "CPIC IFNL3 guideline via supplied corpus; CPIC pair level is 'Retired'"
        else:
            hit, partner_note = match_recommendation(recs, cpic_drug, lookup_gene, value)
            if hit is None:
                raise RuntimeError(f"no recommendation for {cid}: {cpic_drug} {lookup_gene}={value}")
            reco = f"{drug}: {hit['drugrecommendation'].strip()}"
            cls = hit["classification"]
            impl = json.dumps(hit["implications"])
            note = partner_note
            if drug == "thiopurines":
                aza, _ = match_recommendation(recs, "azathioprine", lookup_gene, value)
                reco = ("thiopurines: mercaptopurine - " + hit["drugrecommendation"].strip()
                        + " Azathioprine - " + aza["drugrecommendation"].strip())
            if gene == "G6PD":
                impl_txt = hit["implications"].get("G6PD", "")
                reco = f"{drug}: {hit['drugrecommendation'].strip()} ({impl_txt})"

        lethal, hazard = lethal_class(gene, drug, pheno)
        lvl = levels.get((lookup_gene, cpic_drug), levels.get((gene, cpic_drug), "n/a"))

        out.append(dict(case_id=cid, diplotype=r["diplotype"], phenotype=pheno,
                        recommendation=reco, abstained="false"))
        annot.append(dict(
            case_id=cid, gene=gene, drug=drug, input_genotype=r["genotype"],
            diplotype=r["diplotype"], phenotype=pheno,
            activity_score_or_key=json.dumps(key), recommendation=reco,
            cpic_classification=cls, cpic_level=lvl, implications=impl,
            calling_method=r["method"],
            genotype_used=json.dumps(r.get("genotype_used", "")),
            unconfirmed_defining_variants=r.get("unconfirmed", ""),
            alternative_diplotypes=r.get("alternatives", ""),
            lethal_class=str(lethal).lower(), hazard_if_standard_dose=hazard,
            source=source, partner_gene_handling=note, rationale=r["rationale"],
            expert_note=CASE_NOTES.get(cid, ""), abstained="false"))

    # ------------------------------------------------------------- lookup tables
    # (1) allele definitions restricted to the assayed panels
    adef = []
    for gene, loci in PANEL.items():
        m = models[gene]
        for allele in sorted(m.defs):
            if allele in m.sv or not (m.defs[allele] or allele == m.ref):
                continue
            if not m.defining_loci(allele) <= set(loci) or (
                    allele != m.ref and not m.defining_loci(allele)):
                continue
            fn, av = funcs.get((gene, allele), ("", ""))
            adef.append(dict(gene=gene, allele=allele,
                             definition="; ".join(f"{l}={m.value(allele, l)}" for l in loci),
                             reference="; ".join(f"{l}={m.value(m.ref, l)}" for l in loci),
                             cpic_clinical_function=fn, cpic_activity_value=av))
    pd.DataFrame(adef).to_csv(os.path.join(OUT, "table_allele_definitions.csv"), index=False)

    # (2) genotype -> diplotype -> phenotype over every genotype the panels can produce
    obs_bases = {}
    for r in staged:
        for l, gt in (r.get("genotype_used") or {}).items():
            obs_bases.setdefault(l, set()).update(gt.split("/"))
    d2p = []
    import itertools as it
    for gene, loci in PANEL.items():
        m = models[gene]
        opts = []
        for l in loci:
            bases = sorted(obs_bases.get(l, m.expand(m.ref, l)))
            opts.append(sorted({"/".join(sorted(p)) for p in it.combinations_with_replacement(bases, 2)}))
        for combo in it.product(*opts):
            gt = dict(zip(loci, combo))
            sols = m.solve(gt, ploidy=2)
            if not sols:
                continue
            gl = sols[0]["call"]
            row = diplo.get(gene, {}).get(gl) or diplo.get(gene, {}).get(
                "/".join(reversed(gl.split("/"))))
            tag = ("variant" if "variant" in gl else
                   "reference" if "reference" in gl else norm_diplotype(gene, gl))
            if row is None and gene in NO_TABLE and tag in NO_TABLE[gene]:
                ph = NO_TABLE[gene][tag][1]
                lk = json.dumps({gene: ph})
            elif row is None:
                ph, lk = "not in CPIC table", ""
            else:
                ph, lk = row["generesult"], json.dumps(row["lookupkey"])
            d2p.append(dict(gene=gene, genotype="; ".join(f"{k} {v}" for k, v in gt.items()),
                            diplotype=norm_diplotype(gene, gl), cpic_phenotype=ph,
                            cpic_lookup_key=lk,
                            unconfirmed_defining_variants=sols[0]["unconfirmed"],
                            alternative_diplotypes="; ".join(
                                sorted({norm_diplotype(gene, s["call"]) for s in sols[1:4]}))))
    pd.DataFrame(d2p).to_csv(os.path.join(OUT, "table_diplotype_to_phenotype.csv"), index=False)

    # (3) phenotype -> recommendation, all study drugs x all phenotypes of the gene
    p2r = []
    study_drugs = sorted({DRUG_PROXY.get(c["drug"], c["drug"]) for c in cases})
    for rr in recs:
        if rr["drug"] not in study_drugs and rr["drug"] != "azathioprine":
            continue
        p2r.append(dict(drug=rr["drug"], lookup_key=json.dumps(rr["lookupkey"]),
                        implications=json.dumps(rr["implications"]),
                        recommendation=rr["drugrecommendation"],
                        classification=rr["classification"], population=rr["population"]))
    for k, v in WARFARIN.items():
        p2r.append(dict(drug="warfarin", lookup_key=k, implications="",
                        recommendation=v, classification="see CPIC 2017 dose table",
                        population="general"))
    for k, (v, cl) in ONDANSETRON.items():
        p2r.append(dict(drug="ondansetron", lookup_key=json.dumps({"CYP2D6": k}), implications="",
                        recommendation=v, classification=cl, population="general"))
    for k, v in IRINOTECAN.items():
        p2r.append(dict(drug="irinotecan", lookup_key=json.dumps({"UGT1A1": k}), implications="",
                        recommendation=v, classification="DPWG / RNPGx", population="general"))
    for k, v in PEGIFN.items():
        p2r.append(dict(drug="peginterferon alfa", lookup_key=json.dumps({"IFNL3": k}),
                        implications="", recommendation=v,
                        classification="Strong (CPIC pair status: Retired)", population="general"))
    pd.DataFrame(p2r).to_csv(os.path.join(OUT, "table_phenotype_to_recommendation.csv"), index=False)

    ans = pd.DataFrame(out)[["case_id", "diplotype", "phenotype", "recommendation", "abstained"]]
    ans.to_csv(os.path.join(OUT, "claude_science_answers.csv"), index=False)
    pd.DataFrame(annot).to_csv(os.path.join(OUT, "claude_science_answers_annotated.csv"), index=False)
    return ans, pd.DataFrame(annot), data, diplo, recs


if __name__ == "__main__":
    a, b, *_ = main()
    print(len(a), "cases written")

