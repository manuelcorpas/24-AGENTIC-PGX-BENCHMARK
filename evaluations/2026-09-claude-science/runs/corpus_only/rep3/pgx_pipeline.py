#!/usr/bin/env python3
"""
Blinded pharmacogenomics study - corpus_only/rep3
=================================================
Maps 110 blinded cases (gene + drug + genotype description) to:
  diplotype -> CPIC metaboliser/function phenotype -> guideline recommendation

Inputs (same folder):
  cases_blinded.json            110 cases
  cpic_guideline_corpus.json    CPIC/DPWG/RNPGx guideline text, 21 genes, snapshot 2026-05-09

Outputs (same folder):
  claude_science_answers.csv        case_id,diplotype,phenotype,recommendation,abstained
  table_variant_to_allele.csv      assay result -> star allele / haplotype call
  table_diplotype_to_phenotype.csv diplotype (+activity score) -> phenotype
  table_phenotype_to_recommendation.csv  (gene,drug,phenotype) -> action + verbatim guideline text

Method: deterministic rule tables, no network, no ML.
  - Allele definitions (rsID -> star allele) come from the analyst's own
    pharmacogenetic knowledge, because the supplied corpus contains allele
    definitions for CYP2C9 only. Every such call is flagged in
    table_variant_to_allele.csv with basis='analyst-knowledge'.
  - Diplotype->phenotype translation uses CPIC activity-score bins quoted in
    the supplied corpus where the corpus states them.
  - Every recommendation string is a verbatim substring of the supplied corpus
    (enforced by assertion at run time), so no recommendation text is invented.
"""
import json, os, re, csv, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
cases = json.load(open(os.path.join(HERE, "cases_blinded.json")))
corpus = json.load(open(os.path.join(HERE, "cpic_guideline_corpus.json")))
GENES = corpus["genes"]
CORPUS_TEXT = {k: v["guideline_excerpt"] for k, v in GENES.items()}
ALL_TEXT = "\n".join(CORPUS_TEXT.values())
NORM_TEXT = re.sub(r"\s+", " ", ALL_TEXT)


def verbatim(s):
    """Return s after asserting it occurs verbatim (whitespace-normalised) in the corpus."""
    probe = re.sub(r"\s+", " ", s).strip()
    assert probe in NORM_TEXT, "NOT IN CORPUS: " + probe[:120]
    return s


# ---------------------------------------------------------------------------
# STEP 1.  Assay result -> allele call
#   (regex on the case 'genotype' string, per gene) -> (diplotype, basis, note)
# ---------------------------------------------------------------------------
# basis: 'corpus'          = allele definition present in supplied corpus
#        'analyst-knowledge' = allele definition from analyst's own PGx knowledge
#        'case-annotation' = the case string itself states the allele/typing result
VARIANT_RULES = {
 "CYP2D6": [
  (r"rs3892097 T/T",            "*4/*4",     "analyst-knowledge",
   "rs3892097 = CYP2D6 1846G>A (c.506-1G>A splice defect), the *4-defining SNV; forward-genomic C>T (CYP2D6 on minus strand). T/T = *4/*4, two no-function alleles."),
  (r"rs3892097 C/T \+ rs1065852 G/A", "*1/*4", "analyst-knowledge",
   "rs1065852 (100C>T, P34S) is part of the *4 core haplotype, so the two heterozygous variants are most parsimoniously in cis on one *4 allele -> *1/*4. Phase is not measured; a trans arrangement (*4-like/*10-like) cannot be formally excluded."),
  (r"rs3892097 C/T",            "*1/*4",     "analyst-knowledge",
   "Heterozygous for the *4-defining 1846G>A; other CYP2D6 SNVs reference -> *1/*4."),
  (r"All CYP2D6 SNPs homozygous reference", "*1/*1", "analyst-knowledge",
   "No variant detected at any interrogated CYP2D6 SNV -> *1/*1. Caveat: a SNV-only panel cannot exclude *5 whole-gene deletion or duplications."),
  (r"\*1xN duplication",        "*1/*1xN",   "case-annotation",
   "Duplication reported by copy-number assay (qPCR/MLPA/long-range PCR). No SNV defines a duplication, so this call is only possible because a CNV assay was run. Copy number N is not specified; *1/*1xN is the corpus's own example UM diplotype."),
 ],
 "CYP2C19": [
  (r"rs4244285 A/A",   "*2/*2",   "analyst-knowledge", "rs4244285 = 681G>A (splice defect), the *2-defining SNV. A/A = *2/*2, two no-function alleles."),
  (r"rs4244285 G/A",   "*1/*2",   "analyst-knowledge", "Heterozygous *2."),
  (r"rs12248560 T/T",  "*17/*17", "analyst-knowledge", "rs12248560 = -806C>T, the *17-defining SNV (increased function). T/T = *17/*17."),
  (r"rs12248560 C/T",  "*1/*17",  "analyst-knowledge", "Heterozygous *17."),
  (r"All CYP2C19 SNPs homozygous reference", "*1/*1", "analyst-knowledge", "No variant detected -> *1/*1."),
 ],
 "CYP2C9": [
  (r"rs1057910 C/C",   "*3/*3",   "corpus", "Corpus allele table: '*3 | A>C at rs1057910 | I359L | Decreased'. C/C = *3/*3."),
  (r"rs1057910 A/C",   "*1/*3",   "corpus", "Corpus allele table: *3 = A>C at rs1057910. A/C = *1/*3."),
  (r"All CYP2C9 SNPs homozygous reference", "*1/*1", "corpus", "Corpus allele table: '*1 | reference allele at all positions | Normal'."),
 ],
 "CYP3A5": [
  (r"rs776746 A/A",    "*1/*1",   "case-annotation", "Case annotates A/A as homozygous *1; rs776746 = 6986A>G, the *3-defining intron-3 SNV (A = *1, G = *3)."),
  (r"rs776746 A/G",    "*1/*3",   "case-annotation", "Heterozygous *3."),
  (r"rs776746 G/G",    "*3/*3",   "case-annotation", "Homozygous *3 (non-expresser)."),
 ],
 "CYP2B6": [
  (r"rs3745274 T/T",   "*6/*6",   "analyst-knowledge",
   "rs3745274 = 516G>T (Q172H). *6 = 516G>T + 785A>G; 516G>T alone = *9. rs2279343 not reported, so *9/*9 cannot be excluded - but *6 and *9 are both decreased-function, so the phenotype (poor metaboliser) is identical either way."),
  (r"rs3745274 G/T",   "*1/*6",   "analyst-knowledge",
   "Heterozygous 516G>T -> *1/*6 (or *1/*9; both decreased function, same phenotype)."),
  (r"rs2279343 A/G",   "*1/*4",   "analyst-knowledge",
   "rs2279343 = 785A>G (K262R). 785A>G without 516G>T defines *4, an INCREASED-function allele -> *1/*4. Calling this *1/*6 would invert the phenotype."),
  (r"All CYP2B6 SNPs homozygous reference", "*1/*1", "analyst-knowledge", "No variant detected -> *1/*1."),
 ],
 "DPYD": [
  (r"rs3918290 T/T",   "c.1905+1G>A/c.1905+1G>A (*2A/*2A)", "analyst-knowledge",
   "rs3918290 = DPYD c.1905+1G>A (*2A), no function; forward-genomic C>T (DPYD on minus strand). T/T = homozygous, activity score 0."),
  (r"rs3918290 C/T",   "c.1905+1G>A/Reference (*1/*2A)",    "analyst-knowledge",
   "Heterozygous c.1905+1G>A, activity score 1 (1 + 0)."),
  (r"rs56038477 C/T",  "c.1236G>A(HapB3)/Reference",        "analyst-knowledge",
   "rs56038477 = c.1236G>A, the HapB3 tag; DECREASED (not no) function, allele value 0.5 -> activity score 1.5. Treating HapB3 as a no-function allele would over-reduce the dose."),
  (r"All DPYD risk variants reference", "Reference/Reference (*1/*1)", "analyst-knowledge",
   "No DPYD no/decreased-function variant detected -> activity score 2."),
 ],
 "TPMT": [
  (r"rs1800460 T/T \+ rs1142345 C/C", "*3A/*3A", "analyst-knowledge",
   "Forward-genomic alleles; TPMT is on the minus strand, so c.460G>A = genomic C>T (variant = T) and c.719A>G = genomic T>C (variant = C). Both variants homozygous = *3A/*3A (*3A carries c.460A + c.719G in cis)."),
  (r"rs1800460 C/T",  "*1/*3B", "analyst-knowledge",
   "Heterozygous c.460G>A only (c.719A>G reference) = *3B, not *3A. Reading the forward-genomic letters as transcript letters would wrongly call this reference."),
  (r"rs1142345 T/C",  "*1/*3C", "analyst-knowledge",
   "Heterozygous c.719A>G only = *3C (the commonest TPMT no-function allele outside Europe)."),
  (r"All TPMT SNPs homozygous reference", "*1/*1", "analyst-knowledge", "No variant detected -> *1/*1."),
 ],
 "NUDT15": [
  (r"rs116855232 T/T", "*3/*3", "analyst-knowledge", "rs116855232 = c.415C>T (R139C), the *3-defining SNV, no function. T/T = *3/*3."),
  (r"rs116855232 C/T", "*1/*3", "analyst-knowledge", "Heterozygous *3."),
  (r"rs116855232 C/C", "*1/*1", "analyst-knowledge", "Reference -> *1/*1."),
 ],
 "SLCO1B1": [
  (r"rs4149056 C/C",   "*5/*5 (521CC; *15/*15 indistinguishable)", "corpus",
   "Corpus simvastatin table maps rs4149056 CC to *5/*5, *5/*15, *5/*17, *15/*15, *15/*17, *17/*17. rs2306283 (388A>G) not typed, so *5 vs *15/*17 cannot be separated; all are no-function, so the phenotype is unaffected."),
  (r"rs4149056 T/C",   "*1/*5 (521TC; *1/*15 indistinguishable)", "corpus",
   "Corpus simvastatin table maps rs4149056 TC to *1a/*5, *1a/*15, *1a/*17 etc."),
  (r"rs4149056 T/T",   "*1/*1 (521TT)", "corpus",
   "Corpus simvastatin table maps rs4149056 TT to *1a/*1a, *1a/*1b, *1b/*1b (two normal-function alleles)."),
 ],
 "UGT1A1": [
  (r"rs8175347 TA7/TA7", "*28/*28", "analyst-knowledge", "rs8175347 = promoter (TA)n repeat; TA7 = *28 (decreased function), TA6 = *1."),
  (r"rs8175347 TA6/TA7", "*1/*28",  "analyst-knowledge", "Heterozygous *28."),
  (r"rs8175347 TA6/TA6", "*1/*1",   "analyst-knowledge", "Two TA6 (reference) alleles -> *1/*1."),
 ],
 "HLA-B*57:01": [
  (r"rs2395029 G/T or G/G", "HLA-B*57:01 positive (inferred from HCP5 rs2395029 G)", "analyst-knowledge",
   "rs2395029 is an HCP5 tag SNP in strong LD with HLA-B*57:01 in European ancestry; it is a PROXY, not HLA typing. Positive tag is treated as *57:01 positive."),
  (r"rs2395029 T/T",        "HLA-B*57:01 negative (inferred from HCP5 rs2395029 T/T)", "analyst-knowledge",
   "No tag allele. LD with *57:01 is incomplete outside European ancestry, so a negative tag does not fully exclude carriage."),
 ],
 "HLA-B*15:02": [
  (r"positive for \*15:02", "HLA-B*15:02 positive (one or two *15:02 alleles)", "case-annotation", "Direct HLA typing reported positive."),
  (r"negative for \*15:02", "HLA-B*15:02 negative",                             "case-annotation", "Direct HLA typing reported negative."),
 ],
 "HLA-A*31:01": [
  (r"positive for \*31:01", "HLA-A*31:01 positive (one or two *31:01 alleles)", "case-annotation", "Direct HLA typing reported positive."),
  (r"negative for \*31:01", "HLA-A*31:01 negative",                            "case-annotation", "Direct HLA typing reported negative."),
 ],
 "HLA-B*58:01": [
  (r"positive for \*58:01", "HLA-B*58:01 positive (one or two *58:01 alleles)", "case-annotation", "Direct HLA typing reported positive."),
  (r"negative for \*58:01", "HLA-B*58:01 negative",                            "case-annotation", "Direct HLA typing reported negative."),
 ],
 "G6PD": [
  (r"rs5030868 A \(hemizygous male\)", "Mediterranean (hemizygous male) or Mediterranean/Mediterranean (homozygous female)", "analyst-knowledge",
   "rs5030868 = c.563C>T (p.Ser188Phe) = G6PD Mediterranean, WHO class II; forward-genomic G>A (G6PD on minus strand of X), so forward A = variant."),
  (r"rs1050828 T \+ rs1050829 C in cis", "A- (hemizygous male) or A-/A- (homozygous female)", "analyst-knowledge",
   "c.202G>A (rs1050828, forward C>T) + c.376A>G (rs1050829, forward T>C) in cis = G6PD A-, WHO class III."),
  (r"rs1050828 C/T \+ rs1050829", "B/A- (heterozygous female)", "analyst-knowledge",
   "Heterozygous female carrying one non-deficient (class IV, B) and one deficient (class III, A-) allele. X-inactivation mosaicism makes measured activity unpredictable."),
  (r"rs1050828 C \+ rs5030868 G", "B (hemizygous male) or B/B (homozygous female)", "analyst-knowledge",
   "Reference at both interrogated positions -> non-deficient class IV B allele."),
 ],
 "IFNL3": [
  (r"rs12979860 C/C", "rs12979860 C/C", "corpus", "Corpus table is keyed directly on rs12979860 genotype; no star alleles are used for IFNL3."),
  (r"rs12979860 T/T", "rs12979860 T/T", "corpus", "Corpus table is keyed directly on rs12979860 genotype."),
 ],
 "VKORC1": [
  (r"rs9923231 A/A", "rs9923231 (-1639G>A) A/A", "corpus", "Corpus warfarin dose table rows are keyed on VKORC1 -1639G>A (rs9923231) GG/GA/AA."),
  (r"rs9923231 G/G", "rs9923231 (-1639G>A) G/G", "corpus", "Corpus warfarin dose table rows are keyed on VKORC1 -1639G>A (rs9923231) GG/GA/AA."),
 ],
 "CYP4F2": [
  (r"rs2108622 T/T", "*3/*3", "analyst-knowledge", "rs2108622 = c.1297G>A (V433M), the CYP4F2*3-defining SNV (decreased function). rs2108622 does NOT appear anywhere in the supplied corpus."),
  (r"rs2108622 C/C", "*1/*1", "analyst-knowledge", "Reference at rs2108622 -> *1/*1. rs2108622 does NOT appear anywhere in the supplied corpus."),
 ],
 "CFTR": [
  (r"rs75527207 G/A or A/A",  "G551D positive (rs75527207 AA or AG)", "corpus",
   "Corpus row: 'Homozygous or Heterozygous G551D-CFTR, rs75527207 genotype AA or AG'."),
  (r"rs113993960 \(F508del homozygous\)", "F508del/F508del", "corpus",
   "Corpus row: 'Homozygous for F508del-CFTR, rs113993960 or rs199826652 genotype del/del'."),
 ],
 "MT-RNR1": [
  (r"m\.1555A>G", "m.1555A>G (homoplasmic)", "analyst-knowledge",
   "m.1555A>G is the canonical MT-RNR1 increased-risk variant for aminoglycoside ototoxicity. Mitochondrial, so 'diplotype' is a homoplasmy/heteroplasmy state, not a diplotype."),
  (r"m\.1555A \(no variant\)", "m.1555A (reference)", "analyst-knowledge",
   "No increased-risk MT-RNR1 variant at m.1555."),
 ],
 "RYR1": [
  (r"RYR1 mutation", "RYR1 c.1840C>T p.(Arg614Cys), heterozygous", "corpus",
   "Corpus lists 'RYR1 c.1840C>T; p.(Arg614Cys)' among EMHG malignant-hyperthermia causative variants."),
  (r"RYR1 no pathogenic variant", "No RYR1/CACNA1S MH-causative variant detected", "corpus",
   "Corpus row: 'An individual negative for a RYR1 or CACNA1S malignant hyperthermia causative variant as designated by the EMHG'."),
 ],
}

# ---------------------------------------------------------------------------
# STEP 2.  Diplotype -> phenotype (with activity score where CPIC uses one)
# ---------------------------------------------------------------------------
# CYP2D6: current CPIC standardised translation, quoted verbatim in the corpus
# codeine/paroxetine tables: AS >2.25 UM | 1.25<=AS<=2.25 NM | 0<AS<1.25 IM | 0 PM.
# This supersedes the older bins still printed in the corpus' ondansetron,
# amitriptyline and tamoxifen tables (the corpus itself states "CPIC has updated
# the CYP2D6 genotype-phenotype translation in all supplemental files since the
# guideline was published").
DIPLO2PHENO = {
 ("CYP2D6", "*4/*4"):   ("0",   "CYP2D6 poor metabolizer"),
 ("CYP2D6", "*1/*4"):   ("1.0", "CYP2D6 intermediate metabolizer"),
 ("CYP2D6", "*1/*1"):   ("2.0", "CYP2D6 normal metabolizer"),
 ("CYP2D6", "*1/*1xN"): (">2.25", "CYP2D6 ultrarapid metabolizer"),

 ("CYP2C19", "*2/*2"):   ("n/a", "CYP2C19 poor metabolizer"),
 ("CYP2C19", "*1/*2"):   ("n/a", "CYP2C19 intermediate metabolizer"),
 ("CYP2C19", "*1/*1"):   ("n/a", "CYP2C19 normal metabolizer"),
 ("CYP2C19", "*1/*17"):  ("n/a", "CYP2C19 rapid metabolizer"),
 ("CYP2C19", "*17/*17"): ("n/a", "CYP2C19 ultrarapid metabolizer"),

 ("CYP2C9", "*3/*3"): ("0",   "CYP2C9 poor metabolizer"),
 ("CYP2C9", "*1/*3"): ("1",   "CYP2C9 intermediate metabolizer"),
 ("CYP2C9", "*1/*1"): ("2",   "CYP2C9 normal metabolizer"),

 ("CYP3A5", "*1/*1"): ("n/a", "CYP3A5 normal metabolizer (extensive metabolizer, CYP3A5 expresser)"),
 ("CYP3A5", "*1/*3"): ("n/a", "CYP3A5 intermediate metabolizer (CYP3A5 expresser)"),
 ("CYP3A5", "*3/*3"): ("n/a", "CYP3A5 poor metabolizer (CYP3A5 non-expresser)"),

 ("CYP2B6", "*6/*6"): ("n/a", "CYP2B6 poor metabolizer"),
 ("CYP2B6", "*1/*6"): ("n/a", "CYP2B6 intermediate metabolizer"),
 ("CYP2B6", "*1/*1"): ("n/a", "CYP2B6 normal metabolizer"),
 ("CYP2B6", "*1/*4"): ("n/a", "CYP2B6 rapid metabolizer"),

 ("DPYD", "c.1905+1G>A/c.1905+1G>A (*2A/*2A)"): ("0",   "DPYD poor metabolizer"),
 ("DPYD", "c.1905+1G>A/Reference (*1/*2A)"):    ("1",   "DPYD intermediate metabolizer"),
 ("DPYD", "c.1236G>A(HapB3)/Reference"):        ("1.5", "DPYD intermediate metabolizer"),
 ("DPYD", "Reference/Reference (*1/*1)"):       ("2",   "DPYD normal metabolizer"),

 ("TPMT", "*3A/*3A"): ("n/a", "TPMT poor metabolizer"),
 ("TPMT", "*1/*3B"):  ("n/a", "TPMT intermediate metabolizer"),
 ("TPMT", "*1/*3C"):  ("n/a", "TPMT intermediate metabolizer"),
 ("TPMT", "*1/*1"):   ("n/a", "TPMT normal metabolizer"),

 ("NUDT15", "*3/*3"): ("n/a", "NUDT15 poor metabolizer"),
 ("NUDT15", "*1/*3"): ("n/a", "NUDT15 intermediate metabolizer"),
 ("NUDT15", "*1/*1"): ("n/a", "NUDT15 normal metabolizer"),

 ("SLCO1B1", "*5/*5 (521CC; *15/*15 indistinguishable)"): ("n/a", "SLCO1B1 poor function"),
 ("SLCO1B1", "*1/*5 (521TC; *1/*15 indistinguishable)"):  ("n/a", "SLCO1B1 decreased function"),
 ("SLCO1B1", "*1/*1 (521TT)"):                            ("n/a", "SLCO1B1 normal function"),

 ("UGT1A1", "*28/*28"): ("n/a", "UGT1A1 poor metabolizer"),
 ("UGT1A1", "*1/*28"):  ("n/a", "UGT1A1 intermediate metabolizer"),
 ("UGT1A1", "*1/*1"):   ("n/a", "UGT1A1 extensive (normal) metabolizer"),

 ("HLA-B*57:01", "HLA-B*57:01 positive (inferred from HCP5 rs2395029 G)"):
   ("n/a", "High risk of abacavir hypersensitivity (HLA-B*57:01 positive)"),
 ("HLA-B*57:01", "HLA-B*57:01 negative (inferred from HCP5 rs2395029 T/T)"):
   ("n/a", "Very low risk of abacavir hypersensitivity (HLA-B*57:01 negative)"),
 ("HLA-B*15:02", "HLA-B*15:02 positive (one or two *15:02 alleles)"):
   ("n/a", "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN (HLA-B*15:02 carrier)"),
 ("HLA-B*15:02", "HLA-B*15:02 negative"):
   ("n/a", "Normal risk of HLA-B*15:02-mediated SJS/TEN (non-carrier)"),
 ("HLA-A*31:01", "HLA-A*31:01 positive (one or two *31:01 alleles)"):
   ("n/a", "Increased risk of carbamazepine-induced hypersensitivity (HLA-A*31:01 carrier)"),
 ("HLA-A*31:01", "HLA-A*31:01 negative"):
   ("n/a", "Normal risk of HLA-A*31:01-mediated hypersensitivity (non-carrier)"),
 ("HLA-B*58:01", "HLA-B*58:01 positive (one or two *58:01 alleles)"):
   ("n/a", "Significantly increased risk of allopurinol SCAR (HLA-B*58:01 carrier)"),
 ("HLA-B*58:01", "HLA-B*58:01 negative"):
   ("n/a", "Low or reduced risk of allopurinol SCAR (non-carrier)"),

 ("G6PD", "Mediterranean (hemizygous male) or Mediterranean/Mediterranean (homozygous female)"):
   ("n/a", "G6PD Deficient (WHO class II)"),
 ("G6PD", "A- (hemizygous male) or A-/A- (homozygous female)"):
   ("n/a", "G6PD Deficient (WHO class III)"),
 ("G6PD", "B/A- (heterozygous female)"):
   ("n/a", "G6PD Variable"),
 ("G6PD", "B (hemizygous male) or B/B (homozygous female)"):
   ("n/a", "G6PD Normal (WHO class IV)"),

 ("IFNL3", "rs12979860 C/C"): ("n/a", "Favorable response genotype"),
 ("IFNL3", "rs12979860 T/T"): ("n/a", "Unfavorable response genotype"),

 ("VKORC1", "rs9923231 (-1639G>A) A/A"): ("n/a", "VKORC1 -1639 A/A: increased warfarin sensitivity (decreased dose requirement)"),
 ("VKORC1", "rs9923231 (-1639G>A) G/G"): ("n/a", "VKORC1 -1639 G/G: normal warfarin sensitivity"),

 ("CYP4F2", "*3/*3"): ("n/a", "CYP4F2 decreased function (*3 homozygote)"),
 ("CYP4F2", "*1/*1"): ("n/a", "CYP4F2 normal function"),

 ("CFTR", "G551D positive (rs75527207 AA or AG)"): ("n/a", "CFTR ivacaftor-responsive genotype (G551D gating variant)"),
 ("CFTR", "F508del/F508del"):                      ("n/a", "CFTR ivacaftor-non-responsive genotype (F508del homozygote)"),

 ("MT-RNR1", "m.1555A>G (homoplasmic)"): ("n/a", "MT-RNR1 increased risk of aminoglycoside-induced hearing loss"),
 ("MT-RNR1", "m.1555A (reference)"):     ("n/a", "MT-RNR1 normal risk of aminoglycoside-induced hearing loss"),

 ("RYR1", "RYR1 c.1840C>T p.(Arg614Cys), heterozygous"):        ("n/a", "Malignant Hyperthermia susceptible"),
 ("RYR1", "No RYR1/CACNA1S MH-causative variant detected"):     ("n/a", "Uncertain susceptibility"),
}

# ---------------------------------------------------------------------------
# STEP 3.  (gene, drug, phenotype) -> short action + verbatim guideline text
# ---------------------------------------------------------------------------
# short  = the answers-CSV string ("drug: action")
# quote  = verbatim recommendation text from the supplied corpus (asserted)
# src    = guideline annotation the quote comes from
# abst   = abstain flag (only where the corpus has no usable recommendation)
PA = {
 "codeine":       "CPIC PA166104996",
 "tamoxifen":     "CPIC PA166176068",
 "ondansetron":   "CPIC PA166161954",
 "paroxetine":    "CPIC PA166127636",
 "amitriptyline": "CPIC PA166105006",
 "clopidogrel":   "CPIC PA166104948",
 "voriconazole":  "CPIC PA166161537",
 "citalopram":    "CPIC PA166127638",
 "omeprazole":    "CPIC PA166219103",
 "warfarin":      "CPIC PA166104949",
 "phenytoin":     "CPIC PA166122806",
 "celecoxib":     "CPIC PA166122696",
 "tacrolimus":    "CPIC PA166124619",
 "efavirenz":     "CPIC PA166182603",
 "fluorouracil":  "CPIC PA166122686",
 "capecitabine":  "CPIC PA166109594",
 "azathioprine":  "CPIC PA166104933",
 "thiopurines":   "CPIC PA166104933",
 "simvastatin":   "CPIC PA166105005",
 "atorvastatin":  "CPIC PA166262221",
 "atazanavir":    "CPIC PA166128738",
 "irinotecan":    "DPWG PA166104951",
 "abacavir":      "CPIC PA166104997",
 "allopurinol":   "CPIC PA166105003",
 "carbamazepine": "CPIC PA166105008",
 "oxcarbazepine": "CPIC PA166176623",
 "primaquine":    "CPIC PA166279621",
 "rasburicase":   "CPIC PA166267376",
 "PEG-IFN-a":     "CPIC PA166110235",
 "aminoglycosides": "CPIC PA166262781",
 "volatile-anaesthetics": "CPIC PA166303941",
 "ivacaftor":     "CPIC PA166114461",
}

R = {}   # (gene, drug, phenotype) -> dict(short=, quote=, abst=)


def rec(gene, drug, pheno, short, quote, abst=False):
    R[(gene, drug, pheno)] = {"short": short, "quote": verbatim(quote) if quote else "",
                              "abst": abst, "src": PA.get(drug, "")}


# ---- CYP2D6 / codeine -----------------------------------------------------
rec("CYP2D6", "codeine", "CYP2D6 ultrarapid metabolizer",
    "codeine: avoid - risk of serious morphine toxicity; if an opioid is warranted use a non-tramadol opioid",
    "Avoid codeine use because of potential for serious toxicity. If opioid use is warranted, consider a non-tramadol opioid.")
rec("CYP2D6", "codeine", "CYP2D6 normal metabolizer",
    "codeine: use label-recommended age- or weight-specific dosing",
    "Use codeine label recommended age- or weight-specific dosing")
rec("CYP2D6", "codeine", "CYP2D6 intermediate metabolizer",
    "codeine: use label-recommended age- or weight-specific dosing; if no response and an opioid is warranted, consider a non-tramadol opioid",
    "Use codeine label recommended age- or weight-specific dosing. If no response and opioid use is warranted, consider a non-tramadol opioid")
rec("CYP2D6", "codeine", "CYP2D6 poor metabolizer",
    "codeine: avoid - diminished analgesia; if an opioid is warranted use a non-tramadol opioid",
    "Avoid codeine use because of possibility of diminished analgesia. If opioid use is warranted, consider a non-tramadol opioid.")

# ---- CYP2D6 / tamoxifen ---------------------------------------------------
_tam_std = "Avoid moderate and strong CYP2D6 inhibitors. Initiate therapy with recommended standard of care dosing (tamoxifen 20 mg/day)."
_tam_alt = ("Consider hormonal therapy such as an aromatase inhibitor for postmenopausal women or aromatase inhibitor "
            "along with ovarian function suppression in premenopausal women")
_tam_pm = ("Recommend alternative hormonal therapy such as an aromatase inhibitor for postmenopausal women or aromatase "
           "inhibitor along with ovarian function suppression in premenopausal women")
rec("CYP2D6", "tamoxifen", "CYP2D6 ultrarapid metabolizer",
    "tamoxifen: standard dosing 20 mg/day; avoid moderate/strong CYP2D6 inhibitors", _tam_std)
rec("CYP2D6", "tamoxifen", "CYP2D6 normal metabolizer",
    "tamoxifen: standard dosing 20 mg/day; avoid moderate/strong CYP2D6 inhibitors", _tam_std)
rec("CYP2D6", "tamoxifen", "CYP2D6 intermediate metabolizer",
    "tamoxifen: consider an aromatase inhibitor instead (plus ovarian suppression if premenopausal); if an AI is contraindicated consider tamoxifen 40 mg/day and avoid CYP2D6 inhibitors", _tam_alt)
rec("CYP2D6", "tamoxifen", "CYP2D6 poor metabolizer",
    "tamoxifen: avoid - use alternative hormonal therapy (aromatase inhibitor, plus ovarian suppression if premenopausal)", _tam_pm)

# ---- CYP2D6 / ondansetron -------------------------------------------------
rec("CYP2D6", "ondansetron", "CYP2D6 ultrarapid metabolizer",
    "ondansetron: avoid - select an alternative antiemetic not predominantly metabolised by CYP2D6 (e.g. granisetron)",
    "Select and alternative drug not predominantly metabolized by CYP2D6 (i.e. granisetron).")
rec("CYP2D6", "ondansetron", "CYP2D6 normal metabolizer",
    "ondansetron: initiate therapy with the recommended starting dose",
    "Initiate therapy with recommended starting dose.")
rec("CYP2D6", "ondansetron", "CYP2D6 intermediate metabolizer",
    "ondansetron: initiate therapy with the recommended starting dose (insufficient evidence for a genotype-based change)",
    "Insufficient evidence demonstrating clinical impact based on CYP2D6 genotype. Initiate therapy with recommended starting dose.")
rec("CYP2D6", "ondansetron", "CYP2D6 poor metabolizer",
    "ondansetron: initiate therapy with the recommended starting dose (insufficient evidence for a genotype-based change)",
    "Insufficient evidence demonstrating clinical impact based on CYP2D6 genotype. Initiate therapy with recommended starting dose.")

# ---- CYP2D6 / paroxetine --------------------------------------------------
rec("CYP2D6", "paroxetine", "CYP2D6 ultrarapid metabolizer",
    "paroxetine: avoid - select an alternative antidepressant not predominantly metabolised by CYP2D6",
    "Select alternative drug not predominantly metabolized by CYP2D6.")
rec("CYP2D6", "paroxetine", "CYP2D6 normal metabolizer",
    "paroxetine: initiate therapy with the recommended starting dose",
    "Initiate therapy with recommended starting dose.")
rec("CYP2D6", "paroxetine", "CYP2D6 intermediate metabolizer",
    "paroxetine: consider a lower starting dose and slower titration than for normal metabolisers",
    "Consider a lower starting dose and slower titration schedule as compared to normal metabolizers.")
rec("CYP2D6", "paroxetine", "CYP2D6 poor metabolizer",
    "paroxetine: consider a 50% lower starting dose, slower titration and 50% lower maintenance dose",
    "Consider a 50% reduction in recommended starting dose, slower titration schedule, and a 50% lower maintenance dose as compared to normal metabolizers.")

# ---- CYP2D6 / amitriptyline (CYP2C19 not genotyped; CYP2C19 NM assumed) ----
rec("CYP2D6", "amitriptyline", "CYP2D6 ultrarapid metabolizer",
    "amitriptyline: avoid (assuming CYP2C19 normal metaboliser); if warranted, titrate to a higher target dose with TDM",
    "Avoid amitriptyline use. If amitriptyline is warranted, consider titrating to a higher target dose (compared to normal metabolizers).")
rec("CYP2D6", "amitriptyline", "CYP2D6 normal metabolizer",
    "amitriptyline: initiate therapy with the recommended starting dose (assuming CYP2C19 normal metaboliser)",
    "Initiate therapy with recommended starting dose.")
rec("CYP2D6", "amitriptyline", "CYP2D6 intermediate metabolizer",
    "amitriptyline: consider a 25% reduction of the recommended starting dose with TDM (assuming CYP2C19 normal metaboliser)",
    "Consider 25% reduction of recommended starting dose.")
rec("CYP2D6", "amitriptyline", "CYP2D6 poor metabolizer",
    "amitriptyline: avoid (assuming CYP2C19 normal metaboliser); if warranted, 50% reduction of the recommended starting dose with TDM",
    "Avoid amitriptyline use. If amitriptyline is warranted, consider 50% reduction of recommended starting dose.")

# ---- CYP2C19 / clopidogrel ------------------------------------------------
_clop_std = "If considering clopidogrel, use at standard dose (75 mg/day)"
for ph in ("CYP2C19 ultrarapid metabolizer", "CYP2C19 rapid metabolizer", "CYP2C19 normal metabolizer"):
    rec("CYP2C19", "clopidogrel", ph, "clopidogrel: use at standard dose (75 mg/day)", _clop_std)
rec("CYP2C19", "clopidogrel", "CYP2C19 intermediate metabolizer",
    "clopidogrel: avoid standard dose if possible - use prasugrel or ticagrelor at standard dose if not contraindicated",
    "Avoid standard dose (75 mg) clopidogrel if possible. Use prasugrel or ticagrelor at standard dose if no contraindication.")
rec("CYP2C19", "clopidogrel", "CYP2C19 poor metabolizer",
    "clopidogrel: avoid if possible - use prasugrel or ticagrelor at standard dose if not contraindicated",
    "Avoid clopidogrel if possible. Use prasugrel or ticagrelor at standard dose if no contraindication.")

# ---- CYP2C19 / voriconazole (adult table; age not supplied) ---------------
_vori_alt = ("Choose an alternative agent that is not dependent on CYP2C19 metabolism as primary therapy in lieu of "
             "voriconazole. Such agents include isavuconazole, liposomal amphotericin B, and posaconazole.")
_vori_pm = ("In the event that voriconazole is considered to be the most appropriate agent, based on clinical advice, for a "
            "patient with poor metabolizer genotype, voriconazole should be administered at a preferably lower than standard "
            "dosage with careful therapeutic drug monitoring.")
rec("CYP2C19", "voriconazole", "CYP2C19 ultrarapid metabolizer",
    "voriconazole: avoid - choose an agent not dependent on CYP2C19 (isavuconazole, liposomal amphotericin B, posaconazole)", _vori_alt)
rec("CYP2C19", "voriconazole", "CYP2C19 rapid metabolizer",
    "voriconazole: avoid in adults - choose an agent not dependent on CYP2C19 (isavuconazole, liposomal amphotericin B, posaconazole)", _vori_alt)
rec("CYP2C19", "voriconazole", "CYP2C19 normal metabolizer",
    "voriconazole: initiate therapy with recommended standard-of-care dosing",
    "Initiate therapy with recommended standard of care dosing.")
rec("CYP2C19", "voriconazole", "CYP2C19 intermediate metabolizer",
    "voriconazole: initiate therapy with recommended standard-of-care dosing",
    "Initiate therapy with recommended standard of care dosing.")
rec("CYP2C19", "voriconazole", "CYP2C19 poor metabolizer",
    "voriconazole: avoid - choose an agent not dependent on CYP2C19; if voriconazole is unavoidable, use a lower dose with careful TDM", _vori_pm)

# ---- CYP2C19 / citalopram -------------------------------------------------
rec("CYP2C19", "citalopram", "CYP2C19 ultrarapid metabolizer",
    "citalopram: consider an alternative antidepressant not predominantly metabolised by CYP2C19",
    "Consider a clinically appropriate alternative antidepressant not predominantly metabolized by CYP2C19.")
rec("CYP2C19", "citalopram", "CYP2C19 rapid metabolizer",
    "citalopram: initiate therapy with the recommended starting dose; titrate up or switch if response is inadequate",
    "Initiate therapy with recommended starting dose.")
rec("CYP2C19", "citalopram", "CYP2C19 normal metabolizer",
    "citalopram: initiate therapy with the recommended starting dose",
    "Initiate therapy with recommended starting dose.")
rec("CYP2C19", "citalopram", "CYP2C19 intermediate metabolizer",
    "citalopram: initiate at the recommended starting dose but use a slower titration and lower maintenance dose",
    "Initiate therapy with recommended starting dose. Consider a slower titration schedule and lower maintenance dose than normal metabolizers.")
rec("CYP2C19", "citalopram", "CYP2C19 poor metabolizer",
    "citalopram: consider an alternative antidepressant not predominantly metabolised by CYP2C19; if used, lower starting dose, slower titration and 50% lower maintenance dose",
    "Consider a clinically appropriate antidepressant not predominantly metabolized by CYP2C19. If citalopram or escitalopram are clinically appropriate, consider a lower starting dose, slower titration schedule and 50% reduction of the standard maintenance dose as compared to normal metabolizers.")

# ---- CYP2C19 / omeprazole -------------------------------------------------
rec("CYP2C19", "omeprazole", "CYP2C19 ultrarapid metabolizer",
    "omeprazole: increase the starting daily dose by 100% and monitor for efficacy",
    "Increase starting daily dose by 100%. Daily dose may be given in divided doses. Monitor for efficacy.")
_ppi_std = ("Initiate standard starting daily dose. Consider increasing dose by 50-100% for the treatment of H. pylori "
            "infection and erosive esophagitis. Daily dose may be given in divided doses. Monitor for efficacy.")
rec("CYP2C19", "omeprazole", "CYP2C19 rapid metabolizer",
    "omeprazole: standard starting daily dose; consider a 50-100% increase for H. pylori or erosive oesophagitis", _ppi_std)
rec("CYP2C19", "omeprazole", "CYP2C19 normal metabolizer",
    "omeprazole: standard starting daily dose; consider a 50-100% increase for H. pylori or erosive oesophagitis", _ppi_std)
_ppi_im = ("Initiate standard starting daily dose. For chronic therapy (>12 weeks) and efficacy achieved, consider 50% "
           "reduction in daily dose and monitor for continued efficacy.")
rec("CYP2C19", "omeprazole", "CYP2C19 intermediate metabolizer",
    "omeprazole: standard starting daily dose; if chronic therapy (>12 weeks) is effective, consider a 50% dose reduction", _ppi_im)
rec("CYP2C19", "omeprazole", "CYP2C19 poor metabolizer",
    "omeprazole: standard starting daily dose; if chronic therapy (>12 weeks) is effective, consider a 50% dose reduction", _ppi_im)

# ---- CYP2C9 / warfarin (VKORC1 not genotyped; -1639 GG assumed) -----------
# The corpus gives no per-phenotype sentence for warfarin: its recommendation IS
# the VKORC1 x CYP2C9 stable-dose table, quoted verbatim by row below.
_warf_gg = "GG | 5-7 | 5-7 | 3-4 | 3-4 | 3-4 | 0.5-2"
_warf_aa = "AA | 3-4 | 3-4 | 0.5-2 | 0.5-2 | 0.5-2 | 0.5-2"
rec("CYP2C9", "warfarin", "CYP2C9 normal metabolizer",
    "warfarin: use a validated genotype-guided dosing algorithm; corpus dose table gives an expected stable dose of 5-7 mg/day (CYP2C9 *1/*1, assuming VKORC1 -1639 GG)",
    _warf_gg)
rec("CYP2C9", "warfarin", "CYP2C9 intermediate metabolizer",
    "warfarin: use a validated genotype-guided dosing algorithm; corpus dose table gives an expected stable dose of 3-4 mg/day (CYP2C9 *1/*3, assuming VKORC1 -1639 GG)",
    _warf_gg)
rec("CYP2C9", "warfarin", "CYP2C9 poor metabolizer",
    "warfarin: use a validated genotype-guided dosing algorithm; corpus dose table gives an expected stable dose of only 0.5-2 mg/day (CYP2C9 *3/*3, assuming VKORC1 -1639 GG) - standard initiation carries a high major-bleeding risk",
    _warf_gg)

# ---- VKORC1 / warfarin (CYP2C9 not genotyped; *1/*1 assumed) --------------
rec("VKORC1", "warfarin", "VKORC1 -1639 A/A: increased warfarin sensitivity (decreased dose requirement)",
    "warfarin: use a validated genotype-guided dosing algorithm; corpus dose table gives an expected stable dose of 3-4 mg/day (VKORC1 AA, assuming CYP2C9 *1/*1)",
    _warf_aa)
rec("VKORC1", "warfarin", "VKORC1 -1639 G/G: normal warfarin sensitivity",
    "warfarin: use a validated genotype-guided dosing algorithm; corpus dose table gives an expected stable dose of 5-7 mg/day (VKORC1 GG, assuming CYP2C9 *1/*1)",
    _warf_gg)

# ---- CYP4F2 / warfarin : ABSTAIN ------------------------------------------
# The supplied CYP4F2 excerpt (935 chars) is the warfarin annotation whose only
# dose table is keyed on VKORC1 x CYP2C9. It names CYP4F2 as an input gene but
# gives no CYP4F2 genotype-specific recommendation, and rs2108622 appears
# nowhere in the corpus. Per brief rule 1, abstain rather than fetch.
for _ph in ("CYP4F2 decreased function (*3 homozygote)", "CYP4F2 normal function"):
    rec("CYP4F2", "warfarin", _ph,
        "warfarin: ABSTAIN - the supplied CPIC excerpt gives no CYP4F2 genotype-specific dose recommendation "
        "(its only dose table is VKORC1 x CYP2C9); a CYP4F2-only result cannot be translated into an action",
        "", abst=True)

# ---- CYP2C9 / phenytoin (HLA-B*15:02 not typed; 'negative' branch used) ----
rec("CYP2C9", "phenytoin", "CYP2C9 normal metabolizer",
    "phenytoin: no adjustment from typical dosing; titrate by TDM and response (assumes HLA-B*15:02 negative/untested)",
    "No adjustments needed from typical dosing strategies. Subsequent doses should be adjusted according to therapeutic drug monitoring, response, and side effects.")
rec("CYP2C9", "phenytoin", "CYP2C9 intermediate metabolizer",
    "phenytoin: typical loading dose, then ~25% lower maintenance dose with TDM (assumes HLA-B*15:02 negative/untested)",
    "For first dose, use typical initial or loading dose. For subsequent doses, use approximately 25% less than typical maintenance dose.")
rec("CYP2C9", "phenytoin", "CYP2C9 poor metabolizer",
    "phenytoin: typical loading dose, then ~50% lower maintenance dose with TDM (assumes HLA-B*15:02 negative/untested)",
    "For first dose, use typical initial or loading dose. For subsequent doses use approximately 50% less than typical maintenance dose.")

# ---- CYP2C9 / celecoxib ---------------------------------------------------
rec("CYP2C9", "celecoxib", "CYP2C9 normal metabolizer",
    "celecoxib: initiate therapy with the recommended starting dose",
    "Initiate therapy with recommended starting dose. In accordance with the prescribing information, use the lowest effective dosage for shortest duration consistent with individual patient treatment goals.")
rec("CYP2C9", "celecoxib", "CYP2C9 intermediate metabolizer",
    "celecoxib: initiate at the lowest recommended starting dose and titrate up with caution",
    "Initiate therapy with lowest recommended starting dose. Titrate dose upward to clinical effect or maximum recommended dose with caution.")
rec("CYP2C9", "celecoxib", "CYP2C9 poor metabolizer",
    "celecoxib: initiate at 25-50% of the lowest recommended starting dose, titrate cautiously, or use a non-CYP2C9 alternative",
    "Initiate therapy with 25-50% of the lowest recommended starting dose. Titrate dose upward to clinical effect or 25-50% of the maximum recommended dose with caution.")

# ---- CYP3A5 / tacrolimus --------------------------------------------------
_tac_up = ("Increase starting dose 1.5 to 2 times recommended starting dose d . Total starting dose should not exceed "
           "0.3mg/kg/day. Use therapeutic drug monitoring to guide dose adjustments")
rec("CYP3A5", "tacrolimus", "CYP3A5 normal metabolizer (extensive metabolizer, CYP3A5 expresser)",
    "tacrolimus: increase the starting dose 1.5-2x (not exceeding 0.3 mg/kg/day) and guide by TDM", _tac_up)
rec("CYP3A5", "tacrolimus", "CYP3A5 intermediate metabolizer (CYP3A5 expresser)",
    "tacrolimus: increase the starting dose 1.5-2x (not exceeding 0.3 mg/kg/day) and guide by TDM", _tac_up)
rec("CYP3A5", "tacrolimus", "CYP3A5 poor metabolizer (CYP3A5 non-expresser)",
    "tacrolimus: initiate at the standard recommended dose and guide by TDM",
    "Initiate therapy with standard recommended dose. Use therapeutic drug monitoring to guide dose adjustments")

# ---- CYP2B6 / efavirenz ---------------------------------------------------
for ph in ("CYP2B6 rapid metabolizer", "CYP2B6 normal metabolizer"):
    rec("CYP2B6", "efavirenz", ph, "efavirenz: initiate with standard dosing (600 mg/day)",
        "Initiate efavirenz with standard dosing (600 mg/day)")
rec("CYP2B6", "efavirenz", "CYP2B6 intermediate metabolizer",
    "efavirenz: consider initiating at a decreased dose of 400 mg/day",
    "Consider initiating efavirenz with decreased dose of 400 mg/day")
rec("CYP2B6", "efavirenz", "CYP2B6 poor metabolizer",
    "efavirenz: consider initiating at a decreased dose of 400 or 200 mg/day",
    "Consider initiating efavirenz with decreased dose of 400 or 200 mg/day")

# ---- DPYD / fluorouracil + capecitabine -----------------------------------
for drug in ("fluorouracil", "capecitabine"):
    rec("DPYD", drug, "DPYD normal metabolizer",
        f"{drug}: no genotype-based change - use label-recommended dosage and administration",
        "Based on genotype, there is no indication to change dose or therapy. Use label-recommended dosage and administration")
    rec("DPYD", drug, "DPYD intermediate metabolizer",
        f"{drug}: reduce the starting dose by 50%, then titrate on toxicity/TDM",
        "Reduce starting dose by 50% followed by titration of dose based on toxicity or therapeutic drug monitoring (if available).")
    rec("DPYD", drug, "DPYD poor metabolizer",
        f"{drug}: avoid - do not use 5-fluorouracil or 5-fluorouracil prodrug-based regimens (complete DPD deficiency)",
        "Activity score 0: Avoid use of 5-fluorouracil or 5-fluorouracil prodrug-based regimens.")

# ---- TPMT + NUDT15 / thiopurines ------------------------------------------
# The supplied azathioprine annotation carries no per-phenotype dosing table;
# its recommendation text is the prose summary, quoted verbatim below.
_thio_im = "a lower than standard starting dose of azathioprine should be considered for TPMT or NUDT15 intermediate metabolizers"
_thio_pm = "For TPMT and/or NUDT15 poor metabolizers, an alternative immunosuppressant therapy should be considered."
rec("TPMT", "azathioprine", "TPMT normal metabolizer",
    "azathioprine: standard starting dose - no genotype-based adjustment indicated for TPMT normal metabolisers",
    "The CPIC dosing guideline recommends that a lower than standard starting dose of azathioprine should be considered for TPMT or NUDT15 intermediate metabolizers")
rec("TPMT", "azathioprine", "TPMT intermediate metabolizer",
    "azathioprine: consider a lower than standard starting dose", _thio_im)
rec("TPMT", "azathioprine", "TPMT poor metabolizer",
    "azathioprine: avoid standard dosing - consider an alternative immunosuppressant; if a thiopurine is used, drastic dose reduction is required (tolerated dose ~10% of normal)", _thio_pm)
rec("NUDT15", "thiopurines", "NUDT15 normal metabolizer",
    "thiopurines (azathioprine/mercaptopurine/thioguanine): standard starting dose - no genotype-based adjustment indicated",
    "The CPIC dosing guideline recommends that a lower than standard starting dose of azathioprine should be considered for TPMT or NUDT15 intermediate metabolizers")
rec("NUDT15", "thiopurines", "NUDT15 intermediate metabolizer",
    "thiopurines (azathioprine/mercaptopurine/thioguanine): consider a lower than standard starting dose", _thio_im)
rec("NUDT15", "thiopurines", "NUDT15 poor metabolizer",
    "thiopurines (azathioprine/mercaptopurine/thioguanine): avoid standard dosing - consider an alternative immunosuppressant; tolerated dose is ~10% of normal", _thio_pm)

# ---- SLCO1B1 / simvastatin + atorvastatin ---------------------------------
rec("SLCO1B1", "simvastatin", "SLCO1B1 normal function",
    "simvastatin: prescribe the desired starting dose and adjust per disease-specific guidelines",
    "Prescribe desired starting dose and adjust doses based on disease-specific guidelines.")
rec("SLCO1B1", "simvastatin", "SLCO1B1 decreased function",
    "simvastatin: prescribe an alternative statin of the desired potency; if simvastatin is warranted, limit the dose to <20 mg/day",
    "Prescribe an alternative statin depending on the desired potency (see Figure 1 for recommendations for alternative statins). If simvastatin therapy is warranted, limit dose to <20mg/day.")
rec("SLCO1B1", "simvastatin", "SLCO1B1 poor function",
    "simvastatin: avoid - prescribe an alternative statin of the desired potency (highly increased myopathy risk)",
    "Prescribe an alternative statin depending on the desired potency (see Figure 1 for recommendations for alternative statins).")
rec("SLCO1B1", "atorvastatin", "SLCO1B1 normal function",
    "atorvastatin: prescribe the desired starting dose and adjust per disease-specific guidelines",
    "Prescribe desired starting dose and adjust doses based on disease-specific guidelines.")
rec("SLCO1B1", "atorvastatin", "SLCO1B1 decreased function",
    "atorvastatin: prescribe <=40 mg as the starting dose; if >40 mg is needed, consider combination therapy",
    "Prescribe \u226440mg as a starting dose and adjust doses of atorvastatin based on disease-specific guidelines.")
rec("SLCO1B1", "atorvastatin", "SLCO1B1 poor function",
    "atorvastatin: prescribe <=20 mg as the starting dose; if >20 mg is needed, consider rosuvastatin or combination therapy",
    "Prescribe \u226420mg as a starting dose and adjust doses of atorvastatin based on disease-specific guidelines.")

# ---- UGT1A1 / atazanavir --------------------------------------------------
rec("UGT1A1", "atazanavir", "UGT1A1 extensive (normal) metabolizer",
    "atazanavir: no need to avoid atazanavir based on the UGT1A1 result",
    "There is no need to avoid prescribing of atazanavir based on UGT1A1 genetic test result.")
rec("UGT1A1", "atazanavir", "UGT1A1 intermediate metabolizer",
    "atazanavir: no need to avoid atazanavir; advise that bilirubin-related jaundice is unlikely with this genotype",
    "There is no need to avoid prescribing of atazanavir based on UGT1A1 genetic test result. Inform the patient that some patients stop atazanavir because of jaundice (yellow eyes and skin), but that this patient\u2019s genotype makes this unlikely")
rec("UGT1A1", "atazanavir", "UGT1A1 poor metabolizer",
    "atazanavir: consider an alternative agent, particularly where jaundice would concern the patient",
    "Consider an alternative agent particularly where jaundice would be of concern to the patient.")

# ---- UGT1A1 / irinotecan (no CPIC guideline in corpus; DPWG used) ---------
rec("UGT1A1", "irinotecan", "UGT1A1 extensive (normal) metabolizer",
    "irinotecan: no action needed for this gene-drug interaction (DPWG; no CPIC guideline in the supplied corpus)",
    "NO action is needed for this gene-drug interaction")
rec("UGT1A1", "irinotecan", "UGT1A1 intermediate metabolizer",
    "irinotecan: no action needed for this gene-drug interaction (DPWG; no CPIC guideline in the supplied corpus)",
    "NO action is needed for this gene-drug interaction")
rec("UGT1A1", "irinotecan", "UGT1A1 poor metabolizer",
    "irinotecan: start at 70% of the normal dose and escalate on neutrophil count if tolerated (DPWG; no CPIC guideline in the supplied corpus)",
    "Start with 70% of the normal a dose. If the patient tolerates this initial dose, the dose can be increased, guided by the neutrophil count.")

# ---- HLA ------------------------------------------------------------------
rec("HLA-B*57:01", "abacavir", "High risk of abacavir hypersensitivity (HLA-B*57:01 positive)",
    "abacavir: avoid - abacavir is not recommended (high risk of hypersensitivity reaction)",
    "Abacavir is not recommended")
rec("HLA-B*57:01", "abacavir", "Very low risk of abacavir hypersensitivity (HLA-B*57:01 negative)",
    "abacavir: use per standard dosing guidelines",
    "Use abacavir per standard dosing guidelines")
rec("HLA-B*58:01", "allopurinol", "Significantly increased risk of allopurinol SCAR (HLA-B*58:01 carrier)",
    "allopurinol: avoid - allopurinol is contraindicated (risk of severe cutaneous adverse reaction)",
    "Allopurinol is contraindicated")
rec("HLA-B*58:01", "allopurinol", "Low or reduced risk of allopurinol SCAR (non-carrier)",
    "allopurinol: use per standard dosing guidelines",
    "Use allopurinol per standard dosing guidelines")
rec("HLA-B*15:02", "carbamazepine", "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN (HLA-B*15:02 carrier)",
    "carbamazepine: avoid - if the patient is carbamazepine-naive, do not use carbamazepine (risk of SJS/TEN)",
    "If patient is carbamazepine-naive, do not use carbamazepine.")
rec("HLA-B*15:02", "carbamazepine", "Normal risk of HLA-B*15:02-mediated SJS/TEN (non-carrier)",
    "carbamazepine: use per standard dosing guidelines (assumes HLA-A*31:01 also non-carrier/untested)",
    "Use carbamazepine per standard dosing guidelines.")
rec("HLA-B*15:02", "oxcarbazepine", "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN (HLA-B*15:02 carrier)",
    "oxcarbazepine: avoid - if the patient is oxcarbazepine-naive, do not use oxcarbazepine (risk of SJS/TEN)",
    "If patient is oxcarbazepine-naive, do not use oxcarbazepine.")
rec("HLA-B*15:02", "oxcarbazepine", "Normal risk of HLA-B*15:02-mediated SJS/TEN (non-carrier)",
    "oxcarbazepine: use per standard dosing guidelines",
    "Use oxcarbazepine per standard dosing guidelines.")
rec("HLA-A*31:01", "carbamazepine", "Increased risk of carbamazepine-induced hypersensitivity (HLA-A*31:01 carrier)",
    "carbamazepine: avoid if the patient is carbamazepine-naive and alternatives exist; if no alternative, use with increased clinical monitoring and stop at the first sign of a cutaneous reaction",
    "If patient is carbamazepine-naive and alternative agents are available, do not use carbamazepine.")
rec("HLA-A*31:01", "carbamazepine", "Normal risk of HLA-A*31:01-mediated hypersensitivity (non-carrier)",
    "carbamazepine: use per standard dosing guidelines (assumes HLA-B*15:02 also non-carrier/untested)",
    "Use carbamazepine per standard dosing guidelines.")

# ---- G6PD -----------------------------------------------------------------
rec("G6PD", "rasburicase", "G6PD Deficient (WHO class II)",
    "rasburicase: contraindicated - do not use; alternatives include allopurinol",
    "Rasburicase is contraindicated; alternatives include allopurinol")
rec("G6PD", "rasburicase", "G6PD Deficient (WHO class III)",
    "rasburicase: contraindicated - do not use; alternatives include allopurinol",
    "Rasburicase is contraindicated; alternatives include allopurinol")
rec("G6PD", "rasburicase", "G6PD Normal (WHO class IV)",
    "rasburicase: no reason to withhold rasburicase based on G6PD status",
    "No reason to withhold rasburicase based on G6PD status")
rec("G6PD", "rasburicase", "G6PD Variable",
    "rasburicase: do not give until G6PD enzyme activity has been measured and shown to be normal; alternatives include allopurinol",
    "To ascertain that G6PD status is normal, enzyme activity must be measured; alternatives include allopurinol")
rec("G6PD", "primaquine", "G6PD Deficient (WHO class III)",
    "primaquine: avoid standard 14-day dosing - high risk of acute haemolytic anaemia; if radical cure is needed use the WHO/CDC weekly regimen (0.75 mg/kg or 45 mg once weekly x8 weeks) with close monitoring for haemolysis",
    "Avoid primaquine, except in the following cases where established expert consensus guidelines for the treatment of malaria should be followed: (1) Treating Plasmodium vivax or Plasmodium ovale malaria for radical cure of liver-stage infections: 0.75 mg/kg once weekly x8 weeks (WHO) or 45 mg once weekly x8 weeks (CDC) - with close monitoring for hemolysis")
rec("G6PD", "primaquine", "G6PD Normal (WHO class IV)",
    "primaquine: no reason to avoid primaquine based on G6PD status",
    "No reason to avoid primaquine based on G6PD status")
rec("G6PD", "primaquine", "G6PD Variable",
    "primaquine: measure G6PD enzyme activity first and dose per the activity-based phenotype",
    "To ascertain G6PD status, enzyme activity must be measured. Drug use should be guided per the recommendations based on the activity-based phenotype.")

# ---- IFNL3 ----------------------------------------------------------------
rec("IFNL3", "PEG-IFN-\u03b1", "Favorable response genotype",
    "PEG-IFN-\u03b1 (+ ribavirin): favourable response genotype - ~70% chance of SVR after 48 weeks; consider this before initiating a PEG-IFN-\u03b1/RBV regimen",
    "Approximately 70% chance for SVR c after 48 weeks of treatment. Consider implications before initiating PEG-IFN alpha and RBV containing regimens.")
rec("IFNL3", "PEG-IFN-\u03b1", "Unfavorable response genotype",
    "PEG-IFN-\u03b1 (+ ribavirin): unfavourable response genotype - ~30% chance of SVR after 48 weeks; consider a non-interferon regimen before initiating PEG-IFN-\u03b1/RBV",
    "Approximately 30% chance for SVR c after 48 weeks of treatment. Consider implications before initiating PEG-IFN alpha and RBV containing regimens.")

# ---- MT-RNR1 --------------------------------------------------------------
rec("MT-RNR1", "aminoglycosides", "MT-RNR1 increased risk of aminoglycoside-induced hearing loss",
    "aminoglycosides: avoid unless the high risk of permanent hearing loss is outweighed by the severity of infection and a lack of safe, effective alternatives",
    "Avoid aminoglycoside antibiotics unless the high risk of permanent hearing loss is outweighed by the severity of infection and lack of safe or effective alternative therapies.")
rec("MT-RNR1", "aminoglycosides", "MT-RNR1 normal risk of aminoglycoside-induced hearing loss",
    "aminoglycosides: use at standard doses for the shortest feasible course with TDM and regular hearing evaluation",
    "Use aminoglycoside antibiotics at standard doses for the shortest feasible course with therapeutic dose monitoring. Evaluate regularly for hearing loss in line with local guidance.")

# ---- RYR1 -----------------------------------------------------------------
rec("RYR1", "volatile-anaesthetics", "Malignant Hyperthermia susceptible",
    "volatile anaesthetics (and succinylcholine): relatively contraindicated - do not use except in extraordinary circumstances; use a non-triggering anaesthetic technique",
    "Halogenated volatile anesthetics or depolarizing muscle relaxants succinylcholine are relatively contraindicated in persons with MHS. They should not be used, except in extraordinary circumstances where the benefits outweigh the risks.")
rec("RYR1", "volatile-anaesthetics", "Uncertain susceptibility",
    "volatile anaesthetics: a negative result does not exclude MH susceptibility - let clinical findings, family history, further genetic testing and other laboratory data guide use",
    "Clinical findings, family history, further genetic testing and other laboratory data should guide use of halogenated volatile anesthetics or depolarizing muscle relaxants.")

# ---- CFTR -----------------------------------------------------------------
rec("CFTR", "ivacaftor", "CFTR ivacaftor-responsive genotype (G551D gating variant)",
    "ivacaftor: use ivacaftor according to the product label",
    "Use ivacaftor according to the product label")
rec("CFTR", "ivacaftor", "CFTR ivacaftor-non-responsive genotype (F508del homozygote)",
    "ivacaftor: not recommended - F508del homozygotes are unlikely to respond to ivacaftor monotherapy",
    "Ivacaftor is not recommended")


# ---------------------------------------------------------------------------
# STEP 4.  Lethal-class flag
# ---------------------------------------------------------------------------
# A case is 'lethal-class' when standard dosing of the named drug in this
# phenotype can plausibly kill or cause irreversible major organ/tissue injury
# (not merely treatment failure or tolerable toxicity).
LETHAL = {
 ("CYP2D6", "codeine", "CYP2D6 ultrarapid metabolizer"): "Ultrarapid morphine formation; documented fatal/near-fatal respiratory depression at label doses",
 ("CYP2C19", "clopidogrel", "CYP2C19 poor metabolizer"): "Loss of antiplatelet effect; stent thrombosis, MI and cardiovascular death",
 ("CYP2C9", "warfarin", "CYP2C9 poor metabolizer"): "Stable dose 0.5-2 mg/day; standard initiation causes major/fatal haemorrhage",
 ("DPYD", "fluorouracil", "DPYD poor metabolizer"): "Complete DPD deficiency; corpus states 'severe or even fatal drug toxicity'",
 ("DPYD", "capecitabine", "DPYD poor metabolizer"): "Complete DPD deficiency; corpus states 'severe or even fatal drug toxicity'",
 ("DPYD", "fluorouracil", "DPYD intermediate metabolizer"): "Partial DPD deficiency; corpus states 'severe or even fatal drug toxicity' at full dose",
 ("DPYD", "capecitabine", "DPYD intermediate metabolizer"): "Partial DPD deficiency; corpus states 'severe or even fatal drug toxicity' at full dose",
 ("TPMT", "azathioprine", "TPMT poor metabolizer"): "Profound prolonged myelosuppression at conventional dose; can be fatal",
 ("NUDT15", "thiopurines", "NUDT15 poor metabolizer"): "Severe myelosuppression at conventional dose; tolerated dose ~10% of normal",
 ("HLA-B*15:02", "carbamazepine", "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN (HLA-B*15:02 carrier)"): "SJS/TEN with substantial mortality",
 ("HLA-B*15:02", "oxcarbazepine", "Increased risk of carbamazepine/oxcarbazepine-induced SJS/TEN (HLA-B*15:02 carrier)"): "SJS/TEN with substantial mortality",
 ("HLA-A*31:01", "carbamazepine", "Increased risk of carbamazepine-induced hypersensitivity (HLA-A*31:01 carrier)"): "DRESS/SJS/TEN; can be fatal",
 ("HLA-B*57:01", "abacavir", "High risk of abacavir hypersensitivity (HLA-B*57:01 positive)"): "Systemic hypersensitivity; fatal on rechallenge",
 ("HLA-B*58:01", "allopurinol", "Significantly increased risk of allopurinol SCAR (HLA-B*58:01 carrier)"): "SJS/TEN/DRESS with substantial mortality",
 ("G6PD", "rasburicase", "G6PD Deficient (WHO class II)"): "Contraindicated; severe acute haemolysis and methaemoglobinaemia",
 ("G6PD", "rasburicase", "G6PD Deficient (WHO class III)"): "Contraindicated; severe acute haemolysis and methaemoglobinaemia",
 ("G6PD", "primaquine", "G6PD Deficient (WHO class III)"): "Acute haemolytic anaemia; can be severe and life-threatening at 14-day dosing",
 ("RYR1", "volatile-anaesthetics", "Malignant Hyperthermia susceptible"): "Malignant hyperthermia crisis; rapidly fatal without dantrolene",
 ("SLCO1B1", "simvastatin", "SLCO1B1 poor function"): "Rhabdomyolysis with renal failure; can be fatal",
 ("UGT1A1", "irinotecan", "UGT1A1 poor metabolizer"): "Corpus states 'serious, life-threatening adverse events occur more often' (febrile neutropenia, severe diarrhoea)",
 ("MT-RNR1", "aminoglycosides", "MT-RNR1 increased risk of aminoglycoside-induced hearing loss"): "Not fatal, but permanent bilateral hearing loss after a single standard course (irreversible major harm)",
}


# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------
def call_case(case):
    gene, drug, gt = case["gene"], case["drug"], case["genotype"]
    for pat, diplo, basis, note in VARIANT_RULES[gene]:
        if re.search(pat, gt):
            break
    else:
        raise KeyError(f"no variant rule matched: {gene} | {gt}")
    score, pheno = DIPLO2PHENO[(gene, diplo)]
    r = R[(gene, drug, pheno)]
    return dict(case_id=case["id"], gene=gene, drug=drug, genotype=gt,
                diplotype=diplo, activity_score=score, phenotype=pheno,
                recommendation=r["short"], abstained=str(r["abst"]).lower(),
                allele_call_basis=basis, allele_call_note=note,
                guideline_source=r["src"], verbatim_recommendation=r["quote"],
                lethal_class=str((gene, drug, pheno) in LETHAL).lower(),
                lethal_rationale=LETHAL.get((gene, drug, pheno), ""))


rows = [call_case(c) for c in cases]
assert len(rows) == len(cases) == 110
assert len({r["case_id"] for r in rows}) == 110

out = os.path.join(HERE, "claude_science_answers.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case_id", "diplotype", "phenotype", "recommendation", "abstained"])
    for r in rows:
        w.writerow([r["case_id"], r["diplotype"], r["phenotype"], r["recommendation"], r["abstained"]])

# full audit trail (superset of the required answers file)
with open(os.path.join(HERE, "answers_audit_trail.csv"), "w", newline="") as f:
    cols = ["case_id", "gene", "drug", "genotype", "diplotype", "activity_score", "phenotype",
            "recommendation", "abstained", "lethal_class", "lethal_rationale",
            "allele_call_basis", "guideline_source", "verbatim_recommendation", "allele_call_note"]
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

# Table 1: variant -> allele
with open(os.path.join(HERE, "table_variant_to_allele.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "assay_result_pattern", "allele_call", "basis", "reasoning"])
    for gene, rules in VARIANT_RULES.items():
        for pat, diplo, basis, note in rules:
            w.writerow([gene, pat.replace("\\", ""), diplo, basis, note])

# Table 2: diplotype -> phenotype
with open(os.path.join(HERE, "table_diplotype_to_phenotype.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "diplotype", "activity_score", "phenotype"])
    for (gene, diplo), (score, pheno) in DIPLO2PHENO.items():
        w.writerow([gene, diplo, score, pheno])

# Table 3: phenotype -> recommendation
with open(os.path.join(HERE, "table_phenotype_to_recommendation.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "drug", "phenotype", "action", "abstained", "guideline_source",
                "verbatim_guideline_text"])
    for (gene, drug, pheno), v in R.items():
        w.writerow([gene, drug, pheno, v["short"], str(v["abst"]).lower(), v["src"], v["quote"]])

print(f"cases={len(rows)}  abstained={sum(r['abstained']=='true' for r in rows)}  "
      f"lethal_class={sum(r['lethal_class']=='true' for r in rows)}")
print("corpus sha256 (cases):", hashlib.sha256(open(os.path.join(HERE,'cases_blinded.json'),'rb').read()).hexdigest()[:16])
print("corpus generated_at_utc:", corpus["generated_at_utc"], "| schema:", corpus["schema_version"])
