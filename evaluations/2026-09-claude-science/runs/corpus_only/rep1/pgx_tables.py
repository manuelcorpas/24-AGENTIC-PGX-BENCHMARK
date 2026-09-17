"""
pgx_tables.py -- knowledge tables for the blinded pharmacogenomics study.

Three tables, kept separate so each can be audited on its own:

  ALLELE_DEFS      marker (rsID / typing result / CNV call) -> star allele or
                   haplotype, with allele function. Source column records
                   whether the definition is stated in the supplied CPIC
                   corpus or supplied by the analyst from standard
                   pharmacogenetic allele nomenclature ("reasoning").

  GENOTYPE_RULES   per-gene ordered list of (regex, diplotype, activity_score,
                   phenotype) rules applied to the raw genotype string of a
                   case. Regexes, not string equality, because the blinded
                   cases carry the same genotype in several textual dresses.

  RECOMMENDATIONS  (gene, drug, phenotype) -> recommendation string, guideline
                   strength and the PharmGKB annotation id the wording is
                   taken from.

Nothing here is fetched from the network. Sources are (a) the supplied
cpic_guideline_corpus.json (schema v3.1, snapshot 2026-05-09/10) and
(b) the analyst's own reasoning where the corpus does not contain an
allele-definition table.
"""

CORPUS_PATH = "cpic_guideline_corpus.json"
CASES_PATH = "cases_blinded.json"

# ---------------------------------------------------------------------------
# 1. Allele / marker definitions
# ---------------------------------------------------------------------------
# columns: gene, marker, variant_allele, defines, allele_function, source
ALLELE_DEFS = [
    ("CYP2C9",  "rs1799853",   "T",        "CYP2C9*2",  "Decreased function",
     "corpus: PA166104949 allele table ('*2 | C>T at rs1799853 | R144C | Decreased')"),
    ("CYP2C9",  "rs1057910",   "C",        "CYP2C9*3",  "No/decreased function (activity score 0)",
     "corpus: PA166104949 allele table ('*3 | A>C at rs1057910 | I359L | Decreased'); "
     "AS 0 for *3 inferred from corpus phenotype rows (*1/*3 = AS 1, *3/*3 = AS 0-0.5)"),
    ("CYP2C19", "rs4244285",   "A",        "CYP2C19*2", "No function", "reasoning (c.681G>A, splice defect)"),
    ("CYP2C19", "rs12248560",  "T",        "CYP2C19*17","Increased function", "reasoning (c.-806C>T)"),
    ("CYP2D6",  "rs3892097",   "T",        "CYP2D6*4",  "No function",
     "reasoning (c.506-1G>A / 1846G>A splice defect; CYP2D6 on minus strand, forward-genomic C>T)"),
    ("CYP2D6",  "rs1065852",   "A",        "CYP2D6*10 (also carried in cis on *4)", "Decreased function (*10)",
     "reasoning (c.100C>T, P34S; forward-genomic G>A)"),
    ("CYP2D6",  "CNV *1xN",    "dup",      "CYP2D6*1xN","Increased copies of normal function allele",
     "reasoning (structural duplication; corpus lists *1/*1xN as ultrarapid example diplotype)"),
    ("CYP2B6",  "rs3745274",   "T",        "CYP2B6*9 (or *6 when in cis with rs2279343 G)", "Decreased function",
     "reasoning (c.516G>T, Q172H)"),
    ("CYP2B6",  "rs2279343",   "G",        "CYP2B6*4 (alone) / part of *6",  "Increased function (*4)",
     "reasoning (c.785A>G, K262R)"),
    ("CYP3A5",  "rs776746",    "G",        "CYP3A5*3",  "No function (non-expresser)",
     "reasoning (c.219-237A>G cryptic splice); case text states A = *1"),
    ("CYP4F2",  "rs2108622",   "T",        "CYP4F2*3",  "Not assigned by CPIC",
     "reasoning (c.1297G>A, V433M); no CYP4F2 allele or phenotype table in supplied corpus"),
    ("VKORC1",  "rs9923231",   "A",        "VKORC1 -1639A (haplotype A)", "Increased warfarin sensitivity",
     "corpus: PA166104949 dose table is indexed by rs9923231 GG/GA/AA"),
    ("DPYD",    "rs3918290",   "T",        "DPYD*2A (c.1905+1G>A)", "No function",
     "reasoning; corpus lists c.1905+1G>A as a no-function allele"),
    ("DPYD",    "rs56038477",  "T",        "DPYD HapB3 (c.1236G>A / c.1129-5923C>G)", "Decreased function",
     "reasoning for the rsID; corpus lists c.1129-5923C>G as a decreased function allele"),
    ("TPMT",    "rs1800460",   "T",        "TPMT c.460G>A (*3B; *3A when in cis with rs1142345)", "No function",
     "reasoning (A154T); case text states forward-genomic, TPMT reported on minus strand"),
    ("TPMT",    "rs1142345",   "C",        "TPMT c.719A>G (*3C; *3A when in cis with rs1800460)", "No function",
     "reasoning (Y240C)"),
    ("NUDT15",  "rs116855232", "T",        "NUDT15*3 (c.415C>T, R139C)", "No function", "reasoning"),
    ("SLCO1B1", "rs4149056",   "C",        "SLCO1B1*5 (c.521T>C, V174A; *15/*17 haplotypes)", "No function",
     "corpus: PA166105005 second table maps rs4149056 TT/TC/CC directly onto phenotype"),
    ("UGT1A1",  "rs8175347",   "TA7",      "UGT1A1*28 (promoter (TA)7)", "Decreased function",
     "reasoning; corpus lists *28 as a decreased function allele"),
    ("G6PD",    "rs1050828",   "T",        "G6PD c.202G>A (V68M); with c.376A>G in cis = G6PD A- (class III)",
     "Deficient (class III)", "reasoning; corpus lists A- as a class II-III deficient allele"),
    ("G6PD",    "rs1050829",   "C",        "G6PD c.376A>G (N126D); A (class IV) alone, A- with c.202G>A in cis",
     "Non-deficient alone", "reasoning; corpus lists B/A as a Normal genotype"),
    ("G6PD",    "rs5030868",   "A",        "G6PD Mediterranean (c.563C>T, S188F)", "Deficient (class II)",
     "reasoning; corpus lists Mediterranean as a class II deficient allele"),
    ("HLA-B*57:01", "rs2395029", "G",      "HCP5 tag SNP for HLA-B*57:01", "Risk allele carrier",
     "reasoning (HCP5 rs2395029 is in near-complete LD with HLA-B*57:01)"),
    ("HLA-A*31:01", "HLA typing", "positive", "HLA-A*31:01", "Risk allele carrier", "corpus: PA166105008"),
    ("HLA-B*15:02", "HLA typing", "positive", "HLA-B*15:02", "Risk allele carrier", "corpus: PA166105008"),
    ("HLA-B*58:01", "HLA typing", "positive", "HLA-B*58:01", "Risk allele carrier", "corpus: PA166105003"),
    ("IFNL3",   "rs12979860",  "T",        "IFNL3 (IL28B) unfavourable T allele", "Unfavourable response",
     "corpus: PA166110235 table is indexed by rs12979860 CC vs CT/TT"),
    ("MT-RNR1", "m.1555A>G",   "G",        "MT-RNR1 m.1555A>G", "Increased risk of AIHL",
     "corpus: PA166229081 lists m.1555A>G as an increased-risk variant"),
    ("RYR1",    "c.1840C>T",   "T",        "RYR1 p.Arg614Cys", "MH causative (EMHG)",
     "corpus: PA166303941 lists RYR1 c.1840C>T; p.(Arg614Cys) as an MH causative variant"),
    ("CFTR",    "rs113993960", "del",      "CFTR F508del", "Not ivacaftor-responsive",
     "corpus: PA166114461 ('Homozygous for F508del-CFTR, rs113993960 ... genotype del/del')"),
    ("CFTR",    "rs75527207",  "A",        "CFTR G551D", "Ivacaftor-responsive (gating)",
     "corpus: PA166114461 ('Homozygous or Heterozygous G551D-CFTR, rs75527207 genotype AA or AG')"),
]

# ---------------------------------------------------------------------------
# 2. Genotype-string -> diplotype / activity score / phenotype
# ---------------------------------------------------------------------------
# Each entry: (regex, diplotype, activity_score, phenotype, note)
# Rules are tried in order; the first match wins. activity_score is None where
# the gene has no CPIC activity-score system.
GENOTYPE_RULES = {
    "CYP2D6": [
        (r"\*1xN duplication",
         "*1/*1xN", ">2.25", "CYP2D6 Ultrarapid Metabolizer",
         "Copy-number call of a duplicated normal-function allele; corpus lists *1/*1xN "
         "as an ultrarapid example diplotype. Exact N and the second allele are not "
         "resolvable from a CNV assay, but any *1xN with a normal-function partner gives AS > 2.25."),
        (r"rs3892097\s+T/T",
         "*4/*4", "0", "CYP2D6 Poor Metabolizer",
         "Homozygous rs3892097 variant = two no-function *4 alleles."),
        (r"rs3892097\s+C/T\s*\+\s*rs1065852\s+G/A",
         "*1/*4 (unphased; *4-without-c.100C>T / *10 combination not formally excluded)",
         "1.0", "CYP2D6 Intermediate Metabolizer",
         "c.100C>T (rs1065852) is carried in cis on the *4 haplotype, so a *1/*4 subject is "
         "heterozygous at both markers; *4/*10 would require rs1065852 A/A and is excluded. "
         "The residual unphased alternative (*4 lacking c.100C>T together with *10) would give "
         "AS 0.25, which is also Intermediate Metabolizer, so the phenotype call is unaffected."),
        (r"rs3892097\s+C/T",
         "*1/*4", "1.0", "CYP2D6 Intermediate Metabolizer",
         "One no-function *4 allele plus one normal-function allele; AS 1.0. Assigned by the "
         "current CPIC activity-score translation (codeine/paroxetine tables in the corpus: "
         "IM = 0 < AS < 1.25). The older amitriptyline/ondansetron/tamoxifen tables in the same "
         "corpus place *1/*4 in a Normal-metabolizer band; those annotations state that CPIC has "
         "since updated the CYP2D6 genotype-phenotype translation."),
        (r"All CYP2D6 SNPs homozygous reference",
         "*1/*1", "2.0", "CYP2D6 Normal Metabolizer",
         "No variant detected on the assayed SNP panel. A SNP panel cannot exclude the "
         "gene-deletion (*5) or duplication structural alleles, so *1/*1 is the panel-level call."),
    ],
    "CYP2C19": [
        (r"rs4244285\s+A/A", "*2/*2", None, "CYP2C19 Poor Metabolizer",
         "Homozygous *2 (two no-function alleles)."),
        (r"rs4244285\s+G/A", "*1/*2", None, "CYP2C19 Intermediate Metabolizer",
         "One normal-function plus one no-function allele."),
        (r"rs12248560\s+T/T", "*17/*17", None, "CYP2C19 Ultrarapid Metabolizer",
         "Two increased-function *17 alleles."),
        (r"rs12248560\s+C/T", "*1/*17", None, "CYP2C19 Rapid Metabolizer",
         "One normal-function plus one increased-function allele."),
        (r"All CYP2C19 SNPs homozygous reference", "*1/*1", None, "CYP2C19 Normal Metabolizer",
         "No variant detected on the assayed SNP panel."),
    ],
    "CYP2C9": [
        (r"rs1057910\s+C/C", "*3/*3", "0", "CYP2C9 Poor Metabolizer",
         "Two *3 alleles; corpus phenytoin table gives *3/*3 activity score 0-0.5 = Poor Metabolizer."),
        (r"rs1057910\s+A/C", "*1/*3", "1", "CYP2C9 Intermediate Metabolizer (activity score 1)",
         "One normal-function plus one *3 allele; corpus phenytoin and celecoxib tables both "
         "give *1/*3 = activity score 1 Intermediate Metabolizer."),
        (r"All CYP2C9 SNPs homozygous reference", "*1/*1", "2", "CYP2C9 Normal Metabolizer",
         "No variant detected on the assayed SNP panel."),
    ],
    "CYP2B6": [
        (r"rs3745274\s+T/T",
         "*6/*6 or *9/*9 (c.516G>T homozygous; rs2279343 not reported)", None,
         "CYP2B6 Poor Metabolizer",
         "c.516G>T defines *9 on its own and *6 when in cis with c.785A>G (rs2279343). "
         "rs2279343 was not reported, so *6/*6 and *9/*9 cannot be distinguished; both are "
         "decreased-function alleles, so the phenotype (two decreased-function alleles = Poor "
         "Metabolizer) is the same either way."),
        (r"rs3745274\s+G/T",
         "*1/*6 or *1/*9 (c.516G>T heterozygous; rs2279343 not reported)", None,
         "CYP2B6 Intermediate Metabolizer",
         "One normal-function plus one decreased-function allele either way."),
        (r"rs2279343\s+A/G", "*1/*4", None, "CYP2B6 Rapid Metabolizer",
         "c.785A>G alone defines the increased-function *4 allele; c.516G>T is reference, "
         "which excludes *6. One normal plus one increased-function allele."),
        (r"All CYP2B6 SNPs homozygous reference", "*1/*1", None, "CYP2B6 Normal Metabolizer",
         "No variant detected on the assayed SNP panel."),
    ],
    "CYP3A5": [
        (r"rs776746\s+A/A", "*1/*1", None, "CYP3A5 Extensive Metabolizer (CYP3A5 expresser)",
         "Case text states A = *1; two functional alleles."),
        (r"rs776746\s+A/G", "*1/*3", None, "CYP3A5 Intermediate Metabolizer (CYP3A5 expresser)",
         "One functional plus one non-functional (*3) allele."),
        (r"rs776746\s+G/G", "*3/*3", None, "CYP3A5 Poor Metabolizer (CYP3A5 non-expresser)",
         "Two non-functional *3 alleles."),
    ],
    "CYP4F2": [
        (r"rs2108622\s+T/T", "*3/*3", None, "No CPIC phenotype defined for CYP4F2",
         "CYP4F2*3 homozygous. The supplied warfarin excerpt (PA166104949) contains only the "
         "VKORC1 x CYP2C9 dose table and the CYP2C9 allele table; it carries no CYP4F2 allele "
         "table, phenotype assignment or CYP4F2-specific dosing action."),
        (r"rs2108622\s+C/C", "*1/*1", None, "No CPIC phenotype defined for CYP4F2",
         "CYP4F2 reference. Same corpus limitation as above."),
    ],
    "VKORC1": [
        (r"rs9923231\s+A/A", "-1639A/-1639A (rs9923231 A/A)", None,
         "Increased warfarin sensitivity (low dose requirement)",
         "Read off the corpus warfarin dose table, whose AA row carries the lowest predicted "
         "stable doses at every CYP2C9 genotype."),
        (r"rs9923231\s+G/G", "-1639G/-1639G (rs9923231 G/G)", None,
         "Normal warfarin sensitivity (reference)",
         "GG row of the corpus warfarin dose table."),
    ],
    "DPYD": [
        (r"rs3918290\s+T/T", "c.[1905+1G>A];[1905+1G>A] (*2A/*2A)", "0", "DPYD Poor Metabolizer",
         "Two no-function alleles; corpus lists c.[1905+1G>A];[1905+1G>A] as an activity score 0 example."),
        (r"rs3918290\s+C/T", "c.[1905+1G>A];[=] (*1/*2A)", "1", "DPYD Intermediate Metabolizer",
         "One no-function plus one normal allele; corpus lists c.[1905+1G>A];[=] as an "
         "activity score 1 example."),
        (r"rs56038477\s+C/T", "c.[1129-5923C>G];[=] (HapB3 heterozygous)", "1.5",
         "DPYD Intermediate Metabolizer",
         "rs56038477 (c.1236G>A) is the HapB3 tag marker, in complete LD with the causal "
         "c.1129-5923C>G decreased-function variant; corpus lists c.[1129-5923C>G];[=] as an "
         "intermediate metaboliser example (activity score 1.5)."),
        (r"All DPYD risk variants reference", "c.[=];[=] (*1/*1)", "2", "DPYD Normal Metabolizer",
         "No risk variant detected."),
    ],
    "TPMT": [
        (r"rs1800460\s+T/T\s*\+\s*rs1142345\s+C/C", "*3A/*3A", None, "TPMT Poor Metabolizer",
         "Homozygous variant at both c.460G>A and c.719A>G. The two variants occur in cis on "
         "the *3A haplotype, so the homozygous-both genotype is *3A/*3A (two no-function "
         "alleles). The formally possible *3B/*3C compound configuration is a recognised "
         "rarity and is also a two-no-function-allele Poor Metabolizer."),
        (r"rs1800460\s+C/T", "*1/*3B", None, "TPMT Intermediate Metabolizer",
         "Heterozygous c.460G>A only; one no-function plus one normal-function allele."),
        (r"rs1142345\s+T/C", "*1/*3C", None, "TPMT Intermediate Metabolizer",
         "Heterozygous c.719A>G only; one no-function plus one normal-function allele."),
        (r"All TPMT SNPs homozygous reference", "*1/*1", None, "TPMT Normal Metabolizer",
         "No variant detected on the assayed SNP panel."),
    ],
    "NUDT15": [
        (r"rs116855232\s+T/T", "*3/*3", None, "NUDT15 Poor Metabolizer",
         "Two no-function *3 alleles."),
        (r"rs116855232\s+C/T", "*1/*3", None, "NUDT15 Intermediate Metabolizer",
         "One normal plus one no-function allele."),
        (r"rs116855232\s+C/C", "*1/*1", None, "NUDT15 Normal Metabolizer", "Reference."),
    ],
    "SLCO1B1": [
        (r"rs4149056\s+C/C", "*5/*5 (or *15/*15; rs4149056 C/C)", None, "SLCO1B1 Poor Function",
         "Two no-function c.521C alleles; corpus maps rs4149056 CC directly onto the "
         "low/poor function phenotype."),
        (r"rs4149056\s+T/C", "*1/*5 (or *1/*15; rs4149056 T/C)", None, "SLCO1B1 Decreased Function",
         "One normal plus one no-function allele; corpus maps rs4149056 TC onto the "
         "intermediate/decreased function phenotype."),
        (r"rs4149056\s+T/T", "*1/*1 (rs4149056 T/T)", None, "SLCO1B1 Normal Function",
         "Two normal-function alleles."),
    ],
    "UGT1A1": [
        (r"TA7/TA7", "*28/*28", None, "UGT1A1 Poor Metabolizer",
         "Two decreased-function (TA)7 promoter alleles."),
        (r"TA6/TA7", "*1/*28", None, "UGT1A1 Intermediate Metabolizer",
         "One reference plus one decreased-function allele."),
        (r"TA6/TA6", "*1/*1", None, "UGT1A1 Extensive (Normal) Metabolizer",
         "Two reference-function alleles."),
    ],
    "G6PD": [
        (r"rs5030868 A", "Mediterranean (hemizygous male) or Mediterranean/Mediterranean", None,
         "G6PD Deficient",
         "c.563C>T = G6PD Mediterranean, a class II deficient allele. One X carrying a class II "
         "allele, or two class II alleles, is the corpus Deficient row."),
        (r"rs1050828 C/T\s*\+\s*rs1050829", "B/A- (heterozygous female)", None, "G6PD Variable",
         "One non-deficient class IV allele (B) and one deficient class III allele (A-). The "
         "corpus Variable row is exactly this configuration and notes that X-linked mosaicism "
         "means an enzyme activity test is needed to assign the phenotype."),
        (r"rs1050828 T\s*\+\s*rs1050829 C in cis", "A- (hemizygous male) or A-/A-", None,
         "G6PD Deficient",
         "c.202G>A and c.376A>G in cis = G6PD A-, a class III deficient allele."),
        (r"rs1050828 C\s*\+\s*rs5030868 G", "B (hemizygous male) or B/B", None, "G6PD Normal",
         "No deficiency variant detected; class IV (B) allele(s) only."),
    ],
    "HLA-A*31:01": [
        (r"positive", "HLA-A*31:01 positive (*31:01/*X or *31:01/*31:01)", None,
         "HLA-A*31:01 positive (increased risk of carbamazepine-induced SJS/TEN, DRESS and MPE)",
         "Corpus carrier row."),
        (r"negative", "HLA-A*31:01 negative (*X/*X)", None,
         "HLA-A*31:01 negative (normal or reduced risk)", "Corpus non-carrier row."),
    ],
    "HLA-B*15:02": [
        (r"positive", "HLA-B*15:02 positive (*15:02/*X or *15:02/*15:02)", None,
         "HLA-B*15:02 positive (increased risk of SJS/TEN)", "Corpus carrier row."),
        (r"negative", "HLA-B*15:02 negative (*X/*X)", None,
         "HLA-B*15:02 negative (normal or reduced risk)", "Corpus non-carrier row."),
    ],
    "HLA-B*58:01": [
        (r"positive", "HLA-B*58:01 positive (*58:01/*X or *58:01/*58:01)", None,
         "HLA-B*58:01 positive (significantly increased risk of allopurinol-induced SCAR)",
         "Corpus carrier row."),
        (r"negative", "HLA-B*58:01 negative (*X/*X)", None,
         "HLA-B*58:01 negative (low or reduced risk of allopurinol SCAR)", "Corpus non-carrier row."),
    ],
    "HLA-B*57:01": [
        (r"rs2395029 G/T or G/G", "HLA-B*57:01 positive (*57:01/*X or *57:01/*57:01)", None,
         "HLA-B*57:01 positive (high risk of abacavir hypersensitivity)",
         "rs2395029 (HCP5) G is the tag allele for HLA-B*57:01; a positive tag result is "
         "treated as the corpus carrier row. A tag SNP is not a typing result, so the "
         "diplotype is reported as *57:01-positive rather than as a resolved HLA diplotype."),
        (r"rs2395029 T/T", "HLA-B*57:01 negative (*X/*X)", None,
         "HLA-B*57:01 negative (very low risk of abacavir hypersensitivity)",
         "No HCP5 tag risk allele."),
    ],
    "IFNL3": [
        (r"rs12979860 C/C", "rs12979860 C/C", None, "IFNL3 favourable response genotype",
         "Corpus CC row."),
        (r"rs12979860 T/T", "rs12979860 T/T", None, "IFNL3 unfavourable response genotype",
         "Corpus 'CT or TT' unfavourable row."),
    ],
    "MT-RNR1": [
        (r"m\.1555A>G", "m.1555A>G (homoplasmic)", None,
         "MT-RNR1 increased risk of aminoglycoside-induced hearing loss",
         "Corpus lists m.1555A>G among the increased-risk variants."),
        (r"m\.1555A\b", "m.1555A (reference)", None,
         "MT-RNR1 normal risk of aminoglycoside-induced hearing loss",
         "No detectable increased-risk variant."),
    ],
    "RYR1": [
        (r"RYR1 mutation", "RYR1 c.1840C>T p.(Arg614Cys) heterozygous", None,
         "Malignant Hyperthermia susceptible",
         "Corpus lists RYR1 c.1840C>T; p.(Arg614Cys) as an EMHG malignant-hyperthermia "
         "causative variant."),
        (r"no pathogenic variant", "RYR1 no MH causative variant detected", None,
         "Uncertain susceptibility",
         "Corpus assigns 'Uncertain susceptibility', not 'not susceptible', to a negative "
         "result: about half of contracture-test-confirmed MH has no identified genetic cause."),
    ],
    "CFTR": [
        (r"G551D positive", "G551D/other (rs75527207 A/G or A/A)", None,
         "CFTR gating variant present (ivacaftor-responsive)",
         "Corpus row 'Homozygous or Heterozygous G551D-CFTR, rs75527207 genotype AA or AG'."),
        (r"F508del homozygous", "F508del/F508del", None,
         "CFTR F508del homozygous (not ivacaftor-responsive)",
         "Corpus row 'Homozygous for F508del-CFTR'."),
    ],
}

# ---------------------------------------------------------------------------
# 3. Phenotype x drug -> recommendation
# ---------------------------------------------------------------------------
# key: (gene, drug, phenotype) -> (recommendation, strength, source_pa_id)
# `recommendation` is written as "<drug>: <action>" as required by the brief.
R = {}

def _add(gene, drug, pheno, rec, strength, pa):
    R[(gene, drug, pheno)] = (rec, strength, pa)

# --- CYP2D6 ---------------------------------------------------------------
_add("CYP2D6", "codeine", "CYP2D6 Ultrarapid Metabolizer",
     "codeine: avoid; potential for serious toxicity from increased morphine formation - "
     "if an opioid is warranted use a non-tramadol opioid (e.g. morphine or a non-opioid analgesic)",
     "Strong", "PA166104996")
_add("CYP2D6", "codeine", "CYP2D6 Normal Metabolizer",
     "codeine: use label-recommended age- or weight-specific dosing", "Strong", "PA166104996")
_add("CYP2D6", "codeine", "CYP2D6 Intermediate Metabolizer",
     "codeine: use label-recommended age- or weight-specific dosing; if no response and an "
     "opioid is warranted, consider a non-tramadol opioid", "Moderate", "PA166104996")
_add("CYP2D6", "codeine", "CYP2D6 Poor Metabolizer",
     "codeine: avoid because of diminished analgesia; if an opioid is warranted use a "
     "non-tramadol opioid", "Strong", "PA166104996")

_add("CYP2D6", "tamoxifen", "CYP2D6 Ultrarapid Metabolizer",
     "tamoxifen: initiate standard of care dosing (20 mg/day) and avoid moderate and strong "
     "CYP2D6 inhibitors", "Strong", "PA166176068")
_add("CYP2D6", "tamoxifen", "CYP2D6 Normal Metabolizer",
     "tamoxifen: initiate standard of care dosing (20 mg/day) and avoid moderate and strong "
     "CYP2D6 inhibitors", "Strong", "PA166176068")
_add("CYP2D6", "tamoxifen", "CYP2D6 Intermediate Metabolizer",
     "tamoxifen: consider an aromatase inhibitor (postmenopausal) or an aromatase inhibitor "
     "plus ovarian function suppression (premenopausal); if an aromatase inhibitor is "
     "contraindicated consider tamoxifen 40 mg/day; avoid CYP2D6 inhibitors",
     "Optional (activity score 1.0 row, no *10 allele)", "PA166176068")
_add("CYP2D6", "tamoxifen", "CYP2D6 Poor Metabolizer",
     "tamoxifen: recommend alternative hormonal therapy - an aromatase inhibitor "
     "(postmenopausal) or an aromatase inhibitor plus ovarian function suppression "
     "(premenopausal); tamoxifen 40 mg/day only if an aromatase inhibitor is contraindicated",
     "Strong", "PA166176068")

_add("CYP2D6", "ondansetron", "CYP2D6 Ultrarapid Metabolizer",
     "ondansetron: select an alternative antiemetic not predominantly metabolised by CYP2D6 "
     "(e.g. granisetron)", "Moderate", "PA166161954")
_add("CYP2D6", "ondansetron", "CYP2D6 Normal Metabolizer",
     "ondansetron: initiate therapy with the recommended starting dose", "Strong", "PA166161954")
_add("CYP2D6", "ondansetron", "CYP2D6 Intermediate Metabolizer",
     "ondansetron: initiate therapy with the recommended starting dose (insufficient evidence "
     "of clinical impact for intermediate metabolisers)", "No recommendation", "PA166161954")
_add("CYP2D6", "ondansetron", "CYP2D6 Poor Metabolizer",
     "ondansetron: initiate therapy with the recommended starting dose (insufficient evidence "
     "of clinical impact for poor metabolisers)", "No recommendation", "PA166161954")

_add("CYP2D6", "paroxetine", "CYP2D6 Ultrarapid Metabolizer",
     "paroxetine: select an alternative antidepressant not predominantly metabolised by CYP2D6",
     "Moderate", "PA166127636")
_add("CYP2D6", "paroxetine", "CYP2D6 Normal Metabolizer",
     "paroxetine: initiate therapy with the recommended starting dose", "Strong", "PA166127636")
_add("CYP2D6", "paroxetine", "CYP2D6 Intermediate Metabolizer",
     "paroxetine: consider a lower starting dose and a slower titration schedule than for "
     "normal metabolisers", "Optional", "PA166127636")
_add("CYP2D6", "paroxetine", "CYP2D6 Poor Metabolizer",
     "paroxetine: consider a 50% reduction in starting dose, slower titration and a 50% lower "
     "maintenance dose", "Moderate", "PA166127636")

_add("CYP2D6", "amitriptyline", "CYP2D6 Ultrarapid Metabolizer",
     "amitriptyline: avoid (potential lack of efficacy); consider an alternative drug not "
     "metabolised by CYP2D6 [assumes CYP2C19 normal metaboliser - CYP2C19 not genotyped]",
     "Strong", "PA166105006")
_add("CYP2D6", "amitriptyline", "CYP2D6 Normal Metabolizer",
     "amitriptyline: initiate therapy with the recommended starting dose [assumes CYP2C19 "
     "normal metaboliser - CYP2C19 not genotyped]", "Strong", "PA166105006")
_add("CYP2D6", "amitriptyline", "CYP2D6 Intermediate Metabolizer",
     "amitriptyline: consider a 25% reduction of the recommended starting dose with "
     "therapeutic drug monitoring [assumes CYP2C19 normal metaboliser - CYP2C19 not genotyped]",
     "Moderate", "PA166105006")
_add("CYP2D6", "amitriptyline", "CYP2D6 Poor Metabolizer",
     "amitriptyline: avoid (potential for side effects); if a TCA is warranted consider a 50% "
     "reduction of the recommended starting dose with therapeutic drug monitoring [assumes "
     "CYP2C19 normal metaboliser - CYP2C19 not genotyped]", "Strong", "PA166105006")

# --- CYP2C19 --------------------------------------------------------------
_add("CYP2C19", "clopidogrel", "CYP2C19 Ultrarapid Metabolizer",
     "clopidogrel: use at standard dose (75 mg/day)", "Strong", "PA166104948")
_add("CYP2C19", "clopidogrel", "CYP2C19 Rapid Metabolizer",
     "clopidogrel: use at standard dose (75 mg/day)", "Strong", "PA166104948")
_add("CYP2C19", "clopidogrel", "CYP2C19 Normal Metabolizer",
     "clopidogrel: use at standard dose (75 mg/day)", "Strong", "PA166104948")
_add("CYP2C19", "clopidogrel", "CYP2C19 Intermediate Metabolizer",
     "clopidogrel: avoid standard dose if possible - use prasugrel or ticagrelor at standard "
     "dose if not contraindicated", "Strong", "PA166104948")
_add("CYP2C19", "clopidogrel", "CYP2C19 Poor Metabolizer",
     "clopidogrel: avoid if possible - use prasugrel or ticagrelor at standard dose if not "
     "contraindicated", "Strong", "PA166104948")

_add("CYP2C19", "voriconazole", "CYP2C19 Ultrarapid Metabolizer",
     "voriconazole: choose an alternative agent not dependent on CYP2C19 metabolism "
     "(isavuconazole, liposomal amphotericin B, posaconazole)", "Moderate", "PA166161537")
_add("CYP2C19", "voriconazole", "CYP2C19 Rapid Metabolizer",
     "voriconazole: in adults choose an alternative agent not dependent on CYP2C19 metabolism "
     "(isavuconazole, liposomal amphotericin B, posaconazole); in children initiate standard "
     "dosing with therapeutic drug monitoring", "Moderate", "PA166161537")
_add("CYP2C19", "voriconazole", "CYP2C19 Normal Metabolizer",
     "voriconazole: initiate therapy with standard of care dosing", "Strong", "PA166161537")
_add("CYP2C19", "voriconazole", "CYP2C19 Intermediate Metabolizer",
     "voriconazole: initiate therapy with standard of care dosing", "Moderate", "PA166161537")
_add("CYP2C19", "voriconazole", "CYP2C19 Poor Metabolizer",
     "voriconazole: choose an alternative agent not dependent on CYP2C19 metabolism; if "
     "voriconazole is still the most appropriate agent, use a lower than standard dose with "
     "careful therapeutic drug monitoring", "Moderate", "PA166161537")

_add("CYP2C19", "citalopram", "CYP2C19 Ultrarapid Metabolizer",
     "citalopram: consider an alternative antidepressant not predominantly metabolised by "
     "CYP2C19; if citalopram is used and efficacy is inadequate at standard maintenance "
     "dosing, consider titrating to a higher maintenance dose", "Strong", "PA166127638")
_add("CYP2C19", "citalopram", "CYP2C19 Rapid Metabolizer",
     "citalopram: initiate with the recommended starting dose; if response is inadequate, "
     "consider a higher maintenance dose or an alternative antidepressant not predominantly "
     "metabolised by CYP2C19", "Optional", "PA166127638")
_add("CYP2C19", "citalopram", "CYP2C19 Normal Metabolizer",
     "citalopram: initiate therapy with the recommended starting dose", "Strong", "PA166127638")
_add("CYP2C19", "citalopram", "CYP2C19 Intermediate Metabolizer",
     "citalopram: initiate with the recommended starting dose but consider a slower titration "
     "schedule and a lower maintenance dose than for normal metabolisers", "Moderate",
     "PA166127638")
_add("CYP2C19", "citalopram", "CYP2C19 Poor Metabolizer",
     "citalopram: consider an alternative antidepressant not predominantly metabolised by "
     "CYP2C19; if citalopram is clinically appropriate use a lower starting dose, slower "
     "titration and a 50% reduced maintenance dose (FDA maximum 20 mg/day in CYP2C19 poor "
     "metabolisers because of QT prolongation)", "Strong", "PA166127638")

_add("CYP2C19", "omeprazole", "CYP2C19 Ultrarapid Metabolizer",
     "omeprazole: increase the starting daily dose by 100% (may be divided) and monitor for "
     "efficacy", "Optional", "PA166219103")
_add("CYP2C19", "omeprazole", "CYP2C19 Rapid Metabolizer",
     "omeprazole: initiate the standard starting daily dose; consider a 50-100% dose increase "
     "for H. pylori eradication or erosive oesophagitis and monitor for efficacy", "Moderate",
     "PA166219103")
_add("CYP2C19", "omeprazole", "CYP2C19 Normal Metabolizer",
     "omeprazole: initiate the standard starting daily dose; consider a 50-100% dose increase "
     "for H. pylori eradication or erosive oesophagitis and monitor for efficacy", "Moderate",
     "PA166219103")
_add("CYP2C19", "omeprazole", "CYP2C19 Intermediate Metabolizer",
     "omeprazole: initiate the standard starting daily dose; for chronic therapy (>12 weeks) "
     "with efficacy achieved, consider a 50% dose reduction", "Optional", "PA166219103")
_add("CYP2C19", "omeprazole", "CYP2C19 Poor Metabolizer",
     "omeprazole: initiate the standard starting daily dose; for chronic therapy (>12 weeks) "
     "with efficacy achieved, consider a 50% dose reduction", "Moderate", "PA166219103")

# --- CYP2C9 ---------------------------------------------------------------
_add("CYP2C9", "celecoxib", "CYP2C9 Normal Metabolizer",
     "celecoxib: initiate therapy with the recommended starting dose, at the lowest effective "
     "dosage for the shortest duration", "Strong", "PA166191841")
_add("CYP2C9", "celecoxib", "CYP2C9 Intermediate Metabolizer (activity score 1)",
     "celecoxib: initiate therapy with the lowest recommended starting dose and titrate upward "
     "with caution, monitoring blood pressure and kidney function", "Moderate", "PA166191841")
_add("CYP2C9", "celecoxib", "CYP2C9 Poor Metabolizer",
     "celecoxib: initiate therapy with 25-50% of the lowest recommended starting dose (upward "
     "titration only after steady state, at least 8 days), or alternatively use an NSAID not "
     "metabolised by CYP2C9 (aspirin, ketorolac, naproxen, sulindac)", "Moderate", "PA166191841")

_add("CYP2C9", "phenytoin", "CYP2C9 Normal Metabolizer",
     "phenytoin: no adjustment from typical dosing; adjust by therapeutic drug monitoring "
     "[assumes HLA-B*15:02 negative - not genotyped]", "Strong", "PA166122806")
_add("CYP2C9", "phenytoin", "CYP2C9 Intermediate Metabolizer (activity score 1)",
     "phenytoin: give a typical initial/loading dose, then use approximately 25% less than the "
     "typical maintenance dose with therapeutic drug monitoring [assumes HLA-B*15:02 negative - "
     "not genotyped]", "Moderate", "PA166122806")
_add("CYP2C9", "phenytoin", "CYP2C9 Poor Metabolizer",
     "phenytoin: give a typical initial/loading dose, then use approximately 50% less than the "
     "typical maintenance dose with therapeutic drug monitoring [assumes HLA-B*15:02 negative - "
     "not genotyped]", "Strong", "PA166122806")

_add("CYP2C9", "warfarin", "CYP2C9 Normal Metabolizer",
     "warfarin: dose from the CPIC VKORC1 x CYP2C9 table - 5-7 mg/day for *1/*1 with VKORC1 "
     "-1639 GG or GA, 3-4 mg/day with VKORC1 AA [VKORC1 not genotyped in this case]",
     "Strong", "PA166104949")
_add("CYP2C9", "warfarin", "CYP2C9 Intermediate Metabolizer (activity score 1)",
     "warfarin: dose from the CPIC VKORC1 x CYP2C9 table - 3-4 mg/day for *1/*3 with VKORC1 "
     "-1639 GG or GA, 0.5-2 mg/day with VKORC1 AA [VKORC1 not genotyped in this case]",
     "Strong", "PA166104949")
_add("CYP2C9", "warfarin", "CYP2C9 Poor Metabolizer",
     "warfarin: dose from the CPIC VKORC1 x CYP2C9 table - 0.5-2 mg/day for *3/*3 at every "
     "VKORC1 -1639 genotype [VKORC1 not genotyped in this case]", "Strong", "PA166104949")

# --- VKORC1 ---------------------------------------------------------------
_add("VKORC1", "warfarin", "Increased warfarin sensitivity (low dose requirement)",
     "warfarin: use a reduced starting dose - CPIC table predicts 3-4 mg/day with CYP2C9 "
     "*1/*1 or *1/*2 and 0.5-2 mg/day with any *3-containing genotype [CYP2C9 not genotyped "
     "in this case]", "Strong", "PA166104949")
_add("VKORC1", "warfarin", "Normal warfarin sensitivity (reference)",
     "warfarin: CPIC table predicts 5-7 mg/day with CYP2C9 *1/*1 or *1/*2, 3-4 mg/day with "
     "*1/*3, *2/*2 or *2/*3, and 0.5-2 mg/day with *3/*3 [CYP2C9 not genotyped in this case]",
     "Strong", "PA166104949")

# --- CYP2B6 ---------------------------------------------------------------
_add("CYP2B6", "efavirenz", "CYP2B6 Normal Metabolizer",
     "efavirenz: initiate with standard dosing (600 mg/day)", "Strong", "PA166182603")
_add("CYP2B6", "efavirenz", "CYP2B6 Rapid Metabolizer",
     "efavirenz: initiate with standard dosing (600 mg/day)", "Strong", "PA166182603")
_add("CYP2B6", "efavirenz", "CYP2B6 Intermediate Metabolizer",
     "efavirenz: consider initiating with a decreased dose of 400 mg/day", "Moderate",
     "PA166182603")
_add("CYP2B6", "efavirenz", "CYP2B6 Poor Metabolizer",
     "efavirenz: consider initiating with a decreased dose of 400 or 200 mg/day", "Moderate",
     "PA166182603")

# --- CYP3A5 ---------------------------------------------------------------
_add("CYP3A5", "tacrolimus", "CYP3A5 Extensive Metabolizer (CYP3A5 expresser)",
     "tacrolimus: increase the starting dose to 1.5-2 times the recommended starting dose "
     "(not exceeding 0.3 mg/kg/day) and use therapeutic drug monitoring", "Strong", "PA166124619")
_add("CYP3A5", "tacrolimus", "CYP3A5 Intermediate Metabolizer (CYP3A5 expresser)",
     "tacrolimus: increase the starting dose to 1.5-2 times the recommended starting dose "
     "(not exceeding 0.3 mg/kg/day) and use therapeutic drug monitoring", "Strong", "PA166124619")
_add("CYP3A5", "tacrolimus", "CYP3A5 Poor Metabolizer (CYP3A5 non-expresser)",
     "tacrolimus: initiate therapy with the standard recommended dose and use therapeutic drug "
     "monitoring", "Strong", "PA166124619")

# --- DPYD -----------------------------------------------------------------
for _drug, _pa in (("fluorouracil", "PA166122686"), ("capecitabine", "PA166109594")):
    _add("DPYD", _drug, "DPYD Normal Metabolizer",
         "%s: no genotype-based change - use label-recommended dosage and administration" % _drug,
         "Strong", _pa)
    _add("DPYD", _drug, "DPYD Intermediate Metabolizer",
         "%s: reduce the starting dose by 50%%, then titrate by toxicity or therapeutic drug "
         "monitoring" % _drug,
         "Strong (activity score 1) / Moderate (activity score 1.5)", _pa)
    _add("DPYD", _drug, "DPYD Poor Metabolizer",
         "%s: avoid - use an alternative, non-fluoropyrimidine regimen (activity score 0)" % _drug,
         "Strong", _pa)

# --- TPMT / NUDT15 --------------------------------------------------------
_add("TPMT", "azathioprine", "TPMT Normal Metabolizer",
     "azathioprine: start at the normal standard dose", "Strong", "PA166104933")
_add("TPMT", "azathioprine", "TPMT Intermediate Metabolizer",
     "azathioprine: consider a lower than standard starting dose (further reduction if the "
     "patient is also a NUDT15 intermediate metaboliser - NUDT15 not genotyped)", "Strong",
     "PA166104933")
_add("TPMT", "azathioprine", "TPMT Poor Metabolizer",
     "azathioprine: consider an alternative immunosuppressant (drastic dose reduction is "
     "required if a thiopurine is used)", "Strong", "PA166104933")
_add("NUDT15", "thiopurines", "NUDT15 Normal Metabolizer",
     "thiopurines (azathioprine/mercaptopurine/thioguanine): start at the normal standard dose",
     "Strong", "PA166104933")
_add("NUDT15", "thiopurines", "NUDT15 Intermediate Metabolizer",
     "thiopurines (azathioprine/mercaptopurine/thioguanine): consider a lower than standard "
     "starting dose (further reduction if the patient is also a TPMT intermediate metaboliser - "
     "TPMT not genotyped)", "Strong", "PA166104933")
_add("NUDT15", "thiopurines", "NUDT15 Poor Metabolizer",
     "thiopurines (azathioprine/mercaptopurine/thioguanine): consider an alternative "
     "immunosuppressant", "Strong", "PA166104933")

# --- SLCO1B1 --------------------------------------------------------------
_add("SLCO1B1", "simvastatin", "SLCO1B1 Normal Function",
     "simvastatin: prescribe the desired starting dose and adjust by disease-specific "
     "guidelines", "Strong", "PA166105005")
_add("SLCO1B1", "simvastatin", "SLCO1B1 Decreased Function",
     "simvastatin: prescribe an alternative statin of the desired potency; if simvastatin is "
     "warranted limit the dose to <20 mg/day", "Strong", "PA166105005")
_add("SLCO1B1", "simvastatin", "SLCO1B1 Poor Function",
     "simvastatin: prescribe an alternative statin of the desired potency (highly increased "
     "myopathy risk)", "Strong", "PA166105005")
_add("SLCO1B1", "atorvastatin", "SLCO1B1 Normal Function",
     "atorvastatin: prescribe the desired starting dose and adjust by disease-specific "
     "guidelines", "Strong", "PA166262221")
_add("SLCO1B1", "atorvastatin", "SLCO1B1 Decreased Function",
     "atorvastatin: prescribe <=40 mg as the starting dose, with awareness of increased "
     "myopathy risk especially at 40 mg", "Moderate", "PA166262221")
_add("SLCO1B1", "atorvastatin", "SLCO1B1 Poor Function",
     "atorvastatin: prescribe <=20 mg as the starting dose; if >20 mg is needed consider "
     "rosuvastatin or combination therapy", "Moderate", "PA166262221")

# --- UGT1A1 ---------------------------------------------------------------
_add("UGT1A1", "atazanavir", "UGT1A1 Extensive (Normal) Metabolizer",
     "atazanavir: no need to avoid atazanavir on the basis of UGT1A1 genotype", "Strong",
     "PA166128738")
_add("UGT1A1", "atazanavir", "UGT1A1 Intermediate Metabolizer",
     "atazanavir: no need to avoid atazanavir; inform the patient that jaundice is unlikely "
     "with this genotype", "Strong", "PA166128738")
_add("UGT1A1", "atazanavir", "UGT1A1 Poor Metabolizer",
     "atazanavir: consider an alternative agent, particularly where jaundice would be of "
     "concern to the patient (high likelihood of bilirubin-related discontinuation)", "Strong",
     "PA166128738")
_add("UGT1A1", "irinotecan", "UGT1A1 Extensive (Normal) Metabolizer",
     "irinotecan: no action needed for this gene-drug interaction (standard dosing; *1/*1 is "
     "the only genotype cleared for intensified doses >=240 mg/m2)",
     "DPWG: no action / RNPGx", "PA166104951 + PA166127626 (no CPIC guideline for irinotecan)")
_add("UGT1A1", "irinotecan", "UGT1A1 Intermediate Metabolizer",
     "irinotecan: no action needed for this gene-drug interaction (standard dosing)",
     "DPWG: no action / RNPGx", "PA166104951 + PA166127626 (no CPIC guideline for irinotecan)")
_add("UGT1A1", "irinotecan", "UGT1A1 Poor Metabolizer",
     "irinotecan: reduce the starting dose - begin at 70% of the normal dose and escalate as "
     "tolerated, guided by neutrophil count; intensified doses >=240 mg/m2 are contraindicated",
     "DPWG / RNPGx (life-threatening adverse events)",
     "PA166104951 + PA166127626 (no CPIC guideline for irinotecan)")

# --- G6PD -----------------------------------------------------------------
_add("G6PD", "primaquine", "G6PD Normal",
     "primaquine: no reason to avoid primaquine on the basis of G6PD status", "Strong",
     "PA166279621")
_add("G6PD", "primaquine", "G6PD Deficient",
     "primaquine: avoid, except under expert malaria-treatment consensus regimens - "
     "0.75 mg/kg once weekly for 8 weeks (WHO) or 45 mg once weekly for 8 weeks (CDC) with "
     "close monitoring for haemolysis for P. vivax/P. ovale radical cure, or a single "
     "0.25 mg/kg gametocytocidal dose for P. falciparum", "Strong", "PA166279621")
_add("G6PD", "primaquine", "G6PD Variable",
     "primaquine: measure G6PD enzyme activity before prescribing and dose per the "
     "activity-based phenotype (X-linked mosaicism makes the genotype non-predictive)",
     "Moderate", "PA166279621")
_add("G6PD", "rasburicase", "G6PD Normal",
     "rasburicase: no reason to withhold rasburicase on the basis of G6PD status", "Strong",
     "PA166119846")
_add("G6PD", "rasburicase", "G6PD Deficient",
     "rasburicase: contraindicated - use an alternative such as allopurinol", "Strong",
     "PA166119846")
_add("G6PD", "rasburicase", "G6PD Variable",
     "rasburicase: measure G6PD enzyme activity to confirm normal status before use; "
     "alternatives include allopurinol", "Moderate", "PA166119846")

# --- HLA ------------------------------------------------------------------
_add("HLA-A*31:01", "carbamazepine",
     "HLA-A*31:01 positive (increased risk of carbamazepine-induced SJS/TEN, DRESS and MPE)",
     "carbamazepine: if the patient is carbamazepine-naive and alternatives are available, do "
     "not use carbamazepine; if no alternative exists, use with increased clinical monitoring "
     "and stop at the first sign of a cutaneous reaction", "Strong", "PA166105008")
_add("HLA-A*31:01", "carbamazepine",
     "HLA-A*31:01 negative (normal or reduced risk)",
     "carbamazepine: use per standard dosing guidelines (an HLA-A*31:01 negative result does "
     "not eliminate risk)", "Strong", "PA166105008")
_add("HLA-B*15:02", "carbamazepine",
     "HLA-B*15:02 positive (increased risk of SJS/TEN)",
     "carbamazepine: do not use in a carbamazepine-naive patient; avoid other aromatic "
     "anticonvulsants (phenytoin, oxcarbazepine) as alternatives", "Strong", "PA166105008")
_add("HLA-B*15:02", "carbamazepine",
     "HLA-B*15:02 negative (normal or reduced risk)",
     "carbamazepine: use per standard dosing guidelines (a negative result does not eliminate "
     "risk)", "Strong", "PA166105008")
_add("HLA-B*15:02", "oxcarbazepine",
     "HLA-B*15:02 positive (increased risk of SJS/TEN)",
     "oxcarbazepine: do not use in an oxcarbazepine-naive patient", "Strong", "PA166176623")
_add("HLA-B*15:02", "oxcarbazepine",
     "HLA-B*15:02 negative (normal or reduced risk)",
     "oxcarbazepine: use per standard dosing guidelines", "Strong", "PA166176623")
_add("HLA-B*57:01", "abacavir",
     "HLA-B*57:01 positive (high risk of abacavir hypersensitivity)",
     "abacavir: not recommended - use an alternative antiretroviral (significantly increased "
     "risk of hypersensitivity reaction)", "Strong", "PA166104997")
_add("HLA-B*57:01", "abacavir",
     "HLA-B*57:01 negative (very low risk of abacavir hypersensitivity)",
     "abacavir: use per standard dosing guidelines", "Strong", "PA166104997")
_add("HLA-B*58:01", "allopurinol",
     "HLA-B*58:01 positive (significantly increased risk of allopurinol-induced SCAR)",
     "allopurinol: contraindicated - use an alternative urate-lowering therapy", "Strong",
     "PA166105003")
_add("HLA-B*58:01", "allopurinol",
     "HLA-B*58:01 negative (low or reduced risk of allopurinol SCAR)",
     "allopurinol: use per standard dosing guidelines", "Strong", "PA166105003")

# --- IFNL3 / MT-RNR1 / RYR1 / CFTR ---------------------------------------
_add("IFNL3", "PEG-IFN-α", "IFNL3 favourable response genotype",
     "PEG-IFN-alpha (with ribavirin): favourable response genotype - approximately 70% chance "
     "of SVR after 48 weeks (about 90% with a protease inhibitor); weighs in favour of using a "
     "PEG-IFN-alpha/ribavirin regimen", "Strong", "PA166110235")
_add("IFNL3", "PEG-IFN-α", "IFNL3 unfavourable response genotype",
     "PEG-IFN-alpha (with ribavirin): unfavourable response genotype - approximately 30% "
     "chance of SVR after 48 weeks; consider the implications before initiating a "
     "PEG-IFN-alpha/ribavirin regimen", "Strong", "PA166110235")
_add("MT-RNR1", "aminoglycosides",
     "MT-RNR1 increased risk of aminoglycoside-induced hearing loss",
     "aminoglycosides: avoid unless the high risk of permanent hearing loss is outweighed by "
     "the severity of infection and the lack of a safe or effective alternative; if used, use "
     "the lowest dose and shortest duration with therapeutic drug monitoring and frequent "
     "audiometry", "Strong", "PA166229081")
_add("MT-RNR1", "aminoglycosides",
     "MT-RNR1 normal risk of aminoglycoside-induced hearing loss",
     "aminoglycosides: use standard doses for the shortest feasible course with therapeutic "
     "drug monitoring and regular hearing evaluation", "Strong", "PA166229081")
_add("RYR1", "volatile-anaesthetics", "Malignant Hyperthermia susceptible",
     "volatile anaesthetics (and succinylcholine): relatively contraindicated - use a "
     "non-triggering anaesthetic technique; do not use except in extraordinary circumstances "
     "where the benefits outweigh the risks", "Strong", "PA166303941")
_add("RYR1", "volatile-anaesthetics", "Uncertain susceptibility",
     "volatile anaesthetics (and succinylcholine): a negative result does not exclude "
     "malignant hyperthermia susceptibility - let clinical findings, family history and "
     "further testing guide use", "Strong", "PA166303941")
_add("CFTR", "ivacaftor", "CFTR gating variant present (ivacaftor-responsive)",
     "ivacaftor: use according to the product label", "Strong", "PA166114461")
_add("CFTR", "ivacaftor", "CFTR F508del homozygous (not ivacaftor-responsive)",
     "ivacaftor: not recommended - F508del homozygotes are unlikely to respond to ivacaftor "
     "monotherapy", "Moderate", "PA166114461")

RECOMMENDATIONS = R

# ---------------------------------------------------------------------------
# 4. Lethal-class designation
# ---------------------------------------------------------------------------
# A case is lethal-class when standard (genotype-naive) dosing of the named
# drug in a patient with that phenotype carries a documented risk of death or
# of a life-threatening event, as opposed to reduced efficacy or reversible
# toxicity. `basis` records the corpus wording relied on, or flags a judgement
# that goes beyond the corpus wording.
LETHAL_CLASS = {
    ("DPYD", "DPYD Poor Metabolizer"):
        "corpus: 'Complete DPD deficiency and increased risk for severe or even fatal drug "
        "toxicity'; full-dose fluoropyrimidine in a DPD-null patient is a recognised cause of "
        "treatment-related death",
    ("DPYD", "DPYD Intermediate Metabolizer"):
        "corpus: 'increased risk for severe or even fatal drug toxicity when treated with "
        "fluoropyrimidine drugs'",
    ("TPMT", "TPMT Poor Metabolizer"):
        "standard-dose thiopurine in a TPMT-null patient causes profound, potentially fatal "
        "myelosuppression (judgement: the corpus annotation states only that an alternative "
        "immunosuppressant should be considered)",
    ("NUDT15", "NUDT15 Poor Metabolizer"):
        "standard-dose thiopurine in a NUDT15-null patient causes profound, potentially fatal "
        "myelosuppression; corpus notes tolerated dose is ~10% of normal (judgement)",
    ("CYP2D6", "CYP2D6 Ultrarapid Metabolizer"):
        "codeine only: corpus 'Increased formation of morphine leading to higher risk of "
        "toxicity' - the basis of the paediatric respiratory-depression deaths; not lethal for "
        "tamoxifen/ondansetron/paroxetine/amitriptyline, where the risk is loss of efficacy",
    ("HLA-B*15:02", "HLA-B*15:02 positive (increased risk of SJS/TEN)"):
        "corpus: increased risk of Stevens-Johnson syndrome and toxic epidermal necrolysis, "
        "both of which carry substantial mortality",
    ("HLA-A*31:01",
     "HLA-A*31:01 positive (increased risk of carbamazepine-induced SJS/TEN, DRESS and MPE)"):
        "corpus: increased risk of SJS/TEN and DRESS",
    ("HLA-B*57:01", "HLA-B*57:01 positive (high risk of abacavir hypersensitivity)"):
        "corpus: significantly increased risk of abacavir hypersensitivity, a systemic "
        "reaction that is fatal on rechallenge",
    ("HLA-B*58:01",
     "HLA-B*58:01 positive (significantly increased risk of allopurinol-induced SCAR)"):
        "corpus: significantly increased risk of severe cutaneous adverse reaction "
        "(SJS/TEN/DRESS)",
    ("RYR1", "Malignant Hyperthermia susceptible"):
        "corpus: relatively contraindicated; malignant hyperthermia is an acute, "
        "rapidly fatal intraoperative event if untreated",
    ("G6PD", "G6PD Deficient"):
        "corpus: 'High risk of acute hemolytic anemia'; rasburicase is contraindicated and can "
        "precipitate massive haemolysis and methaemoglobinaemia",
    ("UGT1A1", "UGT1A1 Poor Metabolizer"):
        "irinotecan only: corpus DPWG wording 'Serious, life-threatening adverse events occur "
        "more often in patients with this genetic variation'; not lethal for atazanavir, where "
        "the issue is jaundice and non-adherence",
    ("CYP2C19", "CYP2C19 Poor Metabolizer"):
        "clopidogrel only: corpus 'increased risk for adverse cardiac and cerebrovascular "
        "events' - stent thrombosis and stroke are frequently fatal (judgement: the corpus "
        "does not use the word fatal)",
    ("CYP2C19", "CYP2C19 Intermediate Metabolizer"):
        "clopidogrel only: same reasoning as poor metabolizer, lower magnitude (judgement)",
}
# phenotype/drug pairs where the lethal-class key above must be restricted
LETHAL_DRUG_RESTRICTION = {
    ("CYP2D6", "CYP2D6 Ultrarapid Metabolizer"): {"codeine"},
    ("UGT1A1", "UGT1A1 Poor Metabolizer"): {"irinotecan"},
    ("CYP2C19", "CYP2C19 Poor Metabolizer"): {"clopidogrel"},
    ("CYP2C19", "CYP2C19 Intermediate Metabolizer"): {"clopidogrel"},
}
