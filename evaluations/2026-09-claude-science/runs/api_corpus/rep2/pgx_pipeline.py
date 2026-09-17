#!/usr/bin/env python3
"""
Blinded pharmacogenomics study - genotype -> diplotype -> CPIC phenotype -> CPIC recommendation.

Inputs  (same directory):
    cases_blinded.json            110 blinded cases (id, gene, drug, genotype)
    cpic_guideline_corpus.json    CPIC/PharmGKB guideline text snapshot 2026-05-09
    cpic_api_snapshot/*.json      CPIC data API (api.cpicpgx.org/v1) tables, fetched by fetch_cpic_api.py

Outputs (same directory):
    claude_science_answers.csv        case_id,diplotype,phenotype,recommendation,abstained
    table_diplotype_to_phenotype.csv  diplotype -> phenotype lookup table actually used
    table_phenotype_to_recommendation.csv  phenotype/lookupkey -> CPIC recommendation table used
    table_variant_to_allele.csv       rsID/variant -> star allele definition table used
    case_audit.csv                    full per-case audit trail (evidence, source, lethal class, notes)

Method summary
--------------
1. Each case genotype string is parsed into a structured observation:
   typed positions (rsID -> two observed bases, or one for X/mito hemizygosity), or a
   non-SNV observation (HLA allele status, TA-repeat length, copy-number, mito variant,
   "all reference").
2. For star-allele genes the diplotype is CALLED, not looked up by hand: candidate
   haplotype pairs are enumerated from the CPIC allele definition tables
   (allele_definition x allele_location_value x sequence_location) under the panel
   assumption that untyped positions are reference. A pair is retained only if every
   locus (typed loci, plus any locus where an allele mandates a non-reference base)
   can be reconciled with the observed genotype. Strand is checked against the CPIC
   plus-strand representation and complemented if, and only if, the observed bases do
   not occur at that locus but their complements do (logged per case).
3. Phenotype comes from the CPIC `diplotype` table (diplotype -> generesult + lookupkey).
4. Recommendation comes from the CPIC `recommendation` table, matched on drug +
   gene lookupkey, with companion genes set to CPIC's own "No Result" key where the case
   provides only one gene. Where CPIC publishes no machine-readable recommendation
   (warfarin, ondansetron, peginterferon alfa) the supplied guideline corpus is used.
5. Abstention: any case where the call is not unique at the phenotype level, where the
   input cannot be mapped to a CPIC allele/diplotype, or where no CPIC recommendation
   exists for the gene-drug pair.
"""
import json, os, re, csv, itertools, collections, hashlib, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
SNAP = os.path.join(HERE, "cpic_api_snapshot")

# ---------------------------------------------------------------- load inputs
cases = json.load(open(os.path.join(HERE, "cases_blinded.json")))
corpus = json.load(open(os.path.join(HERE, "cpic_guideline_corpus.json")))
T = {t: json.load(open(os.path.join(SNAP, t + ".json")))
     for t in ["gene", "pair", "drug", "guideline", "allele", "allele_definition",
               "allele_location_value", "sequence_location", "gene_result",
               "gene_result_lookup", "recommendation", "term"]}
DIP = json.load(open(os.path.join(SNAP, "diplotype_by_gene.json")))
API_META = json.load(open(os.path.join(SNAP, "_fetch_meta.json")))

DRUGID2NAME = {d["drugid"]: d["name"] for d in T["drug"]}
NAME2DRUGID = {}
for d in T["drug"]:
    NAME2DRUGID.setdefault(d["name"].lower(), d["drugid"])

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPL = {"A": "T", "T": "A", "C": "G", "G": "C"}

# ------------------------------------------------- allele definition structures
AD = {a["id"]: a for a in T["allele_definition"]}
SL = {s["id"]: s for s in T["sequence_location"]}

def _bases(val):
    """Expand a CPIC variantallele string into a set of allowed values.

    Single IUPAC letters expand to their base set; 'X or Y' alternatives are split;
    indels / repeat tokens are kept as opaque strings and compared literally.
    """
    if val is None:
        return None
    out = set()
    for part in re.split(r"\s+or\s+", val.strip()):
        p = part.strip()
        if len(p) == 1 and p in IUPAC:
            out |= set(IUPAC[p])
        elif p:
            out.add(p)
    return out or None

ALLELE_DEF = collections.defaultdict(dict)      # gene -> allele -> {locid: set(bases)}
for v in T["allele_location_value"]:
    d = AD.get(v["alleledefinitionid"]); s = SL.get(v["locationid"])
    if not d or not s:
        continue
    ALLELE_DEF[d["genesymbol"]].setdefault(d["name"], {})[v["locationid"]] = _bases(v["variantallele"])

RSID2LOC = collections.defaultdict(dict)        # gene -> rsid -> locid
for s in T["sequence_location"]:
    if s.get("dbsnpid"):
        RSID2LOC[s["genesymbol"]][s["dbsnpid"]] = s["id"]

REF_ALLELE = {}                                 # gene -> reference allele name
for gene, alls in ALLELE_DEF.items():
    for cand in ("*1", "Reference", "B (reference)", "rs9923231 reference (C)"):
        if cand in alls:
            REF_ALLELE[gene] = cand
            break

_HGVS_SNV = re.compile(r"[gcmn]\.[\d+\-_*]+([ACGT])>([ACGT])")

def _ref_from_location(loc):
    """Reference base parsed from the CPIC genomic HGVS string of the locus.

    Needed because some CPIC reference alleles (e.g. CYP2C19 *1) carry no positional
    rows at all, so the reference base cannot be read off the reference haplotype.
    """
    s = SL.get(loc) or {}
    for field in ("chromosomelocation", "genelocation", "name"):
        txt = s.get(field) or ""
        m = _HGVS_SNV.search(txt)
        if m:
            return {m.group(1)}
    return None

_REF_CACHE = {}

def ref_base(gene, loc):
    key = (gene, loc)
    if key in _REF_CACHE:
        return _REF_CACHE[key]
    ra = REF_ALLELE.get(gene)
    val = None
    if ra and loc in ALLELE_DEF[gene].get(ra, {}):
        val = ALLELE_DEF[gene][ra][loc]
    if val is None:
        val = _ref_from_location(loc)
    _REF_CACHE[key] = val
    return val

def dip_vocabulary(gene):
    """Allele names that appear in the CPIC diplotype table for this gene."""
    voc = set()
    for row in DIP.get(gene, []):
        if "/" in row["diplotype"]:
            a, b = row["diplotype"].split("/", 1)
            voc.add(a); voc.add(b)
        else:
            voc.add(row["diplotype"])
    return voc

DIP_INDEX = {}                                  # gene -> {frozenset-ish key: row}
for gene, rows in DIP.items():
    idx = {}
    for r in rows:
        d = r["diplotype"]
        idx[d] = r
        if "/" in d:
            a, b = d.split("/", 1)
            idx[b + "/" + a] = r
    DIP_INDEX[gene] = idx

# --------------------------------------------------------- the diplotype caller
def mandatory_variant_loci(gene, allele):
    """Loci where this allele cannot carry the reference base."""
    out = set()
    for loc, bs in ALLELE_DEF[gene][allele].items():
        rb = ref_base(gene, loc)
        if rb is None:
            continue
        if bs and not (bs & rb):
            out.add(loc)
    return out

def call_diplotypes(gene, observed, ploidy=2, tag_snp=False):
    """
    observed: {locid: [base, ...]} for typed loci (1 base = one copy, i.e. hemizygous).

    Strict mode (tag_snp=False): untyped loci are assumed homozygous reference - the
    semantics of a targeted SNV panel that reports only the variants it found.

    Tag-SNP mode (tag_snp=True): untyped loci are treated as unknown, and an allele may
    be called if it is ANCHORED by an observed variant locus (its other defining
    variants being untyped). This is the linkage assumption needed for haplotypes that
    are assayed through a single marker variant (e.g. DPYD HapB3 via rs56038477).

    Returns sorted list of diplotype strings present in the CPIC diplotype table.
    """
    alls = ALLELE_DEF.get(gene, {})
    voc = dip_vocabulary(gene)
    var_loci = set()
    for loc, obs in observed.items():
        rb = ref_base(gene, loc)
        if rb is None or any(b not in rb for b in obs):
            var_loci.add(loc)
    cands = []
    for name in alls:
        if name not in voc:
            continue
        must = mandatory_variant_loci(gene, name)
        if must <= var_loci:
            cands.append(name)
        elif tag_snp and (must & var_loci):
            cands.append(name)

    def ok(loc, a1, a2):
        s1 = alls[a1].get(loc) or ref_base(gene, loc)
        s2 = alls[a2].get(loc) or ref_base(gene, loc)
        if s1 is None or s2 is None:
            return True
        obs = observed.get(loc)
        if obs is None:
            if tag_snp:
                return True                   # untyped -> unknown
            rb = ref_base(gene, loc)          # untyped -> reference / reference
            return bool(s1 & rb) and bool(s2 & rb)
        if len(obs) == 1:                     # single copy typed (hemizygous)
            return bool(s1 & set(obs))
        x, y = obs
        return (x in s1 and y in s2) or (y in s1 and x in s2)

    hits = []
    if ploidy == 1:
        for a in cands:
            loci = set(observed) | mandatory_variant_loci(gene, a)
            if all(ok(l, a, a) for l in loci) and a in DIP_INDEX[gene]:
                hits.append(a)
        return sorted(set(hits))
    for a1, a2 in itertools.combinations_with_replacement(sorted(cands), 2):
        loci = set(observed) | mandatory_variant_loci(gene, a1) | mandatory_variant_loci(gene, a2)
        if all(ok(l, a1, a2) for l in loci):
            for d in (a1 + "/" + a2, a2 + "/" + a1):
                if d in DIP_INDEX[gene]:
                    hits.append(DIP_INDEX[gene][d]["diplotype"])
                    break
    return sorted(set(hits))

def primary_diplotype(gene, dips):
    """Canonical representative of a set of equivalent candidate diplotypes: the one
    built from the fewest non-reference alleles, then the shortest/lowest-numbered."""
    ref = REF_ALLELE.get(gene)
    def rank(d):
        alleles = d.split("/")
        nonref = sum(1 for a in alleles if a != ref)
        nums = []
        for a in alleles:
            m = re.match(r"\*(\d+)", a)
            nums.append(int(m.group(1)) if m else 999)
        return (nonref, sorted(nums), len(d), d)
    return sorted(dips, key=rank)[0]

# ------------------------------------------------------------- genotype parsing
GT_RE = re.compile(r"(rs\d+)\s+([ACGT])(?:\s*/\s*([ACGT]))?(?![\w/])")
TA_RE = re.compile(r"TA(\d)\s*/\s*TA(\d)")

def parse_snv_observations(gene, gtext):
    """
    Returns (observed {locid: [bases]}, notes, unresolved[list of rsids]).
    Strand is complemented only when the observed bases are absent from the CPIC
    representation at that locus but their complements are present.
    """
    observed, notes, unresolved = {}, [], []
    for m in GT_RE.finditer(gtext):
        rsid, a1, a2 = m.group(1), m.group(2), m.group(3)
        obs = [a1] + ([a2] if a2 else [])
        loc = RSID2LOC.get(gene, {}).get(rsid)
        if loc is None:
            unresolved.append(rsid)
            continue
        known = set()
        for name, dd in ALLELE_DEF[gene].items():
            if loc in dd and dd[loc]:
                known |= dd[loc]
        if not set(obs) <= known:
            comp = [COMPL[b] for b in obs]
            if set(comp) <= known:
                notes.append(f"{rsid}: observed {'/'.join(obs)} complemented to {'/'.join(comp)} "
                             f"to match CPIC plus-strand representation")
                obs = comp
            else:
                unresolved.append(f"{rsid}(alleles {'/'.join(obs)} not in CPIC definition)")
                continue
        prev = observed.get(loc)
        observed[loc] = obs if prev is None else prev
    return observed, notes, unresolved

# ----------------------------------------------------- recommendation retrieval
REC = T["recommendation"]
POP_CHOICE = {          # population subsets that need an explicit clinical assumption
    "voriconazole": "adults",
    "phenytoin": "PHT naive",
    "carbamazepine": "CBZ naive",
    "oxcarbazepine": "OXC naive",
    "clopidogrel": "CVI ACS PCI ",
}
CLASS_DRUGS = {         # study drug label -> CPIC drug names in that class
    "thiopurines": ["azathioprine", "mercaptopurine", "thioguanine"],
    "aminoglycosides": ["amikacin", "gentamicin", "kanamycin", "plazomicin",
                        "streptomycin", "tobramycin"],
    "volatile-anaesthetics": ["desflurane", "enflurane", "halothane", "isoflurane",
                              "methoxyflurane", "sevoflurane"],
    "PEG-IFN-α": ["peginterferon alfa-2a", "peginterferon alfa-2b"],
}

def get_recommendation(drug_label, gene_key, lookup_value):
    """
    Match CPIC recommendation rows on drug + {gene: lookup_value}; companion genes in a
    multi-gene key are set to CPIC's "No Result".
    Returns (rows, cpic_drug_names_used, population_used, homogeneous_text_bool).
    """
    names = CLASS_DRUGS.get(drug_label, [drug_label.lower()])
    out, used = [], []
    for n in names:
        did = NAME2DRUGID.get(n)
        if not did:
            continue
        rows = [r for r in REC if r["drugid"] == did]
        keep = []
        for r in rows:
            lk = r.get("lookupkey") or {}
            if lk.get(gene_key) != lookup_value:
                continue
            if any(v != "No Result" for k, v in lk.items() if k != gene_key):
                continue
            keep.append(r)
        if keep:
            used.append(n)
            out += keep
    pop = POP_CHOICE.get(drug_label.lower())
    if pop:
        sel = [r for r in out if r["population"].strip() == pop.strip()]
        out = sel or out
    pops = sorted(set(r["population"] for r in out))
    texts = sorted(set(r["drugrecommendation"].strip() for r in out))
    return out, used, pops, len(texts) <= 1

def squeeze(text, maxlen=600):
    t = " ".join(text.split())
    if len(t) <= maxlen:
        return t
    cut = t[:maxlen]
    sp = cut.rfind(". ")
    return (cut[:sp + 1] if sp > 60 else cut.rstrip() + " ...")

# ------------------------------------------- corpus-sourced recommendations (no API rows)
CORPUS_REC = {}
for gene, entry in corpus["genes"].items():
    CORPUS_REC[gene] = entry
WARFARIN_TABLE = {          # from CPIC warfarin guideline table in the supplied corpus
    ("*1/*1", "GG"): "5-7", ("*1/*1", "GA"): "5-7", ("*1/*1", "AA"): "3-4",
    ("*1/*3", "GG"): "3-4", ("*1/*3", "GA"): "3-4", ("*1/*3", "AA"): "0.5-2",
    ("*3/*3", "GG"): "0.5-2", ("*3/*3", "GA"): "0.5-2", ("*3/*3", "AA"): "0.5-2",
}
# CPIC ondansetron/CYP2D6 recommendations - guideline table in the supplied corpus
# (the CPIC recommendation table in the API has no machine-readable ondansetron rows).
ONDANSETRON_REC = {
    "Ultrarapid Metabolizer": "Select an alternative drug not predominantly metabolized by CYP2D6 "
                              "(e.g. granisetron).",
    "Normal Metabolizer": "Initiate therapy with recommended starting dose.",
    "Intermediate Metabolizer": "Insufficient evidence demonstrating clinical impact based on CYP2D6 "
                                "genotype. Initiate therapy with recommended starting dose.",
    "Poor Metabolizer": "Insufficient evidence demonstrating clinical impact based on CYP2D6 "
                        "genotype. Initiate therapy with recommended starting dose.",
}
# CPIC IFNL3/peginterferon alfa guideline table in the supplied corpus.
IFNL3_REC = {
    "Favorable response genotype":
        "Approximately 70% chance of SVR after 48 weeks of PEG-IFN-alpha + RBV treatment "
        "(~90% with a protease-inhibitor combination); weighs in favour of using PEG-IFN-alpha "
        "containing regimens. Consider implications before initiating therapy.",
    "Unfavorable response genotype":
        "Approximately 30% chance of SVR after 48 weeks of PEG-IFN-alpha + RBV treatment "
        "(~60% with a protease-inhibitor combination); consider implications, including "
        "non-interferon-based alternatives, before initiating therapy.",
}

# ------------------------------------------------------------------ lethal class
# "Lethal class" = standard/label dosing in this genotype is expected to cause severe
# (grade >=3), irreversible or life-threatening harm, or a life-threatening therapeutic
# failure. Encoded as (gene, drug regex, phenotype/lookup regex, mechanism).
LETHAL_RULES = [
    ("CYP2D6",  r"codeine",              r"Ultrarapid",            "morphine overdose: fatal respiratory depression"),
    ("CYP2D6",  r"amitriptyline",        r"Poor Metabolizer",      "tricyclic accumulation: seizures, cardiac arrhythmia"),
    ("CYP2C19", r"clopidogrel",          r"Poor Metabolizer|Intermediate Metabolizer",
                                                                   "loss of antiplatelet effect: stent thrombosis, MI, death (ACS/PCI)"),
    ("CYP2C19", r"voriconazole",         r"Poor Metabolizer",      "voriconazole accumulation: hepatotoxicity, neurotoxicity"),
    ("CYP2C19", r"voriconazole",         r"Ultrarapid|Rapid",      "subtherapeutic exposure: failure of therapy for invasive fungal infection"),
    ("CYP2C9",  r"warfarin",             r"^0\.0$|Poor Metabolizer", "over-anticoagulation: major/fatal haemorrhage"),
    ("CYP2C9",  r"phenytoin",            r"^0\.0$|Poor Metabolizer", "phenytoin accumulation: severe neurotoxicity"),
    ("DPYD",    r"fluorouracil|capecitabine", r"^0\.0$|Poor Metabolizer",
                                                                   "fluoropyrimidine overdose: fatal myelosuppression, mucositis, neurotoxicity"),
    ("DPYD",    r"fluorouracil|capecitabine", r"^0\.5$|^1\.0$|^1\.5$|Intermediate Metabolizer",
                                                                   "partial DPD deficiency: severe (grade >=3), sometimes fatal fluoropyrimidine toxicity at full dose"),
    ("TPMT",    r"azathioprine|thiopurine|mercaptopurine", r"Poor Metabolizer",
                                                                   "thiopurine accumulation: life-threatening myelosuppression"),
    ("NUDT15",  r"azathioprine|thiopurine|mercaptopurine", r"Poor Metabolizer",
                                                                   "thiopurine accumulation: life-threatening myelosuppression"),
    ("SLCO1B1", r"simvastatin",          r"Poor Function",         "statin myopathy / rhabdomyolysis"),
    ("UGT1A1",  r"irinotecan",           r"Poor Metabolizer",      "severe neutropenia and diarrhoea"),
    ("G6PD",    r"rasburicase",          r"Deficient",             "acute haemolysis and methaemoglobinaemia; rasburicase is contraindicated in G6PD deficiency"),
    ("G6PD",    r"primaquine",           r"Deficient",             "acute haemolytic anaemia at standard primaquine dosing"),
    ("G6PD",    r"rasburicase|primaquine", r"Variable",            "possible acute haemolysis - G6PD activity must be measured before use"),
    ("HLA-B",   r"abacavir",             r"\*57:01 positive",      "abacavir hypersensitivity syndrome, fatal on re-challenge"),
    ("HLA-B",   r"carbamazepine|oxcarbazepine", r"\*15:02 positive", "Stevens-Johnson syndrome / toxic epidermal necrolysis"),
    ("HLA-A",   r"carbamazepine",        r"\*31:01 positive",      "DRESS, SJS/TEN, maculopapular exanthema"),
    ("HLA-B",   r"allopurinol",          r"\*58:01 positive",      "allopurinol-induced severe cutaneous adverse reaction (SJS/TEN, DRESS)"),
    ("RYR1",    r"anaesthetic|flurane",  r"Malignant Hyperthermia", "malignant hyperthermia crisis"),
    ("MT-RNR1", r"aminoglycoside|micin|mycin", r"increased risk",  "permanent bilateral sensorineural hearing loss"),
]

def lethal_class(gene_key, drug, *pheno_or_key):
    """Match against every available phenotype descriptor (CPIC phenotype term and the
    lookup key, which is an activity score for CYP2D6/CYP2C9/DPYD)."""
    fields = [str(x) for x in pheno_or_key if x]
    for g, drx, prx, mech in LETHAL_RULES:
        if g != gene_key or not re.search(drx, drug, re.I):
            continue
        if any(re.search(prx, f, re.I) for f in fields):
            return True, mech
    return False, ""

# ------------------------------------------------------------------- resolution
GENE_KEY = {"HLA-B*57:01": "HLA-B", "HLA-B*15:02": "HLA-B", "HLA-B*58:01": "HLA-B",
            "HLA-A*31:01": "HLA-A"}
CNV_CAVEAT = ("SNV panel cannot detect gene deletion (*5) or duplication; untyped positions "
              "assumed reference per CPIC panel convention")

def phenotype_row(gene, diplotype):
    return DIP_INDEX.get(gene, {}).get(diplotype)

def resolve(case):
    gene, drug, gt = case["gene"], case["drug"], case["genotype"]
    gk = GENE_KEY.get(gene, gene)
    r = {"case_id": case["id"], "gene": gene, "gene_key": gk, "drug": drug, "genotype": gt,
         "diplotype": "", "phenotype": "", "lookup_value": "", "recommendation": "",
         "abstained": False, "abstain_reason": "", "rec_source": "", "pheno_source": "",
         "candidates": "", "notes": [], "cpic_level": "", "rec_classification": "",
         "implications": "", "rec_population": ""}
    p = T["pair"]
    lvl = [x for x in p if x["genesymbol"] == gk
           and DRUGID2NAME.get(x["drugid"], "").lower() in
           [d.lower() for d in CLASS_DRUGS.get(drug, [drug])]]
    r["cpic_level"] = "/".join(sorted({x["cpiclevel"] for x in lvl})) if lvl else "not listed"

    # ---- 1. genotype -> diplotype + phenotype
    if gk in ("HLA-B", "HLA-A"):
        allele = gene.split("-", 1)[1].split("*", 1)[1]
        allele = "*" + allele
        if re.search(r"negative|no risk allele|T/T", gt, re.I):
            status = f"{allele} negative"
        elif re.search(r"positive|any positive|G/T|G/G", gt, re.I):
            status = f"{allele} positive"
        else:
            r.update(abstained=True, abstain_reason="HLA allele status not interpretable")
            return r
        if "rs2395029" in gt:
            r["notes"].append("HLA-B*57:01 status inferred from the HCP5 proxy SNP rs2395029 "
                              "(strong linkage but not direct HLA typing; imperfect concordance)")
        r["diplotype"] = f"HLA-{gene.split('*')[1] if False else allele.lstrip('*')} " \
                         f"{'carrier' if 'positive' in status else 'non-carrier'}"
        r["diplotype"] = f"{gk}{allele} {'positive (carrier)' if 'positive' in status else 'negative'}"
        r["phenotype"] = f"{gk}{status}"
        r["lookup_value"] = status
        r["pheno_source"] = "CPIC allele status (recommendation lookupkey vocabulary)"

    elif gene == "MT-RNR1":
        if re.search(r"m\.1555A>G", gt):
            r["diplotype"] = "m.1555A>G (homoplasmic)"
            r["lookup_value"] = "increased risk of aminoglycoside-induced hearing loss"
            r["pheno_source"] = "CPIC diplotype table (m.1555A>G)"
        elif re.search(r"m\.1555A\b|no variant", gt):
            r["diplotype"] = "reference (no MT-RNR1 increased-risk variant detected)"
            r["lookup_value"] = "normal risk of aminoglycoside-induced hearing loss"
            r["pheno_source"] = ("CPIC guideline table: individuals with no detectable MT-RNR1 "
                                 "increased-risk variant are normal risk")
        else:
            r.update(abstained=True, abstain_reason="MT-RNR1 variant not interpretable")
            return r
        r["phenotype"] = "MT-RNR1 " + r["lookup_value"]
        r["notes"].append("only the m.1555 position is reported; other MT-RNR1 risk variants "
                          "(m.1095T>C, m.1494C>T) not assessed")

    elif gene == "RYR1":
        if re.search(r"no pathogenic|no variant", gt, re.I):
            dip = "Reference/Reference"
        elif re.search(r"c\.1840C>T|Arg614Cys", gt):
            dip = "c.1840C>T/Reference"
        else:
            r.update(abstained=True, abstain_reason="RYR1 variant not identified")
            return r
        row = phenotype_row("RYR1", dip)
        r["diplotype"] = dip
        r["phenotype"] = row["generesult"]
        r["lookup_value"] = row["lookupkey"]["RYR1"]
        r["pheno_source"] = "CPIC diplotype table (RYR1)"
        if dip == "Reference/Reference":
            r["notes"].append("CPIC term is 'Uncertain Susceptibility': a negative RYR1/CACNA1S "
                              "result does not exclude malignant hyperthermia susceptibility")

    elif gene == "CFTR":
        if re.search(r"G551D", gt):
            dip = "G551D/ivacaftor non-responsive CFTR sequence"
            r["notes"].append("case states 'G/A or A/A'; heterozygous G551D is the conservative "
                              "call and gives the same CPIC result as G551D/G551D")
        elif re.search(r"F508del", gt):
            dip = "ivacaftor non-responsive CFTR sequence/ivacaftor non-responsive CFTR sequence"
            r["notes"].append("F508del is not in the CPIC ivacaftor-responsive variant list; "
                              "it maps to 'ivacaftor non-responsive CFTR sequence'")
        else:
            r.update(abstained=True, abstain_reason="CFTR variant not interpretable")
            return r
        row = phenotype_row("CFTR", dip)
        r["diplotype"] = dip
        r["phenotype"] = row["generesult"]
        r["lookup_value"] = row["lookupkey"]["CFTR"]
        r["pheno_source"] = "CPIC diplotype table (CFTR)"

    elif gene == "IFNL3":
        if re.search(r"C/C", gt):
            r["diplotype"] = "rs12979860 reference (C)/rs12979860 reference (C)"
            r["phenotype"] = "Favorable response genotype"
        elif re.search(r"T/T", gt):
            r["diplotype"] = "rs12979860 variant (T)/rs12979860 variant (T)"
            r["phenotype"] = "Unfavorable response genotype"
        else:
            r.update(abstained=True, abstain_reason="IFNL3 genotype not interpretable")
            return r
        r["lookup_value"] = r["phenotype"]
        r["pheno_source"] = "CPIC IFNL3/peginterferon guideline table (supplied corpus)"

    elif gene in ("VKORC1", "CYP4F2"):
        obs, notes, unres = parse_snv_observations(gene, gt)
        r["notes"] += notes
        dips = call_diplotypes(gene, obs) if obs else []
        if gene == "VKORC1":
            var = "variant (T)" if re.search(r"A/A", gt) else "reference (C)"
            r["diplotype"] = f"rs9923231 {var}/rs9923231 {var}"
            r["phenotype"] = ("n/a - CPIC assigns no phenotype term to VKORC1; rs9923231 "
                              + ("A/A (-1639 AA) = low warfarin dose requirement"
                                 if var.startswith("variant") else
                                 "G/G (-1639 GG) = reference dose requirement"))
        else:
            allele = "*3" if re.search(r"T/T", gt) else "*1"
            r["diplotype"] = f"{allele}/{allele}"
            r["phenotype"] = ("n/a - CPIC assigns no phenotype term to CYP4F2; "
                              f"rs2108622 {'T/T (*3/*3)' if allele=='*3' else 'C/C (*1/*1)'}")
        r["lookup_value"] = ""
        r["pheno_source"] = ("CPIC allele definitions; CPIC publishes no phenotype table for "
                             "this gene (no rows in gene_result / diplotype)")
        r["candidates"] = ";".join(dips)

    elif gene == "UGT1A1":
        m = TA_RE.search(gt)
        if not m:
            r.update(abstained=True, abstain_reason="UGT1A1 repeat genotype not interpretable")
            return r
        # conventional (TA)n nomenclature -> CPIC allele: (TA)6 = *1, (TA)7 = *28
        conv = {"5": "*36", "6": "*1", "7": "*28", "8": "*37"}
        a1, a2 = conv.get(m.group(1)), conv.get(m.group(2))
        if not a1 or not a2:
            r.update(abstained=True, abstain_reason="unrecognised UGT1A1 TA-repeat length")
            return r
        dip = f"{a1}/{a2}"
        row = phenotype_row("UGT1A1", dip)
        r["diplotype"] = dip
        r["phenotype"] = row["generesult"]
        r["lookup_value"] = row["lookupkey"]["UGT1A1"]
        r["pheno_source"] = "CPIC diplotype table (UGT1A1)"
        r["notes"].append("rs8175347 (merged into rs3064744) TA-repeat length mapped to CPIC "
                          "alleles: (TA)6=*1, (TA)7=*28; the *80 (rs887829) and *6/*27 coding "
                          "variants were not typed")

    elif gene == "CYP2D6" and re.search(r"\*1xN|duplication", gt, re.I):
        cands = ["*1/*1x2", "*1/*1x\u22653"]
        rows = [phenotype_row("CYP2D6", c) for c in cands]
        phenos = {x["generesult"] for x in rows if x}
        r["diplotype"] = "*1/*1xN (N>=2)"
        r["candidates"] = ";".join(cands)
        if len(phenos) != 1:
            r.update(abstained=True, abstain_reason="copy number N unknown and phenotype differs by N")
            return r
        r["phenotype"] = phenos.pop()
        r["lookup_value"] = rows[0]["lookupkey"]["CYP2D6"]      # activity score for N = 2
        r["pheno_source"] = "CPIC diplotype table (CYP2D6 *1/*1x2 and *1/*1x>=3)"
        r["notes"].append("copy-number assay identifies a duplication but not which allele is "
                          "duplicated; with no SNVs reported both alleles are taken as *1, giving "
                          "activity score >=3 and ultrarapid metaboliser for any N>=2")
        r["notes"].append("if the duplicated allele were a no-function allele (e.g. *4xN) the "
                          "phenotype would be entirely different - long-range PCR or phased data "
                          "are needed to exclude this")

    else:
        # ---- generic SNV path: CYP2D6/CYP2C19/CYP2C9/CYP2B6/DPYD/TPMT/NUDT15/SLCO1B1/G6PD
        obs, notes, unres = parse_snv_observations(gene, gt)
        r["notes"] += notes
        if unres:
            r.update(abstained=True,
                     abstain_reason="variant(s) not in the CPIC allele definition table: "
                                    + ", ".join(unres))
            return r
        # a locus reported as a single base carries one copy (X-linked hemizygosity)
        hemizygous = bool(obs) and all(len(v) == 1 for v in obs.values())
        dips = call_diplotypes(gene, obs, ploidy=1 if hemizygous else 2)
        if not dips:
            dips = call_diplotypes(gene, obs, ploidy=1 if hemizygous else 2, tag_snp=True)
            if dips:
                r["notes"].append("no diplotype was consistent under the strict panel assumption; "
                                  "the call uses the tag-variant (linkage) assumption stated in the "
                                  "case, i.e. the observed marker variant implies the rest of the "
                                  "haplotype (its other defining variants were not typed)")
        if hemizygous and re.search(r"homozygous female", gt, re.I):
            dips2 = call_diplotypes(gene, {k: (v * 2)[:2] for k, v in obs.items()}, ploidy=2)
            ph1 = {phenotype_row(gene, d)["generesult"] for d in dips}
            ph2 = {phenotype_row(gene, d)["generesult"] for d in dips2}
            if ph1 == ph2:
                r["notes"].append("hemizygous male and homozygous female readings give the same "
                                  "CPIC phenotype: " + "/".join(sorted(ph1)))
                r["_hemi_single"] = primary_diplotype(gene, dips)
                dips = sorted(set(dips) | set(dips2))
            else:
                r.update(abstained=True, abstain_reason="hemizygous and homozygous readings of the "
                         "same input give different phenotypes")
                return r
        if not dips:
            r.update(abstained=True,
                     abstain_reason="no CPIC diplotype is consistent with the observed genotype "
                                    "under the assumption that untyped positions are reference")
            return r
        rows = [phenotype_row(gene, d) for d in dips]
        phenos = {x["generesult"] for x in rows}
        keys = {json.dumps(x["lookupkey"], sort_keys=True) for x in rows}
        r["candidates"] = ";".join(dips)
        if len(phenos) > 1:
            r.update(abstained=True, diplotype=" or ".join(dips),
                     abstain_reason="genotype is consistent with >1 diplotype with different "
                                    "phenotypes (" + " / ".join(sorted(phenos)) + "); "
                                    "phase/copy number cannot be resolved from the input")
            return r
        r["phenotype"] = phenos.pop()
        r["pheno_source"] = f"CPIC diplotype table ({gene})"
        if r.get("_hemi_single"):
            s = r["_hemi_single"]
            r["diplotype"] = f"{s} (hemizygous male) or {s}/{s} (homozygous female)"
        elif len(dips) == 1:
            r["diplotype"] = dips[0]
        else:
            best = primary_diplotype(gene, dips)
            others = [d for d in dips if d != best]
            show = others[:3] + (["..."] if len(others) > 3 else [])
            r["diplotype"] = f"{best} (equivalent alternatives: {', '.join(show)})"
            r["notes"].append(f"{len(dips)} diplotypes are consistent with the unphased, "
                              "SNV-panel-limited genotype and all give the same CPIC phenotype; "
                              "the call reported is the one with the fewest non-reference alleles")
        if len(keys) > 1:
            r["notes"].append("candidate diplotypes differ in activity score (" +
                              ", ".join(sorted(json.loads(k)[gene] for k in keys)) +
                              "); the recommendation was checked to be identical for all")
        r["lookup_value"] = rows[0]["lookupkey"][gene]
        r["all_lookup_values"] = sorted(json.loads(k)[gene] for k in keys)
        # sensitivity check: would the phenotype change if untyped loci were treated as
        # unknown rather than reference (i.e. if the panel simply did not cover them)?
        relaxed = call_diplotypes(gene, obs, ploidy=1 if hemizygous else 2, tag_snp=True)
        rel_ph = sorted({phenotype_row(gene, d)["generesult"] for d in relaxed}) if relaxed else []
        r["relaxed_phenotypes"] = rel_ph
        if rel_ph and set(rel_ph) != {r["phenotype"]}:
            r["notes"].append("if the untyped positions are treated as unknown rather than "
                              "reference, the genotype would also be compatible with: "
                              + ", ".join(p for p in rel_ph if p != r["phenotype"])
                              + " - the reported call relies on the panel convention that "
                                "unreported positions are reference")
        if gene in ("CYP2D6", "CYP2C19", "CYP2C9", "CYP2B6", "TPMT", "NUDT15", "SLCO1B1"):
            r["notes"].append(CNV_CAVEAT)

    # ---- 2. phenotype/allele status -> recommendation
    lethal, mech = lethal_class(gk, drug, r["phenotype"], r["lookup_value"])
    r["lethal_class"], r["harm_mechanism"] = lethal, mech

    if drug.lower() == "warfarin":
        _warfarin_recommendation(r)
    elif drug.lower() == "ondansetron":
        key = r["phenotype"]
        txt = ONDANSETRON_REC.get(key)
        if txt:
            r["recommendation"] = f"ondansetron: {txt}"
            r["rec_source"] = ("CPIC ondansetron/CYP2D6 guideline table, supplied corpus "
                               "(PharmGKB PA166161954); no machine-readable rows in the CPIC API")
            r["rec_classification"] = "Moderate (UM) / Strong (NM) / No recommendation (IM, PM)"
        else:
            r.update(abstained=True, abstain_reason="no CPIC ondansetron recommendation for "
                                                    f"phenotype {key}")
    elif drug == "PEG-IFN-α":
        txt = IFNL3_REC.get(r["phenotype"])
        r["recommendation"] = f"PEG-IFN-α (peginterferon alfa-2a/2b) + ribavirin: {txt}"
        r["rec_source"] = ("CPIC IFNL3/peginterferon guideline, supplied corpus "
                           "(PharmGKB PA166110235); no rows in the CPIC API recommendation table")
        r["rec_classification"] = "Strong"
        r["notes"].append("the CPIC gene-drug pair level for IFNL3/peginterferon alfa is now "
                          "'Retired' in the CPIC pair table (drug largely superseded by "
                          "direct-acting antivirals), so this is no longer a current Level A pair")
    elif drug.lower() == "irinotecan":
        r["recommendation"] = ("irinotecan: no CPIC recommendation published - CPIC lists "
                               "UGT1A1/irinotecan as a level A pair but has no guideline for it "
                               "(usedforrecommendation = false, no guideline id); only DPWG and "
                               "RNPGx guidance exists, which is not CPIC")
        r["rec_source"] = "CPIC pair table + supplied corpus (DPWG/RNPGx annotations only)"
        r.update(abstained=True,
                 abstain_reason="no CPIC Level A recommendation text exists for UGT1A1-irinotecan; "
                                "phenotype is reported but the recommendation is abstained")
    else:
        rows, used, pops, homo = get_recommendation(drug, gk, r["lookup_value"])
        if not rows:
            r.update(abstained=True,
                     abstain_reason=f"no CPIC recommendation row for {gk}={r['lookup_value']} "
                                    f"and {drug}")
        else:
            if not homo:
                texts = sorted({x["drugrecommendation"].strip() for x in rows})
                r["recommendation"] = f"{drug}: " + " | ".join(squeeze(t, 200) for t in texts)
                r["notes"].append("CPIC gives different text for the drugs/populations in this "
                                  "class; all are reported")
            else:
                r["recommendation"] = f"{drug}: {squeeze(rows[0]['drugrecommendation'])}"
            # if several diplotypes with different activity scores survived, confirm identical text
            if r.get("all_lookup_values") and len(r["all_lookup_values"]) > 1:
                alts = set()
                for lv in r["all_lookup_values"]:
                    rws, _, _, _ = get_recommendation(drug, gk, lv)
                    alts |= {x["drugrecommendation"].strip() for x in rws}
                if len(alts) > 1:
                    r.update(abstained=True,
                             abstain_reason="candidate diplotypes share a phenotype but CPIC gives "
                                            "different recommendations for their activity scores")
            r["rec_source"] = ("CPIC data API recommendation table (guideline "
                               + "/".join(sorted({str(x['guidelineid']) for x in rows}))
                               + "; drugs matched: " + ", ".join(used) + ")")
            r["rec_classification"] = "/".join(sorted({str(x["classification"]) for x in rows}))
            r["rec_population"] = "/".join(pops)
            imps = sorted({(x.get("implications") or {}).get(gk, "") for x in rows})
            r["implications"] = squeeze(" ".join(i for i in imps if i), 400)
            if len(pops) == 1 and pops[0] not in ("general",):
                r["notes"].append(f"CPIC recommendation is population-specific; the '{pops[0]}' "
                                  "row was used because the case gives no clinical context")
    return r

def _warfarin_recommendation(r):
    """CPIC warfarin dosing from the guideline table in the supplied corpus."""
    gene = r["gene"]
    src = ("CPIC warfarin guideline (CYP2C9, VKORC1, CYP4F2) dose table in the supplied corpus "
           "(PharmGKB PA166104949). CPIC states warfarin 'does not follow simple diplotype to "
           "phenotype translation' and publishes no machine-readable recommendation rows.")
    r["rec_source"] = src
    r["rec_classification"] = "see guideline (algorithm-based dosing)"
    if gene == "CYP2C9":
        dip = r["diplotype"].split(" ")[0]
        doses = {v: WARFARIN_TABLE[(dip, v)] for v in ("GG", "GA", "AA")
                 if (dip, v) in WARFARIN_TABLE}
        if not doses:
            r.update(abstained=True, abstain_reason="CYP2C9 diplotype not in the CPIC warfarin "
                                                    "dose table (table covers *1/*2/*3 only)")
            return
        uniq = sorted(set(doses.values()))
        if len(uniq) == 1:
            r["recommendation"] = (f"warfarin: use a reduced starting dose of {uniq[0]} mg/day "
                                   f"(CPIC table, all VKORC1 -1639 genotypes) with INR-guided "
                                   f"titration; consider an alternative anticoagulant if dosing "
                                   f"cannot be managed")
        else:
            r["recommendation"] = ("warfarin: dose depends on VKORC1 -1639 (rs9923231), which was "
                                   "not typed - " + ", ".join(f"{k} {v} mg/day" for k, v in doses.items())
                                   + "; use the CPIC/pharmacogenetic dosing algorithm with INR monitoring")
            r["notes"].append("VKORC1 genotype absent: the CPIC dose depends on it, so the "
                              "recommendation is conditional")
    elif gene == "VKORC1":
        v = "AA" if "variant" in r["diplotype"] else "GG"
        doses = {dip: WARFARIN_TABLE[(dip, v)] for dip in ("*1/*1", "*1/*3", "*3/*3")}
        r["recommendation"] = (f"warfarin: with VKORC1 -1639 {v}, CPIC dose is "
                               f"{doses['*1/*1']} mg/day if CYP2C9 *1/*1, "
                               f"{doses['*1/*3']} mg/day if CYP2C9 *1/*3 and "
                               f"{doses['*3/*3']} mg/day if CYP2C9 *3/*3; CYP2C9 was not typed, so "
                               f"use the CPIC dosing algorithm with INR monitoring")
        r["notes"].append("CYP2C9 genotype absent: the CPIC dose depends on it, so the "
                          "recommendation is conditional")
    else:   # CYP4F2
        r["recommendation"] = ("warfarin: no CPIC dose recommendation can be derived from CYP4F2 "
                               "alone - the CPIC dose table is indexed by VKORC1 x CYP2C9 and "
                               "CYP4F2 is only a minor dose modifier")
        r.update(abstained=True,
                 abstain_reason="CPIC warfarin dosing is algorithm-based and indexed on VKORC1 and "
                                "CYP2C9; CYP4F2 genotype alone does not yield a dose recommendation "
                                "(no CYP4F2 phenotype table and no CYP4F2 recommendation rows in "
                                "the CPIC API; the supplied corpus excerpt contains no CYP4F2 "
                                "dosing statement)")

# ==================================================================== outputs
def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

def main():
    results = [resolve(c) for c in cases]

    # 1. answers
    write_csv(os.path.join(HERE, "claude_science_answers.csv"),
              ["case_id", "diplotype", "phenotype", "recommendation", "abstained"],
              [[r["case_id"], r["diplotype"], r["phenotype"], r["recommendation"],
                "true" if r["abstained"] else "false"] for r in results])

    # 2a. audit trail
    write_csv(os.path.join(HERE, "case_audit.csv"),
              ["case_id", "gene", "drug", "genotype", "diplotype", "candidate_diplotypes",
               "phenotype", "cpic_lookup_value", "recommendation", "abstained", "abstain_reason",
               "lethal_class", "harm_mechanism", "phenotype_if_untyped_loci_unknown",
               "cpic_pair_level", "rec_classification",
               "rec_population", "cpic_implication", "phenotype_source", "recommendation_source",
               "notes"],
              [[r["case_id"], r["gene"], r["drug"], r["genotype"], r["diplotype"], r["candidates"],
                r["phenotype"], r["lookup_value"], r["recommendation"],
                "true" if r["abstained"] else "false", r["abstain_reason"],
                "true" if r.get("lethal_class") else "false", r.get("harm_mechanism", ""),
                "; ".join(r.get("relaxed_phenotypes") or []),
                r["cpic_level"], r["rec_classification"], r["rec_population"], r["implications"],
                r["pheno_source"], r["rec_source"], " | ".join(r["notes"])] for r in results])

    # 2b. diplotype -> phenotype table (all diplotypes used, plus the full table for the
    #     genes whose tables are small enough to be auditable)
    used_dips = collections.defaultdict(set)
    for r, c in zip(results, cases):
        for d in (r["candidates"].split(";") if r["candidates"] else [r["diplotype"]]):
            if d:
                used_dips[r["gene"]].add(d)
    rows = []
    for gene in sorted(DIP):
        small = len(DIP[gene]) <= 100
        for row in DIP[gene]:
            in_use = row["diplotype"] in used_dips.get(gene, set())
            if not (small or in_use):
                continue
            rows.append([gene, row["diplotype"], row["generesult"],
                         row["totalactivityscore"], row["function1"], row["function2"],
                         json.dumps(row["lookupkey"]), "yes" if in_use else "no"])
    write_csv(os.path.join(HERE, "table_diplotype_to_phenotype.csv"),
              ["gene", "diplotype", "cpic_phenotype", "total_activity_score", "function_allele_1",
               "function_allele_2", "cpic_lookupkey", "used_by_a_case"], rows)

    # 2c. phenotype -> recommendation table for every study gene-drug pair
    rows = []
    study_pairs = sorted({(c["gene"], c["drug"]) for c in cases})
    for gene, drug in study_pairs:
        gk = GENE_KEY.get(gene, gene)
        names = CLASS_DRUGS.get(drug, [drug.lower()])
        seen = set()
        for n in names:
            did = NAME2DRUGID.get(n)
            for rec in [x for x in REC if x["drugid"] == did]:
                lk = rec.get("lookupkey") or {}
                if gk not in lk:
                    continue
                if any(v != "No Result" for k, v in lk.items() if k != gk):
                    continue
                key = (gene, drug, n, rec["population"], lk.get(gk))
                if key in seen:
                    continue
                seen.add(key)
                rows.append([gene, drug, n, lk.get(gk), rec["population"],
                             (rec.get("phenotypes") or {}).get(gk, ""),
                             " ".join((rec["drugrecommendation"] or "").split()),
                             rec["classification"],
                             " ".join(((rec.get("implications") or {}).get(gk, "")).split()),
                             "CPIC API recommendation table"])
        if gene == "CYP2D6" and drug == "ondansetron":
            for ph, txt in ONDANSETRON_REC.items():
                rows.append([gene, drug, "ondansetron", ph, "general", ph, txt, "see guideline", "",
                             "CPIC guideline table, supplied corpus (PA166161954)"])
        if drug == "PEG-IFN-α":
            for ph, txt in IFNL3_REC.items():
                rows.append([gene, drug, "peginterferon alfa-2a/2b", ph, "general", ph, txt,
                             "Strong", "", "CPIC guideline table, supplied corpus (PA166110235)"])
        if drug == "warfarin":
            for (dip, vk), dose in sorted(WARFARIN_TABLE.items()):
                rows.append([gene, drug, "warfarin", f"CYP2C9 {dip} + VKORC1 -1639 {vk}", "general",
                             "", f"{dose} mg/day starting dose (CPIC dose table)",
                             "see guideline", "",
                             "CPIC guideline table, supplied corpus (PA166104949)"])
    write_csv(os.path.join(HERE, "table_phenotype_to_recommendation.csv"),
              ["gene", "study_drug_label", "cpic_drug", "cpic_lookup_value", "population",
               "cpic_phenotype", "cpic_recommendation", "classification", "cpic_implication",
               "source"], rows)

    # 2d. variant -> allele definition table for the variants named in the cases
    rsids = sorted({m.group(0) for c in cases for m in re.finditer(r"rs\d+", c["genotype"])})
    rows = []
    for c_gene in sorted({c["gene"] for c in cases}):
        gene = c_gene
        for rsid in rsids:
            loc = RSID2LOC.get(gene, {}).get(rsid)
            if loc is None:
                continue
            s = SL[loc]
            for name, dd in sorted(ALLELE_DEF[gene].items()):
                if loc in dd and dd[loc]:
                    am = next((a for a in T["allele"]
                               if a["genesymbol"] == gene and a["name"] == name), {})
                    rows.append([gene, rsid, s["name"], s["chromosomelocation"], name,
                                 "".join(sorted(dd[loc])),
                                 am.get("clinicalfunctionalstatus") or am.get("functionalstatus") or "",
                                 am.get("activityvalue") or "",
                                 "yes" if name in dip_vocabulary(gene) else "no"])
    write_csv(os.path.join(HERE, "table_variant_to_allele.csv"),
              ["gene", "rsid", "cpic_variant_name", "genomic_location", "allele",
               "allele_base_plus_strand", "cpic_allele_function", "cpic_activity_value",
               "in_cpic_diplotype_table"], rows)

    # console summary
    n_abs = sum(r["abstained"] for r in results)
    n_leth = sum(bool(r.get("lethal_class")) for r in results)
    print(f"cases={len(results)} abstained={n_abs} lethal_class={n_leth}")
    print("abstentions:")
    for r in results:
        if r["abstained"]:
            print(f"  {r['case_id']} {r['gene']:12s} {r['drug']:16s} {r['abstain_reason'][:110]}")
    return results

if __name__ == "__main__":
    RESULTS = main()
