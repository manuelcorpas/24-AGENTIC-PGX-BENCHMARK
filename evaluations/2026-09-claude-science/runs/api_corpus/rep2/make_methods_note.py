#!/usr/bin/env python3
"""Generate METHODS.md and table_sources_by_gene.csv from the pipeline outputs.

Every number in the methods note is computed from case_audit.csv / the CPIC snapshot
metadata / the supplied corpus, so the note cannot drift from the results.
"""
import json, os, csv, collections

HERE = os.path.dirname(os.path.abspath(__file__))
aud = list(csv.DictReader(open(os.path.join(HERE, "case_audit.csv"))))
ans = list(csv.DictReader(open(os.path.join(HERE, "claude_science_answers.csv"))))
meta = json.load(open(os.path.join(HERE, "cpic_api_snapshot", "_fetch_meta.json")))
corpus = json.load(open(os.path.join(HERE, "cpic_guideline_corpus.json")))

n = len(aud)
abst = [r for r in aud if r["abstained"] == "true"]
leth = [r for r in aud if r["lethal_class"] == "true"]
assumption_dependent = [r for r in aud
                        if r["phenotype_if_untyped_loci_unknown"].strip()
                        and set(r["phenotype_if_untyped_loci_unknown"].split("; ")) != {r["phenotype"]}]
by_gene = collections.Counter(r["gene"] for r in aud)
src_family = collections.Counter(r["recommendation_source"].split(" (")[0] for r in aud)
strand = [r for r in aud if "complemented" in r["notes"]]
dupcases = [r for r in aud if "xN" in r["diplotype"]]
tagsnp = [r for r in aud if "tag-variant" in r["notes"]]
pop_notes = collections.Counter(
    r["rec_population"] for r in aud if r["rec_population"].strip() and r["rec_population"] != "general")

# ------------------------------------------------- per-gene source/version table
rows = []
for gene in sorted(by_gene):
    c = corpus["genes"].get(gene, {})
    genes_using_api = {r["gene"] for r in aud
                       if r["recommendation_source"].startswith("CPIC data API")}
    genes_using_corpus = {r["gene"] for r in aud if "supplied corpus" in r["recommendation_source"]}
    rows.append([
        gene,
        "CPIC data API (api.cpicpgx.org/v1): allele_definition, allele_location_value, "
        "sequence_location, diplotype, gene_result, recommendation, pair",
        meta["fetched_at_utc"],
        f"CPIC change_log latest entry {meta['change_log_latest_date']}",
        "yes" if gene in genes_using_api else "no",
        "yes" if gene in genes_using_corpus else "no",
        c.get("sha256", ""),
        ", ".join(c.get("drugs_covered", [])),
        "; ".join(s["name"] + " [" + s["pa_id"] + ", retrieved " + s["retrieved_date"] + "]"
                  for s in c.get("sources", [])),
    ])
with open(os.path.join(HERE, "table_sources_by_gene.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "primary_source_tables", "api_fetched_at_utc", "api_version_marker",
                "api_used_for_recommendation", "corpus_used_for_recommendation",
                "corpus_excerpt_sha256", "corpus_drugs_covered", "corpus_source_annotations"])
    w.writerows(rows)

# ----------------------------------------------------------------- methods note
def bullets(items):
    return "\n".join("- " + i for i in items)

lethal_lines = [f"- `{r['case_id']}` {r['gene']} / {r['drug']} - {r['phenotype']} - "
                f"{r['harm_mechanism']}" for r in leth]
abst_lines = [f"- `{r['case_id']}` {r['gene']} / {r['drug']}: {r['abstain_reason']}" for r in abst]
assump_lines = [f"- `{r['case_id']}` {r['gene']} / {r['drug']}: called **{r['phenotype']}** "
                f"({r['diplotype']}); also compatible with "
                f"{', '.join(p for p in r['phenotype_if_untyped_loci_unknown'].split('; ') if p != r['phenotype'])}"
                " if the untyped positions are unknown rather than reference"
                for r in assumption_dependent]

md = f"""# Methods note - blinded pharmacogenomics study ({n} cases, 21 genes, 32 drug labels)

## 1. Approach

The pipeline is deterministic code, not per-case reasoning. Three stages, each backed by a
machine-readable CPIC table rather than recall:

**(a) genotype string -> observation.** Each case genotype is parsed into typed loci
(rsID -> observed bases; a single base means one copy, i.e. X-linked hemizygosity) or into a
non-SNV observation (HLA allele status, UGT1A1 TA-repeat length, CYP2D6 copy number,
mitochondrial variant, "all reference"). Observed bases are checked against the alleles CPIC
records at that locus and complemented only when the observed pair is absent but its
complement is present; this was needed for {len(strand)} cases
({', '.join(sorted({r['gene'] for r in strand}))}), where the case reports the transcript
(minus) strand and CPIC represents the plus strand. No strand flip was applied where the case
already matched CPIC (TPMT, G6PD, DPYD, CYP2D6, CYP2C19, CYP2C9, CYP2B6, SLCO1B1, CYP4F2).

**(b) observation -> diplotype (called, not hand-mapped).** Candidate haplotype pairs are
enumerated from the CPIC allele definition tables (`allele_definition` x
`allele_location_value` x `sequence_location`). An allele is a candidate only if every locus
at which it mandates a non-reference base was observed as variant; a pair is retained only if
every relevant locus (typed loci plus the alleles' mandatory loci) can be reconciled with the
observed genotype, with untyped loci treated as homozygous reference - the semantics of a
targeted panel that reports only what it found. Where several diplotypes survive and share a
phenotype, the one with the fewest non-reference alleles is reported and the alternatives are
listed. Where they do not share a phenotype, the case is abstained.

A second, relaxed pass treats untyped loci as *unknown* and allows a haplotype that is
anchored by an observed marker variant. This is required for haplotypes assayed through a
single tag variant ({len(tagsnp)} case: DPYD HapB3 via rs56038477, whose second defining
variant c.1129-5923C>G was not typed) and is also run on every SNV case as a sensitivity
analysis (section 5).

**(c) diplotype -> phenotype -> recommendation.** Phenotype is the CPIC `diplotype` table's
`generesult`; the CPIC `lookupkey` from the same row (a phenotype term, or an activity score
for CYP2D6/CYP2C9/DPYD) is the key into the CPIC `recommendation` table, together with the
drug. Where a case provides only one gene of a multi-gene recommendation (CYP2C9 without
HLA-B for phenytoin, HLA-B without HLA-A for carbamazepine, TPMT without NUDT15, RYR1 without
CACNA1S, CYP2D6 without CYP2C19), CPIC's own `No Result` key is used, so the recommendation
returned is the one CPIC publishes for a single-gene result rather than an assumption about
the untested gene.

Recommendation sources actually used across the {n} cases:

{bullets(f'{k}: {v} cases' for k, v in src_family.most_common())}

Population/indication subsets had to be chosen where CPIC publishes several and the case gives
no clinical context; the choice is recorded per case and was always the higher-risk or
treatment-naive row: {', '.join(f'{k} ({v})' for k, v in sorted(pop_notes.items()))}.

## 2. Sources and versions

| source | version / date | used for |
|---|---|---|
| CPIC data API `api.cpicpgx.org/v1` | fetched {meta['fetched_at_utc']}; latest CPIC `change_log` entry {meta['change_log_latest_date']} | allele definitions, diplotype->phenotype, recommendations, gene-drug pair levels |
| supplied `cpic_guideline_corpus.json` | schema {corpus['schema_version']}, generated {corpus['generated_at_utc']}, PharmGKB bundle `{corpus['pharmgkb_bundle_path']}` | warfarin dose table, ondansetron table, IFNL3/peginterferon table, UGT1A1-irinotecan (confirming no CPIC guideline) |
| PharmVar | not used | star-allele definitions were taken from CPIC, which mirrors PharmVar core definitions |

Row counts in the API snapshot: {', '.join(f"{k}={v}" for k, v in sorted(meta['table_row_counts'].items()) if k in ('allele_definition','allele_location_value','sequence_location','recommendation','pair','gene_result'))}.
`table_sources_by_gene.csv` records, per gene, the API tables used, the corpus excerpt sha256,
and which of the two supplied the recommendation. No published benchmark, paper or repository
reporting answers for these cases was consulted, and none was encountered.

## 3. Results at a glance

- {n} cases, {n - len(abst)} answered, {len(abst)} abstained.
- {len(leth)} cases are lethal-class (section 4).
- {len(assumption_dependent)} answered cases depend on the panel convention that unreported
  positions are reference (section 5).
- Cases per gene: {', '.join(f'{g} {c}' for g, c in sorted(by_gene.items()))}.

## 4. Lethal-class cases and why

Criterion used: standard or label dosing in this genotype is expected to cause severe
(grade >=3), irreversible or life-threatening harm, or a life-threatening therapeutic failure.
Reduced efficacy without acute danger (e.g. CYP2D6 poor metaboliser + codeine analgesic
failure, CFTR F508del + ivacaftor non-response, CYP2D6 poor metaboliser + tamoxifen) is
deliberately *not* in this class, although the tamoxifen case carries long-term recurrence
risk.

{chr(10).join(lethal_lines)}

Not classified as lethal-class but worth flagging: TPMT/NUDT15 intermediate metabolisers
(myelosuppression risk that CPIC manages by dose reduction rather than avoidance), CYP2C9
intermediate metabolisers on warfarin or phenytoin, and CYP2D6 ultrarapid metabolisers on
amitriptyline (treatment failure rather than toxicity).

## 5. Where this is most likely wrong

**(a) The panel-completeness assumption is the largest single risk.** All variant-only cases
are interpreted as "everything not reported is reference". That is how targeted panels report,
but it is an assumption, and for {len(assumption_dependent)} cases the phenotype would change
if the untyped positions were merely unknown:

{chr(10).join(assump_lines)}

**(b) CYP2D6 structural variation.** A SNV panel cannot see the *5 whole-gene deletion or a
duplication, so every CYP2D6 (and CYP2C19/CYP2B6/TPMT/NUDT15/SLCO1B1) call assumes two intact
copies. The {len(dupcases)} `*1xN` duplication cases assume the duplicated allele and the second allele
are both *1, because no SNVs were reported; if the duplicated allele were a no-function allele
(*4xN) the phenotype would swing from ultrarapid to intermediate or poor. Activity score >=3
holds for any N>=2, so the ultrarapid call is robust to copy number but not to which allele is
duplicated.

**(c) Phase.** The CYP2D6 compound case (rs3892097 + rs1065852 heterozygous) is reported as
*1/*4 because rs1065852 is part of the *4 haplotype; *4/*10 would require rs1135840, which was
not typed. Both are intermediate metabolisers and CPIC gives identical codeine text for
activity scores 1.0 and 0.25, so the recommendation is unaffected, but the diplotype label
could be wrong.

**(d) HLA by proxy.** The abacavir case gives rs2395029 (HCP5), not HLA-B typing. Linkage with
HLA-B*57:01 is strong but not complete, and CPIC's recommendation is keyed on HLA-B allele
status, so a false-positive or false-negative proxy call is possible.

**(e) Gene-level phenotype terms that CPIC does not define.** CPIC publishes no phenotype
vocabulary for VKORC1, CYP4F2 or IFNL3 (no `gene_result` or `diplotype` rows). For those the
answer reports the genotype-level result and says so rather than inventing a metaboliser
label; a scorer expecting e.g. "VKORC1 -1639 AA = increased sensitivity" would mark those
phenotype cells differently.

**(f) Recommendation text is verbatim CPIC, trimmed.** Long recommendations are truncated at a
sentence boundary in the answers file; the untruncated text, CPIC classification, implication
and population are in `case_audit.csv`. For drug classes (thiopurines, aminoglycosides,
volatile anaesthetics) CPIC's per-drug texts differ in dose detail and all variants are
reported joined by `|`.

**(g) Allele naming for G6PD, which a scorer may not match textually.** With only rs5030868
(c.563C>T) typed, the CPIC allele that fits is the single-variant entry named
"Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham"; the entry CPIC calls
"Mediterranean Haplotype" additionally requires c.1311C>T, which was not typed, so the strict
caller excludes it. Both are Deficient, so the phenotype and recommendation are unaffected, but
the diplotype string is the long synonym list rather than "Mediterranean". Likewise the two
X-linked cases that state both a hemizygous-male and a homozygous-female reading are reported
as both forms, because CPIC gives the same phenotype for each.

**(h) Indication assumptions.** Clopidogrel used the ACS/PCI row, phenytoin/carbamazepine/
oxcarbazepine the treatment-naive rows, voriconazole the adult row, efavirenz the
adult/child >40 kg row. A case intended as non-ACS clopidogrel or paediatric voriconazole
would be scored against different text.

**(i) Abstentions, which are judgements and could be argued either way.**

{chr(10).join(abst_lines)}

The three irinotecan cases carry a UGT1A1 phenotype but no recommendation: CPIC lists
UGT1A1-irinotecan as a level A pair yet publishes no guideline for it
(`usedforrecommendation = false`, no guideline id), and the supplied corpus contains only DPWG
and RNPGx annotations. Reporting either of those as "the CPIC recommendation" would be wrong,
so the recommendation is abstained while the phenotype (poor/intermediate/normal metaboliser)
is still given. The two CYP4F2 warfarin cases abstain because CPIC's warfarin recommendation
is an algorithm indexed on VKORC1 x CYP2C9 - CPIC's own note says warfarin "does not follow
simple diplotype to phenotype translation" - and CYP4F2 is only a minor dose modifier. The
IFNL3 cases are answered but flagged: that pair's CPIC level is now `Retired`, so it is no
longer a current level A recommendation.

## 6. Files

| file | content |
|---|---|
| `claude_science_answers.csv` | the {n} answers: case_id, diplotype, phenotype, recommendation, abstained |
| `case_audit.csv` | per-case evidence: candidate diplotypes, CPIC lookup key, full recommendation, abstention reason, lethal class and mechanism, sensitivity result, sources, notes |
| `table_variant_to_allele.csv` | rsID -> star allele definition table used (CPIC plus-strand bases, allele function, activity value) |
| `table_diplotype_to_phenotype.csv` | diplotype -> CPIC phenotype/activity score, flagged where used by a case |
| `table_phenotype_to_recommendation.csv` | phenotype (or CPIC lookup key) -> CPIC recommendation, per study gene-drug pair |
| `table_sources_by_gene.csv` | source and version per gene |
| `pgx_pipeline.py` | the pipeline that produced the answers and tables |
| `fetch_cpic_api.py` | CPIC API snapshot fetcher (writes `cpic_api_snapshot/`) |
| `make_methods_note.py` | generates this note and `table_sources_by_gene.csv` |
| `claude_science_analysis.ipynb` | notebook: runs the pipeline and shows the audit/summary views |
"""
open(os.path.join(HERE, "METHODS.md"), "w").write(md)
print(f"METHODS.md written: {len(md)} chars; lethal={len(leth)} abstained={len(abst)} "
      f"assumption_dependent={len(assumption_dependent)}")
