# Methods note - blinded pharmacogenomics study (110 cases, 21 genes, 32 drug labels)

## 1. Approach

The pipeline is deterministic code, not per-case reasoning. Three stages, each backed by a
machine-readable CPIC table rather than recall:

**(a) genotype string -> observation.** Each case genotype is parsed into typed loci
(rsID -> observed bases; a single base means one copy, i.e. X-linked hemizygosity) or into a
non-SNV observation (HLA allele status, UGT1A1 TA-repeat length, CYP2D6 copy number,
mitochondrial variant, "all reference"). Observed bases are checked against the alleles CPIC
records at that locus and complemented only when the observed pair is absent but its
complement is present; this was needed for 5 cases
(CYP3A5, VKORC1), where the case reports the transcript
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
single tag variant (1 case: DPYD HapB3 via rs56038477, whose second defining
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

Recommendation sources actually used across the 110 cases:

- CPIC data API recommendation table: 95 cases
- CPIC warfarin guideline: 7 cases
- CPIC ondansetron/CYP2D6 guideline table, supplied corpus: 3 cases
- CPIC pair table + supplied corpus: 3 cases
- CPIC IFNL3/peginterferon guideline, supplied corpus: 2 cases

Population/indication subsets had to be chosen where CPIC publishes several and the case gives
no clinical context; the choice is recorded per case and was always the higher-risk or
treatment-naive row: CBZ naive (4), CVI ACS PCI  (5), OXC naive (2), PHT naive (3), adults (5), child >40kg_adult (4).

## 2. Sources and versions

| source | version / date | used for |
|---|---|---|
| CPIC data API `api.cpicpgx.org/v1` | fetched 2026-09-16T14:27:50.593227+00:00; latest CPIC `change_log` entry 2026-08-03 | allele definitions, diplotype->phenotype, recommendations, gene-drug pair levels |
| supplied `cpic_guideline_corpus.json` | schema v3.1, generated 2026-05-10T04:40:15Z, PharmGKB bundle `SPECS/cpic_rag_cache/_pharmgkb_bundle/20260509T203526Z` | warfarin dose table, ondansetron table, IFNL3/peginterferon table, UGT1A1-irinotecan (confirming no CPIC guideline) |
| PharmVar | not used | star-allele definitions were taken from CPIC, which mirrors PharmVar core definitions |

Row counts in the API snapshot: allele_definition=1354, allele_location_value=3160, gene_result=101, pair=635, recommendation=2115, sequence_location=1210.
`table_sources_by_gene.csv` records, per gene, the API tables used, the corpus excerpt sha256,
and which of the two supplied the recommendation. No published benchmark, paper or repository
reporting answers for these cases was consulted, and none was encountered.

## 3. Results at a glance

- 110 cases, 105 answered, 5 abstained.
- 29 cases are lethal-class (section 4).
- 15 answered cases depend on the panel convention that unreported
  positions are reference (section 5).
- Cases per gene: CFTR 2, CYP2B6 4, CYP2C19 20, CYP2C9 9, CYP2D6 20, CYP3A5 3, CYP4F2 2, DPYD 7, G6PD 6, HLA-A*31:01 2, HLA-B*15:02 4, HLA-B*57:01 2, HLA-B*58:01 2, IFNL3 2, MT-RNR1 2, NUDT15 3, RYR1 2, SLCO1B1 6, TPMT 4, UGT1A1 6, VKORC1 2.

## 4. Lethal-class cases and why

Criterion used: standard or label dosing in this genotype is expected to cause severe
(grade >=3), irreversible or life-threatening harm, or a life-threatening therapeutic failure.
Reduced efficacy without acute danger (e.g. CYP2D6 poor metaboliser + codeine analgesic
failure, CFTR F508del + ivacaftor non-response, CYP2D6 poor metaboliser + tamoxifen) is
deliberately *not* in this class, although the tamoxifen case carries long-term recurrence
risk.

- `74597c5503ec` CYP2D6 / codeine - Ultrarapid Metabolizer - morphine overdose: fatal respiratory depression
- `d456dcf2b837` CYP2D6 / amitriptyline - Poor Metabolizer - tricyclic accumulation: seizures, cardiac arrhythmia
- `ca9e14feca53` CYP2C19 / clopidogrel - Poor Metabolizer - loss of antiplatelet effect: stent thrombosis, MI, death (ACS/PCI)
- `97d4e4559a13` CYP2C19 / clopidogrel - Intermediate Metabolizer - loss of antiplatelet effect: stent thrombosis, MI, death (ACS/PCI)
- `00b16db7406d` CYP2C19 / voriconazole - Poor Metabolizer - voriconazole accumulation: hepatotoxicity, neurotoxicity
- `1fb04aa9e16c` CYP2C19 / voriconazole - Rapid Metabolizer - subtherapeutic exposure: failure of therapy for invasive fungal infection
- `0c21c969a5c4` CYP2C19 / voriconazole - Ultrarapid Metabolizer - subtherapeutic exposure: failure of therapy for invasive fungal infection
- `97c3baebed4e` CYP2C9 / warfarin - Poor Metabolizer - over-anticoagulation: major/fatal haemorrhage
- `e85ecfa887cb` CYP2C9 / phenytoin - Poor Metabolizer - phenytoin accumulation: severe neurotoxicity
- `b94a43e574db` DPYD / fluorouracil - Poor Metabolizer - fluoropyrimidine overdose: fatal myelosuppression, mucositis, neurotoxicity
- `c57caa44a38a` DPYD / fluorouracil - Intermediate Metabolizer - partial DPD deficiency: severe (grade >=3), sometimes fatal fluoropyrimidine toxicity at full dose
- `6b8d2db91001` DPYD / fluorouracil - Intermediate Metabolizer - partial DPD deficiency: severe (grade >=3), sometimes fatal fluoropyrimidine toxicity at full dose
- `c6549d637d49` DPYD / capecitabine - Poor Metabolizer - fluoropyrimidine overdose: fatal myelosuppression, mucositis, neurotoxicity
- `ae4a7215ae64` DPYD / capecitabine - Intermediate Metabolizer - partial DPD deficiency: severe (grade >=3), sometimes fatal fluoropyrimidine toxicity at full dose
- `08e18e809fcb` TPMT / azathioprine - Poor Metabolizer - thiopurine accumulation: life-threatening myelosuppression
- `f3dd19d8a337` NUDT15 / thiopurines - Poor Metabolizer - thiopurine accumulation: life-threatening myelosuppression
- `d1a12c4a16c4` SLCO1B1 / simvastatin - Poor Function - statin myopathy / rhabdomyolysis
- `4f19dd3f5283` UGT1A1 / irinotecan - Poor Metabolizer - severe neutropenia and diarrhoea
- `e70379d3dbac` HLA-B*57:01 / abacavir - HLA-B*57:01 positive - abacavir hypersensitivity syndrome, fatal on re-challenge
- `7181a24a65d5` HLA-B*15:02 / carbamazepine - HLA-B*15:02 positive - Stevens-Johnson syndrome / toxic epidermal necrolysis
- `569a842e1a9a` HLA-B*15:02 / oxcarbazepine - HLA-B*15:02 positive - Stevens-Johnson syndrome / toxic epidermal necrolysis
- `85976e86289a` HLA-A*31:01 / carbamazepine - HLA-A*31:01 positive - DRESS, SJS/TEN, maculopapular exanthema
- `ebf82c3bcfdd` HLA-B*58:01 / allopurinol - HLA-B*58:01 positive - allopurinol-induced severe cutaneous adverse reaction (SJS/TEN, DRESS)
- `d49e9f88d9c2` G6PD / rasburicase - Deficient - acute haemolysis and methaemoglobinaemia; rasburicase is contraindicated in G6PD deficiency
- `3786ab9a929d` G6PD / rasburicase - Deficient - acute haemolysis and methaemoglobinaemia; rasburicase is contraindicated in G6PD deficiency
- `2a8544a44d5f` G6PD / primaquine - Deficient - acute haemolytic anaemia at standard primaquine dosing
- `7b560d699828` G6PD / primaquine - Variable - possible acute haemolysis - G6PD activity must be measured before use
- `369c2bcff713` MT-RNR1 / aminoglycosides - MT-RNR1 increased risk of aminoglycoside-induced hearing loss - permanent bilateral sensorineural hearing loss
- `1229a4123378` RYR1 / volatile-anaesthetics - Malignant Hyperthermia Susceptibility - malignant hyperthermia crisis

Not classified as lethal-class but worth flagging: TPMT/NUDT15 intermediate metabolisers
(myelosuppression risk that CPIC manages by dose reduction rather than avoidance), CYP2C9
intermediate metabolisers on warfarin or phenytoin, and CYP2D6 ultrarapid metabolisers on
amitriptyline (treatment failure rather than toxicity).

## 5. Where this is most likely wrong

**(a) The panel-completeness assumption is the largest single risk.** All variant-only cases
are interpreted as "everything not reported is reference". That is how targeted panels report,
but it is an assumption, and for 15 cases the phenotype would change
if the untyped positions were merely unknown:

- `3d0e97b86ad0` CYP2D6 / codeine: called **Intermediate Metabolizer** (*1/*4); also compatible with Indeterminate, Poor Metabolizer if the untyped positions are unknown rather than reference
- `97c3baebed4e` CYP2C9 / warfarin: called **Poor Metabolizer** (*3/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `9b5e3093ab14` CYP2C9 / warfarin: called **Intermediate Metabolizer** (*1/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `e85ecfa887cb` CYP2C9 / phenytoin: called **Poor Metabolizer** (*3/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `c9726ef5c899` CYP2C9 / phenytoin: called **Intermediate Metabolizer** (*1/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `b0db8db6df72` CYP2C9 / celecoxib: called **Poor Metabolizer** (*3/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `6b95f7ee2785` CYP2C9 / celecoxib: called **Intermediate Metabolizer** (*1/*3); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `240c62812cbc` CYP2B6 / efavirenz: called **Rapid Metabolizer** (*1/*4); also compatible with Intermediate Metabolizer if the untyped positions are unknown rather than reference
- `d1a12c4a16c4` SLCO1B1 / simvastatin: called **Poor Function** (*5/*5); also compatible with Indeterminate, Possible Decreased Function if the untyped positions are unknown rather than reference
- `51a382a7ca5c` SLCO1B1 / simvastatin: called **Decreased Function** (*1/*5); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `9a62df2e13c7` SLCO1B1 / atorvastatin: called **Poor Function** (*5/*5); also compatible with Indeterminate, Possible Decreased Function if the untyped positions are unknown rather than reference
- `3480326d71f3` SLCO1B1 / atorvastatin: called **Decreased Function** (*1/*5); also compatible with Indeterminate if the untyped positions are unknown rather than reference
- `3786ab9a929d` G6PD / rasburicase: called **Deficient** (A- 202A_376G (hemizygous male) or A- 202A_376G/A- 202A_376G (homozygous female)); also compatible with Deficient with CNSHA if the untyped positions are unknown rather than reference
- `2a8544a44d5f` G6PD / primaquine: called **Deficient** (A- 202A_376G); also compatible with Deficient with CNSHA if the untyped positions are unknown rather than reference
- `7b560d699828` G6PD / primaquine: called **Variable** (A- 202A_376G/B (reference) (equivalent alternatives: A/Asahi)); also compatible with Deficient if the untyped positions are unknown rather than reference

**(b) CYP2D6 structural variation.** A SNV panel cannot see the *5 whole-gene deletion or a
duplication, so every CYP2D6 (and CYP2C19/CYP2B6/TPMT/NUDT15/SLCO1B1) call assumes two intact
copies. The 5 `*1xN` duplication cases assume the duplicated allele and the second allele
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

- `4f19dd3f5283` UGT1A1 / irinotecan: no CPIC Level A recommendation text exists for UGT1A1-irinotecan; phenotype is reported but the recommendation is abstained
- `912145392abc` UGT1A1 / irinotecan: no CPIC Level A recommendation text exists for UGT1A1-irinotecan; phenotype is reported but the recommendation is abstained
- `6eeeaffed889` UGT1A1 / irinotecan: no CPIC Level A recommendation text exists for UGT1A1-irinotecan; phenotype is reported but the recommendation is abstained
- `e55c1b0dd1a0` CYP4F2 / warfarin: CPIC warfarin dosing is algorithm-based and indexed on VKORC1 and CYP2C9; CYP4F2 genotype alone does not yield a dose recommendation (no CYP4F2 phenotype table and no CYP4F2 recommendation rows in the CPIC API; the supplied corpus excerpt contains no CYP4F2 dosing statement)
- `807ece818520` CYP4F2 / warfarin: CPIC warfarin dosing is algorithm-based and indexed on VKORC1 and CYP2C9; CYP4F2 genotype alone does not yield a dose recommendation (no CYP4F2 phenotype table and no CYP4F2 recommendation rows in the CPIC API; the supplied corpus excerpt contains no CYP4F2 dosing statement)

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
| `claude_science_answers.csv` | the 110 answers: case_id, diplotype, phenotype, recommendation, abstained |
| `case_audit.csv` | per-case evidence: candidate diplotypes, CPIC lookup key, full recommendation, abstention reason, lethal class and mechanism, sensitivity result, sources, notes |
| `table_variant_to_allele.csv` | rsID -> star allele definition table used (CPIC plus-strand bases, allele function, activity value) |
| `table_diplotype_to_phenotype.csv` | diplotype -> CPIC phenotype/activity score, flagged where used by a case |
| `table_phenotype_to_recommendation.csv` | phenotype (or CPIC lookup key) -> CPIC recommendation, per study gene-drug pair |
| `table_sources_by_gene.csv` | source and version per gene |
| `pgx_pipeline.py` | the pipeline that produced the answers and tables |
| `fetch_cpic_api.py` | CPIC API snapshot fetcher (writes `cpic_api_snapshot/`) |
| `make_methods_note.py` | generates this note and `table_sources_by_gene.csv` |
| `claude_science_analysis.ipynb` | notebook: runs the pipeline and shows the audit/summary views |
