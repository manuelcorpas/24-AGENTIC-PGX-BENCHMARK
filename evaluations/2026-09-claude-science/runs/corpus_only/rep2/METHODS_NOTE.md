# Methods note - blinded pharmacogenomics study

Run on 110 blinded cases (`cases_blinded.json`) against the supplied CPIC
guideline corpus (`cpic_guideline_corpus.json`). No network source was used at any point; no
published benchmark, paper or repository reporting answers to these cases was consulted or
encountered. No network access was requested or granted in this session; the only permission
requested was read/write access to the local folder holding the two input files.

## 1. Approach

A deterministic rule-based pipeline (`pgx_pipeline.py`), not case-by-case free-text reasoning.
Every case goes through the same four stages, and the intermediate state of each stage is
written out for audit (`case_audit_trail.csv`).

1. **Parse.** The free-text genotype description is reduced either to per-locus variant
   dosages (0, 1 or 2 copies of the forward-genomic variant base; 1 for a hemizygous X call),
   or to a non-SNV descriptor: copy-number call, HLA typing call, mitochondrial variant,
   TA-repeat length, or "all reference".
2. **Assemble a diplotype.** Per-gene resolvers combine locus dosages into a star-allele (or
   haplotype) pair using the allele definitions in `allele_definitions.csv`, honouring the cis
   relationships that make a naive per-locus count wrong (below).
3. **Translate to a phenotype.** CPIC activity scores where the gene uses them (CYP2D6,
   CYP2C9, DPYD), allele-function pairs otherwise. 60 distinct gene/diplotype
   combinations were reached; all are listed in `diplotype_to_phenotype.csv`.
4. **Look up the recommendation.** `phenotype_to_recommendation.csv` is a table transcribed
   from the corpus guideline text, keyed on gene, drug and phenotype. Each row records which
   corpus source (PharmGKB annotation ID and issuing body) it came from.

## 2. Sources and versions

| item | value |
|---|---|
| corpus schema | v3.1 |
| corpus generated (UTC) | 2026-05-10T04:40:15Z |
| guideline snapshot date stated in brief | 2026-05-09 |
| genes in corpus | 21 |
| `cases_blinded.json` sha256 | `bcda9e79e870a15ed378c4f9d7c97e42c52db4984599c4015dd1ee73b9b6ae1c` |
| `cpic_guideline_corpus.json` sha256 | `f5ef71428715cacaeaed5f593c35fe04796732923e685b5366d725ec8f13b232` |

All recommendation text derives from the corpus. Every gene-drug pair in the case set is
covered by the corpus, with one qualification: **UGT1A1-irinotecan has no CPIC guideline in
the corpus** (only DPWG PA166104951 and RNPGx PA166127626 annotations). Those three
recommendations are therefore not CPIC Level A, and each one says so in its text.

Allele *definitions* - which rsID defines which star allele, and which forward-genomic base
is the variant - are not in the corpus. They come from standard pharmacogene nomenclature and
are tabulated explicitly in `allele_definitions.csv` so they can be checked independently.

## 3. Interpretation decisions that change answers

These are the points where a defensible alternative reading gives a different result. Each is
implemented in one place in the code and is easy to flip.

- **Strand.** Six genes in this case set sit on the minus strand (CYP2D6, TPMT, DPYD, G6PD,
  CYP4F2, and VKORC1 by cDNA convention), and the cases give forward-genomic bases, as they
  state. Reading `rs1142345 T/C` or `rs1800460 C/T` as cDNA rather than forward-genomic would
  invert the TPMT calls. Forward-genomic reference/variant bases are recorded per rsID.
- **CYP2D6 activity score 1.0 is an intermediate metabolizer.** The corpus codeine and
  paroxetine tables print the CPIC standardised bands (UM >2.25, NM 1.25-2.25, IM 0 < AS
  < 1.25, PM 0), whereas the corpus ondansetron, amitriptyline and tamoxifen tables still
  print the superseded band that put AS 1.0 in "normal/extensive metabolizer". The
  standardised bands were applied to all five CYP2D6 drugs, because each of those three
  annotations carries the note that CPIC has updated the genotype-phenotype translation in all
  supplemental files since publication. This affects 5 cases
  (b07cff8c99c0, 3d0e97b86ad0, 9d75134030af, f73ec8cf3d6f, abf51805207b); under the old translation they would be normal
  metabolizers on standard dosing.
- **cis relationships.** `rs3892097 C/T + rs1065852 G/A` is one *4 allele, not two reduced
  alleles: c.100C>T is carried on the *4 haplotype, and *4/*10 is excluded because it would
  make rs1065852 homozygous. Likewise TPMT c.460G>A + c.719A>G in cis is a single *3A, and
  G6PD c.202G>A + c.376A>G in cis is a single A- allele. A per-locus additive score would call
  all three poor metabolizers.
- **CYP2D6 *1xN duplications** (5 cases) were called *1/*1xN, activity score 3,
  ultrarapid. This assumes the non-duplicated allele is normal function and that N = 2; a
  copy-number assay does not phase the duplication to a star allele. If the second allele were
  in fact no-function, the same patient would be a normal metabolizer on standard dosing. The
  assumption is recorded per case in the audit trail.
- **CYP2B6 c.516G>T without c.785A>G** was reported as *6 by panel convention; strictly it is
  *9. Both are decreased function, so the phenotype and recommendation are unchanged.
- **SLCO1B1** rs4149056 alone cannot separate *5 from *15 (which also carries c.388A>G); *15
  is reported by convention. Function assignment is identical.
- **Single-gene warfarin cases.** The corpus dose table is indexed jointly on VKORC1 -1639 and
  CYP2C9, and the CYP2C9-only and VKORC1-only cases supply just one axis. Rather than assume
  the missing genotype, the recommendation gives the dose for each possible value of the other
  gene and states that it was not supplied. Only CYP2C9 *3/*3 is unconditional (0.5-2 mg/day
  at every VKORC1 genotype).
- **G6PD ploidy.** G6PD is X-linked. Where a case does not fix sex, the diplotype is given in
  both forms (hemizygous male / homozygous female); the predicted phenotype is the same.
  The one heterozygous-female case is "Variable", which is a distinct CPIC phenotype requiring
  an enzyme activity measurement, not a deficiency call.
- **RYR1 negative** is reported as CPIC's "uncertain susceptibility", not as normal: roughly
  half of contracture-test-confirmed malignant hyperthermia has no identified variant.
- **MT-RNR1 negative** rests on position m.1555 only; m.1095T>C and m.1494C>T were not
  reported and also confer increased risk.

## 4. Abstentions (2 of 110)

| case_id | gene | drug | diplotype | phenotype |
|---|---|---|---|---|
| e55c1b0dd1a0 | CYP4F2 | warfarin | *3/*3 | Decreased CYP4F2 function (*3/*3) |
| 807ece818520 | CYP4F2 | warfarin | *1/*1 | Normal CYP4F2 function |

Reason: the supplied CPIC warfarin excerpt indexes its dose table on VKORC1 -1639 and CYP2C9
only. It contains no CYP4F2 row and no CYP4F2-specific dosing statement, so a CYP4F2-only
genotype yields no corpus-derivable recommendation. The diplotype and the direction of the
functional effect are still reported; the recommendation is withheld. (The published CPIC
warfarin guideline does carry a CYP4F2 statement, but it is not in this corpus snapshot and
the brief forbids fetching it.)

No other case was abstained on. In particular, the deliberately awkward inputs - copy-number
duplications with no defining SNV, compound cis genotypes, the heterozygous G6PD female, the
HLA tag SNP rather than direct typing - are all resolvable with a stated assumption, and a
stated assumption is more useful here than an abstention. Every such assumption is in the
`note` column of `case_audit_trail.csv`.

## 5. Lethal-class cases (29 of 110)

"Lethal class" means standard dosing in this genotype carries a documented risk of a fatal or
immediately life-threatening event, so the genotype-guided action is the difference between
routine care and a death. This is the set where an error in this pipeline would matter most.

| case_id | gene | drug | phenotype | why |
|---|---|---|---|---|
| 74597c5503ec | CYP2D6 | codeine | Ultrarapid metabolizer (activity score 3) | Ultrarapid morphine formation; fatal respiratory depression reported in children and in breastfed infants of UM mothers. |
| ca9e14feca53 | CYP2C19 | clopidogrel | Poor metabolizer | Loss of antiplatelet effect; stent thrombosis and fatal myocardial infarction. |
| 97d4e4559a13 | CYP2C19 | clopidogrel | Intermediate metabolizer | Reduced active metabolite; increased risk of stent thrombosis and fatal cardiovascular events. |
| 1fb04aa9e16c | CYP2C19 | voriconazole | Rapid metabolizer | Subtherapeutic exposure; failure of therapy for invasive fungal infection. |
| 0c21c969a5c4 | CYP2C19 | voriconazole | Ultrarapid metabolizer | Subtherapeutic exposure; failure of therapy for invasive fungal infection. |
| 97c3baebed4e | CYP2C9 | warfarin | Poor metabolizer (activity score 0) | Standard 5 mg induction in *3/*3 causes marked over-anticoagulation and major haemorrhage. |
| 9b5e3093ab14 | CYP2C9 | warfarin | Intermediate metabolizer (activity score 1) | Standard induction over-anticoagulates *1/*3; major bleeding risk. |
| b94a43e574db | DPYD | fluorouracil | Poor metabolizer (activity score 0) | Complete DPD deficiency; standard 5-FU is reliably fatal without dose avoidance. |
| c57caa44a38a | DPYD | fluorouracil | Intermediate metabolizer (activity score 1) | Partial DPD deficiency; severe or fatal myelosuppression and mucositis at full dose. |
| 6b8d2db91001 | DPYD | fluorouracil | Intermediate metabolizer (activity score 1.5) | Partial DPD deficiency; severe or fatal myelosuppression and mucositis at full dose. |
| c6549d637d49 | DPYD | capecitabine | Poor metabolizer (activity score 0) | Complete DPD deficiency; standard capecitabine is reliably fatal without dose avoidance. |
| ae4a7215ae64 | DPYD | capecitabine | Intermediate metabolizer (activity score 1) | Partial DPD deficiency; severe or fatal toxicity at full dose. |
| 08e18e809fcb | TPMT | azathioprine | Poor metabolizer | Life-threatening myelosuppression at conventional thiopurine doses. |
| d00994ea4b9e | TPMT | azathioprine | Intermediate metabolizer | Substantial risk of severe myelosuppression at conventional doses. |
| 13d1cd0cfa0f | TPMT | azathioprine | Intermediate metabolizer | Substantial risk of severe myelosuppression at conventional doses. |
| f3dd19d8a337 | NUDT15 | thiopurines | Poor metabolizer | Severe, life-threatening myelosuppression; tolerated dose is roughly 10% of normal. |
| 753560fff22e | NUDT15 | thiopurines | Intermediate metabolizer | Increased risk of severe myelosuppression at conventional doses. |
| d1a12c4a16c4 | SLCO1B1 | simvastatin | Poor function | Markedly increased simvastatin acid exposure; rhabdomyolysis can be fatal. |
| 4f19dd3f5283 | UGT1A1 | irinotecan | Poor metabolizer | Life-threatening neutropenia and diarrhoea at full dose (DPWG lists death among the clinical effects). |
| e70379d3dbac | HLA-B*57:01 | abacavir | High risk of abacavir hypersensitivity | Abacavir hypersensitivity syndrome; fatal on re-challenge. |
| 7181a24a65d5 | HLA-B*15:02 | carbamazepine | Carrier of HLA-B*15:02: increased risk of severe cutaneous adverse reaction | SJS/TEN with mortality up to ~30% for TEN. |
| 569a842e1a9a | HLA-B*15:02 | oxcarbazepine | Carrier of HLA-B*15:02: increased risk of severe cutaneous adverse reaction | SJS/TEN with substantial mortality. |
| 85976e86289a | HLA-A*31:01 | carbamazepine | Carrier of HLA-A*31:01: increased risk of severe cutaneous adverse reaction | SJS/TEN and DRESS, both potentially fatal. |
| ebf82c3bcfdd | HLA-B*58:01 | allopurinol | Significantly increased risk of allopurinol SCAR | Allopurinol-induced SCAR (DRESS/SJS/TEN) with reported mortality up to ~25%. |
| d49e9f88d9c2 | G6PD | rasburicase | Deficient (WHO class II-III) | Acute intravascular haemolysis and methaemoglobinaemia; rasburicase is contraindicated. |
| 3786ab9a929d | G6PD | rasburicase | Deficient (WHO class II-III) | Acute intravascular haemolysis and methaemoglobinaemia; rasburicase is contraindicated. |
| 2a8544a44d5f | G6PD | primaquine | Deficient (WHO class II-III) | Acute haemolytic anaemia, which can be severe enough to be fatal at standard 14-day dosing. |
| c1e4c84b49dd | VKORC1 | warfarin | Increased warfarin sensitivity (low dose requirement) | Low dose requirement; standard induction causes over-anticoagulation and major haemorrhage. |
| 1229a4123378 | RYR1 | volatile-anaesthetics | Malignant hyperthermia susceptible | Malignant hyperthermia; untreated mortality approaches 70%. |

A second tier (10 cases) carries serious or irreversible harm that is
not usually fatal, and is flagged `severe` rather than `lethal` in
`case_audit_trail.csv`: aminoglycoside-induced permanent hearing loss in the MT-RNR1 carrier,
statin myopathy short of rhabdomyolysis, phenytoin and voriconazole accumulation, citalopram
QT risk, tricyclic accumulation, tamoxifen under-treatment, ineffective ivacaftor in an
F508del homozygote, and the one G6PD variable-phenotype case. (The recommendation lookup
carries a `severe` row for a G6PD variable phenotype on rasburicase as well, but no case in
this set has that combination.)

| case_id | gene | drug | phenotype | why |
|---|---|---|---|---|
| f9e4a13551fd | CYP2D6 | tamoxifen | Poor metabolizer (activity score 0) | Reduced endoxifen; higher breast-cancer recurrence risk (harm from under-treatment). |
| d456dcf2b837 | CYP2D6 | amitriptyline | Poor metabolizer (activity score 0) | TCA accumulation with cardiac conduction and anticholinergic toxicity. |
| 00b16db7406d | CYP2C19 | voriconazole | Poor metabolizer | Supratherapeutic troughs with hepatotoxicity and neurotoxicity. |
| 7b7e74f7aa47 | CYP2C19 | citalopram | Poor metabolizer | Higher exposure with QT prolongation risk; the FDA caps the dose at 20 mg/day. |
| e85ecfa887cb | CYP2C9 | phenytoin | Poor metabolizer (activity score 0) | Phenytoin accumulation causing ataxia, encephalopathy and nystagmus; rarely fatal. |
| 51a382a7ca5c | SLCO1B1 | simvastatin | Decreased function | Increased myopathy risk. |
| 9a62df2e13c7 | SLCO1B1 | atorvastatin | Poor function | Increased atorvastatin exposure and myopathy risk. |
| 7b560d699828 | G6PD | primaquine | Variable | X-inactivation mosaicism means a heterozygous female can be functionally deficient; standard 14-day primaquine could cause significant haemolysis if enzyme activity is not measured first. |
| 18325551d245 | CFTR | ivacaftor | Unlikely to respond to ivacaftor (F508del homozygous) | No benefit; harm is the opportunity cost of ineffective therapy. |
| 369c2bcff713 | MT-RNR1 | aminoglycosides | Increased risk of aminoglycoside-induced hearing loss | Permanent bilateral sensorineural hearing loss after even a single dose; irreversible but not fatal. |

Note the asymmetry the pipeline has to get right: for CYP2D6 and codeine the ultrarapid
metabolizer is the lethal case (opioid toxicity) while the poor metabolizer is only an
efficacy failure, whereas for CYP2C19 and clopidogrel it is the poor metabolizer that is
lethal (thrombosis). A pipeline that conflates "poor metabolizer" with "danger" gets one of
those two backwards.

## 6. Where this is most likely wrong

In descending order of my own expectation of error:

1. **The CYP2D6 activity-score-1.0 band (5 cases).** If the intended answer set was
   generated from the older tables that are still printed verbatim in the corpus for
   ondansetron, amitriptyline and tamoxifen, these should be normal metabolizers on standard
   dosing rather than intermediate metabolizers. I judged the standardised translation
   correct, but this is a genuine fork and it is the largest single block of cases at risk.
2. **The *1xN duplication cases (5 cases).** Ultrarapid is the clinically safe call
   and the conventional one, but it is an inference from an unphased copy-number result, not a
   measurement. If the intended answer is an abstention, these are wrong - and the codeine one
   is in the lethal class.
3. **UGT1A1-irinotecan (3 cases).** The corpus has no CPIC guideline, so the action comes from
   DPWG. If the expected answer is an abstention (or a CPIC recommendation from outside this
   snapshot), these differ.
4. **CYP4F2 (2 cases).** I abstained on a gene-drug pair that CPIC does address in the full
   guideline. If the grader expects the published "consider a modestly higher warfarin dose in
   *3 carriers of European ancestry", an abstention scores as a miss.
5. **Recommendation wording granularity.** The corpus frequently carries two or three tables
   per gene-drug pair from different guideline vintages, and sometimes separate adult and
   paediatric recommendations that genuinely disagree - CYP2C19 rapid metabolizers on
   voriconazole get "choose an alternative agent" as adults and "standard dose plus therapeutic
   drug monitoring" as children. Since no case states an age, adult recommendations were used
   and the paediatric divergence is named in the text. Any string-matching comparison against
   a single expected phrasing will penalise this even where the clinical action agrees.
6. **Diplotype nomenclature.** `*1/*15` versus `*1/*5`, `*6/*6` versus `*9/*9`,
   `c.[1905+1G>A];[=]` versus `*2A/*1`, and the ploidy-explicit G6PD strings are all
   defensible renderings of the same call. The phenotype column is the more robust comparison
   point.
7. **Single-gene warfarin cases (5 cases).** Reporting a conditional dose range is right but
   unusual; a grader expecting one number will not match it.

What I am confident about: the HLA, MT-RNR1, RYR1, CFTR, IFNL3, G6PD, DPYD, TPMT, NUDT15 and
SLCO1B1 calls, and the direction of every lethal-class recommendation. The CYP3A5 direction in
particular is a common inversion - the non-expresser (*3/*3) takes the standard tacrolimus
dose and the expresser needs 1.5-2x - and it is right in this output.

## 7. Files

| file | contents |
|---|---|
| `claude_science_answers.csv` | the 110 answers: case_id, diplotype, phenotype, recommendation, abstained |
| `case_audit_trail.csv` | per-case intermediate state: parsed kind, diplotype, activity score, phenotype key, harm class and reason, and the assumption note |
| `allele_definitions.csv` | rsID to star allele, cDNA change, forward-genomic reference and variant bases, strand, function, CPIC activity value |
| `diplotype_to_phenotype.csv` | every gene/diplotype reached, with activity score and phenotype |
| `phenotype_to_recommendation.csv` | the full recommendation lookup, with harm class and corpus source per row |
| `pgx_pipeline.py` | the pipeline; `python pgx_pipeline.py <dir>` regenerates every CSV |
| `pgx_analysis.ipynb` | notebook that runs the pipeline and summarises the results |
| `build_methods_note.py` | regenerates this note from the outputs |
| `run_provenance.json` | corpus version, input checksums, and summary counts |
