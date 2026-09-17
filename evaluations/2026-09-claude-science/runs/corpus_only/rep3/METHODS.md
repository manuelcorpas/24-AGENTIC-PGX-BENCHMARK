# Methods note — blinded pharmacogenomics study (`corpus_only/rep3`)

110 blinded cases, 21 genes, 32 drugs, 35 gene–drug pairs. All 110 cases were answered;
2 carry an abstention on the recommendation field. No network access was used or
requested at any point, and no published benchmark, paper or repository reporting
answers to these cases was searched for or consulted.

---

## 1. Approach

A deterministic three-table rule pipeline (`pgx_pipeline.py`), no statistical model and
no free-text generation in the answer path:

1. **Assay result → allele call** (`table_variant_to_allele.csv`, 61 rules).
   Each case `genotype` string is matched by a per-gene regex to an explicit allele/haplotype
   call. Every rule carries a `basis` field recording where the allele definition came from:
   `corpus` (14 rules), `case-annotation` (10 rules, where the case string itself states the
   allele or the typing result), or `analyst-knowledge` (37 rules). Across the 110 cases the
   basis split is 23 / 16 / 71.
2. **Diplotype → phenotype** (`table_diplotype_to_phenotype.csv`, 60 rules), with the CPIC
   activity score recorded where CPIC uses one (CYP2D6, CYP2C9, DPYD).
3. **(gene, drug, phenotype) → recommendation** (`table_phenotype_to_recommendation.csv`,
   109 rules). Each rule holds a short drug-prefixed action (what goes in the answers file)
   **and** the verbatim guideline text it was derived from. A run-time assertion checks that
   every verbatim string is a whitespace-normalised substring of the supplied corpus, so the
   pipeline cannot run if any recommendation text was misremembered or invented.

`claude_science_answers.csv` is the required five-column output.
`answers_audit_trail.csv` is the same 110 rows plus activity score, allele-call basis and
reasoning, guideline source, verbatim guideline text, and the lethal-class flag and rationale.
`pgx_analysis.ipynb` drives `pgx_pipeline.py` and displays all three tables and the summaries.

### Why part of the mapping cannot be corpus-grounded

The corpus is guideline *recommendation* text, not allele-definition data. Of the 24 rsIDs
used by the cases, only 6 appear anywhere in it (rs1057910, rs4149056, rs9923231, rs12979860,
rs75527207, rs113993960), and the only allele-definition table present is CYP2C9's
(`*3 | A>C at rs1057910 | I359L | Decreased`). The rsID → star-allele step therefore comes
from the analyst's own pharmacogenetic knowledge for 71 of 110 cases. This is disclosed
per rule rather than hidden: `table_variant_to_allele.csv` marks every such call
`analyst-knowledge` and states the reasoning. The corpus *is* sufficient for the
phenotype → recommendation step for 33 of 35 gene–drug pairs.

## 2. Sources and versions

- `cpic_guideline_corpus.json`, schema `v3.1`, `generated_at_utc` 2026-05-10T04:40:15Z,
  PharmGKB bundle `20260509T203526Z` (snapshot 2026-05-09). 21 genes, 157 kB of text.
- `cases_blinded.json`, sha256 prefix `bcda9e79e870a15e`, 110 cases.
- Guideline body per recommendation is recorded in the `guideline_source` column
  (PharmGKB PA accession + issuing body). 33 of 35 gene–drug pairs are CPIC.
  **UGT1A1/irinotecan is not**: the corpus carries only DPWG (PA166104951) and RNPGx
  (PA166127626) annotations for it, no CPIC guideline. Those three answers therefore give
  the DPWG recommendation, explicitly labelled as DPWG in the answer string.
- Where the corpus prints both a current and a superseded table for the same gene–drug pair
  (CYP2D6/codeine, CYP2D6/paroxetine, CYP2C19/citalopram, CYP2C19/clopidogrel,
  SLCO1B1/simvastatin), the current table was used.

### Assumptions forced by single-gene input

Several CPIC recommendations are multi-gene, but each case supplies one gene. The
assumption is stated inside the answer string in every affected case:

| Cases | Missing input | Assumption |
|---|---|---|
| CYP2C9/warfarin (3) | VKORC1 | VKORC1 −1639 GG; dose read from the corpus VKORC1 × CYP2C9 table row `GG` |
| VKORC1/warfarin (2) | CYP2C9 | CYP2C9 \*1/\*1; dose read from table row `AA` / `GG` |
| CYP2C9/phenytoin (3) | HLA-B\*15:02 | HLA-B\*15:02 negative/untested (the guideline's other branch says "do not use phenytoin" irrespective of CYP2C9) |
| CYP2D6/amitriptyline (4) | CYP2C19 | CYP2C19 normal metaboliser row of the CYP2C19 × CYP2D6 matrix |
| HLA-B\*15:02-negative carbamazepine (1) | HLA-A\*31:01 | also non-carrier/untested |
| HLA-A\*31:01 carbamazepine (2) | HLA-B\*15:02 | not reported |
| CYP2C19/voriconazole (5) | age | adult table (the corpus gives separate adult and paediatric recommendations that differ for rapid metabolisers) |

### CYP2D6 activity-score translation — a deliberate, consequential choice

The corpus contains **two incompatible CYP2D6 activity-score binnings**. The updated
codeine and paroxetine tables use AS > 2.25 UM / 1.25 ≤ AS ≤ 2.25 NM / 0 < AS < 1.25 IM /
0 PM, whereas the older ondansetron, amitriptyline and tamoxifen tables put AS 1.0–2.0 in
"normal metabolizer". The corpus itself states, for tamoxifen and ondansetron, that "CPIC
has updated the CYP2D6 genotype-phenotype translation in all supplemental files since the
guideline was published". I applied the current translation uniformly, so **\*1/\*4
(AS 1.0) is an intermediate metaboliser for all five CYP2D6 drugs**. This changes the
phenotype label on all 5 `*1/*4` cases relative to the older binning, and the *recommendation*
on 2 of them (`f73ec8cf3d6f` paroxetine: lower starting dose + slower titration instead of
standard dose; `abf51805207b` amitriptyline: 25% dose reduction instead of standard dose).
For the other 3 (`b07cff8c99c0` and `3d0e97b86ad0` codeine, `9d75134030af` tamoxifen) the
action is the same under either binning — and the corpus tamoxifen table itself labels
AS 1 as "Normal metabolizer or Intermediate metabolizer (controversy remains)". This is the
single largest source of systematic disagreement if the answer key was built from the
printed per-drug tables instead of the current translation.

## 3. Lethal-class cases (22 of 110)

Criterion: standard (label) dosing of the named drug, in the called phenotype, can
plausibly kill the patient or cause irreversible major organ/tissue injury — not merely
treatment failure or tolerable toxicity. Per-case rationale is in
`answers_audit_trail.csv` (`lethal_class`, `lethal_rationale`).

| Case | Gene / drug | Phenotype | Mechanism of lethal harm |
|---|---|---|---|
| 74597c5503ec | CYP2D6 / codeine | ultrarapid metaboliser | Runaway morphine formation; documented fatal and near-fatal respiratory depression at label doses |
| b94a43e574db, c6549d637d49 | DPYD / fluorouracil, capecitabine | poor metaboliser (AS 0) | Complete DPD deficiency; corpus language is "severe or even fatal drug toxicity" |
| c57caa44a38a, 6b8d2db91001, ae4a7215ae64 | DPYD / fluorouracil, capecitabine | intermediate metaboliser (AS 1 / 1.5) | Partial DPD deficiency; the same corpus sentence applies, and full-dose 5-FU in heterozygotes has caused grade 4–5 toxicity |
| 08e18e809fcb | TPMT / azathioprine | poor metaboliser (\*3A/\*3A) | Profound prolonged myelosuppression at conventional dose |
| f3dd19d8a337 | NUDT15 / thiopurines | poor metaboliser (\*3/\*3) | Severe myelosuppression; tolerated dose ≈10% of normal |
| 1229a4123378 | RYR1 / volatile anaesthetics | MH susceptible | Malignant hyperthermia crisis; rapidly fatal without dantrolene |
| d49e9f88d9c2, 3786ab9a929d | G6PD / rasburicase | deficient (class II, class III) | Contraindicated; severe acute haemolysis and methaemoglobinaemia |
| 2a8544a44d5f | G6PD / primaquine | deficient (class III) | Acute haemolytic anaemia on 14-day dosing |
| 7181a24a65d5, 569a842e1a9a | HLA-B\*15:02 / carbamazepine, oxcarbazepine | carrier | SJS/TEN with substantial mortality |
| 85976e86289a | HLA-A\*31:01 / carbamazepine | carrier | DRESS / SJS / TEN |
| ebf82c3bcfdd | HLA-B\*58:01 / allopurinol | carrier | SJS/TEN/DRESS with substantial mortality |
| e70379d3dbac | HLA-B\*57:01 / abacavir | positive | Systemic hypersensitivity; fatal on rechallenge |
| ca9e14feca53 | CYP2C19 / clopidogrel | poor metaboliser | Loss of antiplatelet effect; stent thrombosis, MI, cardiovascular death |
| 97c3baebed4e | CYP2C9 / warfarin | poor metaboliser (\*3/\*3) | Stable dose 0.5–2 mg/day; standard initiation causes major or fatal haemorrhage |
| d1a12c4a16c4 | SLCO1B1 / simvastatin | poor function (521 CC) | Rhabdomyolysis with renal failure |
| 4f19dd3f5283 | UGT1A1 / irinotecan | poor metaboliser (\*28/\*28) | Corpus: "serious, life-threatening adverse events occur more often" — febrile neutropenia, severe diarrhoea |
| 369c2bcff713 | MT-RNR1 / aminoglycosides | increased risk (m.1555A>G) | Not fatal, but permanent bilateral deafness after a single standard course — included as irreversible major harm |

Borderline calls I made explicitly: DPYD **heterozygotes** are counted as lethal-class
because the corpus applies the "severe or even fatal" language to the intermediate
metaboliser row as well; CYP2D6 **poor** metabolisers on codeine are *not* (the harm is
analgesic failure); CYP2C9 poor metabolisers on phenytoin and CYP2B6 poor metabolisers on
efavirenz are *not* (dose-dependent toxicity, reversible); tamoxifen poor metabolisers are
*not* (recurrence risk, not acute harm).

## 4. Abstentions (2 of 110)

`e55c1b0dd1a0` and `807ece818520` — **CYP4F2 / warfarin** (rs2108622 T/T and C/C).
The diplotype (\*3/\*3, \*1/\*1) and function phenotype are reported, but the
recommendation field abstains. Reason: the supplied CYP4F2 excerpt is the warfarin
annotation whose only dose table is keyed on VKORC1 × CYP2C9; it names CYP4F2 as an input
gene but contains **no CYP4F2 genotype-specific dose modifier**, and `rs2108622` does not
appear anywhere in the 157 kB corpus. A CYP4F2-only result cannot be translated into an
action from the supplied text, and brief rule 1 directs abstention over fetching. (I do
know the published CPIC position — a small upward dose adjustment for \*3 carriers of
non-African ancestry — but that figure is not in the corpus, so it is not asserted here.)

Cases I considered abstaining on and answered anyway, with the caveat carried in the
answer or the audit trail: the 4 CYP2D6 `*1xN` duplications (copy number N not given, but
a CNV assay result *was* given and the corpus lists `*1/*1xN` as its own example UM
diplotype); the CYP2D6 compound case (phase not measured); the G6PD heterozygous female
("Variable" is itself the correct CPIC phenotype, and "measure enzyme activity" is the
correct CPIC action); the two HCP5-tag abacavir cases.

## 5. Where I am most likely wrong

Ordered by my estimate of probability × consequence.

1. **CYP2D6 AS 1.0 binning (5 cases, 2 of them on the recommendation).** If the key follows
   each drug's printed table rather than the current CPIC translation, `*1/*4` should read
   "normal metaboliser" for amitriptyline and paroxetine (and "NM or IM" for tamoxifen),
   and the amitriptyline recommendation would change from "25% dose reduction" to "standard
   starting dose". This is my most likely systematic error.
2. **CYP2D6 `*1xN` (4 cases).** I called `*1/*1xN` → ultrarapid. If the intended answer is
   an abstention (copy number and the non-duplicated allele are genuinely unspecified), or
   `*1xN/*1xN`, these 4 are wrong. Note the failure is asymmetric: calling UM triggers
   "avoid codeine", which is the safe direction.
3. **CYP2D6 compound case `3d0e97b86ad0`** (rs3892097 C/T + rs1065852 G/A). I called
   `*1/*4` on cis-parsimony, since 100C>T is part of the `*4` core haplotype. Unphased data
   formally permits a trans arrangement (a `*4`-like / `*10`-like diplotype, AS 0.25),
   which would make this a poor-to-low intermediate metaboliser and change the codeine
   action from "use label dosing" to "avoid".
4. **CYP2B6 `*6` vs `*9`, and `*4` vs `*6` (2 cases).** For `rs3745274 T/T` I reported
   `*6/*6`; with rs2279343 reference it is strictly `*9/*9`. The phenotype is identical
   (both decreased function), so only the diplotype string is at risk. Conversely
   `rs2279343 A/G` alone is `*1/*4` — an *increased*-function allele, rapid metaboliser,
   standard 600 mg/day. If the key treated 785A>G as a `*6` marker it would call
   intermediate metaboliser and a reduced dose; I believe my call is right but this is the
   case where a plausible alternative reading inverts the direction of the answer.
5. **TPMT `*3B` vs `*3A` (`d00994ea4b9e`).** rs1800460 heterozygous alone is `*3B`; a key
   that assumes the commoner `*3A` haplotype would disagree on the diplotype string.
   Phenotype (intermediate metaboliser) is unaffected.
6. **SLCO1B1 diplotype strings.** Without rs2306283 (388A>G), `*5` cannot be separated from
   `*15`/`*17`. I reported `*5/*5` and `*1/*5` with the ambiguity noted; a key using
   `*15/*15` / `*1/*15` would differ on the string only, not the function phenotype.
7. **Warfarin dose numbers (5 cases).** These depend on the assumed second gene. If the key
   expects only "use a genotype-guided algorithm" without numbers, my answers are still
   compatible; if it expects a different second-gene assumption, the mg/day figures differ.
8. **UGT1A1/irinotecan source (3 cases).** If the key expects a CPIC recommendation, no
   such text exists in this corpus and my DPWG-labelled answer will not match. An
   abstention was the alternative; I judged a labelled real guideline more useful than a
   blank.
9. **HCP5 proxy for HLA-B\*57:01 (2 cases).** rs2395029 is a tag SNP, not HLA typing. LD
   with `*57:01` is strong in European ancestry and weaker elsewhere, so the negative case
   (`f17cf45007e4`) carries a residual false-negative risk that my "use abacavir per
   standard dosing" answer does not convey as strongly as a real typing result would.
10. **Phenotype nomenclature.** I used the CPIC terms as printed in this corpus (e.g.
    "CYP3A5 poor metabolizer (CYP3A5 non-expresser)", "SLCO1B1 poor function", "G6PD
    Variable", "Uncertain susceptibility"). Exact-string scoring against differently
    worded labels would penalise answers that are conceptually correct.

Two classes of error the design deliberately avoids: no recommendation text is invented
(the verbatim assertion enforces this), and no case is silently dropped or guessed — an
unmatched genotype string raises rather than falling through to a default.

## 6. Files

| File | Contents |
|---|---|
| `claude_science_answers.csv` | required output: `case_id,diplotype,phenotype,recommendation,abstained`, 110 rows |
| `answers_audit_trail.csv` | the same rows plus activity score, allele-call basis and reasoning, guideline source, verbatim guideline text, lethal-class flag and rationale |
| `table_variant_to_allele.csv` | assay result → allele call, with basis and reasoning (61 rules) |
| `table_diplotype_to_phenotype.csv` | diplotype → phenotype, with activity score (60 rules) |
| `table_phenotype_to_recommendation.csv` | (gene, drug, phenotype) → action + verbatim guideline text + source (109 rules) |
| `pgx_pipeline.py` | the pipeline; running it regenerates every CSV above |
| `pgx_analysis.ipynb` | executed notebook: corpus inventory, the three tables, coverage, abstentions, lethal-class list |
