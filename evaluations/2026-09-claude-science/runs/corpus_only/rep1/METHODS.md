# Methods note — blinded pharmacogenomics study (corpus-only run)

Model: **Claude Opus 5** (`claude-opus-5`), running as Claude Science.
Date of run: 2026-09-15. Working directory: `RUNS/corpus_only/rep1`.

## 1. Approach

The task is a genotype → diplotype → phenotype → recommendation chain, so I built it
as three explicit, separately auditable lookup tables and a deterministic resolver,
rather than answering each case as free text. Every case therefore has a traceable
path from its genotype string to a named table row, and re-running the code
reproduces the answers exactly.

| Stage | Artefact | What it does |
|---|---|---|
| marker → allele | `table_allele_definitions.csv` | rsID / HLA typing result / CNV call → star allele or haplotype, with allele function and the source of the definition |
| genotype → diplotype + phenotype | `table_diplotype_to_phenotype.csv` | per-gene ordered regex rules over the raw case genotype string → diplotype, activity score, CPIC phenotype |
| phenotype → recommendation | `table_phenotype_to_recommendation.csv` | (gene, drug, phenotype) → recommendation wording, guideline strength, PharmGKB annotation id |

Code: `pgx_tables.py` (the tables), `build_answers.py` (loader, resolver, writers),
`analysis.ipynb` (notebook that runs the pipeline and reproduces the summary
figures quoted below). Outputs: `claude_science_answers.csv` (the required
deliverable) and `case_rationale.csv` (one row per case with activity score,
guideline strength, source annotation id, lethal-class flag and the reasoning note
for the diplotype call).

Regex matching rather than string equality was necessary because the same genotype
appears in several textual dresses across the 110 cases (`rs12248560 C/T` and
`rs12248560 C/T (heterozygous)`; two different phrasings of the same CYP2D6
duplication call).

Two automated input checks run before resolution: all 110 cases must have a
non-empty guideline excerpt for their gene, and the named drug must appear in that
gene's `drugs_covered` list. Both passed for all 110 cases.

## 2. Sources and versions

Only the two supplied local files were used:

| File | sha256 |
|---|---|
| `cases_blinded.json` | `bcda9e79e870a15ed378c4f9d7c97e42c52db4984599c4015dd1ee73b9b6ae1c` |
| `cpic_guideline_corpus.json` | `f5ef71428715cacaeaed5f593c35fe04796732923e685b5366d725ec8f13b232` |

Corpus metadata: `schema_version` v3.1, `generated_at_utc` 2026-05-10T04:40:15Z,
PharmGKB bundle `SPECS/cpic_rag_cache/_pharmgkb_bundle/20260509T203526Z`, 21 genes,
per-annotation `retrieved_date` 2026-05-10. Annotation sources are CPIC for every
gene–drug pair in the case set except UGT1A1/irinotecan, which the corpus carries as
DPWG (PA166104951) and RNPGx (PA166127626) annotations — there is no CPIC irinotecan
guideline in the corpus. The NUDT15/TPMT annotation (PA166104933) cites the CPIC
2025 thiopurine update (PMID 41618934, doi 10.1002/cpt.70209).

**No network source was used.** No web search, no API call, no literature fetch,
no lookup of any published benchmark, paper or repository relating to these cases;
nothing of that kind was encountered. The session made one local-filesystem access
request (to read `BRIEF.md` and the two inputs in the granted folder); no network
access request was offered or made.

Where the corpus does **not** contain the information needed for a step, that is
recorded in the `source` column of `table_allele_definitions.csv`. The corpus
contains an explicit allele-definition table for only one gene (CYP2C9: `*2` =
rs1799853 C>T, `*3` = rs1057910 A>C), and direct rsID-indexed phenotype rows for
two more (SLCO1B1 rs4149056 TT/TC/CC; VKORC1 rs9923231 GG/GA/AA in the warfarin
dose table; IFNL3 rs12979860 CC vs CT/TT). Every other marker → star-allele
definition is my own reasoning from standard pharmacogenetic allele nomenclature,
as Rule 3 permits. This is the main place where the answers depend on knowledge
outside the supplied text, and it is where I expect errors to live (§5).

### Reproducing the run

```
python build_answers.py        # writes the answers CSV and all three tables
```

`analysis.ipynb` ships with executed outputs. Its cells were run in-process against
the same module (`exec` in a shared namespace, stdout and the final expression
captured into each cell's outputs) because a Jupyter kernel cannot bind a local
socket in this sandbox. The outputs in the notebook are therefore genuine results
of running the code, produced without a kernel round-trip; re-running the notebook
under a normal Jupyter installation gives the same figures.

The notebook also carries the verification assertions: case ids identical to the
input and in the same order, 110 rows, no nulls or duplicates, no lethal-class case
abstained, and no lethal-class case left on plain standard dosing (0 of 20).

## 3. Results

- 110 cases, 110 rows, no duplicate or missing ids, ids emitted exactly as given.
- **108 answered, 2 abstentions** (both CYP4F2/warfarin).
- 20 cases classified lethal-class.

Abstentions — `e55c1b0dd1a0` (CYP4F2 `rs2108622 T/T`) and `807ece818520`
(CYP4F2 `rs2108622 C/C`). The diplotype is confidently assignable (`*3/*3`,
`*1/*1`), but the supplied warfarin excerpt (PA166104949) contains only the
VKORC1 × CYP2C9 predicted-dose table and the CYP2C9 allele table. It carries no
CYP4F2 allele table, no CYP4F2 phenotype assignment and no CYP4F2-specific dosing
action. Per Rule 1 I abstained rather than supplying the CYP4F2 dose modifier from
memory. The diplotype and the reason are still reported in
`case_rationale.csv`.

## 4. Lethal-class cases

I define lethal-class as: standard, genotype-naive dosing of the named drug in a
patient with that phenotype carries a documented risk of **death or a
life-threatening event** — as distinct from reduced efficacy (CYP2D6 poor
metaboliser on codeine, CFTR F508del on ivacaftor, IFNL3 unfavourable genotype) or
reversible / non-fatal toxicity (MT-RNR1 — permanent hearing loss is grave but not
lethal; SLCO1B1 poor function — myopathy; UGT1A1 poor metaboliser on atazanavir —
jaundice and non-adherence).

20 of 110 cases qualify:

| n | Gene / phenotype | Drug | Basis |
|---|---|---|---|
| 2 | DPYD poor metaboliser (AS 0) | fluorouracil, capecitabine | corpus: "complete DPD deficiency and increased risk for severe or even fatal drug toxicity" |
| 3 | DPYD intermediate metaboliser (AS 1, 1.5) | fluorouracil ×2, capecitabine | corpus: "increased risk for severe or even fatal drug toxicity" |
| 1 | TPMT poor metaboliser (*3A/*3A) | azathioprine | profound, potentially fatal myelosuppression at standard dose (judgement — see below) |
| 1 | NUDT15 poor metaboliser (*3/*3) | thiopurines | same; corpus notes tolerated dose is ~10% of normal (judgement) |
| 1 | CYP2D6 ultrarapid metaboliser | codeine | corpus: "increased formation of morphine leading to higher risk of toxicity" — the mechanism behind the paediatric respiratory-depression deaths |
| 2 | HLA-B*15:02 positive | carbamazepine, oxcarbazepine | corpus: increased risk of SJS/TEN |
| 1 | HLA-A*31:01 positive | carbamazepine | corpus: increased risk of SJS/TEN and DRESS |
| 1 | HLA-B*57:01 positive | abacavir | corpus: significantly increased hypersensitivity risk; fatal on rechallenge |
| 1 | HLA-B*58:01 positive | allopurinol | corpus: significantly increased risk of SCAR |
| 1 | RYR1 MH-susceptible | volatile anaesthetics | malignant hyperthermia is acute and rapidly fatal untreated; corpus: relatively contraindicated |
| 3 | G6PD deficient (A−, Mediterranean) | rasburicase ×2, primaquine | corpus: "high risk of acute hemolytic anemia"; rasburicase contraindicated |
| 1 | UGT1A1 poor metaboliser (*28/*28) | irinotecan | corpus DPWG: "serious, life-threatening adverse events occur more often" |
| 2 | CYP2C19 poor and intermediate metaboliser | clopidogrel | corpus: "increased risk for adverse cardiac and cerebrovascular events" — stent thrombosis and stroke are frequently fatal (judgement) |

Four of these designations go beyond the corpus wording and are flagged as
judgement in the `lethal_basis` column: TPMT and NUDT15 poor metabolisers (the
corpus annotation says only that an alternative immunosuppressant should be
considered), and the two clopidogrel cases (the corpus does not use the word
"fatal"). Including DPYD intermediate metabolisers is a deliberately inclusive
reading — the corpus applies the phrase "severe or even fatal" to the
intermediate row as well as the poor row, and full-dose fluoropyrimidine in an
AS-1 patient has caused deaths.

Every lethal-class case received a protective recommendation: avoid or
dose-reduce, never standard dosing.

## 5. Where I am most likely wrong

Ordered by expected contribution to error.

1. **Strand and orientation of the variant calls (largest risk, and it lands on the
   lethal-class genes).** Several case strings advertise a forward-genomic
   orientation for a gene transcribed on the minus strand — CYP2D6 rs3892097,
   TPMT rs1800460/rs1142345, G6PD rs1050828/rs1050829/rs5030868. I read all of
   them as "variant allele = the non-reference allele named in the string", so
   `rs3892097 T/T` → `*4/*4` poor metaboliser and `rs1800460 T/T + rs1142345 C/C`
   → `*3A/*3A` poor metaboliser. If any of these orientations is inverted
   relative to what the case author intended, the phenotype flips between normal
   and poor and the recommendation inverts with it. I consider my reading correct,
   but the corpus contains no allele-definition table to check it against, so
   this is unverifiable from the supplied material.

2. **CYP2D6 activity score 1.0 → intermediate or normal metaboliser (3–5 cases).**
   The supplied corpus is internally inconsistent here. The codeine and paroxetine
   tables use the current translation (intermediate metaboliser = 0 < AS < 1.25,
   so `*1/*4` is IM); the amitriptyline, ondansetron and tamoxifen tables list
   `*1/*4` explicitly as a *normal* metaboliser example. Those three annotations
   each state that CPIC has since updated the CYP2D6 genotype–phenotype
   translation, so I applied the updated translation (IM) uniformly to all five
   CYP2D6 drugs. If the reference answers follow each drug's table as printed, the
   tamoxifen, ondansetron and amitriptyline `*1/*4` cases are wrong on phenotype.
   The recommendation consequence is small for ondansetron (standard dose either
   way) and amitriptyline (standard dose vs 25% reduction) but material for
   tamoxifen: standard tamoxifen 20 mg/day versus switching to an aromatase
   inhibitor.

3. **Gene-siloed warfarin cases (7 cases).** The CPIC warfarin recommendation is a
   joint CYP2C9 × VKORC1 dose table, but the cases supply one locus at a time
   (3 CYP2C9, 2 VKORC1, 2 CYP4F2). I abstained on CYP4F2 and gave *conditional*
   doses for the CYP2C9 and VKORC1 cases, stating the missing locus in the
   recommendation text. A grader expecting a single dose number, or expecting the
   unstated locus to be assumed reference, will score these differently. The
   table's ancestry-specific guidance and rs12777823 are not in the excerpt at all.

4. **Diplotype strings where the panel cannot resolve phase or copy number.**
   `All <gene> SNPs homozygous reference` cannot exclude the CYP2D6 `*5`
   whole-gene deletion or a duplication, so my `*1/*1` calls there are
   panel-level, not true diplotypes. CYP2B6 `rs3745274 T/T` is `*6/*6` or `*9/*9`
   depending on unreported rs2279343 — I report both (phenotype is poor
   metaboliser either way, but an exact-string grader will mark the diplotype
   wrong). TPMT homozygous-at-both is `*3A/*3A` on the standard cis
   configuration, with the rare `*3B`/`*3C` trans alternative not formally
   excluded. The CYP2D6 "compound" case (`rs3892097 C/T + rs1065852 G/A`) I call
   `*1/*4`, because c.100C>T rides in cis on `*4` so a `*1/*4` subject is
   heterozygous at both markers, and `*4/*10` would require rs1065852 A/A; if the
   case author intended `*4/*10`, my diplotype is wrong while the phenotype
   (intermediate metaboliser) still holds.

5. **UGT1A1/irinotecan is not CPIC in this corpus (3 cases).** The brief asks for
   the "CPIC Level A" recommendation; the corpus carries only DPWG and RNPGx
   annotations for irinotecan. I answered from DPWG (70% starting dose for
   `*28/*28`) and flagged the source rather than abstaining. If the intended
   answer is "no CPIC recommendation", these should have been abstentions.

6. **Genes whose corpus excerpt has no diplotype table at all (7 cases).** The
   TPMT/NUDT15 annotation (PA166104933) is a change log plus a summary paragraph;
   it has no phenotype-indexed table. I took the recommendations from its summary
   sentences (lower than standard starting dose for intermediate metabolisers;
   consider an alternative immunosuppressant for poor metabolisers) and supplied
   the diplotype → phenotype step entirely from my own reasoning. Defensible,
   but it is the thinnest textual grounding in the set, and the specific dose
   percentages from the 2025 update are absent.

7. **G6PD heterozygous female → "Variable" (1 case).** This follows the corpus
   Variable row exactly, and the resulting action is "measure enzyme activity
   first". If the reference instead resolves her to Deficient or Normal, both the
   phenotype and the recommendation differ. The diplotype `B/A−` also assumes the
   untyped second allele is reference.

8. **Negative-result phenotype wording (2 cases).** RYR1-negative is "Uncertain
   susceptibility" in CPIC, not "not susceptible", and CYP3A5 phenotype labels use
   the corpus's expresser/non-expresser wording. A reference using different
   labels would mismatch on the phenotype string while agreeing on the action.

9. **Abstention calibration.** I abstained on 2 of 110. The defensible range here
   is roughly 2–9: adding the five single-locus warfarin cases and the three
   irinotecan cases would all be arguable under Rule 1. I chose to answer where a
   phenotype-indexed recommendation existed somewhere in the supplied text for
   that gene–drug pair, and to abstain only where the supplied text contains no
   recommendation for the gene at all. Judged the other way, I am under-abstaining
   rather than over-abstaining.

10. **Phenotype label strings vs exact-match grading.** Several phenotype labels
    carry qualifiers taken from the corpus tables ("CYP2C9 Intermediate
    Metabolizer (activity score 1)", "CYP3A5 Extensive Metabolizer (CYP3A5
    expresser)"). These are the corpus's own wording, but exact string matching
    against a shorter reference label would count them as misses even though the
    substance agrees.
