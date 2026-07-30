# Revision results (CELL-GENOMICS-D-26-00551)

> **CORRECTED 2026-07-27.** Every number in this file has been regenerated from
> post-fix data. The earlier version of this file carried a pre-fix N1 table, a
> pre-fix N6 variance table, a three-cohort N5 result, and a section arguing that
> execution costs accuracy when rules are supplied. All four are superseded and
> have been replaced in place. The corrections, with their numeric effect, are
> listed in `CORRECTIONS.md`.
>
> Authoritative sources for everything below:
>
> - `data/v3_five_cell_live_report.txt` (matched factorial, both scorers)
> - `data/v3_five_cell_live_stats.txt` (intervals and variance decomposition)
> - `data/v3_ancestry_four_cohorts.txt` (ancestry, four cohorts)
> - `real-genome-arm/n0/n0_result.json` (end-to-end mapping concordance)
>
> These are the numbers written into manuscript v30. See `CHANGES-v29-to-v30.md`
> in the manuscript folder.

All numbers below are live runs from this repository, scored by
`code/61-rescore-matched.py` under both the baseline scorer and the frozen
clinical-equivalence scorer (SCORING-PREREG.md). Total API spend: about $75.

## Roster: three of nine models could not be re-run

`claude-opus-4-20250514` and `claude-sonnet-4-20250514` return 404, retired by
the vendor between submission and revision; `mistral-large-2411` was withdrawn.
Substituted with the same-tier current models. **A benchmark whose correctness
depends on specific model versions decays as vendors retire them; correctness
executed in versioned code does not.** Mistral is reported separately: its tier
returns HTTP 429 under concurrency (447 of 462 calls), which is a rate limit and
not a model failure, and it is re-run paced.

## N1: matched factorial, five cells, eight models common to all cells

Identical patient text, the same named target drug, one output schema, one
scorer. n = 2,640 per cell.

Baseline scorer. A1 is phenotype accuracy, A2 is recommendation accuracy.

| cell | knowledge | mechanism | phenotype | recommendation | lethal-class | errors | coverage |
|---|---|---|---|---|---|---|---|
| free_generation | none | model | 0.7436 | 0.6241 | 0.8393 | 54 | 0.9977 |
| rag_generation | prose | model | 0.8197 | 0.5508 | 0.5923 | 137 | 0.9917 |
| rag_execution | prose | code | 0.9652 | 0.9705 | 0.9315 | 23 | 0.9973 |
| skill_execution | rules | code | 0.9674 | 0.9729 | 0.9256 | 25 | 0.9996 |
| skill_generation | rules | model | 0.9644 | 0.9705 | 0.9226 | 26 | 0.9981 |

Under the frozen clinical-equivalence scorer, only the two generation-from-prose
cells move (free_generation A1 0.7686, rag_generation A1 0.8845); every execution
cell is identical under both scorers, which is what determinism predicts.

Four findings, in order of how much they matter.

**1. Knowledge representation dominates mechanism.** Validated structured rules
reach 0.96 to 0.97 whatever the mechanism. Retrieved prose reaches 0.55 generated
and 0.97 executed. Free recall reaches 0.62. No mechanism change approaches the
margin that the knowledge representation buys. This is the paper's real result.

**2. Execution rescues weak knowledge and matches strong knowledge.** With prose,
moving the decision from generation to execution lifts recommendation accuracy
0.5508 to 0.9705 and cuts lethal-class errors from 137 to 23. With rules already
supplied, execution and generation are level: 25 lethal errors against 26,
recommendation 0.9729 against 0.9705. Execution is a safety net for weak
knowledge and costs nothing when the knowledge is already strong.

**3. `rag_execution` equals `skill_generation`.** Extracting from a guideline and
executing the result lands in the same place as hand-authored rule tables
(0.9705 recommendation in both cells; 23 lethal errors against 26). Hand-authored
rule tables are not a precondition for the result. This is new, deployable, and
absent from the submitted paper.

**4. Retrieval alone degrades safety relative to no knowledge at all.**
Free-prompted makes 54 lethal-class errors; retrieval-augmented makes 137.
Handing a model a guideline it must reason over more than doubles its
lethal-class errors. N4 shows this is not an indexing artefact.

## N4: drug-keyed retrieval changes nothing

Same six models, gene-keyed against drug-keyed indexing:

| cell | gene-keyed lethal errors | drug-keyed lethal errors |
|---|---|---|
| rag_generation | 98 | 96 |
| rag_execution | 22 | 20 |

The retrieval safety deficit is **not** an indexing artefact. The 16b prototype
suggested drug-keyed chunking eliminated drug substitution; under the matched
harness, where the prompt already names the target drug, better indexing buys
nothing. Clean negative result.

## N6: variance decomposition

| cell | case variance | model variance |
|---|---|---|
| free_generation | 0.082605 | 0.004778 |
| rag_generation | 0.090589 | 0.002251 |
| rag_execution | 0.008855 | 0.000859 |
| skill_execution | 0.011123 | 0.000551 |
| skill_generation | 0.009685 | 0.000610 |

Between-model variance falls by roughly an order of magnitude once the knowledge
is structured, and falls again under execution: 0.004778 free, 0.002251 retrieved
prose, then 0.000551 to 0.000859 in the three cells that use rules or execute.
Case variance also collapses under execution, from about 0.09 in the two
generation-from-prose cells to about 0.009 to 0.011.

Both readings should be stated together. The empirical claim is that model
dependence shrinks as knowledge becomes structured and as the decision moves into
code. The structural claim is that an executed mapping is deterministic given a
valid input call, so the residual model dependence sits entirely in the calling
step. Pair them deliberately: an observation over eight models has a shelf life,
a determinism guarantee does not.

Model-clustered intervals are reported alongside case-clustered ones throughout,
which is R2.5's point. For rag_execution the model-clustered interval is
(0.9439, 0.9845) against a case-clustered (0.9466, 0.9807), and the two-way
clustered interval is (0.9356, 0.9886). Model clustering does widen the interval
in every cell, and the two-way interval is wider still, so the reviewer's
correction stands and is applied; the widening is not the large effect an earlier
pre-fix run suggested.

## N2: PharmCAT

Concordance 83/83 where both engines emit; PharmCAT abstains on 4; 23 of 110
clinically real states cannot be expressed in its input vocabulary at all.

## N5: ancestry

Four cohorts. 373 distinct states, 107 observed in more than one cohort.

| cohort | states | cov(state) | cov(carrier) | standardised | own-states |
|---|---|---|---|---|---|
| 1000G_IBS (n=93) | 137 | 0.1606 | 0.3852 | 0.3627 | 0.0588 |
| CorpasFamily (n=5) | 40 | 0.2750 | 0.3176 | 0.3911 | 0.0000 |
| Peru | 73 | 0.1918 | 0.4190 | 0.3442 | 0.1429 |
| UGR | 289 | 0.0657 | 0.2715 | 0.3216 | 0.0290 |

Raw state-coverage spread 0.2093 collapses to 0.0695 after direct standardisation
to a common state distribution. Standardisation compares cohorts only on states
they share, so it is blind to states a single cohort carries, which is exactly
where the untypeable rare-allele tail sits; on cohort-specific states the spread
is 0.1429. Report both columns. Reporting only the standardised figure would
announce that the disparity had vanished when the analysis had merely stopped
looking at it.

Platform caveat: the family cohort is 23andMe array data and the IBS cohort is
sequencing, so any family-versus-IBS difference could be assay rather than
population.

**Outstanding:** GeT-RM consensus diplotypes for the calling step. The evaluator
(`real-genome-arm/scripts/08_caller_truth_eval.py`) is built and tested, and
refuses to run on an empty truth set; the data requires a documented manual
download (see `real-genome-arm/getrm/README.md`).

## What the paper should claim

Not "versioned skills are more accurate". The defensible claim, supported by
every number above:

> Correctness comes from validated, versioned, structured knowledge, not from
> model recall and not from retrieval over prose. Executing that knowledge in
> code does not by itself buy accuracy where the knowledge is already structured;
> what it buys is determinism, bounded failure, exact error localisation and
> invariance across models. Where the knowledge is weak, execution also buys a
> large safety improvement: given a retrieved guideline, execution lifts
> recommendation accuracy from 0.5508 to 0.9705 and cuts lethal-class errors from
> 137 to 23. Under execution the residual error localises entirely to input
> interpretation, which is separately measurable and separately fixable with a
> validated caller.

## Mistral: a rate limit, demonstrated as such

Same model, same prompts, same day, varying only the request cadence:

| cadence | error rate |
|---|---|
| full concurrency | 96.8% (447 of 462) |
| 3 s between calls | 19.6% |
| 8 s between calls | 2.4% |

Every other model ran at 0 errors throughout. A deterministic per-model failure
that disappears when you slow down is a rate limit, not a model failure, and
scoring those empty responses as format failures would have produced a false
claim about the model's ability to follow an output schema. This is the second
time this trap has appeared in this project; it is now guarded in code
(`--pace`, plus 429 backoff inside `run_one`) rather than in memory.

Headline numbers are reported on the eight-model common set, which is what the
plan of revision already committed to; Mistral is a separately reported
sensitivity check at the paced cadence.

## What was executed, and what is externally blocked

Executed on real data in this session: N1 (five cells, eight models, 13,199
scored evaluations), N2, N4 (5,280 evaluations), N5 ancestry, N6, N7, N9.
Total spend about $75 against a $300 ceiling.

Externally blocked, not incomplete work: the GeT-RM consensus diplotypes needed
for calling-step validation. The CDC pages are JavaScript shells, the NCBI FTP
path 404s, and the consensus tables are journal supplementary material, so
acquisition is a documented manual download (`real-genome-arm/getrm/README.md`).
The evaluator is built and covered by six tests, and refuses to run on an empty
or malformed truth set rather than reporting a vacuous 100%. No truth set was
improvised: fabricating one in a paper about trustworthy pipelines would be the
worst possible failure.

## Error localisation, measured on live data

Treating the model as the caller and the benchmark ground truth as external
truth, across both execution cells, eight models, three replicates
(`real-genome-arm/scripts/08_caller_truth_eval.py`, 5,273 calls):

    call concordance            0.9579   (222 errors)
    end-to-end phenotype        0.965 (rag_execution), 0.967 (skill_execution)

End-to-end accuracy tracks call accuracy to within a point, the residue being
distinct diplotypes that map to the same phenotype. Under execution the pipeline
has no independent failure mode: fix the caller and you fix the pipeline.

Weakest genes for the model-as-caller, i.e. where a deterministic caller buys
the most: TPMT 0.811, CYP2B6 0.875, CYP2C9 0.903, SLCO1B1 0.931, DPYD 0.935.

This is the calling-step validation the plan promised, run against curated
ground truth. GeT-RM remains required for validating the DETERMINISTIC caller
(PyPGx) on real genomes, and that acquisition is still a manual step.
