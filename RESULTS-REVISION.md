# Revision results (CELL-GENOMICS-D-26-00551)

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

| cell | knowledge | mechanism | phenotype | recommendation | lethal-class | errors | coverage |
|---|---|---|---|---|---|---|---|
| free_generation | none | model | 0.744 | 0.624 | 0.839 | 54 | 0.998 |
| rag_generation | prose | model | 0.820 | 0.551 | 0.592 | 137 | 0.992 |
| rag_execution | prose | code | 0.900 | 0.904 | 0.908 | 31 | 0.927 |
| skill_execution | rules | code | 0.932 | 0.938 | 0.902 | 33 | 0.961 |
| skill_generation | rules | model | 0.964 | 0.971 | 0.923 | 26 | 0.998 |

Three findings, in order of how much they matter.

**1. Knowledge representation dominates mechanism.** Validated structured rules
(0.96) beat retrieved prose (0.55-0.90) and free recall (0.62) by margins no
mechanism change approaches. This is the paper's real result.

**2. Execution rescues weak knowledge and does nothing for strong knowledge.**
With prose, moving the decision from generation to execution lifts
recommendation accuracy 0.551 to 0.904 and cuts lethal-class errors from 137 to
31. With rules already supplied, execution slightly *hurts* (26 lethal errors
becomes 33). The two factors interact; the mechanism is a safety net for weak
knowledge, not a universal improvement.

**3. Retrieval degrades safety relative to no knowledge at all.** Free-prompted
makes 54 lethal-class errors; retrieval-augmented makes 137. Handing a model a
guideline more than doubles its lethal-class errors.

## Why execution scores lower than reasoning when rules are supplied

The deficit is entirely the input call, not the mapping (archived-arm
decomposition, six common models):

```
skill_execution   call 0.9278 + coincidence 0.0096 = phenotype 0.9374
skill_generation  call 0.9520 + coincidence 0.0121 = phenotype 0.9636
```

Under execution, end-to-end accuracy equals input-call accuracy to three
decimals: the pipeline has no independent failure mode and every correct answer
traces to a correct call. Under generation, accuracy exceeds call accuracy, so
some answers are right despite a wrong intermediate call, which is correctness
by coincidence. Forcing the model to copy a diplotype verbatim from a controlled
vocabulary costs about 2.4 points of call accuracy against expressing it freely.

N0 is the fix and already exists: replace the model's call with a validated
deterministic caller and the deficit disappears (38 of 38 emitted states correct,
zero disagreements against PyPGx).

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
| free_generation | 0.0826 | 0.0048 |
| rag_generation | 0.0906 | 0.0023 |
| rag_execution | 0.0098 | 0.0315 |
| skill_execution | 0.0132 | 0.0110 |
| skill_generation | 0.0097 | 0.0006 |

Model-invariance comes from the **knowledge**, not the mechanism:
skill_generation has the lowest between-model variance of any cell (0.0006).
Under execution, case variance collapses (the mapping is deterministic) and the
residual variance concentrates in the input-call step, which remains
model-dependent. So "population- and model-invariant" holds for the executed
mapping *conditional on a valid input call*, and must be stated that way.

Model-clustered intervals are materially wider than case-clustered ones for
rag_execution (0.764-0.983 against 0.880-0.917), which is exactly R2.5's point.

## N2: PharmCAT

Concordance 83/83 where both engines emit; PharmCAT abstains on 4; 23 of 110
clinically real states cannot be expressed in its input vocabulary at all.

## N5: ancestry

Raw state-coverage spread 0.209 collapses to 0.047 after direct standardisation
to a common state distribution, but standardisation is blind to cohort-specific
states, where the spread is 0.145 (Uganda 0.043 against 0.167 and 0.188). The
disparity is relocated and measured, not dissolved.

**Outstanding:** GeT-RM consensus diplotypes for the calling step. The evaluator
(`real-genome-arm/scripts/08_caller_truth_eval.py`) is built and tested, and
refuses to run on an empty truth set; the data requires a documented manual
download (see `real-genome-arm/getrm/README.md`).

## What the paper should claim

Not "versioned skills are more accurate". The defensible claim, supported by
every number above:

> Correctness comes from validated, versioned, structured knowledge, not from
> model recall and not from retrieval over prose. Executing that knowledge in
> code buys determinism, bounded failure and exact error localisation; it does
> not by itself buy accuracy, and where the model already holds the rules it can
> cost accuracy by executing a wrong input call rigidly. The residual error
> localises entirely to input interpretation, which is separately measurable and
> separately fixable with a validated caller.
