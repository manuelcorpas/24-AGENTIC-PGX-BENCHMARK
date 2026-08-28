# Post hoc exploratory analyses — NOT part of the frozen scoring pipeline

**Status: exploratory. Not clinically validated. No claim in the manuscript rests
on anything in this directory.**

These scripts were written *after* the primary results were seen, while auditing
the lethal-class avoidance endpoint in `code/61-rescore-matched.py`. They are
retained for transparency: deleting analysis code after seeing its output would
reduce, not improve, the record. They are kept here, outside `code/`, so that
they cannot be mistaken for the pre-specified pipeline.

## What is here

| file | purpose |
|---|---|
| `94-directive-strength.py` | classifies a response on an ordinal directive ladder and separates prohibitions gated on information the case does not supply |
| `95-directive-analysis.py` | applies the above across the five configurations |
| `96-figure-s10-directive-ladder.py` | draft figure, not included in the submission |
| `test_directive_strength.py` | unit tests for the classifier |

## Why it is not an endpoint

The construct is automated, rule-based, and was developed after the results were
seen. It has not been validated against blinded clinical adjudication, which is
the only thing that could establish whether a conditional prohibition ought to
count as an avoidance. Pre-specification protects against selective analysis; it
does not establish construct validity, and neither does this code.

The manuscript reports it, if at all, as a hypothesis for future work.

## What IS in the manuscript

- the pre-specified lexical avoidance endpoint (`61-rescore-matched.py`), reported for transparency
- a bounded vocabulary sensitivity analysis that additionally accepts `do not use`
- the 108 responses scored correct at the recommendation endpoint and incorrect at the avoidance endpoint, reported as differential measurement
