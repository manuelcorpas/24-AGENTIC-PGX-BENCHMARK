# Corrections

A public log of errors found in this project's own evaluation harness, with the
numeric effect of each and the guard that now prevents recurrence.

The manuscript this repository supports argues that agent-mediated pipelines fail
in ways that look like success. Evaluation harnesses have the same disease. Each
entry below produced a confident, plausible, wrong result that would have
survived peer review. All four were found in a single session on 2026-07-26. They
are recorded here rather than quietly fixed, because a paper arguing for
auditability that does not audit itself in public has not made its own case.

Nothing in this log affects the numbers supplied to the editor in the plan of
revision on 2026-07-26. Those numbers come from `real-genome-arm/n0/n0_result.json`,
a deterministic caller compared against a deterministic skill with no model in the
loop, and were re-derived and confirmed on 2026-07-27.

## C1. A parser bug discarded 289 correct calls and inverted a conclusion

**Found:** 2026-07-26. **Fixed:** commit `89186f1`.

`_pgx_rules.py` required a `DIPLOTYPE:` label to extract a called diplotype.
GPT-4.1 frequently answered in the correct controlled vocabulary without the
label, for example `*1/*4`. The parser discarded those responses as unparseable.

**Effect.** GPT-4.1's `rag_execution` accuracy read 0.436 against a true 0.955.
Across the run, 289 correct calls were discarded this way. Because GPT-4.1 sat in
the eight-model common set, the
error propagated into the aggregate and reversed a headline conclusion. The
pre-fix analysis reported that, with rules already supplied, execution slightly
hurts: 26 lethal-class errors becoming 33. Post-fix the two are level, at 25
against 26. The paper's Discussion had been written around the pre-fix reading.

It also inflated an apparent variance result. The pre-fix model-clustered interval
for `rag_execution` was (0.764, 0.983) against a case-clustered (0.880, 0.917),
which looked like a large model-dependence effect. Post-fix it is (0.9439, 0.9845)
against (0.9466, 0.9807), and between-model variance is 0.000859, not 0.0315.

**Guard.** `parse_field` gained an `allow_bare` mode, enabled for execution cells
only, and only for a terse single-line reply carrying no field labels, so a value
can never be salvaged out of prose that answered a different question. Format
compliance is still measured separately by `is_format_compliant()`, so relaxing
the parser does not hide a genuine format failure. Unparseable responses are
counted and reported as `parse_fail` per cell rather than silently dropped. The
re-parse used stored raw text and made no new API calls.

## C2. A scoring rule required the drug name and penalised correct terse answers

**Found:** 2026-07-26.

The recommendation metric required the target drug name to appear in the answer
text. Execution cells emit canonical rule text that always contains the drug name.
Generation cells frequently answered correctly and tersely without repeating it,
and were scored wrong.

**Effect.** The measured advantage of execution over generation was inflated. The
bias ran in the direction of the paper's own hypothesis, which is the worst
direction for a bias to run.

**Guard.** Scoring is on the recommendation content. The frozen
clinical-equivalence scorer and the baseline scorer are both reported for every
cell in `data/v3_five_cell_live_report.txt`, and the scorer tables are SHA-256
frozen (`SCORING-PREREG.md`, `bfab39f4...b626d94`); `--frozen` refuses to score if
they move.

## C3. HTTP 429 rate limiting read as a model's failure to follow a schema

**Found:** 2026-07-26. This is the second occurrence in this project.

Mistral's tier returned empty responses under concurrency. The harness had no
"choices" key to parse and recorded the result as a format failure. The pattern
was deterministic per model and per call index, which is what a real model
limitation looks like.

**Effect.** Left unexamined, the paper would have claimed that a frontier model
cannot follow an output schema. Same model, same prompts, same day, varying only
request cadence:

| cadence | error rate |
|---|---|
| full concurrency | 96.8% (447 of 462) |
| 3 s between calls | 19.6% |
| 8 s between calls | 2.4% |

Every other model ran at zero errors throughout. A deterministic failure that
disappears when you slow down is a rate limit.

**Guard.** `--pace` plus 429 backoff inside `run_one`, in code rather than in
anyone's memory. Headline numbers are reported on the eight-model common set;
Mistral is a separately reported paced sensitivity check. The exclusion criterion
was registered in `SCORING-PREREG.md` section 6 before it was applied, on grounds
independent of the model's scores.

## C4. An empty `AD` field after merge looked like eight uncallable genes

**Found:** 2026-07-26.

After a VCF merge step the `AD` field was `.` for a subset of records. PyPGx
raised on `int(".")` and the affected genes produced no calls.

**Effect.** Eight of seventeen PyPGx genes appeared uncallable, which would have
been reported as a limitation of the deterministic caller rather than as a
pipeline defect.

**Guard.** A depth-strip step in the cohort pipeline, documented in the IBS
pipeline notes.

## C5. A regular expression silently dropped every HLA allele from a truth set

**Found:** 2026-07-27, before any number was reported.

`11_ingest_getrm_consolidated.py` accepted a star allele as `\*[\w.]+`, which
does not admit a colon. GeT-RM writes HLA alleles as `*57:01:01`, so every HLA
column failed the pattern and was discarded as unparseable.

**Effect.** HLA-A and HLA-B were absent from the caller-truth set entirely: 2,828
genotypes instead of 3,554. Those are the loci carrying four of the benchmark's
fourteen lethal-class cases. The ingester reported success. The genes where a
wrong call matters most were the ones going missing, and nothing in the output
said so.

**Guard.** The allele pattern admits colons after a star. Tests cover HLA alleles
both as single calls and as diplotype pairs. The ingester prints every column
that yielded no usable call rather than dropping it in silence.

## C6. SNP genotypes were admitted as star-allele diplotypes

**Found:** 2026-07-27, before any number was reported.

The GeT-RM consolidated table interleaves rsID-specific columns among the
star-allele columns. Their values (`G/A`, `C/C`) match the *shape* of a diplotype
without being one, and the parser accepted them.

**Effect.** Had they entered the truth set, they would have been compared against
PyPGx star-allele calls, mismatched every time, and understated caller
concordance. The bias ran against the deterministic caller, which is the
direction that would have made this paper's own architecture look worse, but a
wrong number in a favourable direction is still a wrong number.

**Guard.** A call is admitted only if it is expressed in an allele vocabulary:
something containing a star or an HLA-style colon, and never a bare nucleotide.
Eleven excluded columns are named in the ingester's output.

## C7. Two errors in the comparison built to answer a referee's point

**Found:** 2026-07-27, before any number was reported.

`15_agent_vs_deterministic.py` compared the executed skill against the caller's
phenotype. Two defects in the comparison itself:

The phenotype tier matcher did not strip hyphens, so `Ultra-rapid Metaboliser`
missed the `ultrarapid` needle, fell through to `rapid`, and was scored as
disagreeing with `Ultrarapid Metabolizer`. A spelling difference reported as a
clinical one.

PyPGx returns `Indeterminate` for CYP4F2, whose alleles it does not assign a
function to. Those states carry no reference answer, and scoring against them
counted eight *reference* abstentions as skill errors.

**Effect.** The deterministic arm read 0.82 to 0.91 accuracy among emitted
answers instead of 1.00, contradicting N0's 38 of 38. Reported unexamined, the
paper would have understated its own deterministic component using a comparison
written to defend it.

**Guard.** Hyphens and spaces are stripped before tier matching.
`is_reference_abstention()` excludes reference-silent states from both arms, so
they face one denominator, which is what N0 did and what this comparison should
have done from the start.

## C8. A metric-comparability error, caught in adversarial self-review

**Found:** 2026-07-27, before publication. **Not committed to any released number.**

The four-cohort real-genome rerun was compared against the submitted 72/51/40
figures and reported as showing that the ancestry ordering collapsed, because the
Latin American cohort scored above the population-matched Iberian one.

The two numbers were different metrics. The submitted figure counts an
indeterminate answer on a determinate state as a failure. The rerun figure was
accuracy among answers the model was willing to give, which excludes abstentions.
The Latin American cohort abstains on 27% of determinate states against 14% for
the Iberian cohort, so excluding abstentions inflated it.

Recomputed on the submitted metric: family 0.694, Iberian 0.606, Latin American
0.534, East African 0.372. The ordering survives; what changes is that the
European anchor drops when the five-member family is replaced by a population
sample, halving the Europe-to-Latin-America gap from 0.16 to 0.07.

Effect had it not been caught: a manuscript and a response letter both claiming a
published finding had been overturned, on the strength of a denominator change.

Guard: `real-genome-arm/scripts/15_agent_vs_deterministic.py` now reports both
denominators, and any comparison against a previously published figure must state
which metric it uses.

## What this log is for

Four of these eight biased results toward this project's own hypothesis; one
biased against it. None was found by peer review. All were found by re-deriving a
number that looked too clean, and three of the last four were found only because
a referee's scrutiny prompted a re-reading of work already believed finished. The general lesson is stated in the manuscript: an evaluation harness
is itself an agent-mediated pipeline, and it fails silently in the same way.

Corrections to this log are welcome as issues on this repository.
