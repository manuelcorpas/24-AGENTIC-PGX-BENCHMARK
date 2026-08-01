# Corrections

A public log of errors found in this project's own evaluation harness, with the
numeric effect of each and the guard that now prevents recurrence.

The manuscript this repository supports argues that agent-mediated pipelines fail
in ways that look like success. Evaluation harnesses have the same disease. Each
entry below produced, or would have produced, a confident, plausible, wrong
result that would have survived peer review. They are recorded here rather than
quietly fixed, because a paper arguing for auditability that does not audit
itself in public has not made its own case.

C1 to C8 were found on 2026-07-26 and 2026-07-27, after the results they affect
had been written up. C9 to C13 are different in one respect worth stating: they
were caught on 2026-07-30 while building the input-normalisation experiment, in
new code, before any number from it reached a manuscript. They are logged on the same footing
anyway. A harness artefact caught early is the same defect as one caught late,
and a log that only records the embarrassing ones would misrepresent the rate.

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

## C9. A variant cap withheld the defining variant, then scored the result as model failure

**Found:** 2026-07-30, before any result existed. **File:** `code/74-input-normalisation.py`.

The input-normalisation experiment shows a model the non-reference genotypes a
sample carries across a gene region and asks for the star-allele diplotype. The
first implementation capped the list at 60 variants. Observed counts run to 681,
with a median of 73.

**Effect had it not been caught.** 318 of 527 prompts were truncated. Truncation
is not random with respect to the answer: it drops variants by position, so for
any sample whose allele-defining variant sits late in the region, the model was
asked to identify a haplotype from input that no longer contained the evidence
for it. Every resulting failure would have been reported as a model unable to
normalise input. This is C5 in a new costume, a harness limit presented as a
finding.

**Guard.** The cap is 800, above the observed maximum, and the build step reports
the truncation count so a reader can see it is zero rather than take it on trust.

## C10. The harness asserted a completeness claim about its own input

**Found:** 2026-07-30. **File:** `code/74-input-normalisation.py`.

The prose rendering closed every prompt with "No other non-reference genotype was
observed in the region examined." That sentence was emitted unconditionally,
including when the variant list had been truncated, where it is false.

**Effect had it not been caught.** The harness would have told the model that an
incomplete list was complete, and then scored the model on the answer. A model
reasoning correctly from a false premise supplied by the evaluator would have
been recorded as wrong. The paper's own argument is that fabricated inputs
produce confident wrong outputs downstream; the evaluator was doing the
fabricating.

**Guard.** `render()` takes a `truncated` flag and every rendering discloses
truncation explicitly. A test asserts the completeness sentence appears only when
the list is complete.

## C11. An output-token cap read as a model declining to answer

**Found:** 2026-07-30. **File:** `code/74-input-normalisation.py`.

The experiment initially reused the matched factorial's model clients, which cap
output at 320 tokens. That is correct for the factorial's one-line task. Input
normalisation makes a model reason about which variants define which haplotype,
and the answer arrives after that reasoning.

**Effect had it not been caught.** Claude was cut off mid-reasoning, before the
`DIPLOTYPE:` line. o3 was worse: it spends its budget on reasoning tokens and
returned empty content on **100 per cent** of calls at both 320 and 2,000 tokens.
The parser, correctly, found no diplotype. Absorbed into the abstention rate,
this would have been published as "o3 abstains on every case", a clean and
completely false finding about a frontier model. Coverage and abstention are the
headline metrics of this experiment, so the artefact would have landed directly
on the number that matters.

**Guard.** `classify()` separates `call`, `abstain` and `truncated_output` using
the output-token count against the cap, so a response cut off by the budget can
never be counted as a refusal. Truncated responses are reported as their own
category and excluded from the scored denominator. The budget was raised to
6,000 tokens, and the units still truncated at that budget were re-run at 14,000
rather than left as an unexplained exclusion.

## C12. HTTP 429 read as a model declining to answer, for the second time

**Found:** 2026-07-30. **File:** `code/74-input-normalisation.py`.

C3 in this log records rate limiting being read as a model's failure to follow a
schema. C11, written the same afternoon, added a `classify()` function to stop a
truncated response being counted as a refusal. That function checked the output
token count and the response text. It did not check whether the call had failed.

A failed request leaves empty text and zero tokens, which on the text alone is
indistinguishable from a terse refusal. So every errored call was recorded as an
abstention.

**Effect.** In the o3 arm, 132 of 450 calls returned HTTP 429 against an
exhausted quota and were written out as abstentions. The first reading of that
run was 254 abstentions, 102 calls, 94 truncated. The true composition is 123
abstentions, 102 calls, 93 truncated and 132 errors. Abstention is the headline
metric of this experiment, and it was overstated by a factor of two. A retry
batch of 94 calls compounded it: every one failed on the same exhausted quota and
was recorded as 94 clean abstentions at a total cost of $0.00, a figure that
should have been read as the tell it was.

**Guard.** `classify()` now takes the error and returns `error` before any other
test, overriding even text that parses as a valid call, because a partial
response banked alongside an error is not a result. Errors are counted and
reported as their own category and never enter the scored denominator. Both o3
result files were re-classified and the affected arm is being re-run.

**Why this one matters most.** It is the fourth appearance of this class in the
log (C3, C5, C9, C11) and the second appearance of rate limiting specifically. It
was introduced inside the guard written to prevent the previous instance, by
someone who had just written that guard, on the same day, having explicitly
enumerated the failure mode in this file. The manuscript claims that knowing a
pipeline can fail silently does not stop it failing silently, and that only a
mechanical check does. This entry is that claim happening to the authors.

## C13. Provider failures left in the denominator, inflating abstention

**Found:** 2026-07-30. **File:** `code/75-evaluate-input-normalisation.py`.

C12 stopped a failed call being recorded as an abstention in the raw rows. The
evaluator then made the same mistake one layer up. It excluded `truncated_output`
from the scored denominator but left `error` rows in it, counted as "did not
emit". A request that failed at the provider tells us nothing about what the
model would have answered, so it belongs out of the denominator exactly as a
truncated response does.

**Effect.** Arm 1 abstention read 0.758 against a true 0.711, and overall coverage
read 0.242 against a true 0.289. Per model the distortion was uneven, because the
errors were not: o3's coverage read 0.286 against a true 0.453. The wrong figures
reached a draft of the response letter and a draft of the manuscript Results
paragraph before the fabrication firewall caught them, and it caught them only
because the number was registered in `72-validate-v30-numbers.py` and recomputed
from raw rows, which disagreed with the summary report.

**Worse, a test asserted the wrong behaviour.** `test_errored_rows_are_counted_as
_abstentions_not_dropped` was written deliberately, to stop errors being silently
dropped. The intent was right and the implementation encoded the opposite error:
not dropping them silently became keeping them in the scored set. A green test
suite therefore certified the bug.

**Guard.** `evaluate()` removes both `truncated_output` and `error` from the
scored set and reports each in its own field. The test now asserts an error row
leaves the denominator, and a second test asserts that a batch of 94 failed calls
yields no abstention rate at all rather than 94 abstentions. Accuracy was
unaffected throughout, because it was already conditioned on emission.

## C14. A gene prefix without a separator, read as an abstention

**Found:** 2026-07-31. **File:** `code/74-input-normalisation.py`.

The parser strips a gene symbol the model may put in front of its answer. It did
so with `^[A-Z0-9]+[- ]+(?=\*)`, which requires a space or a hyphen after the
symbol. `DIPLOTYPE: CYP2D6 *4/*4` therefore parsed and `DIPLOTYPE: CYP2D6*4/*4`
did not, and the second is the standard PharmVar rendering. A well-formed,
frequently correct diplotype was recorded as the model declining to answer.

**Effect.** 796 calls across the deposited rows, 744 of them in the scored set.
Arm 1 coverage read 0.289 against a true 0.542, and abstention 0.711 against a
true 0.458. The distortion was concentrated in the model that most often names
the gene: Claude Opus 4.5 read 0.469 coverage against a true 0.925, and on the
variant-call rendering 0.558 against 0.934.

**It inverted a characterisation, not only a number.** The manuscript reported
that a model without a validated allele-definition table abstains on most real
inputs, which is a safe failure. The corrected rows say it answers nearly
everything and is wrong on about six of every ten answers, which is an unsafe
one. The finding the experiment was built to test survives and is sharper: the
definition artefact barely changes whether the model answers (0.934 to 0.973)
and transforms whether the answer is right (0.396 to 0.967), and against
external GeT-RM consensus the artefact-supplied model still reaches the
deterministic caller (0.784 against 0.774) while two other models do not.

**This is C1 again.** C1 was a parser that discarded 289 correct calls and
inverted a conclusion. The same class recurred here, in the same project, after
that entry was written, in the metric the experiment exists to measure. Parser
strictness was chosen deliberately, because a lenient parser would move the
headline; the lesson is that strictness must be tested against the formats models
actually emit, not only against the ones the author imagined.

**Guard.** The gene the row asked about is now passed to the parser and stripped
case-insensitively with an optional separator; failing that, a symbol shaped like
a gene is stripped. A leading word that is not a gene symbol is still not
stripped, so a hedge (`unknown*1/*2`) remains an abstention rather than becoming
C1 in reverse. Five tests cover the separatorless form, the hyphenated form, a
lowercase symbol, a prefix naming a different gene, and three hedges that must
not parse. `code/79-reparse-normalisation.py` re-derives `call` and `status` from
the stored response text with no new API calls, and reports any call the fix
removes rather than adds; it removed none.

## C15. Presentation wrappers read as abstentions

**Found:** 2026-08-01. **File:** `code/74-input-normalisation.py`.

The input-normalisation parser required `DIPLOTYPE:` to begin a line. Eighty-two
responses in the seven-model analysis instead placed the same explicit marker
and a complete schema-valid diplotype inside Markdown bold. Those calls were
recorded as abstentions. A Gemini pilot used the same marked value inside a
LaTeX box, exposing that any model-specific exclusion would repeat the error at
the model-selection layer.

**Effect.** Re-parsing the stored seven-model rows recovers 87 marked calls and
removes two provisional calls that were followed by an explicit final ABSTAIN;
it changes no raw response text. The corrected per-model values are regenerated by
`code/83-freeze-seven-model-normalisation.py`; the manuscript, Figure 8 and Table
S9 are generated only from that frozen file. No Gemini performance estimate is
reported because a complete confirmatory Gemini run was not collected.

**Guard.** The parser now treats presentation wrappers as irrelevant but still
validates the value immediately following the explicit marker. It does not
search the surrounding reasoning for star alleles, does not strip a hedge before
the value, and does not partially salvage unsupported tandem expressions. If a
response emits the marker more than once, the last marked value is the final
answer. Tests cover Markdown, LaTeX, multiple markers, trailing explanations,
hedges and tandem values. The same rule is applied to every stored response,
irrespective of model.

## What this log is for

Four of the first eight biased results toward this project's own hypothesis; one
biased against it. None was found by peer review. All eight were found by
re-deriving a number that looked too clean, and three of the last four were found
only because a referee's scrutiny prompted a re-reading of work already believed
finished.

C9 to C13 were found a different way, and the difference is the useful part.
Three were caught by writing the test before the code and asking what the harness
would have to get wrong for a clean result to be false. C12 was caught by a cost
of $0.00 on 94 calls that had supposedly produced answers.

C9, C11 and C12 are all the same defect as C5: a limit of the evaluation
apparatus arriving at the analyst's desk wearing the costume of a finding about a
model. C12 is the sharpest case, because it was introduced inside the guard
written to prevent C11, on the same day, by the same hands, in a file that already
documented the failure mode twice. That is the strongest evidence in this log for
the manuscript's actual claim. Knowing the failure mode does not protect you from
it. Only the mechanical guard does, and only for the exact case it checks.

The general lesson is stated in the manuscript: an evaluation harness is itself
an agent-mediated pipeline, and it fails silently in the same way.

Corrections to this log are welcome as issues on this repository.
