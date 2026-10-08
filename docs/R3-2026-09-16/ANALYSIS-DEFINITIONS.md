# Analysis definitions

The main comparison retains all 110 cases and historical references. Five configurations, eight model identifiers and three replicates yield 2,640 attempts per configuration. There are 13,199 returned records. The absent authored-rule-execution record is filled with an explicit zero-score attempt. Empty input files, duplicate evaluation keys, duplicate cases and unexpected design keys are rejected. Missing, unparsed and abstaining responses score zero on the two concordance endpoints; coverage requires a parsed, non-abstaining response. Parse failures include the absent record and can overlap abstentions.

The frozen baseline scorer remains `code/61-rescore-matched.py`, using `10-rescore-v3.py` and the original rule engine. The phenotype-equivalence layer remains `10b-rescore-v3-clinical-equivalence.py`, with registered pattern fingerprint `bfab39f4ebc5bd53aacf6b3805478af29eb0cb06042151fefaa8e3d85b626d94`. It was authored after inspection of outputs; it is frozen, not prospectively preregistered. Baseline and equivalence scores are reported together. No new RYR1 rubric is introduced.

## Three different RYR1 recommendation analyses

| Analysis | Arithmetic | Percentage |
|---|---|---:|
| Original reference, all 110 cases | 2,566.5 / 2,640 | 97.22% |
| Exclude the negative RYR1 case, 109 cases | (2,566.5 - 24) / 2,616 | 97.19% |
| Retain 110 cases, assign its 24 recommendations zero credit | (2,566.5 - 24) / 2,640 | 96.31% |

The last row is a reproducible arithmetic sensitivity calculation, not corrected-reference output from the frozen scorer. The corrected uncertain-susceptibility phenotype is inadequately handled by that scorer. All original prompts, reference files, outputs and scorer code are preserved.

Exclusion changes the largest matched phenotype/recommendation endpoint by 0.7542743953 percentage points, reported as 0.75 at two decimal places; the largest authored-configuration change is 0.0323185988 points, less than 0.04. Exclusion uses 2,616 attempts per configuration and leaves every high-severity denominator at 336 because the excluded case is a negative control.

## Wording sensitivity

High-severity membership is the original authors' `LETHAL` marker in reference text, covering 14 cases and 336 attempts per configuration. It is an operational benchmark flag, not a standard CPIC severity classification or adjudicated harm endpoint. The original rule recognises `avoid` or the `contraindicat` prefix and rejects drug substitutions. The post hoc widened rule also recognises `do not use`. Frozen phenotype equivalence is unrelated to this sensitivity and remains unchanged.

V4 original/widened failures are 46/43, 137/75, 22/22, 24/24 and 23/23 in main-table order. The v3 free/RAG counts are 54/137 original and 50/75 widened. Widening never increases failures. The retrieval-minus-free difference uses paired case-cluster resampling of the 14 distinct cases, 10,000 resamples, seed 20260727: original 27.1 points (95% interval 8.0 to 47.0), widened 9.5 points (-2.7 to 25.3). These intervals describe lexical classifications, not clinical harm.

## Candidate screen

For all 110 cases, select baseline phenotype means >=0.95 in each authored configuration and <=0.10 in RAG generation. Eight candidates meet the post hoc criterion. Review is limited to these candidates. One confirmed reference error among eight selected cases is not an estimated error rate for all 110 cases. Shared errors can be missed when retrieval agrees or the thresholds are not crossed.

## Table 2 and secondary comparisons

Table 2 preserves its existing evaluator (`68-evaluate-extracted-rules.py`) and all 2,639 available authored-rule-execution calls. The authored table covers 2,638/2,639 calls (99.96%); phenotype concordance is 2,552/2,638 among covered calls (96.7%) and 2,552/2,639 across all available calls (96.7%). Its normalisation/evaluation function differs from the matched benchmark scorer. It is not a recommendation metric and must not be substituted for Table 1's attempted denominator.

All eight extracted table rows are reevaluated on those calls. Their conditional phenotype range is 72.8% to 96.2%; the two broadest tables have 73.8% conditional concordance. The legacy caller evaluator on the 5,272 usable V4 execution calls gives 225 discordances and 95.73% concordance. This is a separate input-call endpoint, conditional on the original authored vocabulary.

Figure 3B retains its historical seven-model comparison, which was frozen before curated V4 results were available. Figure 3A contains the eight-model real-genome panel. New V4 benchmark data do now exist; the manuscript no longer claims otherwise. No extra clinical or model evaluation was run for this revision.
