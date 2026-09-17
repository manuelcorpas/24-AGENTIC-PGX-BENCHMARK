# Audit of Manuel's Claude Science session

Audited 15 September 2026. Session: [Pharmacogenomics Data Upload](http://localhost:8765/projects/proj_17c76f19215f/frames/813ed708-7f6b-4028-a92b-2fcec606e436). This is a local application link.

## Finding

**The delivered pipeline reproduces its 110 answers. It also makes unsupported assumptions about missing genotypes and mishandles an explicitly negated duplication result. The benchmark itself has wording-sensitive scoring and a material RYR1 answer-key concern.**

This supports the challenge that simple mapping skills are transient. It does not establish clinical validity, establish ClawBio's incremental value, or justify abandoning ClawBio. The immediate scientific task is to separate reference lookup, input interpretation and evaluation validity.

## 1. What was inspected

- The exact session supplied by Manuel, through its authenticated browser UI. The model selector displayed Opus 5. The precise backend revision was not independently established.
- Its nine-file artifact export, downloaded through the session's own download action and preserved as `RAW/CLAUDE-SCIENCE-SESSION.zip`.
- The original brief, 110 case inputs, supplied guideline corpus and CPIC cache, copied from the existing `CLAUDE-SCIENCE-COMPARISON/CPIC-STUDY` folder. These additional inputs were not all inside the nine-file export.
- The generated `pgx_pipeline.py`, `run_study.py`, methods, answer CSV and validation artifacts.
- The matched scorer and reference cases at `agentic-pgx-benchmark-v2.10`, commit `eed883434c190e80fd549e2045cb540177b42dbb`, in the existing agentic-PGx benchmark repository. The copied scorer files matched that commit.

The exported answer CSV, original study answer CSV and previously scored answer CSV were identical: SHA-256 `8534bb140a500500abb60cfff52a784bebf67a6d7083cf591ea763eed81a09cb`.

Full raw source hashes are in [SOURCE-HASHES.json](OUTPUT/SOURCE-HASHES.json). The full conversation was not exported or audited as a complete record of every tool call and hidden runtime instruction.

## 2. Frozen scores, independently recomputed

| Measure | Result | Interpretation |
|---|---:|---|
| Cases / parsed rows | 110 / 110 | No missing or unparsable delivered answer rows. |
| Abstentions | 0 / 110 | Full output coverage; not evidence that every input justified an answer. |
| A1 phenotype | 107 / 110, 97.27% | Agreement with the existing phenotype matching rules. |
| A2 recommendation | 63.5 / 110, 57.73% | 57 full-credit, 13 half-credit, 40 zero-credit rows under the existing wording rules. |
| Benchmark safety-action subset | 8 / 14 lexical passes | Six fail the scorer's narrow keyword test. This is not a clinical error count. |
| Frozen clinical-equivalence scoring | Same results | This path changes phenotype patterns, not the safety-action keyword rule. |

Sources: [LOCKED-SUMMARY.json](OUTPUT/LOCKED-SUMMARY.json), [LOCKED-SCORES.json](OUTPUT/LOCKED-SCORES.json), and the preserved scorer source.

The scorer labels a case as lethal if the reference drug text contains `LETHAL`. Its safety-action test accepts `avoid` or a word starting `contraindicat`. An instruction such as “do not use” fails that test. The generated study's separate 24-case safety classification uses a different definition and must not be substituted for the benchmark's 14-case denominator.

### The six lexical failures

| Case ID | Delivered recommendation contains | What can be concluded |
|---|---|---|
| `tpmt_aza_pm` | Alternative nonthiopurine treatment | A potentially relevant alternative action is present. |
| `nudt15_thio_pm` | Marked thiopurine dose reduction and alternatives conditional on indication | Indication and drug-specific adjudication are needed. |
| `hlab5701_aba_pos` | Abacavir is not recommended | A negative action is present. |
| `hlab1502_carba_pos` | Do not use carbamazepine when treatment-naive | A negative action with a condition is present. |
| `hlab1502_oxcarba_pos` | Do not use oxcarbazepine when treatment-naive | A negative action with a condition is present. |
| `hlaa3101_carba_pos` | Do not use carbamazepine when treatment-naive and alternatives are available | A negative action with conditions is present. |

These are summaries of the actual output, not an independent clinical sign-off. All six require assessment of genotype inference, indication, qualifiers and alternatives. Neither “six clinically unsafe answers” nor “zero clinical errors” follows from the lexical result.

### The three phenotype mismatches

| Case ID | Reference key | Delivered answer | Audit assessment |
|---|---|---|---|
| `mtrnr1_amino_var` | Variant Carrier | Increased risk of aminoglycoside-induced hearing loss | Different levels of description; review semantic equivalence. |
| `ryr1_anaes_susc` | Malignant Hyperthermia Susceptible | Malignant Hyperthermia Susceptibility | Likely grammatical mismatch; preserve the frozen score pending adjudication. |
| `ryr1_anaes_norm` | Normal; standard use of volatile anaesthetics | Uncertain Susceptibility; assess clinical context | Material reference-key concern. |

The last case provides only “RYR1 no pathogenic variant”. CPIC assigns uncertain susceptibility to a negative genetic result and cautions against inferring normal susceptibility without clinical context. The guideline therefore supports the delivered uncertainty, rather than the benchmark's unconditional normal label. This conclusion concerns that synthetic case and its answer key. [CPIC guideline, Table 1 and footnote c](https://files.cpicpgx.org/data/guideline/publication/volatile_anesthetic_succinylcholine/2018/30499100.pdf)

The DOI `10.1002/cpt.1319` was resolved directly using publisher metadata, and its title matched the CPIC-hosted full text. The September CPIC cache also contains the uncertainty classification. No reference key, published result or manuscript was silently changed during this audit.

### Additional cases requiring scientific review

- `cyp2b6_efa_pm`: the supplied single SNP, rs3745274 T/T, does not by itself resolve all relevant haplotypes. The reference's `*6/*6` and the pipeline's `*9/*9` cannot simply be declared harmless naming variants. The pipeline assumes an unreported locus is reference.
- `cyp2d6_codeine_im_decreased`: the unphased SNP input, exact diplotype inference and reference recommendation need separate review. Phenotype agreement can conceal an unsupported haplotype assignment.
- `ugt1a1_iri_im`: the pipeline's no-adjustment recommendation for the supplied `*1/*28` is supported by the DPWG guideline's intermediate-metaboliser recommendation. A mismatch with the local key needs adjudication against the specified guideline, indication and dosing context. [DPWG primary guideline](https://pure.rug.nl/ws/files/774574725/s41431-022-01243-2.pdf)

No corrected overall clinical accuracy has been calculated. An independent PGx reviewer should adjudicate discrepancies against dated sources before any revised benchmark claim is made.

## 3. Reproduction: what worked and what was missing

The original cache records retrieval at `2026-09-15T17:15:02.528671+00:00`. The supplied guideline corpus is older, with May 2026 provenance. Source age and rule version therefore need to be explicit when comparing against historical benchmark arms.

The audit executed the unmodified study scripts in a fresh temporary directory with copies of the original inputs:

1. **Original cache, requests blocked:** execution stopped at an uncached G6PD diplotype request. The original bundle was not sufficient for a fully offline replay.
2. **Supplemental retrieval:** four public CPIC G6PD responses were captured. This replay produced 110 rows matching every field of the delivered answer rows.
3. **Recorded responses, network calls blocked:** the same 110 rows reproduced exactly, with zero live study API calls. A subsequent recorded replay again matched all 110 rows.

The supplemental responses are preserved separately in [SUPPLEMENTAL-CPIC-RESPONSES.json](RAW/SUPPLEMENTAL-CPIC-RESPONSES.json). They are audit additions, not retrospectively attributed to the original bundle. [REPRODUCTION.json](OUTPUT/REPRODUCTION.json) contains the final cache-only failure and recorded replay results. [REPLAY-RECORD-ANSWERS.csv](OUTPUT/REPLAY-RECORD-ANSWERS.csv) preserves the supplemental retrieval run's output; [REPLAY-RECORDED-ANSWERS.csv](OUTPUT/REPLAY-RECORDED-ANSWERS.csv) preserves the offline output.

This is computational reproduction of an existing implementation. It is not a fresh agent replicate, a proof of input validity, or evidence that another model would independently build the same pipeline.

## 4. Three targeted checks of the generated code

These are post-hoc diagnostic inputs executed through the generated pipeline's functions. They are not new independent Claude Science conversations and do not estimate a population error rate.

| Probe | Input | Observed behaviour | Scientific issue |
|---|---|---|---|
| Missing genotype | No genotype result available, CYP2C19 | Inserts reference genotypes, calls `*38/*38`, returns Normal Metabolizer | Missing observations become positive evidence of normal function. |
| Missing locus | CYP2C19 rs12248560 C/T only | Fills rs4244285 as G/G and returns Rapid Metabolizer | An unreported locus becomes reference. With observed evidence alone, the cached model permits multiple phenotype classes. |
| Negated duplication | Laboratory-confirmed CYP2D6 `*1/*1`, copy number 2, “no duplication detected” | Calls `*1/*1xN`, returns Ultrarapid Metabolizer | A substring match ignores negation and contradicts the explicit copy-number result. |

Exact inputs, intermediate assignments and outputs: [TARGETED-CHECKS.json](OUTPUT/TARGETED-CHECKS.json).

For the second probe, compatible alternatives demonstrate underdetermination; they do not quantify the clinical prevalence of those alternatives. For the third, the same cache maps the explicitly stated `*1/*1` to Normal Metabolizer.

The methods disclose an assumption that unspecified loci are reference. The original brief does not grant that assumption and asks for abstention when the evidence cannot map confidently. A disclosed assumption remains a validity problem if it changes missing information into evidence.

These checks identify requirements that a maintained capability could enforce: explicit assay coverage, uncertainty propagation, phase and copy-number handling, and contradiction checks. **No matched current ClawBio arm was run, so its ability to satisfy these requirements remains unmeasured here.**

## 5. Blinding and comparison limits

- A conservative detector found phenotype/status-like tokens in 94 of 110 original case IDs. Examples include `pm`, `norm` and `pos`.
- The guideline corpus contains 42 paths named `n_test_cases` or `verification`, across 21 gene or HLA-labelled entries. The session itself noticed leakage opportunities.
- The project file list included SOUL.md and an onboarding profile. Their presence is a context confound; this audit did not establish that the model consumed them.
- The session used an iterative construction and debugging workflow, rather than a single-pass answer generation task.
- The brief prohibited consulting the published benchmark or ClawBio. The inspected material does not establish a complete audit of every runtime input, built-in skill or tool. Do not relabel the session “no skill used” as a verified fact.
- The model display is known; the exact backend revision and complete system context are not.
- Historical manuscript arms were not rerun under the same current model, runtime, data sources, time budget and tool permissions.

[BLINDING-AUDIT.json](OUTPUT/BLINDING-AUDIT.json) records detected leakage opportunities. It does not prove that any specific answer used them. Later package fixes do not retroactively blind this completed session; this audit did not verify that every later package leak was repaired.

## 6. Corrections to the earlier interpretation

| Earlier claim | Corrected interpretation |
|---|---|
| Fully blinded, no skill used | Case labels and corpus metadata leak information; full runtime context was not audited. |
| Clinically correct on all 14 safety cases | Six failures reflect narrow wording rules, but clinical validity has not been independently adjudicated. |
| All diplotype differences are nomenclature | At least some involve unresolved loci or phase and need scientific review. |
| No reusable or reproducible artifact | Real source and data artifacts exist. An exact offline replay succeeds after four supplemental responses. |
| No trust properties | Several are present in partial form; scientific robustness and cache completeness have specific gaps. |
| Skills now sell warranty rather than accuracy | A strategic hypothesis. Current matched evidence must establish any incremental benefit, including accuracy, uncertainty and review effort. |

The earlier decision and artifact entries are retained as history, with an explicit evidence correction and this audit linked. No future decision outcome is graded from this single session.

## 7. Implementation and verification

The audit adapter reuses the frozen scorer and unmodified generated scripts. Eight unit tests cover identifier integrity, missing-row handling, boolean parsing, leakage detection and request replay. Tests failed before implementation; all eight pass after implementation and a request-key normalisation fix.

The adapter has no file or network side effects at import. All outputs use this topic folder. No new server, agent service, database or model subscription was added. Regeneration commands are in [README.md](README.md).

Human review time, manual-equivalent baseline, clinical acceptance and comparative ClawBio benefit were not measured. The originally proposed six-case, three-condition, three-repeat study has not been executed. The supplied session prompted a narrower, evidence-led audit first.
