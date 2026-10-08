# Third revision decision record

Prepared 16 September 2026. Scope: implement the agreed third-revision plan and prepare a local, reviewable package. No journal submission, GitHub push or Zenodo publication is part of this execution.

## Evidence consulted

- Submitted manuscript v90, supplement v48, response v68, cover v12 and declaration v32, preserved byte-for-byte in the confidential bundle's RAW directory.
- The editor's September 2026 decision and Reviewer 2's three comments, retained outside the public-data staging directory.
- Historical tag `agentic-pgx-benchmark-v2.10`, original cases, rule engine, matched runner, baseline and registered equivalence scorers, canonical V4 raw responses and merged panels. Principal input/code hashes are in DATA/source-hashes.json; the staged package has a complete file manifest.
- Existing extraction, caller, cohort, input-normalisation, bootstrap and mixed-model evaluators. All new numerical results derive from stored observations; no new model calls or simulated scientific data were introduced.
- Primary guidelines for RYR1 (10.1002/cpt.1319), CYP2D6 consensus (10.1111/cts.12692), abacavir (10.1038/clpt.2011.355), carbamazepine (10.1002/cpt.1004), allopurinol (10.1038/clpt.2012.209) and MT-RNR1 (10.1002/cpt.2309). Candidate dispositions identify the applicable evidence. All manuscript reference DOIs were resolved directly; medRxiv v2 was checked against the primary version API.

## Decisions and alternatives

The original 110-case comparison remains primary because removing a discovered shared error from the headline would hide an important limitation of reference concordance. The complete 109-case exclusion is supplementary. Assigning zero credit to 24 recommendations is a separate arithmetic sensitivity, not a corrected-reference scorer result. A bespoke RYR1 scoring rule was rejected because it would alter the frozen measurement after seeing the finding.

The negative RYR1 reference is corrected only in a separately versioned evaluator file. Historical prompts, cases, responses and scorers remain intact. Candidate review found one reference error, one superseded CYP2D6 interpretation in retrieval and six terminology mismatches. It does not independently adjudicate the remaining cases, and one error among selected candidates is not a benchmark-wide error-rate estimate.

Both wording definitions appear in the main table and figure. The widened rule is post hoc and its interval includes zero. Claims of a general retrieval safety disadvantage or an information-without-action mechanism based only on the strict rule were removed. High-severity membership is the authors' original operational flag, which includes hearing loss; it is not a standard CPIC severity field.

Table 1 uses attempted evaluations, including the absent record. Table 2 retains the existing evaluator and 2,639 available calls so the original extraction comparison remains interpretable. Its metric is phenotype concordance. Different denominator and normalisation conventions are disclosed rather than silently harmonised.

Opaque identifiers are provided only with model-input fields. References and the historical mapping are separate evaluator artefacts. All 550 historical case/configuration prompts remain byte-identical under identifier replacement because the runner never interpolated the identifier. This does not make rule-supplied configurations blind to their supplied knowledge.

Figure 3B's historical seven-model comparison is retained and its freeze timing explained. The now-available V4 curated panel is used for the matched benchmark and recomputed caller comparison. The non-converged RAG-generation REML fit is labelled diagnostic; its non-convergence is not hidden or interpreted as evidence of a fitted variance mechanism.

## Actions and verification

- Implemented full-grid accounting, input rejection, baseline reproduction, both wording sensitivities, exclusion, zero-credit arithmetic, candidate screen, Table 2 reevaluation, caller reevaluation and prompt identity checks with red/green tests.
- Regenerated Figure 2 and relevant supplementary figures from derived observations. Blue/orange, hatching, direct values, non-overlapping labels, final-size rendering, grayscale and three colour-vision simulations were checked. Visual review prompted explicit tests for a title overlap and a label covering a data point.
- Built clean v91/v49/v69/v13/v33 documents and native tracked-change versions. A regression test ensures deleted and inserted figures retain their respective historical and revised bytes. Accepting changes reproduces clean document XML.
- The document verifier compares 793 table cells against generated outputs, preserves 37 protected paragraphs including run formatting, and checks eight replacement image hashes. Original authorship, author contributions, disclosures, reference titles, historical quotations, prompts and existing DOI identities are retained.
- STAR Methods, Key Resources Table, author contributions, declarations, access restrictions and data/code availability were reviewed. Availability text distinguishes the public historical release from locally staged R3 additions. The AI-assistance disclosure accurately includes Codex's code and editing assistance.
- See REVIEW/pytest-final.txt, numerical-claim-ledger.md and the confidential bundle's document-verification.json, reference-resolution.json and final verification record. Publication and submission remain Manuel's actions.

## Residual limits

The authored rule/reference dependence remains a study limitation, not something the correction removes. Frozen equivalence misses MT-RNR1 and the corrected RYR1 label. Clinical candidate review is evidence-led but is not independent blinded adjudication. The review package is reproducible locally, but a new public version and final author approval are still required before submission.
