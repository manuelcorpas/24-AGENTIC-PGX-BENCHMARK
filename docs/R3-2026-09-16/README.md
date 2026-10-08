# Third revision review materials

Prepared 16 September 2026 from preserved V4 Flash benchmark outputs. This folder is a local revision, not a published release. The historical code tag `agentic-pgx-benchmark-v2.10` and Zenodo version DOI `10.5281/zenodo.22141182` remain unchanged.

The revised study measures phenotype concordance and recommendation concordance with authored references on text-formatted cases. The executable rules and references share an authored source. A post hoc screen exposed a shared RYR1 reference error. These materials do not establish end-to-end clinical validation.

## Regeneration

From the repository or the staged deposit's repository root, with Python 3.11 or later and the existing analysis dependencies (`numpy`, `matplotlib`, `scipy`, `statsmodels`, `Pillow`, `pytest`, `python-docx`, `lxml`):

```sh
python docs/R3-2026-09-16/SOURCE/01-regenerate.py
python docs/R3-2026-09-16/SOURCE/02-build-audit.py
python -m pytest tests/test_third_revision.py tests/test_scorer_freeze.py tests/test_rescore_matched.py tests/test_extracted_rules_eval.py tests/test_hierarchical_stats.py tests/test_matched_factorial.py tests/test_caller_truth_eval.py tests/test_figure_r1_layout.py tests/test_r3_supplementary_figures.py tests/test_r3_documents.py tests/test_r3_regeneration.py tests/test_r3_audit.py
```

No model calls or network access are required. The generator first reproduces the frozen returned-record results, then constructs the attempted grid and derives the new analyses. The existing REML fit for RAG generation does not converge; its estimates are retained as diagnostics and labelled accordingly. The bootstrap seeds and original evaluators are unchanged.

The confidential document bundle is separate from the public-data staging area. Its `RAW/` folder holds unchanged v90/v48/v68/v12/v32 sources. With that folder available, rebuild clean documents and native tracked-change versions with:

```sh
python docs/R3-2026-09-16/SOURCE/00-build-documents.py --package /path/to/revision-folder
```

The document builder requires the R3 analyses and figures first. Redlines mark changed paragraphs and table rows, and accept to the clean document XML. They are not word-level edits. Redlines retain historical image bytes in deletions and replacement image bytes in insertions. Figure exports and the preserved original documents are also supplied for comparison. Manuscripts, correspondence and original confidential source documents must not be included in a public data deposit by default.

## Files and interpretation

| File | Purpose |
|---|---|
| `DATA/analysis.json` | Frozen reproduction, full-grid results, both wording rules and intervals, exclusion, zero-credit sensitivity, Table 2, eight candidates with all observed responses |
| `DATA/attempted-scored-rows.json` | All 13,200 attempts; the one missing response is explicit and scores zero |
| `DATA/hierarchical-stats.json`, `DATA/crossed-mixed-model.json` | Recomputed supplementary Tables S4 and S5 |
| `DATA/caller-comparison.json` | Existing caller evaluator applied to the V4 panel, 5,272 usable calls |
| `DATA/model-inputs-v1.json` | Opaque identifiers with gene, drug and genotype only |
| `DATA/evaluator-references-original.json` | Original reference fields, keyed by opaque identifier; evaluator only |
| `DATA/evaluator-references-corrected-v1.json` | Separately versioned RYR1 correction, keyed by opaque identifier; evaluator only |
| `DATA/test-cases-corrected-v1.json` | Compatibility copy of the full historical case schema, with only RYR1 phenotype/recommendation changed; **not model input** |
| `DATA/identifier-map.json` | Historical-to-opaque identifier map; evaluator only |
| `DATA/prompt-identity-proof.json` | SHA-256 hashes and equality checks for all 550 case/configuration prompts |
| `DATA/source-hashes.json` | Hashes of principal source data and analysis code |
| `DATA/candidate-dispositions.json`, `CANDIDATE-AUDIT.md` | Source-linked candidate review, examples and scoring limitations |
| `ANALYSIS-DEFINITIONS.md` | Denominators, arithmetic and post hoc definitions |
| `FIGURES/` | Derived Figure 2 and corrected supplementary figure labels/data, PDF/PNG/TIFF, plus grayscale Figure 2 |
| `DECISION-RECORD.md` | Evidence, choices, tradeoffs, changes and verification |

The historical identifiers often reveal the expected phenotype. Do not pass evaluator references, the identifier map or the compatibility case file to a model. The opaque-input change leaves historical instantiated prompts byte-identical because the original runner did not interpolate the case identifier. This does not make authored-rule conditions blind to their supplied rules; that knowledge is part of the experimental configuration.

Publication of a new code/data version and journal submission remain author actions. Draft notice text is supplied in `DRAFT-README-NOTICE.md` and `DRAFT-ZENODO-NOTICE.md`; neither has been posted.
