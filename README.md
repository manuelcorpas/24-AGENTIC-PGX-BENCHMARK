# 24-AGENTIC-PGX-BENCHMARK

Reproducibility package for:

**Trustworthy agentic genomics through versioned skill libraries: deterministic, auditable pharmacogenomics across nine models**
Corpas, Iacoangeli, Bourdenx, Skene, Aldraimli, Fatumo, Guio (2026). Submitted to *Cell Genomics* (Article).
Companion to the Perspective *Agentic Genomics: From Pipeline Automation to Autonomous Validation* (CELL-GENOMICS-D-26-00316, under review).

## What this is

A controlled three-arm benchmark of nine frontier large language models on pharmacogenomic interpretation, testing whether the source of correctness in agent output is the executed specification or the model. Nine models x 110 CPIC Level A cases x 3 populations x 3 conditions x 3 replicates = 26,730 evaluations.

Conditions: `no_spec` (free-prompted), `cpic_rag` (retrieval-augmented from the CPIC corpus), `with_spec` (a versioned SKILL.md-format specification executed as a contract). Headline findings: retrieval augmentation raises phenotype accuracy yet increases lethal-class safety errors (270 to 414) and introduces cross-drug substitution (eliminated by (gene, drug)-keyed indexing, 0/378); specification-constrained execution is deterministic, population-invariant and model-invariant; an adversarial test corrupts the specification in both directions and 90/90 responses execute the corrupted contract, none reverting to the correct answer.

## Repository layout

```
code/      analysis pipeline (run, rescore, RAG, adversarial, chunking, figures)
specs/     benchmark inputs: test_cases_v3.json, cpic_rag_corpus_v3.json, concordance_spec.md
figures/   publication figures (PNG 300 dpi + TIFF 600 dpi), Figures 1-6
data/      see data/README.md - the dataset is archived on Zenodo, not in git
```

Key scripts (`code/`): `02-run-benchmark-v3.py` (no_spec/with_spec), `15-build-cpic-rag-corpus.py` + `16-run-rag-condition.py` (cpic_rag), `16b-rag-genedrug-chunking.py` ((gene,drug) chunking control), `14-adversarial-scrambled-spec.py` + `14b-adversarial-reverse-spec.py` (bidirectional adversarial), `10-rescore-v3.py` + `10b-rescore-v3-clinical-equivalence.py` (scoring), `18-three-arm-analysis.py` (analysis), `30-` to `36-` (figures), `33-classify-a2-regressions.py` (drug-substitution classification).

## Data

The primary dataset (raw model responses, the locked rescored three-arm dataset, adversarial and chunking results, and figure input CSVs) is archived on Zenodo: **DOI: 10.5281/zenodo.20567743**. Download and unpack it into `data/` to reproduce the analyses and figures. See `data/README.md` for the expected files.

## Reproducing the analysis

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
cp .env.example .env   # add your own ANTHROPIC_API_KEY / OPENAI_API_KEY / DEEPSEEK_API_KEY / GEMINI_API_KEY / MISTRAL_API_KEY
# download the Zenodo dataset into data/, then:
python code/18-three-arm-analysis.py          # regenerate the analysis summary
python code/30-figure1-three-arm-aggregate.py # regenerate Figure 1 (and 31- .. 36- for the rest)
```

All API keys are read from the environment; no credentials are stored in this repository. Re-running the benchmark from scratch (the `*-run-*` scripts) issues live API calls and incurs cost; the locked outputs on Zenodo let you reproduce every figure and statistic without re-querying the models.

## Citation

If you use this benchmark, please cite the paper above and this repository.

## Declaration of interests

M.C. is the founder of ClawBio, whose SKILL.md specification format is evaluated in this study. This repository is an independent reproducibility package; it does not depend on the ClawBio product.

## License

MIT (see `LICENSE`).
