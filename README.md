# 24-AGENTIC-PGX-BENCHMARK

Reproducibility package for:

**Trustworthy agentic genomics through versioned skill libraries**
Corpas, Iacoangeli, Bourdenx, Aldraimli, Jabalameli, Skene, Fatumo, Guio (2026). Under revision at *Cell Genomics* (Research Article; CELL-GENOMICS-D-26-00551).
Companion to the Perspective *Agentic Genomics: From Pipeline Automation to Autonomous Validation* (CELL-GENOMICS-D-26-00316, under review).

## What this is

A controlled **five-condition constraint gradient** benchmarking nine frontier large
language models on pharmacogenomic interpretation, to locate where in an agentic
pipeline correctness must reside. The gradient progressively moves correctness from
the stochastic model into an executed, versioned skill:

1. **free-prompted** (`no_spec`): model answers unaided
2. **retrieval-augmented** (`cpic_rag`): model answers over retrieved CPIC guideline text
3. **skill-reasoning**: model applies a versioned skill's rules and generates the answer
4. **skill-execution**: model supplies a structured input; the validated skill computes the answer in code
5. **answer-supplied positive control** (`with_spec` in earlier dataset releases): parser/infrastructure control

**Evaluation grid.** Every condition is evaluated on the identical grid: 110 CPIC
Level A cases x nine models x three ancestry framings x three replicates =
**8,910 evaluations per condition, 44,550 in total**. The three core conditions
(free-prompted, retrieval-augmented, answer-supplied) contribute 26,730; the two skill
conditions (skill-reasoning, skill-execution) contribute a further **17,820**.

> **Reproducibility note (2026 revision).** An earlier public state of this repository
> described only the three core conditions (26,730 evaluations). The skill-reasoning and
> skill-execution runners and their 17,820 raw evaluations, referenced in the manuscript,
> are now included here (`code/40-`..`code/46-`, `code/57-numbers-engine.py`) and the raw
> skill-arm evaluations are deposited on Zenodo (see **Data**). See
> `code/MODEL-VERSIONS.md` for the model-identifier audit.

## Repository layout

```
code/      analysis pipeline: runners (all 5 conditions), rescorers, RAG,
           adversarial, chunking, five-condition aggregation, numbers engine, figures
specs/     benchmark inputs: test_cases_v3.json, cpic_rag_corpus_v3.json, concordance_spec.md
figures/   publication figures (PNG 300 dpi + TIFF 600 dpi)
data/      see data/README.md - datasets are archived on Zenodo, not in git
real-genome-arm/  three-population real-diplotype pipeline (PyPGx GRCh37)
```

### Key scripts (`code/`)

Core three conditions and scoring:
- `02-run-benchmark-v3.py`: free-prompted (`no_spec`) and answer-supplied (`with_spec`)
- `15-build-cpic-rag-corpus.py` + `16-run-rag-condition.py`: retrieval-augmented (`cpic_rag`)
- `16b-rag-genedrug-chunking.py`: (gene, drug)-keyed chunking control
- `14-adversarial-scrambled-spec.py` + `14b-adversarial-reverse-spec.py`: bidirectional adversarial
- `10-rescore-v3.py` + `10b-rescore-v3-clinical-equivalence.py`: scoring (baseline + clinical-equivalence layer)

Skill conditions (skill-reasoning, skill-execution): **added in the revision**:
- `41-armA9-armBv2.py`: nine-model skill-reasoning (Arm A) and skill-execution (Arm B, controlled vocabulary); defines the model adapters, skill rules and prompt templates reused downstream
- `43-armAB-fullgrid.py`: the definitive skill-arm run on the identical 9-model x 3-pop x 3-rep x 110-case grid = **17,820** calls
- `40-`, `42-`, `42b-`, `41b-`, `41c-`: skill-arm prototypes, population sweep and model re-run helpers

Aggregation, validation, figures:
- `44-aggregate-five-conditions.py`: merges the three core conditions with the two skill conditions into the 44,550-evaluation five-condition dataset
- `56-validate-manuscript-numbers.py`: recompute-everything firewall against the manuscript numbers
- `57-numbers-engine.py`: headline statistics with **case-level cluster-bootstrap 95% confidence intervals** (source of the reported intervals)
- `30-` .. `36-`, `46-figure-study-design.py`: figures

## Data

Datasets are archived on Zenodo (not stored in git; see `data/README.md`):
- **Curated benchmark dataset** (no individual genotypes): DOI 10.5281/zenodo.20567742
- **Raw rescored three-arm + adversarial datasets**: DOI 10.5281/zenodo.21526742
- **Raw skill-arm evaluations (17,820)**: `v3_armAB_fullgrid.json` / `.jsonl`: deposited in Zenodo v1.2.0, DOI 10.5281/zenodo.21526742 (concept DOI 10.5281/zenodo.20567742 resolves to latest).

## Reproducing the analysis

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
cp .env.example .env   # add your own API keys
# download the Zenodo datasets into data/ (or RESULTS/ as referenced in script headers), then:
python code/44-aggregate-five-conditions.py   # build the five-condition dataset
python code/57-numbers-engine.py              # regenerate headline numbers + cluster-bootstrap CIs
python code/56-validate-manuscript-numbers.py # verify every manuscript number
```

Re-running the benchmark from scratch (the `*-run-*` and `4*-arm*` scripts) issues live
API calls and incurs cost; the locked outputs on Zenodo let you reproduce every figure
and statistic without re-querying the models.

## Citation

If you use this benchmark, please cite the paper above and this repository.

## Declaration of interests

M.C. is the founder of ClawBio, whose SKILL.md specification format is evaluated in this
study. This repository is an independent reproducibility package; it does not depend on
the ClawBio product.

## License

MIT (see `LICENSE`).
