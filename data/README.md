# Data

The benchmark datasets are archived on Zenodo (not stored in git; excluded via `.gitignore`).
Download and unpack them into this `data/` directory to reproduce the analyses and figures.
Every script reads and writes here; no renaming or relocation step is required.

## Core three conditions + adversarial (DOI 10.5281/zenodo.21526742)

- `v3_raw_rescored_three_arm.json`: locked rescored three-arm dataset (26,730 rows); input for the core-condition analysis and figures
- `v3_adversarial_scrambled.json`: forward adversarial experiment (lethal -> safe corruption)
- `v3_adversarial_reverse.json`: reverse adversarial experiment (safe -> dangerous corruption)
- `v3_rag_genedrug_chunking.json`: (gene, drug)-keyed chunking control (378 cells)
- `v3_three_arm_a2_regression_classified.csv`: drug-substitution classification
- `v3_three_arm_per_case_a1.csv`, `v3_three_arm_lethal_a3_errors.csv`: figure inputs

## Curated benchmark (DOI 10.5281/zenodo.20567742)

- Curated benchmark dataset with no individual-level genotypes.

## Skill conditions: raw evaluations (deposited with the 2026 revision)

The 17,820 raw skill-arm evaluations referenced by the manuscript (skill-reasoning +
skill-execution, 8,910 each) are produced by `code/43-armAB-fullgrid.py` and deposited as:

- `v3_armAB_fullgrid.json`: 17,820-row list (Counter: A_reasoning 8,910, B_execution 8,910)
- `v3_armAB_fullgrid.jsonl`: same, one evaluation per line (crash-safe checkpoint format)
- supporting precursors: `v3_armA9_armBv2.json`, `v3_armA9_armBv2_POP.json`

`code/44-aggregate-five-conditions.py` merges the three-arm dataset with
`v3_armAB_fullgrid.json` into the 44,550-evaluation five-condition dataset;
`code/57-numbers-engine.py` computes the headline numbers and case-level
cluster-bootstrap confidence intervals from those two inputs.

> These skill-arm files were absent from the earlier public deposit; they are added in the
> revision so that all five conditions and all 44,550 evaluations are reproducible. A new
> Zenodo version carrying these files is minted on deposit and its DOI recorded here.
