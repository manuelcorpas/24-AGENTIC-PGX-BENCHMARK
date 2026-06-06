# Data

The benchmark dataset is archived on Zenodo (not stored in git):

**DOI: 10.5281/zenodo.20567743**

Download the archive and unpack its contents into this `data/` directory. Expected files:

- `v3_raw_rescored_three_arm.json` - locked rescored three-arm dataset (26,730 rows); primary input for analysis and figures
- `v3_adversarial_scrambled.json` - forward adversarial experiment (lethal -> safe corruption)
- `v3_adversarial_reverse.json` - reverse adversarial experiment (safe -> dangerous corruption)
- `v3_rag_genedrug_chunking.json` - (gene, drug)-keyed chunking control
- `v3_three_arm_a2_regression_classified.csv` - drug-substitution classification
- `v3_three_arm_per_case_a1.csv`, `v3_three_arm_lethal_a3_errors.csv` - figure inputs

The figure scripts in `../code/` expect these files at the paths referenced in their headers. The dataset is excluded from git via `.gitignore`.
