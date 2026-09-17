You are running an independent pharmacogenomics study. Work end to end and keep every step reproducible.

## Question

Given a patient genotype description and a named drug, what is the patient's diplotype, CPIC metaboliser (or function) phenotype, and CPIC Level A drug recommendation? How accurately can this be done, and where do errors occur, especially in cases where standard dosing can cause severe or life-threatening harm?

## Inputs

- `cases_blinded.json`: 110 cases, each with `id`, `gene`, `drug`, `genotype`. No answers are included.
- `cpic_guideline_corpus.json`: CPIC guideline text for the 21 genes (snapshot of 2026-05-09). You may also use CPIC, PharmGKB or PharmVar directly, but say which source and version you used for each gene.

## Rules

1. Do not look for or use any published benchmark, paper or repository that reports answers to these cases (including "agentic-pgx-benchmark", ClawBio or Zenodo records by Corpas). If you come across one, stop and say so.
2. Build whatever you think is best: lookup tables, code, or reasoning. Record which method you used.
3. Abstain when the input cannot be mapped with confidence, and explain why. An abstention is better than a guess.

## Required outputs

1. `claude_science_answers.csv` with exactly these columns, one row per case id:
   `case_id,diplotype,phenotype,recommendation,abstained`
   `recommendation` should name the drug and the action, e.g. "codeine: avoid; use a non-CYP2D6 opioid".
   `abstained` is true or false.
2. Any diplotype-to-phenotype and phenotype-to-recommendation tables you built, as CSV.
3. A short methods note: approach, sources and versions, which cases you consider lethal-class and why, and your own estimate of where you are most likely wrong.
4. The code and notebook that produced 1 and 2.
