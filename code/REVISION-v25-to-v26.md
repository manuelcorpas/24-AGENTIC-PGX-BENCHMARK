# Cell Genomics manuscript: v25 to v26 change log

Date: 2026-06-11. Built from the hostile review (`DOCS/BIORXIV-SUBMISSION/REVIEW-hostile-Cell-Genomics-Manuscript-v25.md`). Every new number is bound to the released per-cell data through `numbers_engine.py`, which first reproduces the v25 headline exactly (reproduction gate) before recomputing. v25 files are untouched; v26 are copies.

Files changed:
- `DOCS/BIORXIV-SUBMISSION/Cell-Genomics-Manuscript-v26.docx` (22 edits)
- `DOCS/BIORXIV-SUBMISSION/Cell-Genomics-Supplementary-v26.docx` (nomenclature note + Table S5 + Table S6)

Reproducibility artefacts (this folder, `REVISION-V26/`):
- `numbers_engine.py` -> `v26_numbers.json`, `v26_numbers.txt` (clustered stats + counts; gate passes)
- `mechanism_finding.txt` (HLA vs non-HLA lethal-class breakdown)
- `apply_manuscript_edits.py`, `apply_supplementary_edits.py` (the asserted run-level docx edits)

## New author (your request)
M.Reza Jabalameli added after Mahmoud Aldraimli, with two new affiliations:
13 Broad Institute of Harvard & MIT, Cambridge, MA, USA; 14 Analytical and Translational Genetic Unit, Massachusetts General Hospital, Boston, MA, USA. Added to the author line, the affiliation list, and Author contributions (placeholder role "Writing - review and editing"). New affiliations appended as 13/14 to avoid renumbering the existing 1-12. See AUTHOR TODO.

## M1 Clustered statistics (review's top statistical concern)
The 44,550 cells come from only 110 distinct cases (and 14 distinct lethal-class cases) re-counted across 9 models, 3 population framings of identical genotypes and 3 replicates. v25 reported naive binomial Wilson CIs and two-proportion z-tests that treat these as independent. v26 reports case-level cluster-bootstrap CIs (10,000 resamples) and a model-level sensitivity. Design effects are ~23-43x, so the naive intervals understate uncertainty roughly five-fold. The gradient still separates cleanly.

| Condition | v25 CI (naive Wilson) | v26 CI (case-clustered) |
|-----------|----------------------|--------------------------|
| Free-prompted 80.6% | 79.8-81.4 | 76.1-84.5 |
| Retrieval 89.5% | 88.8-90.1 | 85.1-93.4 |
| Skill-reasoning 95.5% | 95.1-96.0 | 93.1-97.6 |
| Skill-execution 93.3% | 92.7-93.8 | 90.5-95.8 |
| Answer-supplied control 100% | 100-100 | 100-100 |

`Statistical reporting` (STAR Methods) rewritten accordingly; CIs updated at each first mention in Results.

## Safety claim reframed: mechanism-stratified (your approved decision)
v25 claimed retrieval "paradoxically increased lethal-class errors" with `two-proportion z = 6.1, P < 0.001`. That p-value is a pseudoreplication artefact. Clustered over the 14 lethal-class cases the aggregate +12 pp change (24.6% -> 36.6%) is **not significant** (cluster-bootstrap 95% CI -7 to +33 pp; P = 0.25; Wilcoxon signed-rank P = 0.81; 6 cases worse, 8 better). The aggregate is the net of two opposing, mechanism-specific effects, which your own Figure S3 already documents:

- HLA risk-allele loci (5 cases): 13.4% -> 66.1%, **+52.7 pp** (information-without-action).
- Non-HLA metaboliser/other (9 cases): 31.2% -> 20.2%, **-10.9 pp** (retrieval helps).

Changed: Summary, Highlight 2, the retrieval Results paragraph, Figure 3 caption, and STAR Methods now state the mechanism-specific result and remove `z = 6.1 / P < 0.001`. New Table S5 gives the locus-class breakdown. Highlight 1 rewritten (M3): "Executing validated logic, not model reasoning or retrieval, makes the clinical mapping deterministic and auditable" (the measured data show skill-reasoning edges skill-execution on accuracy and lethal-class rate, so execution's case is guarantee/locality, not a lower measured error).

## M4 Reproducibility cluster
- **Dangling reference resolved:** new **Table S6** gives the all-nine-model skill-arm aggregate. Including Mistral on a responding-cells basis leaves the headline unchanged (95.6% reasoning, 93.3% execution); counting its rate-limit errors as wrong lowers it to 85.3% / 83.3% (an API-throughput artefact). This both supplies the promised table and defends the eight-model exclusion as conservative.
- **Mistral 3.7% vs 18% reconciled:** STAR Methods and Limitations now state that 3.7% is the usable (non-error) rate in the primary run (the rest were silent rate-limit errors, no `choices` field), that on responding cells Mistral matched the other models (~97%), and that 18% is the empty rate of the separate paced rerun. Confirmed from the data: the A1=0 Mistral cells contain the literal `<error: 'choices'>`.
- **Skill-arm lethal-class counts added:** "(15.3%, 154 of 1,008 ... versus 8.5%, 86 of 1,008)".
- **"44,550 scored evaluations" -> "44,550 evaluations"** in the Summary (not all attempted cells were scored).

## M5 Construct/COI and population wording
- Significance: "verified across European, Latin American and East African origin individuals" -> "...population framings of identical genotypes", since the curated arm relabels identical genotypes (it tests label-invariance, not cross-population performance on real people). The clinical-grade definition remains stated inline in the Significance.
- Declaration of interests left as the v25 disclosure (already full and explicit).

## M6 Scope/fit calibration
Summary: "the first large-scale controlled evaluation" -> "a large-scale controlled evaluation"; closing sentence "a generalisable, auditable architecture ... at scale" -> "an auditable architecture ... in a codifiable clinical domain, together with the principles by which it should generalise". Title left unchanged (see AUTHOR TODO for the optional scope edit).

## M2 Closed-loop framing
The "Toward 100%" target is now explicitly labelled a prediction that follows by construction from the exact executed mapping, "not a result we run here", with the end-to-end demonstration named as the principal Outlook validation step. The deterministic caller already emits a definite call on a measurable, ancestry-varying fraction (illustratively, on the Peru/AMR cohort the caller is determinate on 70% of observed states and indeterminate on 30%); the full three-cohort closed-loop is the remaining experiment (AUTHOR TODO: EUR/AFR per-state data needed).

## M7 Ethics and data governance (sourced from your papers)
The Ethics paragraph (STAR Methods) was expanded per-cohort with provenance extracted verbatim from the source papers (see `ETHICS-SOURCES.md`):
- European (Corpas family): written consent of all family members for publication of genotype and phenotype; data under a CC0 public-domain dedication (Glusman et al. 2012).
- Latin American (PGP): Declaration of Helsinki; reviewed and approved by the Research and Ethics Committee of the Instituto Nacional de Salud del Peru, authorizations OI-003-11 and OI-087-13 (Guio et al. 2025).
- East African (UGR): General Population Cohort under MRC/UVRI and LSHTM Uganda Research Unit, governed by UVRISEC and UNCST approvals; controlled-access via the UGR Data Access Committee (lead contact S.F., a co-author). Only the study-specific UGR DAC approval number remains to be inserted (see AUTHOR TODO).

## Citations
- Ref 8 (Henricks 2018): the in-text "severe or fatal toxicity in approximately 0.1 to 0.5% of carriers (8)" is not supported by that paper (it reports grade >=3 toxicity, not a 0.1-0.5% carrier fatality rate). Reworded to "severe, sometimes fatal, toxicity, with grade 3 or higher toxicity in a high proportion of DPYD variant carriers (8)".
- Ref 20 (Guio 2025): the "736 individuals across seven subpopulations" is now flagged as an analysed subset of the project's larger sample (the cited paper reports ~1,149 samples).
- Ref 14 (Kara 2026): no change. Earlier suspicion about the `10.64898` DOI prefix was wrong; openRxiv migrated medRxiv DOIs off `10.1101` in late 2025 and it resolves. (Optional: cite the v2 record, 7 authors.)
- Real-genome AMR-vs-AFR `z = 3.5, P < 0.001`: changed to descriptive, with a clustering caveat (denominators are states x models, not independent individuals).

## Figures
Assessed: no regeneration required. The per-model dot-whisker figures (2, 3) show per-model 95% CIs and a point-estimate aggregate tick; the clustering change concerns aggregate inference, now in the text/Methods/Table S5. Figure 3 caption updated to the mechanism-aware framing. The mechanism is shown numerically in the new Table S5 and visually in the existing Figures S3 (information-without-action) and S4. Optional enhancement: add a Figure 3 panel plotting HLA vs non-HLA lethal-class rate by condition (data in Table S5).

## AUTHOR TODO (items only you can complete)
1. **UGR DAC approval number (M7):** the Ethics paragraph is now fully sourced (Corpas CC0 consent; PGP authorizations OI-003-11 and OI-087-13; UGR under UVRISEC and UNCST, accessed via the UGR Data Access Committee). The only field not recorded in any paper is this study's specific UGR DAC approval reference; S.F. (UGR DAC lead contact, co-author) can supply it. Provenance in `ETHICS-SOURCES.md`.
2. **Reference 1 DOI:** post the companion Perspective to bioRxiv so the citation resolves; the clinical-grade definition is already stated inline.
3. **M.R.J. contributions:** confirm/expand the CRediT roles (currently placeholder "Writing - review and editing").
4. **Title (optional, M6):** decide whether to scope "versioned skill libraries" to the single skill / single domain actually evaluated.
5. **EUR/AFR closed-loop (optional, M2):** stage the Corpas-family and Uganda per-state caller calls to run the full three-cohort deterministic-caller closed-loop (the strongest version of the experiment).
6. **PDFs:** export v26 manuscript and supplementary to PDF from Word (no headless converter on this machine), then refresh `Figure*.tiff` if the optional Figure 3 panel is added.
