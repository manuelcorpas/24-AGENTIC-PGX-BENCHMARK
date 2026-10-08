# Post hoc candidate audit, 16 September 2026

Screen: all 110 original-reference cases. Select baseline phenotype means >=0.95 in each of RAG-assisted execution, authored-rule generation and authored-rule execution, and <=0.10 in RAG generation. Eight cases meet the rule. This is a full-case screen followed by candidate review, not independent clinical adjudication of the benchmark. Dispositions concern the observed phenotype discrepancy; they are not blanket approval of every recommendation.

Outcome: one confirmed reference error (RYR1), one outdated retrieval interpretation (CYP2D6 activity score 1.0), and six terminology mismatches. No scoring patterns were changed. All 120 observed responses per candidate, including 24 from each configuration, are retained in DATA/analysis.json. The following examples are the fixed GPT-5.2 replicate 0 pair; complete frequency counts prevent those examples being mistaken for exhaustive response coverage.

Evidence: primary guideline pages were read in the browser, with DOI resolution attempted directly. The CYP2D6 consensus was checked at its Wiley DOI page after resolver access failed. The source-linked dispositions below are reasoned candidate review. RYR1, CYP2D6-amitriptyline and HLA-B*57:01-positive were also hand-checked against the preserved raw responses and references. No new model calls were made.

The terminology cases have frozen equivalence totals 24, 24, 23, 22, 19 and 0 of 24 in the order reported in the manuscript. MT-RNR1 receives zero because the frozen equivalence layer does not recognise that risk phrasing. A terminology disposition does not imply that every response parses or earns equivalence credit.

## Case: cyp2d6_amitrip_im

Gene/drug: CYP2D6 / amitriptyline. Input: rs3892097 C/T. Authored diplotype: *1/*4. Original phenotype: Intermediate Metaboliser.

Disposition: **Outdated retrieval interpretation**. The CPIC/DPWG consensus assigns activity score 1.0 to intermediate metaboliser. Most retrieved responses use the superseded normal/extensive band. Retain the phenotype reference. This disposition does not independently adjudicate the complete drug recommendation.

Source: https://doi.org/10.1111/cts.12692

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.7292 |
| rag_generation | 0.0417 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 1/24 (0.04).

Observed RAG phenotype strings and frequencies:

- 1/24: CYP2D6 Normal metabolizer (activity score ~1.0; per CPIC amitriptyline translation, AS 1.0–2.0 = normal)
- 1/24: CYP2D6 Normal metabolizer (activity score ~1.0; CPIC “normal or intermediate—controversy remains”)
- 5/24: CYP2D6 Normal metabolizer
- 5/24: Normal metabolizer
- 4/24: CYP2D6 normal metabolizer (activity score 1.0)
- 1/24: CYP2D6 normal metabolizer
- 1/24: CYP2D6 Normal Metabolizer (Activity Score 1.0 - one normal function allele with value 1.0 + one no function allele with value 0 = 1.0, which falls within the normal metabolizer range of 1.0-2.0)
- 3/24: [empty]
- 1/24: ** CYP2D6 Normal Metabolizer (activity score 1.0-2.0)
- 1/24: **
- 1/24: CYP2D6 intermediate metabolizer (activity score 1.0)

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: CYP2D6 *1/*4  
PHENOTYPE: CYP2D6 Normal metabolizer (activity score ~1.0; per CPIC amitriptyline translation, AS 1.0–2.0 = normal)  
DRUG: Initiate amitriptyline with the recommended starting dose.  
HAZARD: No increased pharmacogenomic risk expected (normal CYP2D6-mediated clearance of TCAs).
```
Parsed or executed phenotype: CYP2D6 Normal metabolizer (activity score ~1.0; per CPIC amitriptyline translation, AS 1.0–2.0 = normal)
Parsed or executed recommendation: Initiate amitriptyline with the recommended starting dose.

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: *1/*4
```
Parsed or executed phenotype: Intermediate Metaboliser
Parsed or executed recommendation: amitriptyline: 25% dose reduction; monitor

## Case: hlab5701_aba_pos

Gene/drug: HLA-B*57:01 / abacavir. Input: rs2395029 G/T or G/G (any positive). Authored diplotype: HLA-B*57:01 positive. Original phenotype: Positive.

Disposition: **Terminology mismatch**. HLA-B*57:01 carriage corresponds to high abacavir hypersensitivity risk. The case reference uses Positive, while retrieved responses describe risk. Retain reference; report baseline and frozen equivalence together.

Source: https://doi.org/10.1038/clpt.2011.355

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.2500 |
| rag_generation | 0.0000 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 24/24 (1.00).

Observed RAG phenotype strings and frequencies:

- 18/24: High risk of hypersensitivity
- 6/24: High risk of abacavir hypersensitivity

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*57:01/*X  
PHENOTYPE: High risk of hypersensitivity  
DRUG: Abacavir is not recommended  
HAZARD: Significantly increased risk of abacavir hypersensitivity reaction
```
Parsed or executed phenotype: High risk of hypersensitivity
Parsed or executed recommendation: Abacavir is not recommended

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*57:01 positive
```
Parsed or executed phenotype: Positive
Parsed or executed recommendation: abacavir: AVOID (potentially LETHAL hypersensitivity reaction with rechallenge)

## Case: hlab5701_aba_neg

Gene/drug: HLA-B*57:01 / abacavir. Input: rs2395029 T/T (no risk allele). Authored diplotype: HLA-B*57:01 negative. Original phenotype: Negative.

Disposition: **Terminology mismatch**. Absence of HLA-B*57:01 corresponds to very low abacavir hypersensitivity risk. The case reference uses Negative. Retain reference; the frozen equivalence layer does not recognise every risk formulation.

Source: https://doi.org/10.1038/clpt.2011.355

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.2500 |
| rag_generation | 0.0417 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 19/24 (0.79).

Observed RAG phenotype strings and frequencies:

- 18/24: Very low risk of hypersensitivity
- 5/24: Very low risk of abacavir hypersensitivity
- 1/24: Very low risk of abacavir hypersensitivity (HLA-B*57:01 negative)

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: *X/*X (HLA-B*57:01 negative)
PHENOTYPE: Very low risk of hypersensitivity
DRUG: Use abacavir per standard dosing guidelines
HAZARD: Low or reduced risk of abacavir hypersensitivity reaction
```
Parsed or executed phenotype: Very low risk of hypersensitivity
Parsed or executed recommendation: Use abacavir per standard dosing guidelines

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*57:01 negative
```
Parsed or executed phenotype: Negative
Parsed or executed recommendation: abacavir: standard dosing

## Case: hlaa3101_carba_neg

Gene/drug: HLA-A*31:01 / carbamazepine. Input: HLA-A typing negative for *31:01. Authored diplotype: HLA-A*31:01 negative. Original phenotype: Negative.

Disposition: **Terminology mismatch**. The HLA-A*31:01-negative row describes normal or reduced risk of carbamazepine-associated cutaneous reactions. The reference uses Negative. This does not establish absence of risk from other alleles. Retain reference and report the frozen equivalence result.

Source: https://doi.org/10.1002/cpt.1004

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.2083 |
| rag_generation | 0.0833 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 24/24 (1.00).

Observed RAG phenotype strings and frequencies:

- 1/24: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS, and MPE (with respect to HLA-A*31:01)
- 9/24: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS and MPE
- 4/24: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS, and MPE
- 2/24: Normal/reduced risk for carbamazepine-induced SJS/TEN, DRESS, MPE
- 1/24: Normal/reduced risk of carbamazepine-induced SJS/TEN, DRESS and MPE
- 1/24: Normal – non‐carrier of HLA-A*31:01, normal risk of SJS/TEN, DRESS and MPE
- 1/24: Normal risk for carbamazepine-induced SJS/TEN, DRESS and MPE
- 1/24: Non-carrier (normal risk)
- 1/24: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS and MPE.
- 1/24: Homozygous for alleles other than HLA-A*31:01
- 1/24: Homozygous for alleles other than HLA-A*31:01 (Non-carrier)
- 1/24: Homozygous for alleles other than HLA-A*31:01; Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS and MPE

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: *X/*X (HLA-A*31:01 non-carrier; *31:01 negative)
PHENOTYPE: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS, and MPE (with respect to HLA-A*31:01)
DRUG: Use carbamazepine per standard dosing guidelines
HAZARD: No increased HLA-A*31:01–associated risk of severe cutaneous adverse reactions (SJS/TEN, DRESS, MPE)
```
Parsed or executed phenotype: Normal or reduced risk of carbamazepine-induced SJS/TEN, DRESS, and MPE (with respect to HLA-A*31:01)
Parsed or executed recommendation: Use carbamazepine per standard dosing guidelines

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-A*31:01 negative
```
Parsed or executed phenotype: Negative
Parsed or executed recommendation: carbamazepine: standard dosing

## Case: hlab5801_allo_pos

Gene/drug: HLA-B*58:01 / allopurinol. Input: HLA-B typing positive for *58:01. Authored diplotype: HLA-B*58:01 positive. Original phenotype: Positive.

Disposition: **Terminology mismatch**. HLA-B*58:01 positivity is associated with substantially increased risk of allopurinol-associated severe cutaneous reactions. Retrieval describes risk rather than the reference shorthand Positive. Retain reference and both scores.

Source: https://doi.org/10.1038/clpt.2012.209

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.2500 |
| rag_generation | 0.0000 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 23/24 (0.96).

Observed RAG phenotype strings and frequencies:

- 15/24: Significantly increased risk of allopurinol SCAR
- 5/24: Significantly increased risk of allopurinol-induced SCAR
- 1/24: Significantly increased risk of allopurinol‐induced SCAR
- 3/24: Significantly increased risk of allopurinol-induced severe cutaneous adverse reactions (SCAR)

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*58:01/*X  
PHENOTYPE: Significantly increased risk of allopurinol SCAR  
DRUG: Allopurinol is contraindicated (avoid use; select an alternative)  
HAZARD: Significantly increased risk of severe cutaneous adverse reactions (SCAR) including SJS/TEN
```
Parsed or executed phenotype: Significantly increased risk of allopurinol SCAR
Parsed or executed recommendation: Allopurinol is contraindicated (avoid use; select an alternative)

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*58:01 positive
```
Parsed or executed phenotype: Positive
Parsed or executed recommendation: allopurinol: AVOID (potentially LETHAL SJS/TEN/DRESS risk)

## Case: hlab5801_allo_neg

Gene/drug: HLA-B*58:01 / allopurinol. Input: HLA-B typing negative for *58:01. Authored diplotype: HLA-B*58:01 negative. Original phenotype: Negative.

Disposition: **Terminology mismatch**. HLA-B*58:01 negativity is described as low or reduced risk of severe cutaneous reactions, rather than the reference shorthand Negative. A negative test does not eliminate risk. Retain reference and both scores.

Source: https://doi.org/10.1038/clpt.2012.209

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.2083 |
| rag_generation | 0.0000 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 22/24 (0.92).

Observed RAG phenotype strings and frequencies:

- 18/24: Low or reduced risk of allopurinol SCAR
- 2/24: Low or reduced risk of allopurinol-induced SCAR
- 2/24: Low or reduced risk of allopurinol‐induced SCAR
- 2/24: Low or reduced risk of allopurinol-induced severe cutaneous adverse reactions (SCAR)

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: *X/*X (no *58:01 detected)  
PHENOTYPE: Low or reduced risk of allopurinol SCAR  
DRUG: Use allopurinol per standard dosing guidelines  
HAZARD: Low or reduced risk of allopurinol-induced severe cutaneous adverse reactions (SCAR)
```
Parsed or executed phenotype: Low or reduced risk of allopurinol SCAR
Parsed or executed recommendation: Use allopurinol per standard dosing guidelines

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: HLA-B*58:01 negative
```
Parsed or executed phenotype: Negative
Parsed or executed recommendation: allopurinol: standard dosing

## Case: mtrnr1_amino_var

Gene/drug: MT-RNR1 / aminoglycosides. Input: m.1555A>G (homoplasmic mitochondrial variant). Authored diplotype: m.1555A>G. Original phenotype: Variant Carrier.

Disposition: **Terminology mismatch**. The MT-RNR1 m.1555A>G variant confers increased risk of aminoglycoside-associated hearing loss. Retrieval uses that risk phenotype and the reference uses Variant Carrier. The frozen equivalence layer gives no credit; retain and disclose both results rather than extend the scorer.

Source: https://doi.org/10.1002/cpt.2309

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.0000 |
| rag_generation | 0.0000 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 0/24 (0.00).

Observed RAG phenotype strings and frequencies:

- 18/24: MT-RNR1 increased risk of aminoglycoside-induced hearing loss
- 6/24: Increased risk of aminoglycoside-induced hearing loss

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: m.1555A>G (homoplasmic)
PHENOTYPE: MT-RNR1 increased risk of aminoglycoside-induced hearing loss
DRUG: Avoid aminoglycoside antibiotics unless infection severity and lack of safe/effective alternatives outweigh the high risk; if unavoidable, use lowest dose/shortest duration with TDM and frequent hearing assessments.
HAZARD: Very high risk of permanent aminoglycoside-induced hearing loss (ototoxicity)
```
Parsed or executed phenotype: MT-RNR1 increased risk of aminoglycoside-induced hearing loss
Parsed or executed recommendation: Avoid aminoglycoside antibiotics unless infection severity and lack of safe/effective alternatives outweigh the high risk; if unavoidable, use lowest dose/shortest duration with TDM and frequent hearing assessments.

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: m.1555A>G
```
Parsed or executed phenotype: Variant Carrier
Parsed or executed recommendation: aminoglycosides: AVOID (potentially LETHAL ototoxicity — irreversible deafness)

## Case: ryr1_anaes_norm

Gene/drug: RYR1 / volatile-anaesthetics. Input: RYR1 no pathogenic variant. Authored diplotype: RYR1 wildtype. Original phenotype: Normal.

Disposition: **Confirmed reference error**. CPIC Tables 1 and 2 classify a variant-negative result as uncertain susceptibility; clinical findings and family history must inform anaesthetic use. Normal and unconditional standard use are not supported. Correct only the separately versioned reference, preserving historical scoring.

Source: https://doi.org/10.1002/cpt.1319

| Configuration | Baseline phenotype mean |
|---|---:|
| free_generation | 0.3542 |
| rag_generation | 0.0000 |
| rag_execution | 1.0000 |
| skill_generation | 1.0000 |
| skill_execution | 1.0000 |

Frozen RAG equivalence: 0/24 (0.00).

Observed RAG phenotype strings and frequencies:

- 1/24: Uncertain susceptibility (malignant hyperthermia risk not excluded)
- 1/24: Uncertain malignant hyperthermia susceptibility
- 14/24: Uncertain susceptibility
- 3/24: Uncertain malignant-hyperthermia susceptibility
- 2/24: Uncertain susceptibility to malignant hyperthermia
- 1/24: Uncertain susceptibility to malignant hyperthermia (no pathogenic RYR1 variant detected)
- 1/24: Uncertain susceptibility to malignant hyperthermia (no pathogenic variant detected)
- 1/24: Malignant Hyperthermia Uncertain Susceptibility

Example rag_generation, GPT-5.2, replicate 0:

```text
DIPLOTYPE: RYR1 negative for known malignant hyperthermia–causative variant (no pathogenic variant detected)
PHENOTYPE: Uncertain susceptibility (malignant hyperthermia risk not excluded)
DRUG: Use clinical findings, family history, and/or additional testing to guide use of halogenated volatile anesthetics; consider non-triggering alternatives if concern remains
HAZARD: Potential risk of malignant hyperthermia with volatile anesthetics (risk not eliminated by a negative result)
```
Parsed or executed phenotype: Uncertain susceptibility (malignant hyperthermia risk not excluded)
Parsed or executed recommendation: Use clinical findings, family history, and/or additional testing to guide use of halogenated volatile anesthetics; consider non-triggering alternatives if concern remains

Example skill_execution, GPT-5.2, replicate 0:

```text
DIPLOTYPE: RYR1 wildtype
```
Parsed or executed phenotype: Normal
Parsed or executed recommendation: volatile-anaesthetics: standard use
