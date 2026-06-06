# ClawBio Pharmacogenomics Skill Specification (excerpt)

## Task
Given the following genotype data for a single individual, determine:
1. The star allele diplotype for CYP2D6
2. The metaboliser phenotype (Poor/Intermediate/Normal/Ultrarapid)
3. The CPIC recommendation for codeine

## Genotype data
The individual has the following genotypes at CYP2D6 key positions:
- rs3892097 (CYP2D6*4 defining variant): GT = G/A (heterozygous)
- rs1065852 (CYP2D6*10 defining variant): GT = C/C (homozygous reference)
- rs16947 (CYP2D6*2 defining variant): GT = G/A (heterozygous)
- rs1135840: GT = C/G (heterozygous)
- All other CYP2D6 star allele defining positions: homozygous reference

## CPIC Guideline Reference (version 2024-01)
- CYP2D6*1: no variant positions, wild-type (normal function)
- CYP2D6*2: rs16947 G>A (normal function)
- CYP2D6*4: rs3892097 G>A (no function)
- CYP2D6*10: rs1065852 C>T (decreased function)

## Diplotype calling rules
1. Identify which star alleles are present based on defining variants
2. Assign diplotype (two alleles per individual)
3. Map diplotype to activity score: normal=1, decreased=0.5, no function=0
4. Map activity score to phenotype:
   - Activity score 0: Poor Metaboliser
   - Activity score 0.5: Intermediate Metaboliser  
   - Activity score 1.0: Intermediate Metaboliser
   - Activity score 1.5: Normal Metaboliser
   - Activity score 2.0: Normal Metaboliser
   - Activity score >2.0: Ultrarapid Metaboliser

## Codeine CPIC recommendation
- Poor Metaboliser: Avoid codeine. Use non-tramadol, non-codeine analgesic.
- Intermediate Metaboliser: Use codeine with caution at lowest effective dose.
- Normal Metaboliser: Use codeine per standard dosing.
- Ultrarapid Metaboliser: Avoid codeine due to toxicity risk.

## Required output format
Report exactly:
DIPLOTYPE: [allele1]/[allele2]
ACTIVITY_SCORE: [number]
PHENOTYPE: [phenotype]
CODEINE_RECOMMENDATION: [recommendation]
