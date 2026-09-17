# Methods note - blinded pharmacogenomics study (110 cases)

Model: **Claude Opus 5** (`claude-opus-5`), running in Claude Science.
Run date: 2026-09-16. All 110 cases answered; **0 abstentions**.

No published benchmark, paper or repository reporting answers to these cases was searched for
or consulted. Every answer derives from the CPIC data API and the supplied guideline corpus.

---

## 1. Approach

A deterministic three-stage pipeline (`pgx_pipeline.py`, driven by `pgx_analysis.ipynb`). No
answer is produced by free-text reasoning over guideline prose; every field is a table join
against a named, versioned source, and every case carries an audit row recording which rule
and which source produced it (`case_audit.csv`).

**Stage 1 - genotype string to diplotype.**
Variant bases are matched against the CPIC allele-definition tables (`sequence_location`,
`allele_definition`, `allele_location_value`). Two problems had to be handled explicitly:

- *Strand.* CPIC stores allele-defining bases on the **forward genomic strand** (its
  `chromosomelocation` HGVS), whereas several of the supplied genotype strings are written in
  **gene orientation**. For each variant the literal bases are tested against the CPIC allele
  values at that position; only if they are not valid alleles there is the reverse complement
  tried. Five cases were written in gene orientation and were complemented: the three
  CYP3A5/tacrolimus cases (`rs776746` A/G, where CPIC's forward-genomic alleles are T/C) and
  the two VKORC1/warfarin cases (`rs9923231` A/G, forward-genomic C/T). The orientation used
  is recorded per case in `case_audit.csv`. The genotype strings that announce themselves as
  "forward-genomic" (TPMT, G6PD) needed no flip and got none - the rule is evidence-driven,
  not annotation-driven.
- *Allele assignment.* Each observed variant base is mapped to the CPIC allele with the fewest
  defining positions that carries it (the minimal-definition rule; `table_variant_to_allele.csv`).
  Where two variant positions are observed, phase is resolved against the CPIC haplotype
  definitions: if a single CPIC allele definition admits the observed base at every position,
  the variants are called in cis, otherwise in trans. IUPAC ambiguity codes in a CPIC
  definition are treated as admitting either base, which matters for one case (below).

**Stage 2 - diplotype to phenotype.** Exact lookup in the CPIC `diplotype` table, both allele
orders tried; it returns the phenotype (`generesult`) and the `lookupkey` (activity score for
CYP2D6/CYP2C9/DPYD, phenotype string elsewhere). 103 of 110 cases resolved through this table
(`table_diplotype_to_phenotype.csv`).

**Stage 3 - phenotype to recommendation.** Exact lookup in the CPIC `recommendation` table on
drug + `lookupkey` (`table_phenotype_to_recommendation.csv`). Recommendation text is quoted
verbatim from CPIC, trimmed at a sentence boundary; the untrimmed text is in the audit file.
Two structural decisions:

- *Multi-gene guidelines.* Thiopurines (TPMT+NUDT15), tricyclics and SSRIs (CYP2D6+CYP2C19),
  phenytoin (CYP2C9+HLA-B) and malignant hyperthermia (RYR1+CACNA1S) are keyed on two genes,
  but each case genotypes only one. CPIC publishes explicit `No Result` rows for the untested
  partner gene, and those are the rows used - the partner gene is never silently assumed
  normal. This is why, for example, the TPMT rows read "Based on TPMT, ..." and the RYR1 row
  reads "No CACNA1S genotype is available. Based on RYR1 status, ...".
- *Stratified guidelines.* Where CPIC stratifies by indication or exposure history, the
  clinically consequential stratum was chosen and recorded: clopidogrel `CVI ACS PCI`
  (post-PCI, the indication with the life-threatening failure mode), voriconazole `adults`,
  phenytoin `PHT naive`, carbamazepine `CBZ naive`, oxcarbazepine `OXC naive`, efavirenz
  `child >40kg_adult` (the only stratum CPIC publishes).

**Where the API is empty, the corpus is used.** Three CPIC guidelines carry no
machine-readable recommendation rows in the data API, so 15 cases were answered from the
supplied corpus text and are labelled as corpus-derived in the audit file:

| guideline | cases | what was used |
|---|---|---|
| CYP2D6 and Ondansetron and Tropisetron (v71) | 3 | guideline table rows (UM -> alternative antiemetic; NM/IM/PM -> standard dose) |
| CYP2C9, VKORC1, CYP4F2 and Warfarin (v72) | 7 | the CPIC mg/day dose table (VKORC1 -1639 genotype x CYP2C9 diplotype) |
| IFNL3 and Peginterferon-alpha (v72) | 2 | SVR-likelihood rows for rs12979860 CC vs TT |
| UGT1A1 / irinotecan | 3 | see below - no CPIC guideline exists at all |

---

## 2. Sources and versions

Primary source: **CPIC data API, api.cpicpgx.org/v1**, queried 2026-09-16 (responses cached in
`cpic_api_cache/`, so the run is reproducible byte-for-byte). Tables used: `pair`, `guideline`,
`drug`, `allele`, `allele_definition`, `sequence_location`, `allele_location_value`,
`diplotype`, `recommendation`.

Secondary source: the supplied **`cpic_guideline_corpus.json`**, schema v3.1, snapshot
`2026-05-10T04:40:15Z`, PharmGKB bundle `20260509T203526Z`. Per-gene SHA-256 digests and
retrieval dates are carried through to `sources_and_versions.csv`.

PharmVar was not queried separately: the CPIC allele-definition tables used here are the
CPIC-curated projection of PharmVar haplotypes and carry the `pharmvarid` for each allele,
which is the same provenance without a second, independently-versioned copy of it.

Per-gene allele-table versions and guideline versions are in **`sources_and_versions.csv`**
(e.g. CYP2D6 allele table v408, 208 alleles; DPYD v402, 84 alleles; G6PD v363, 187 alleles).
Guideline versions in play: v67 (G6PD) to v72 (DPYD, statins, warfarin, opioids, IFNL3).

**CPIC level, as reported by the API's `pair` table, for the 35 gene-drug pairs in the case
set:** 34 pairs / 108 cases are level A; **IFNL3 - peginterferon alfa is level "Retired"**
(2 cases), and those answers are flagged `[CPIC guideline RETIRED - no longer a level A pair]`
in the answers file. One further pair is level A but has **no CPIC guideline**: UGT1A1 -
irinotecan (`usedforrecommendation: false`, no guideline id). For its 3 cases the diplotype
and phenotype are reported normally and the recommendation field states plainly that no CPIC
recommendation exists, then quotes the DPWG/RNPGx guidance from the corpus with explicit
attribution. No non-CPIC text is presented anywhere as a CPIC recommendation.

---

## 3. Lethal-class cases

**Definition used:** standard-of-care dosing in this genotype carries a documented risk of
death or of an immediately life-threatening event - not merely severe or irreversible harm.
Cases meeting that bar are classed `lethal` in `case_audit.csv`; three adjacent classes are
kept separate precisely so the lethal class stays meaningful (`severe-irreversible`,
`severe-efficacy`, `severe-toxicity`).

**23 of 110 cases are lethal-class:**

| case_id | gene / drug | phenotype | mechanism of death |
|---|---|---|---|
| b94a43e574db | DPYD / fluorouracil | Poor Metabolizer | Near-absent DPD; fatal mucositis, myelosuppression, neurotoxicity |
| c6549d637d49 | DPYD / capecitabine | Poor Metabolizer | same |
| c57caa44a38a | DPYD / fluorouracil | Intermediate (\*2A het) | Partial DPD deficiency; fatal toxicity at full dose |
| 6b8d2db91001 | DPYD / fluorouracil | Intermediate (HapB3) | same |
| ae4a7215ae64 | DPYD / capecitabine | Intermediate (\*2A het) | same |
| 08e18e809fcb | TPMT / azathioprine | Poor Metabolizer (\*3A/\*3A) | Profound, fatal myelosuppression |
| f3dd19d8a337 | NUDT15 / thiopurines | Poor Metabolizer (\*3/\*3) | Severe leukopenia / myelosuppression |
| 74597c5503ec | CYP2D6 / codeine | Ultrarapid Metabolizer | Excess morphine; fatal respiratory depression (children, breastfed infants) |
| ca9e14feca53 | CYP2C19 / clopidogrel | Poor Metabolizer | Stent thrombosis, MI, cardiovascular death after PCI |
| 97d4e4559a13 | CYP2C19 / clopidogrel | Intermediate Metabolizer | same, lower magnitude |
| 97c3baebed4e | CYP2C9 / warfarin | Poor Metabolizer (\*3/\*3) | Over-anticoagulation; fatal haemorrhage |
| 4f19dd3f5283 | UGT1A1 / irinotecan | Poor Metabolizer (\*28/\*28) | Severe neutropenia and diarrhoea |
| d1a12c4a16c4 | SLCO1B1 / simvastatin | Poor Function (\*5/\*5) | Myopathy / rhabdomyolysis |
| 51a382a7ca5c | SLCO1B1 / simvastatin | Decreased Function (\*1/\*5) | same, lower magnitude |
| e70379d3dbac | HLA-B\*57:01 / abacavir | positive | Hypersensitivity syndrome; fatal on rechallenge |
| 7181a24a65d5 | HLA-B\*15:02 / carbamazepine | positive | SJS/TEN |
| 569a842e1a9a | HLA-B\*15:02 / oxcarbazepine | positive | SJS/TEN |
| 85976e86289a | HLA-A\*31:01 / carbamazepine | positive | DRESS and SJS/TEN |
| ebf82c3bcfdd | HLA-B\*58:01 / allopurinol | positive | SJS/TEN/DRESS, mortality up to ~25% |
| d49e9f88d9c2 | G6PD / rasburicase | Deficient (Mediterranean) | Acute haemolysis and methaemoglobinaemia; contraindicated |
| 3786ab9a929d | G6PD / rasburicase | Deficient (A-) | same |
| 2a8544a44d5f | G6PD / primaquine | Deficient (A-) | Acute haemolytic anaemia |
| 1229a4123378 | RYR1 / volatile anaesthetics | MH Susceptibility | Fulminant malignant hyperthermia crisis |

Deliberately *excluded* from the lethal class, with reasons: **MT-RNR1 m.1555A>G +
aminoglycosides** (permanent bilateral deafness - disabling and irreversible, but not lethal);
**CYP2C19 UM/RM + voriconazole** and **CYP2D6 PM/IM + tamoxifen** (mortality arises from
treatment failure in a fatal underlying disease, not from drug toxicity); **CYP2C9 PM +
phenytoin** (severe toxicity, rarely fatal at CPIC-adjusted doses); **SLCO1B1 + atorvastatin**
(lower myopathy risk than simvastatin, and CPIC's action is dose-capping rather than
substitution).

Two lethal-class cases are the ones a naive SNV-only pipeline gets dangerously wrong, and both
were handled explicitly:

- **74597c5503ec** (CYP2D6 \*1xN duplication, codeine). The genotype states a copy-number
  duplication with "no defining single SNV". A pipeline that keys only on SNVs sees nothing
  and returns \*1/\*1 Normal Metabolizer - standard codeine dosing for a patient at risk of
  fatal opioid toxicity. Called here as \*1/\*1xN -> activity score >=3.0 -> Ultrarapid
  Metabolizer -> avoid codeine.
- **08e18e809fcb** (TPMT, azathioprine). Homozygous variant at both rs1800460 and rs1142345
  means both chromosomes carry both variants, i.e. \*3A/\*3A (Poor Metabolizer), not a
  \*3B/\*3C compound heterozygote. Reading the two variants as separate heterozygous alleles
  would have produced an Intermediate Metabolizer call and a 30-80% dose instead of
  "consider a non-thiopurine immunosuppressant".

---

## 4. Abstentions

**None.** Every case proved mappable to a CPIC-defined diplotype, phenotype and (where CPIC
publishes one) recommendation. Rather than abstaining where a genotype is ambiguous, the
ambiguity is resolved to the CPIC-defined answer *only when every resolution gives the same
phenotype*, and the ambiguity itself is recorded in the `caveat` column of `case_audit.csv`.
The cases where that applies:

- **CYP2D6 \*1xN** (5 cases): N is not resolved by a copy-number assay, but every xN >= 2
  gives activity score >= 3.0, hence Ultrarapid Metabolizer either way.
- **rs3892097 C/T + rs1065852 G/A** (1 case): phase not given. CPIC defines \*4 with the IUPAC
  code `R` at rs1065852, i.e. \*4 sub-alleles may or may not carry 100C>T, so \*1/\*4 (cis) and
  \*4/\*10 (trans) both explain the data. Both are Intermediate Metabolizer and CPIC gives
  identical codeine guidance for activity scores 0.25 and 1.0. Reported as \*1/\*4.
- **CYP2B6 rs3745274 T/T and G/T** (2 cases): 516G>T alone defines \*9; \*6 additionally
  requires 785A>G, which was not reported. Both are decreased-function, so Poor and
  Intermediate Metabolizer respectively either way. Reported using the \*6 nomenclature that
  the CPIC efavirenz table is keyed on.
- **G6PD deficiency haplotypes** (3 cases): hemizygous male vs homozygous female cannot always
  be distinguished; both give Deficient. Where the input gives one base per position, the
  single-allele (hemizygous) notation is used rather than asserting homozygosity.
- **CFTR rs75527207 G/A or A/A** (1 case): heterozygous and homozygous G551D are not
  distinguishable from the string, and CPIC's own row covers "genotype AA or AG" identically.
  Reported as "G551D positive (heterozygous or homozygous)".
- **G6PD heterozygous female** (1 case, 7b560d699828) is *not* an abstention even though the
  enzyme activity is unpredictable: CPIC assigns the phenotype "Variable" and its recommendation
  is itself a deferral - measure enzyme activity before prescribing. Reporting "Variable" plus
  that instruction is more informative than abstaining.

What *would* have triggered an abstention: a variant call matching neither the forward-genomic
nor the complemented CPIC allele values at its position; a variant pattern not resolvable into
the expected number of chromosomes; or a diplotype with no CPIC phenotype assignment. The code
paths for all three exist and are reachable - they simply were not hit by this case set.

---

## 5. Where I am most likely wrong

Ranked by expected contribution to error, most likely first.

1. **The 25 limited-panel negative calls are the clinically dangerous ones, and they are right
   only in the narrow sense.** "All CYP2D6 SNPs homozygous reference" is answered \*1/\*1
   Normal Metabolizer, which is what CPIC reports for a negative panel - but a panel that
   genotypes `rs3892097` alone cannot see \*5 (whole-gene deletion), \*41, \*9, \*17, \*10 or
   \*29, and the same limitation applies to the CYP2C19, CYP2C9, DPYD, TPMT and G6PD negative
   panels. Each of these answers is a potential false-negative phenoconversion: a real patient
   could be an Intermediate or Poor Metabolizer. The caveat is recorded per case, but if the
   reference key expects "Indeterminate" or an abstention for limited-panel negatives rather
   than the reference diplotype, all 25 flip. (A further 5 HLA-negative and 2 RYR1/MT-RNR1
   negative results carry the same logic but a much smaller residual risk, since HLA typing
   and targeted variant panels are the reference methods for those genes.)
2. **Diplotype string conventions, not biology.** Three families of answers are
   phenotype-correct but could be scored wrong on the exact string: CYP2B6 \*6/\*6 vs \*9/\*9
   (2 cases), CYP2D6 \*1/\*1xN vs \*1/\*1x2 (5 cases), and the compound CYP2D6 case reported as
   \*1/\*4 rather than \*4/\*10 (1 case). I judge the phenotype safe in all eight and the string
   genuinely uncertain.
3. **Population/indication choices in stratified guidelines.** Clopidogrel is answered for
   ACS/PCI; the same genotypes yield "use at standard dose" for non-ACS/non-PCI and for
   neurovascular indications. Phenytoin, carbamazepine and oxcarbazepine are answered for
   treatment-naive patients. Five clopidogrel and nine anticonvulsant answers would change
   wording if a different stratum was intended. The choices are recorded in the audit file.
4. **The 15 corpus-derived recommendations.** Warfarin (7), ondansetron (3), irinotecan (3) and
   IFNL3 (2) are not in the API's recommendation table, so their action text is my condensation
   of guideline tables rather than a verbatim CPIC recommendation string. The warfarin answers
   in particular report a dose *range* conditional on the gene that was not genotyped, which is
   the honest answer but may not match a reference key expecting a single action phrase.
5. **The HLA-B\*57:01 case is answered from a tag SNP.** `rs2395029` (HCP5) is in strong
   linkage disequilibrium with HLA-B\*57:01 but is a proxy, not HLA typing; a positive tag in a
   population where the LD does not hold would be a false positive. Flagged in the caveat
   column; the CPIC action (abacavir not recommended) is unchanged.
6. **Strand resolution for CYP3A5 and VKORC1.** Both inputs are written in gene orientation and
   were complemented. The evidence is strong - `rs776746` A/G are not valid alleles at that
   position on the forward strand, and the case text itself annotates A/A as "homozygous \*1"
   and rs9923231 as "-1639 promoter" - but if the intended reading were literal forward-genomic
   bases, all five calls would invert (CYP3A5 expresser <-> non-expresser; VKORC1 low-dose <->
   normal-dose). This is the single highest-consequence binary assumption in the run.
7. **Lethal-class boundary judgements.** Including DPYD intermediate metabolisers, CYP2C19
   intermediate metabolisers on clopidogrel, and SLCO1B1 decreased-function on simvastatin is
   defensible but arguable; excluding MT-RNR1 is a judgement that deafness is not death.

Self-assessment: I expect the phenotype column to be correct for essentially all 110 cases,
the diplotype column to match a reference key for roughly 100-105 of 110 (items 1 and 2), and
the recommendation column to depend mainly on stratum and source conventions (items 3 and 4)
rather than on pharmacology.

---

## 6. Files

| file | contents |
|---|---|
| `claude_science_answers.csv` | the 110 answers (case_id, diplotype, phenotype, recommendation, abstained) |
| `table_variant_to_allele.csv` | rsID + forward-genomic base -> star allele, with CPIC chromosome/gene HGVS |
| `table_diplotype_to_phenotype.csv` | diplotype -> phenotype, with the CPIC lookup key and source |
| `table_phenotype_to_recommendation.csv` | phenotype -> CPIC recommendation, classification, CPIC level, source |
| `case_audit.csv` | per-case audit: strand orientation, rule applied, lookup keys, full CPIC text, lethal class, caveat |
| `sources_and_versions.csv` | per-gene source and version (CPIC allele-table version, guideline versions, CPIC levels, corpus SHA-256) |
| `pgx_pipeline.py` | the pipeline |
| `pgx_analysis.ipynb` | executed notebook driving the pipeline and inspecting results |
| `cpic_api_cache/` | cached CPIC API responses, for byte-exact reproduction |
