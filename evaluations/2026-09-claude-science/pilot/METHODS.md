# Methods note — blinded CPIC pharmacogenomics study, 110 cases

Run date 2026-09-15. All outputs in this directory were produced by `run_study.py`
(`python run_study.py`), which imports `pgx_pipeline.py`. `study_notebook.ipynb` reproduces
the run and the summary counts.

## 1. Approach

The task decomposes into three lookups, and I built each one as a table rather than answering
case by case from memory:

    genotype ──(A) constraint solving──▶ diplotype ──(B) CPIC table──▶ phenotype ──(C) CPIC table──▶ recommendation

**(A) Diplotype calling.** Star alleles were called by constraint solving over the CPIC
allele-definition table, the same principle PharmCAT uses. For a gene with assayed loci *L* and
an observed unordered genotype at each locus, the solver enumerates every pair of
non-structural CPIC alleles whose definitions can reproduce the observation at every locus in
*L*, then ranks the survivors by

1. number of **unconfirmed defining variants** — defining positions of the candidate allele that
   were not assayed and therefore cannot be checked (prefer 0);
2. number of non-reference alleles (prefer the reference allele);
3. total number of defining variants.

Two details mattered. CPIC allele definitions use IUPAC ambiguity codes at positions where
sub-alleles disagree (CYP2D6\*4 carries `R` at rs1065852), so matching and the
"definitely differs from reference" test both operate on expanded base sets, and a position whose
ambiguity code still permits the reference base is not counted as a defining variant. CPIC also
lists hybrid haplotypes with no positional definition at all (CYP2D6 `*13+*1`); these were
excluded, because they match any observation vacuously.

After this ranking, **every one of the 83 solver-called cases resolved to a unique
most-parsimonious diplotype** — no case had a second solution with the same assumption burden
and the same number of variant alleles.

**(B) and (C)** are pure lookups against CPIC's own tables — no interpretation of guideline prose.
Phenotype comes from the CPIC `diplotype` table (`diplotype → generesult` plus the `lookupkey`
CPIC uses internally, which is an activity score for CYP2D6/CYP2C9/DPYD and a phenotype string
elsewhere). The recommendation comes from the CPIC `recommendation` table matched on that
`lookupkey`.

**Multi-gene guidelines.** Four guidelines in this case set are keyed on two genes
(amitriptyline: CYP2D6 + CYP2C19; phenytoin: CYP2C9 + HLA-B; carbamazepine: HLA-B + HLA-A;
volatile anaesthetics: RYR1 + CACNA1S; thiopurines: TPMT + NUDT15) while each case types only
one of them. Rather than assume the untyped gene is normal, I used CPIC's own `No Result` key,
which exists for exactly this situation. The `ryr1_anaes_norm` recommendation consequently opens
with "No CACNA1S genotype is available" — that is CPIC's text, not mine.

**Non-SNV inputs** (27 cases: 10 HLA, 6 UGT1A1 repeat, 5 CYP2D6 duplication, 2 each for CFTR,
MT-RNR1 and RYR1) were handled by explicit documented rules: CYP2D6 copy number, the UGT1A1
TATA-box (TA)n repeat
((TA)6 = \*1, (TA)7 = \*28), HLA serotype status, the m.1555A>G mitochondrial variant, RYR1
MH-causative variant status, and CFTR G551D / F508del. Each rule and its reasoning is in the
`rationale` column of `claude_science_answers_annotated.csv`.

### Assumptions, stated explicitly

- **A1 — unmentioned panel loci are reference.** Case descriptions name only non-reference loci.
  A targeted panel reports every locus it assays, so a locus covered by the panel and not flagged
  is homozygous reference. I therefore defined each gene's assayed panel as the union of loci
  mentioned across that gene's cases and filled the unmentioned ones with reference. This is
  load-bearing: without it, `cyp2c19_*_rm` (rs12248560 C/T alone) would be compatible with
  \*2/\*17 (Intermediate Metabolizer) as well as \*1/\*17 (Rapid Metabolizer), which changes the
  clopidogrel recommendation.
- **A2 — two genes are reported in transcript orientation.** The CYP3A5 and VKORC1 cases give
  alleles that are the exact complement of the dbSNP forward-strand alleles CPIC uses
  (`rs776746 A/A = *1/*1`, i.e. 6986A>G orientation; `rs9923231 A = -1639A`, i.e. the
  guideline's G>A orientation). Both genes are on the minus strand. The case text for TPMT and
  G6PD says "forward-genomic" and those alleles matched CPIC directly with no flip.
- **A3 — G6PD hemizygous reading.** Three G6PD cases give a male hemizygous genotype and then its
  female homozygous equivalent. I took the hemizygous reading (stated first, unambiguous); both
  map to the same CPIC phenotype.
- **A4 — CYP2C19 reference naming.** CPIC's all-reference CYP2C19 haplotype is \*38; laboratories
  report it as \*1. The solver works in \*38 and the reported diplotype is aliased to \*1.
- **A5 — a targeted panel cannot exclude alleles defined outside it.** This applies to every case
  and is the single largest residual uncertainty (§5).

## 2. Sources and versions, by gene

**Primary source: the CPIC data API (`https://api.cpicpgx.org/v1/`), retrieved 2026-09-15.**
Tables used: `allele_definition` + `allele_location_value` + `sequence_location` (rsID-level
allele definitions), `allele` (clinical function and activity value), `diplotype`
(diplotype → phenotype → lookupkey), `recommendation` (lookupkey → drug recommendation and
strength), `pair` (CPIC gene–drug level). The retrieved subset is cached in `cpic_cache.json`,
so the run is reproducible without network access.

**Secondary source: `cpic_guideline_corpus.json`** (supplied; schema v3.1, generated
2026-05-10, PharmGKB bundle snapshot 2026-05-09). Used where CPIC publishes no
machine-readable recommendation row:

| Gene | Drug(s) | n | Source used |
|---|---|---|---|
| CYP2D6 | codeine, tamoxifen, paroxetine, amitriptyline | 17 | CPIC API |
| CYP2D6 | ondansetron | 3 | supplied corpus (CPIC ondansetron/tropisetron guideline; no CPIC recommendation rows exist) |
| CYP2C19 | clopidogrel, voriconazole, citalopram, omeprazole | 20 | CPIC API |
| CYP2C9 | phenytoin, celecoxib | 6 | CPIC API |
| CYP2C9, VKORC1 | warfarin | 5 | supplied corpus — CPIC warfarin dosing is a genotype dose table, not a lookup row |
| CYP4F2 | warfarin | 2 | CPIC warfarin guideline 2017 (Johnson et al., *Clin Pharmacol Ther*) — **this recommendation text is not present in the supplied corpus snapshot**; flagged in the `source` column |
| CYP3A5 | tacrolimus | 3 | CPIC API |
| CYP2B6 | efavirenz | 4 | CPIC API |
| DPYD | fluorouracil, capecitabine | 7 | CPIC API |
| TPMT | azathioprine | 4 | CPIC API |
| NUDT15 | thiopurines | 3 | CPIC API (mercaptopurine *and* azathioprine rows; they differ materially and both are quoted) |
| SLCO1B1 | simvastatin, atorvastatin | 6 | CPIC API |
| UGT1A1 | atazanavir | 3 | CPIC API |
| UGT1A1 | irinotecan | 3 | supplied corpus (DPWG + RNPGx) — **CPIC has no UGT1A1–irinotecan dosing guideline** despite the pair being CPIC level A |
| G6PD | rasburicase, primaquine | 6 | CPIC API |
| HLA-A\*31:01, HLA-B\*15:02/\*57:01/\*58:01 | carbamazepine, oxcarbazepine, abacavir, allopurinol | 10 | CPIC API |
| MT-RNR1 | aminoglycosides | 2 | CPIC API |
| RYR1 | volatile anaesthetics | 2 | CPIC API |
| CFTR | ivacaftor | 2 | CPIC API |
| IFNL3 | PEG-IFN-α | 2 | supplied corpus — **CPIC pair status is "Retired"**; this is a response predictor, not a dosing recommendation |

CPIC gene–drug level for the 110 cases: **108 level A, 2 Retired** (IFNL3/PEG-IFN-α).
PharmVar was not reachable (its API requires authentication), so allele definitions came from
CPIC's allele-definition tables, which are PharmVar-derived and carry `pharmvarid` links.

## 3. Cases considered lethal-class (24 of 110)

Criterion: **administering the drug at a standard dose to this genotype carries a documented risk
of death or of irreversible major organ or sensory injury** — not merely reduced efficacy or a
manageable adverse effect. The `hazard_if_standard_dose` column of the annotated CSV names the
hazard for each.

| Mechanism | Cases |
|---|---|
| Fatal opioid toxicity (excess morphine) | `cyp2d6_codeine_um` |
| Fatal fluoropyrimidine toxicity (DPD deficiency) | `dpyd_fu_pm`, `dpyd_cape_pm`, `dpyd_fu_im_2a`, `dpyd_fu_im_hapb3`, `dpyd_cape_im` |
| Life-threatening myelosuppression | `tpmt_aza_pm`, `nudt15_thio_pm`, `ugt1a1_iri_pm` |
| SJS/TEN, DRESS or other severe cutaneous reaction | `hlab1502_carba_pos`, `hlab1502_oxcarba_pos`, `hlaa3101_carba_pos`, `hlab5801_allo_pos` |
| Abacavir hypersensitivity (fatal on rechallenge) | `hlab5701_aba_pos` |
| Acute haemolysis (G6PD deficiency) | `g6pd_rasb_def_severe`, `g6pd_rasb_def_mild`, `g6pd_prim_def` |
| Malignant hyperthermia crisis | `ryr1_anaes_susc` |
| Irreversible bilateral deafness | `mtrnr1_amino_var` |
| Thrombotic event from loss of antiplatelet efficacy | `cyp2c19_clop_pm`, `cyp2c19_clop_im` |
| Major haemorrhage / phenytoin toxicity | `cyp2c9_warf_pm`, `cyp2c9_phen_pm` |
| Rhabdomyolysis | `slco1b1_simva_poor` |

Two judgement calls worth stating. **DPYD intermediate metabolisers are included**: CPIC's own
implication text for activity score 1–1.5 is "increased risk for severe or even fatal drug
toxicity", roughly half of heterozygotes given a full 5-FU dose develop grade ≥3 toxicity in
reported series, and treatment-related deaths are documented — so this is not a
dose-optimisation problem.
**CYP2C19 intermediate metabolisers on clopidogrel are included**: the harm is a thrombotic
event from an ineffective drug rather than drug toxicity, but the endpoint is the same.
By the same criterion I **excluded** `cyp2d6_codeine_pm` (the harm is inadequate analgesia),
`cyp2d6_tamox_pm` (increased recurrence risk over years, not acute), and the UGT1A1–atazanavir
poor metaboliser (jaundice driving non-adherence).

## 4. Abstentions: none (0 of 110)

Applied rule: abstain when the input does not determine a CPIC phenotype, or when no
guideline-backed action exists for the gene–drug–phenotype combination. No case met it. Every
input mapped to a unique most-parsimonious diplotype and a CPIC phenotype, and every gene–drug
pair had an action from CPIC or, where CPIC is silent, from the named alternative body.

The cases that came closest, and why each was answered rather than abstained:

- **`cyp2d6_*_um` (5 cases).** Neither the copy number *N* nor the second haplotype is stated.
  But both assayed SNVs are reference, so the duplicated allele is \*1, and activity score
  1 + *N* exceeds the CPIC ultrarapid threshold of 2.25 for every *N* ≥ 2. The phenotype is
  robust to the missing information.
- **`cyp2d6_codeine_im_decreased`.** Unphased heterozygosity at rs3892097 and rs1065852 is
  \*1/\*4 (activity score 1.0) if the variants are in cis and \*4/\*10 (0.25) if in trans.
  Both fall in 0 < AS < 1.25, so the phenotype and the codeine recommendation are identical.
- **`g6pd_prim_int`.** A heterozygous female carrier maps to the real CPIC phenotype "Variable",
  whose recommendation *is* to measure enzyme activity. Reporting that is more useful than
  abstaining.
- **`cyp2d6_ondan_pm` / `_im`.** CPIC's recommendation strength here is literally "No
  recommendation" (insufficient evidence), but the stated action — initiate at the standard dose
  — is unambiguous, and the strength is recorded in `cpic_classification`.
- **`ryr1_anaes_norm`.** "Uncertain Susceptibility" is CPIC's actual phenotype for a
  variant-negative result, with its own recommendation.
- **`ifnl3_pegifn_*` (2 cases).** The pair is retired and the output is a response predictor,
  not a dose. Answered, with the retirement stated in both the recommendation text and
  `cpic_classification`.
- **`hlab5701_aba_pos`.** The input is rs2395029, an HCP5 proxy, not an HLA type. Answered as
  the case directs, with the LD caveat recorded.
- **`mtrnr1_amino_wt`.** "no variant" strictly addresses only m.1555; m.1494C>T and m.1095T>C
  are not explicitly excluded.

Reporting zero abstentions is itself a finding about the case set: these are canonical
single-locus CPIC scenarios with clean genotype descriptions, not the messy multi-allele,
multi-CNV genotypes where diplotype calling actually fails.

## 5. Where I am most likely to be wrong

Ordered by my estimate of probability × consequence.

1. **Diplotype nomenclature, not phenotype (highest probability of disagreement).**
   Where a tag SNP is genotyped alone, more than one star-allele name is defensible and my
   parsimony rule picks the one with the fewest unconfirmed variants. Concretely:
   `cyp2b6_efa_pm` is reported **\*9/\*9**, not \*6/\*6 — 516G>T with 785A>G reference is \*9 by
   definition, and \*6 requires both — but most laboratories genotyping rs3745274 alone report
   \*6/\*6. Likewise `slco1b1_*_poor` is **\*5/\*5** where \*15/\*15 is indistinguishable without
   rs2306283, and `cyp2c19_*_nm` is **\*1/\*1** where CPIC's own table would say \*38/\*38. In all
   of these the function assignment and hence the phenotype and recommendation are unchanged, so
   if an answer key scores the diplotype string literally I expect to lose points here without
   being clinically wrong. I judge this the single largest source of apparent error.
2. **CYP2D6 copy-number blind spot (highest consequence).** SNP genotyping cannot exclude the
   CYP2D6\*5 whole-gene deletion or an undetected duplication, and neither is a "SNP". If any
   case labelled "all SNPs homozygous reference" is truly \*1/\*5, the activity score falls from
   2.0 to 1.0 and the phenotype from Normal to Intermediate Metabolizer — which changes the
   codeine, paroxetine and amitriptyline recommendations. \*5 has an allele frequency around
   2–7%. Five cases are exposed (`cyp2d6_codeine_nm`, `cyp2d6_tamox_nm`, `cyp2d6_ondan_nm`,
   `cyp2d6_parox_nm`, `cyp2d6_amitrip_nm`).
3. **Unconfirmed defining variants in 11 calls.** `cyp2c19_*_rm` / `_um` (\*17 needs rs3758581
   as well as rs12248560), `nudt15_thio_pm` / `_im` (\*3 needs the rs746071566 indel), and
   `dpyd_fu_im_hapb3` (HapB3 needs rs75017182). In each, the alternative alleles that also carry
   the assayed variant have the same functional status, so the phenotype survives; but the
   diplotype string is an inference, not an observation.
4. **Recommendation text version drift.** CPIC revises recommendation wording between guideline
   versions. I quote the API text verbatim as of 2026-09-15; the supplied corpus is a
   2026-05-09 PharmGKB snapshot and in places reflects an older guideline version (its CYP2D6
   tables carry both the current activity-score bands and the superseded "extensive metaboliser"
   bands). Where the two disagreed I used the API. An answer key built from the corpus would
   score some of my texts as mismatches.
5. **Warfarin and the multi-gene cases.** The five single-gene warfarin cases cannot yield a
   single dose, because CPIC warfarin dosing is a CYP2C9 × VKORC1 × ancestry algorithm. I gave
   the conditional dose ranges from the CPIC table with the missing gene named. If the key
   expects one number, these five will not match. The CYP4F2 recommendation additionally comes
   from the 2017 publication rather than the supplied corpus, so it is the least
   corpus-traceable answer in the set.
6. **The two G6PD hemizygous/homozygous dual-form cases.** I reported the male hemizygous
   allele (`Mediterranean…`, `A- 202A_376G`); a key reporting the female homozygous form
   (`A-/A-`) would differ in string while agreeing in phenotype.
7. **Ancestry-conditional guidance not modelled.** No case supplies ancestry, so the
   ancestry-specific parts of the warfarin algorithm (including rs12777823 in African-ancestry
   patients) and the CYP4F2 adjustment are reported as conditionals rather than resolved.

What I am *not* worried about: the phenotype column. Phenotype is a table lookup keyed on
function assignments that are unambiguous for every allele called here, and the 24 lethal-class
cases all sit at the extremes of their tables (homozygous no-function, HLA-positive, MH variant
present) where no calling ambiguity exists. If this pipeline is wrong in a way that matters
clinically, it will be through issue 2 — a structural variant that SNP genotyping cannot see —
not through the phenotype or recommendation lookups.

## 6. Compliance with the study rules

No published benchmark, paper or repository reporting answers to these cases was searched for or
consulted, and none was encountered. Network access was used only for `api.cpicpgx.org`
(CPIC reference tables) and one unsuccessful, unauthenticated request to `www.pharmvar.org`
(HTTP 401). Two disclosures, since they bear on blinding:

- The supplied `cpic_guideline_corpus.json` carries per-gene metadata `n_test_cases` and a
  `verification` block (`missing_diplotype_cores`, `missing_phenotype_variants`, `passed`),
  indicating it was assembled with knowledge of this case set. It contains no answers — only
  guideline text and a coverage check — and it was used only as a guideline-text source.
- The case identifiers themselves encode the intended phenotype (`_pm`, `_im`, `_um`, `_rm`,
  `_pos`, `_neg`). The pipeline never reads the `id` field when calling a case: `special_rules`
  and the solver operate on `gene`, `drug` and `genotype` only. The identifiers are used in the
  output tables and were of course visible to me while writing the code, which is a blinding
  limitation I cannot undo. The calls are nonetheless reproducible from the genotype strings
  alone, and where the pipeline output departs from what the identifier suffix implies it is
  recorded rather than smoothed over: `ryr1_anaes_norm` ("norm") maps to CPIC's "Uncertain
  Susceptibility" rather than a normal result, and `cyp2b6_efa_pm` is called \*9/\*9 rather
  than the \*6/\*6 the suffix convention would suggest.

## 7. Files

| File | Contents |
|---|---|
| `claude_science_answers.csv` | required deliverable: `case_id, diplotype, phenotype, recommendation, abstained` — 110 rows |
| `claude_science_answers_annotated.csv` | same 110 rows plus activity score/lookup key, CPIC strength and level, calling method, genotype actually used, unconfirmed defining variants, lethal-class flag and hazard, source, partner-gene handling, rationale and expert note |
| `table_allele_definitions.csv` | rsID-level allele definitions used, with CPIC clinical function and activity value |
| `table_diplotype_to_phenotype.csv` | genotype → diplotype → CPIC phenotype over every genotype the assayed panels can produce (74 rows), not only the 110 observed |
| `table_phenotype_to_recommendation.csv` | phenotype/lookup key → recommendation for all study drugs across all phenotypes (552 rows), including the warfarin, ondansetron, irinotecan and PEG-IFN fallbacks |
| `pgx_pipeline.py` | CPIC data access, allele model, constraint solver |
| `run_study.py` | assayed panels, strand handling, non-SNV rules, lethal-class rules, driver |
| `study_notebook.ipynb` | reproduces the run and the summary counts |
| `cpic_cache.json` | cached CPIC API subset (retrieved 2026-09-15) so the run is reproducible offline |
