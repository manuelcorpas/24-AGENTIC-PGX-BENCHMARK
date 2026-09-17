# Methods note — blinded pharmacogenomics study (110 cases, 21 genes)

Run date 2026-09-15. Model: **Claude Opus 5** (model id `claude-opus-5`).

## 1. Approach

The study is a three-stage deterministic pipeline: **reported genotype → diplotype →
CPIC phenotype → CPIC recommendation**. No stage uses free-text reasoning at run time; every
assignment is a table lookup against CPIC's own machine-readable data, and every lookup is
recorded per case in `case_audit.csv`.

**Stage 1 — genotype to diplotype.** For the 59 cases that report a single rsID genotype the
star allele is derived from CPIC's allele-definition tables rather than from prior knowledge:

1. `sequence_location` gives the GRCh38 forward-strand HGVS for each interrogated position
   (e.g. rs3892097 → `g.42128945C>T`), from which the reference and variant base are parsed.
2. The genotype letters as written in the case are matched against that forward-strand
   alphabet. If they do not fit, they are complemented and re-matched, and the orientation
   actually used is recorded in `case_audit.strand_orientation`. None of the 24 rsIDs in this
   case set is a palindromic (A/T or C/G) site, so the strand is always recoverable; the
   pipeline raises rather than guessing if it is not. Five cases required complementation
   (three CYP3A5 rs776746, two VKORC1 rs9923231 — both reported in transcript orientation);
   the remainder, including the cases explicitly labelled "forward-genomic", were already on
   the forward strand.
3. A reference base maps to the gene's reference allele. A variant base maps to the allele
   that CPIC defines **minimally** by that variant — the candidate allele with the fewest
   variant positions relative to the reference haplotype. This yields `*4` for rs3892097 T,
   `*2` for rs4244285 A, `*17` for rs12248560 T, `*3` for rs1057910 C, `*9` for rs3745274 T,
   `*4` for rs2279343 G, `*3B` for rs1800460 T, `*3C` for rs1142345 C, `*5` for rs4149056 C,
   and so on, without any of these being written into the code.

The 51 cases that cannot be keyed on a single SNV go through explicit handlers, each of which
writes its reasoning into `case_audit.basis`: reference-everywhere statements (16), CYP2D6
copy-number duplications (5), two-SNV phasing for CYP2D6 (1) and TPMT (1), the UGT1A1 (TA)n
repeat (6), HLA typing results (8) and the HLA-B\*57:01 tag SNP (2), G6PD haplotypes (6), CFTR
(2), MT-RNR1 (2) and RYR1 (2).

**Stage 2 — diplotype to phenotype.** Lookup in CPIC's `diplotype` table (111,047 rows for the
21 genes), matching CPIC's own diplotype string in either allele order. MT-RNR1 reference
resolves through the `allele` table. CYP4F2, VKORC1 and IFNL3 have no CPIC diplotype-phenotype
table at all — CPIC assigns them no metaboliser or function phenotype — so for those six cases
the phenotype field carries the genotype-level interpretation from the guideline text and says
explicitly that no CPIC phenotype term exists.

**Stage 3 — phenotype to recommendation.** Lookup in CPIC's `recommendation` table on the gene
lookup key. 95 of 110 cases resolve to exactly one CPIC row; none resolved to more than one.
The remaining 15 fall outside CPIC's machine-readable recommendations and are taken from the
supplied corpus with the source named in the answer text and in
`table_phenotype_to_recommendation.csv`:

| cases | why no API row | source used |
|---|---|---|
| 3 ondansetron | guideline exists, recommendation table not populated in the API | supplied corpus, CPIC PA166161954 |
| 7 warfarin (CYP2C9, VKORC1, CYP4F2) | CPIC warfarin guidance is a dose table/algorithm, not a phenotype lookup | supplied corpus, CPIC PA166104949 |
| 3 irinotecan | `guidelineid` is null — CPIC rates UGT1A1–irinotecan level A with testing recommended but has published no guideline | supplied corpus, DPWG PA166104951 + RNPGx PA166127626 |
| 2 PEG-IFN-α | CPIC has **retired** the IFNL3–peginterferon pair (`cpiclevel = "Retired"`) | supplied corpus, CPIC PA166110235 |

Where CPIC keys a recommendation on more than one gene but the case reports only one, the
co-tested gene is set to its reference/normal value and the assumption is recorded on the case:
CYP2C19 Normal Metabolizer for amitriptyline, TPMT or NUDT15 Normal Metabolizer for thiopurines,
HLA-B\*15:02 negative for phenytoin, the other HLA locus negative for carbamazepine, CACNA1S
Uncertain Susceptibility for volatile anaesthetics. Where CPIC gives several populations the
choice is fixed and documented in code: voriconazole `adults`, clopidogrel `CVI ACS PCI`,
phenytoin `PHT naive`, carbamazepine `CBZ naive`, oxcarbazepine `OXC naive`.

## 2. Sources and versions

Per-gene provenance is in `source_versions.csv` (supplied-corpus SHA-256, CPIC `gene` table
version, CPIC lookup method, guidelines used, and the CPIC level of the specific case drugs).

- **CPIC data API**, `https://api.cpicpgx.org/v1/`, accessed 2026-09-15. Tables cached verbatim
  in `cpic_api_raw/`: `gene` (132), `pair` (635), `drug` (324), `guideline` (29),
  `recommendation` (2,115), `diplotype` (111,047), `allele` (1,317), `allele_definition`
  (1,293), `sequence_location` (1,162), `allele_location_value` (2,988). Individual CPIC
  guideline versions used are listed per gene in `source_versions.csv` (e.g. *CYP2D6 and
  Tamoxifen* v71, *DPYD and Fluoropyrimidines* v72, *TPMT, NUDT15 and Thiopurines* v71,
  *SLCO1B1, ABCG2, CYP2C9, and Statins* v72, *G6PD* v67).
- **Supplied corpus** `cpic_guideline_corpus.json`, schema v3.1, generated 2026-05-10T04:40:15Z
  from a PharmGKB bundle retrieved 2026-05-09, 21 genes. Used for the four recommendation
  classes above and to cross-check phenotype terminology.
- PharmVar and PharmGKB were not queried separately; CPIC's allele-definition tables are
  PharmVar-derived and were sufficient.
- Every case drug–gene pair is **CPIC level A** except IFNL3–peginterferon (Retired). Confirmed
  from the `pair` table, not assumed.

No published benchmark, paper or repository reporting answers to these cases was searched for or
consulted, and none was encountered.

## 3. Lethal-class cases

22 of 110 cases are flagged `lethal_class = true` in `case_audit.csv`: cases where
genotype-blind standard dosing carries a documented risk of death or irreversible injury, as
opposed to therapeutic failure or reversible dose-related toxicity.

| case_id | gene | drug | phenotype | why |
|---|---|---|---|---|
| 74597c5503ec | CYP2D6 | codeine | Ultrarapid Metabolizer | ultrarapid morphine formation; fatal respiratory depression, the deaths that prompted the label change |
| b94a43e574db | DPYD | fluorouracil | Poor Metabolizer | complete DPD deficiency; full dose is lethal |
| c6549d637d49 | DPYD | capecitabine | Poor Metabolizer | as above |
| c57caa44a38a, 6b8d2db91001 | DPYD | fluorouracil | Intermediate Metabolizer | partial deficiency; grade 4–5 toxicity at full dose |
| ae4a7215ae64 | DPYD | capecitabine | Intermediate Metabolizer | as above |
| 08e18e809fcb | TPMT | azathioprine | Poor Metabolizer | life-threatening myelosuppression at conventional dose |
| f3dd19d8a337 | NUDT15 | thiopurines | Poor Metabolizer | as above |
| d49e9f88d9c2, 3786ab9a929d | G6PD | rasburicase | Deficient | acute haemolysis and methaemoglobinaemia; rasburicase is contraindicated |
| 2a8544a44d5f | G6PD | primaquine | Deficient | acute haemolytic anaemia |
| 7b560d699828 | G6PD | primaquine | Variable | heterozygous female; skewed X-inactivation means genotype cannot exclude severe haemolysis |
| e70379d3dbac | HLA-B\*57:01 | abacavir | positive | abacavir hypersensitivity syndrome; rechallenge is fatal |
| 7181a24a65d5 | HLA-B\*15:02 | carbamazepine | positive | SJS/TEN |
| 569a842e1a9a | HLA-B\*15:02 | oxcarbazepine | positive | SJS/TEN |
| 85976e86289a | HLA-A\*31:01 | carbamazepine | positive | DRESS and SJS/TEN |
| ebf82c3bcfdd | HLA-B\*58:01 | allopurinol | positive | SJS/TEN and DRESS; allopurinol is contraindicated |
| 369c2bcff713 | MT-RNR1 | aminoglycosides | increased risk | permanent bilateral deafness at any dose |
| 1229a4123378 | RYR1 | volatile-anaesthetics | MH Susceptibility | fulminant malignant hyperthermia |
| ca9e14feca53 | CYP2C19 | clopidogrel | Poor Metabolizer | failure of clopidogrel activation; stent thrombosis and MI |
| 97c3baebed4e | CYP2C9 | warfarin | Poor Metabolizer | severe over-anticoagulation and major bleeding |
| e85ecfa887cb | CYP2C9 | phenytoin | Poor Metabolizer | toxic concentrations and encephalopathy |

Deliberately **not** flagged, with reasons: TPMT/NUDT15 intermediate metabolisers (dose-limiting
but usually reversible cytopenia, unlike DPYD where grade 5 events are documented at full dose);
CYP2C9 poor metabolisers on celecoxib (serious GI/renal toxicity, not typically immediately
life-threatening); CYP2D6 poor metabolisers on codeine (analgesic failure, not toxicity);
SLCO1B1 poor function on simvastatin (rhabdomyolysis is the recognised risk but at the
guideline's dose ceilings it is uncommon and reversible); RYR1 Uncertain Susceptibility
(879d83c0792d) — CPIC is explicit that a negative RYR1 result does not exclude MH, and roughly
half of MH-susceptible individuals carry no identifiable RYR1 variant, but flagging every
genotype-negative patient would empty the category of meaning. This last case is the one where
a "normal" answer is most likely to be clinically misread, which is why the phenotype is
reported as CPIC's *Uncertain Susceptibility* and not as normal.

## 4. Abstentions

**No case was abstained** (`abstained = false` for all 110). All 110 resolved to a CPIC
diplotype and phenotype, and to a recommendation from either CPIC's recommendation table or the
supplied guideline text. Rather than abstaining, cases with residual uncertainty carry the
uncertainty explicitly: 86 cases have a recorded assumption and 18 have a recorded ambiguity in
`case_audit.csv`.

The five CYP2D6 copy-number cases (74597c5503ec, 61741d10f238, 01166e2ea9ad, 9132fcb83a8f,
ae968d6ae64c) came closest to the abstention threshold, because a copy-number assay cannot phase
the duplication to an allele and the cases give no SNV result. I chose to answer them as
`*1/*1x2` → Ultrarapid Metabolizer on three grounds: the case names the duplicated allele as
`*1`; no variant allele is reported, so the non-duplicated chromosome is `*1` by the same
convention used for the "all SNPs homozygous reference" cases; and for codeine specifically the
error is asymmetric — calling Ultrarapid when the truth were Normal withholds a drug, whereas
abstaining when the truth is Ultrarapid withholds the single most consequential warning in the
set. Both assumptions and the alternative interpretation (`*4/*1x2` → activity score 2.0 →
Normal Metabolizer) are recorded on every one of those cases.

## 5. Where I am most likely wrong

Ordered by my own estimate of probability × consequence.

1. **CYP2D6 \*1xN, 5 cases.** Highest-risk call in the study. If the intended answer is an
   abstention, or if the duplicated allele is meant to be unresolved, all five are wrong. The
   phenotype is robust to N but not to the identity of the second chromosome.
2. **Single-gene answers for a multi-gene drug, 7 warfarin cases.** CPIC warfarin dosing needs
   CYP2C9 **and** VKORC1 **and** continental ancestry. Only CYP2C9 \*3/\*3 is determinate on one
   gene (0.5–2 mg/day for every VKORC1 genotype). The other six answers are conditional ranges;
   if the expected answer is a single dose or an abstention, they will not match. The same
   single-gene assumption affects amitriptyline (CYP2C19 assumed normal), phenytoin (HLA-B
   assumed negative) and carbamazepine (second HLA locus assumed negative).
3. **Recommendation wording, all 110.** The answers file carries CPIC's recommendation text
   condensed to ≤200 characters, not a canonical short action string. A grader matching on
   phrases such as "avoid" or "use alternative" will score these differently from a grader
   matching exact text. The verbatim text is in `table_phenotype_to_recommendation.csv`.
4. **Ondansetron phenotype bins, 3 cases.** The published ondansetron table predates CPIC's 2019
   CYP2D6 translation (it calls activity score 1.0 a Normal Metabolizer and 0.5 an Intermediate
   Metabolizer). I used the current translation from the API. No ondansetron case here sits at
   activity score 1.0, so the calls are unaffected, but a grader working from the published table
   would disagree on bin labels.
5. **Non-CPIC recommendation sources, 5 cases.** Irinotecan (3, no CPIC guideline exists) and
   PEG-IFN-α (2, pair retired). If the expected answer is an abstention because there is no CPIC
   Level A recommendation, these are wrong; I judged reporting the DPWG/RNPGx and retired-CPIC
   guidance with the source named to be more useful than a blank.
6. **Allele naming conventions.** CYP2C19's zero-variant haplotype is `*38` in CPIC's definition
   table, not `*1` (`*1` carries 991A>G, which was not interrogated); I report the conventional
   `*1/*2` and `*1/*17`, which is what a clinical laboratory reports but not the strictest reading
   of the data. Similarly SLCO1B1 rs4149056 C/C is reported as `*5/*5` though `*15/*15` is
   indistinguishable without rs2306283 (both Poor Function); CYP2B6 rs3745274 T/T is reported as
   `*9/*9` though `*6/*6` is indistinguishable without rs2279343 (both Poor Metabolizer); TPMT
   double-homozygote is `*3A/*3A` though `*3B/*3C` shares the phenotype; G6PD hemizygotes are
   reported with CPIC's single-allele string.
7. **HLA-B\*57:01 by tag SNP, 2 cases.** rs2395029 is an HCP5 proxy, not an HLA type. Linkage is
   strong in European ancestry and weaker elsewhere; a positive tag with a negative true type
   would wrongly withhold abacavir.
8. **Panel-negative cases, 20 "reference" cases.** "All SNPs homozygous reference" is only as
   good as the panel. A rare no-function allele off-panel would make a reported Normal
   Metabolizer wrong in the dangerous direction. This is a property of the input, not of the
   mapping, but it is the failure mode I would most want flagged in a real clinical report.

## 6. Files

| file | content |
|---|---|
| `claude_science_answers.csv` | the 110 answers: `case_id,diplotype,phenotype,recommendation,abstained` |
| `table_diplotype_to_phenotype.csv` | 60 diplotype → CPIC phenotype rows actually used, with lookup key and source |
| `table_phenotype_to_recommendation.csv` | 106 phenotype → recommendation rows actually used, verbatim text, CPIC classification, source, lethal flag |
| `case_audit.csv` | per-case provenance: mapping method, strand orientation, forward-genomic HGVS, basis, assumptions, ambiguities, sources, lethal flag |
| `source_versions.csv` | per-gene source and version record |
| `fetch_cpic_api.py` | downloads and caches the CPIC API tables |
| `pgx_pipeline.py` | the pipeline that produces every CSV above |
| `pgx_study.ipynb` | notebook walking through the analysis and summarising results |
| `cpic_api_raw/` | verbatim CPIC API responses used |

Reproduce with `python fetch_cpic_api.py && python pgx_pipeline.py` (Python 3.11, standard
library only).

*This is research and informational content, not clinical advice. Prescribing decisions must be
made by a qualified clinician with access to the full patient context and a validated
laboratory report.*
