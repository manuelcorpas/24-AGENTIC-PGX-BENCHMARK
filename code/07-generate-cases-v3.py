#!/usr/bin/env python3
"""
Generate ~120 CPIC-grounded pharmacogenomics test cases for the BiB benchmark v3.

Source: CPIC level-A guidelines (https://cpicpgx.org/guidelines/) — the gene,
diplotype-tier, drug, and recommendation entries below are encoded from the
published CPIC guideline tables. Population frequency notes for AMR/AFR are
either retained from v1 (curated by H. Guio and S. Fatumo) for the original
8 genes, or sourced from CPIC frequency tables and flagged with [CHECK]
where curated review is still pending.

Outputs:
  ../SPECS/test_cases_v3.json      — full case list, runner-ready
  ../SPECS/spot_check_sample.md    — 10 representative cases for manual review
  ../SPECS/generation_log.txt      — case count by gene/tier
"""
from __future__ import annotations
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SPECS = BASE / "SPECS"
SPECS.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# CPIC level-A gene-drug guideline table.
# Each entry encodes: gene, drug(s), per-tier (diplotype, defining variant
# rsID + genotype, activity-score-or-status, recommendation), and per-population
# frequency notes.
# All recommendations are paraphrased from the relevant CPIC guideline.
# Lethal-class cases are flagged with "(potentially LETHAL)" in gt_drug so the
# scorer counts A3 errors.
# ---------------------------------------------------------------------------

GUIDELINES = [
    # ============================== CYP2D6 ==============================
    {
        "gene": "CYP2D6",
        "drug": "codeine",
        "tiers": [
            {"id": "cyp2d6_codeine_pm", "diplotype": "*4/*4", "rsid": "rs3892097 T/T (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Codeine: AVOID (no analgesic benefit; use non-opioid or non-CYP2D6 opioid)"},
            {"id": "cyp2d6_codeine_im", "diplotype": "*1/*4", "rsid": "rs3892097 C/T (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Codeine: standard starting dose; monitor for reduced analgesia"},
            {"id": "cyp2d6_codeine_nm", "diplotype": "*1/*1", "rsid": "All CYP2D6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Codeine: standard dosing"},
            {"id": "cyp2d6_codeine_um", "diplotype": "*1/*1xN (gene duplication)",
             "rsid": "CYP2D6 *1xN duplication detected by copy-number assay (qPCR/MLPA/long-range PCR); duplications are structural variants with no defining single SNV",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "≥3.0",
             "recommendation": "Codeine: AVOID (potentially LETHAL morphine accumulation; respiratory depression risk)"},
            {"id": "cyp2d6_codeine_im_decreased", "diplotype": "*4/*10",
             "rsid": "rs3892097 C/T + rs1065852 G/A (compound)",
             "phenotype": "Intermediate Metaboliser", "activity": "0.5",
             "recommendation": "Codeine: reduce dose or alternative; *10 reduced function compounding *4 null"},
        ],
        "pop_note": {
            "EUR": "CYP2D6*4 ~20% in Europeans; PM ~7%. *1xN (UM) ~1-3%.",
            "AMR": "CYP2D6*4 ~10% in admixed Latin Americans; *10 and *17 more common in indigenous populations.",
            "AFR": "CYP2D6*4 ~6% in Africans; *17 (reduced function) ~20-35%, often missed by EUR-centric panels.",
        },
    },
    {
        "gene": "CYP2D6",
        "drug": "tamoxifen",
        "tiers": [
            {"id": "cyp2d6_tamox_pm", "diplotype": "*4/*4", "rsid": "rs3892097 T/T (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Tamoxifen: alternative endocrine therapy (e.g. aromatase inhibitor)"},
            {"id": "cyp2d6_tamox_im", "diplotype": "*1/*4", "rsid": "rs3892097 C/T (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Tamoxifen: consider higher dose or alternative; reduced endoxifen formation"},
            {"id": "cyp2d6_tamox_nm", "diplotype": "*1/*1", "rsid": "All CYP2D6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Tamoxifen: standard dosing"},
            {"id": "cyp2d6_tamox_um", "diplotype": "*1/*1xN", "rsid": "CYP2D6 *1xN duplication detected by copy-number assay (qPCR/MLPA/long-range PCR); structural variant, no defining single SNV",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "≥3.0",
             "recommendation": "Tamoxifen: standard dosing (UM is acceptable for tamoxifen efficacy)"},
        ],
        "pop_note": {
            "EUR": "CYP2D6*4 ~20% in Europeans; PM ~7% — relevant for breast cancer endocrine therapy.",
            "AMR": "CYP2D6*4 ~10% in admixed Latin Americans; tamoxifen response variable.",
            "AFR": "CYP2D6*17 ~20-35% in Africans causes intermediate metaboliser status not detected by *4-only panels.",
        },
    },
    {
        "gene": "CYP2D6",
        "drug": "ondansetron",
        "tiers": [
            {"id": "cyp2d6_ondan_pm", "diplotype": "*4/*4", "rsid": "rs3892097 T/T",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Ondansetron: standard dosing (PM acceptable for antiemetic efficacy)"},
            {"id": "cyp2d6_ondan_nm", "diplotype": "*1/*1", "rsid": "All CYP2D6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Ondansetron: standard dosing"},
            {"id": "cyp2d6_ondan_um", "diplotype": "*1/*1xN", "rsid": "CYP2D6 *1xN duplication detected by copy-number assay (qPCR/MLPA/long-range PCR); structural variant, no defining single SNV",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "≥3.0",
             "recommendation": "Ondansetron: alternative antiemetic (e.g. granisetron); reduced efficacy in UM"},
        ],
        "pop_note": {
            "EUR": "CYP2D6*1xN (UM) ~1-3% in Europeans.",
            "AMR": "CYP2D6 UM frequency variable in admixed Latin American populations.",
            "AFR": "CYP2D6 UM (*1xN, *2xN) up to 28% in some North/East African populations.",
        },
    },
    {
        "gene": "CYP2D6",
        "drug": "paroxetine",
        "tiers": [
            {"id": "cyp2d6_parox_pm", "diplotype": "*4/*4", "rsid": "rs3892097 T/T",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Paroxetine: alternative SSRI not metabolised by CYP2D6"},
            {"id": "cyp2d6_parox_im", "diplotype": "*1/*4", "rsid": "rs3892097 C/T",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Paroxetine: standard starting dose; monitor for adverse effects"},
            {"id": "cyp2d6_parox_nm", "diplotype": "*1/*1", "rsid": "All CYP2D6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Paroxetine: standard dosing"},
            {"id": "cyp2d6_parox_um", "diplotype": "*1/*1xN", "rsid": "CYP2D6 *1xN duplication detected by copy-number assay (qPCR/MLPA/long-range PCR); structural variant, no defining single SNV",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "≥3.0",
             "recommendation": "Paroxetine: alternative SSRI; reduced efficacy due to rapid clearance"},
        ],
        "pop_note": {
            "EUR": "CYP2D6 PM ~7% in Europeans relevant for SSRI choice.",
            "AMR": "CYP2D6 phenotype distribution variable in admixed populations; testing recommended.",
            "AFR": "CYP2D6*17 + UM duplications affect SSRI dosing; population data sparse.",
        },
    },
    {
        "gene": "CYP2D6",
        "drug": "amitriptyline",
        "tiers": [
            {"id": "cyp2d6_amitrip_pm", "diplotype": "*4/*4", "rsid": "rs3892097 T/T",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Amitriptyline: avoid use; consider alternative not metabolised by CYP2D6"},
            {"id": "cyp2d6_amitrip_im", "diplotype": "*1/*4", "rsid": "rs3892097 C/T",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Amitriptyline: 25% dose reduction; monitor"},
            {"id": "cyp2d6_amitrip_nm", "diplotype": "*1/*1", "rsid": "All CYP2D6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Amitriptyline: standard dosing"},
            {"id": "cyp2d6_amitrip_um", "diplotype": "*1/*1xN", "rsid": "CYP2D6 *1xN duplication detected by copy-number assay (qPCR/MLPA/long-range PCR); structural variant, no defining single SNV",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "≥3.0",
             "recommendation": "Amitriptyline: avoid use; alternative TCA not affected by CYP2D6"},
        ],
        "pop_note": {
            "EUR": "TCA dosing affected by CYP2D6*4 (PM) and *1xN (UM) — common testing in psychiatry.",
            "AMR": "CYP2D6 phenotype variable in admixed populations.",
            "AFR": "CYP2D6*17 prevalent — affects TCA dosing through reduced 2D6 activity.",
        },
    },

    # ============================== CYP2C19 ==============================
    {
        "gene": "CYP2C19",
        "drug": "clopidogrel",
        "tiers": [
            {"id": "cyp2c19_clop_pm", "diplotype": "*2/*2", "rsid": "rs4244285 A/A (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Clopidogrel: AVOID (alternative antiplatelet — prasugrel or ticagrelor; reduced active metabolite, increased thrombotic risk)"},
            {"id": "cyp2c19_clop_im", "diplotype": "*1/*2", "rsid": "rs4244285 G/A (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Clopidogrel: alternative antiplatelet recommended; reduced active metabolite"},
            {"id": "cyp2c19_clop_nm", "diplotype": "*1/*1", "rsid": "All CYP2C19 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Clopidogrel: standard dosing"},
            {"id": "cyp2c19_clop_rm", "diplotype": "*1/*17", "rsid": "rs12248560 C/T (heterozygous)",
             "phenotype": "Rapid Metaboliser", "activity": "1.5",
             "recommendation": "Clopidogrel: standard dosing"},
            {"id": "cyp2c19_clop_um", "diplotype": "*17/*17", "rsid": "rs12248560 T/T (homozygous)",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "1.5+",
             "recommendation": "Clopidogrel: standard dosing (UM may have increased bleeding risk; monitor)"},
        ],
        "pop_note": {
            "EUR": "CYP2C19*2 ~15% in EUR; PM ~2%. *17 ~21%.",
            "AMR": "CYP2C19*2 ~10-15% in admixed Latin Americans.",
            "AFR": "CYP2C19*2 ~15-18% in African populations; PM ~3-5%. *17 ~16-26%.",
        },
    },
    {
        "gene": "CYP2C19",
        "drug": "voriconazole",
        "tiers": [
            {"id": "cyp2c19_vori_pm", "diplotype": "*2/*2", "rsid": "rs4244285 A/A",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Voriconazole: alternative antifungal (toxic accumulation expected)"},
            {"id": "cyp2c19_vori_im", "diplotype": "*1/*2", "rsid": "rs4244285 G/A",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Voriconazole: standard dosing with therapeutic drug monitoring"},
            {"id": "cyp2c19_vori_nm", "diplotype": "*1/*1", "rsid": "All CYP2C19 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Voriconazole: standard dosing"},
            {"id": "cyp2c19_vori_rm", "diplotype": "*1/*17", "rsid": "rs12248560 C/T",
             "phenotype": "Rapid Metaboliser", "activity": "1.5",
             "recommendation": "Voriconazole: alternative antifungal (subtherapeutic levels likely)"},
            {"id": "cyp2c19_vori_um", "diplotype": "*17/*17", "rsid": "rs12248560 T/T",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "1.5+",
             "recommendation": "Voriconazole: alternative antifungal (subtherapeutic levels — treatment failure risk)"},
        ],
        "pop_note": {
            "EUR": "CYP2C19*17 ~21% — RM/UM common, voriconazole TDM recommended.",
            "AMR": "CYP2C19*17 ~15% in admixed Latin Americans.",
            "AFR": "CYP2C19*17 ~16-26% in African populations affecting voriconazole levels.",
        },
    },
    {
        "gene": "CYP2C19",
        "drug": "citalopram",
        "tiers": [
            {"id": "cyp2c19_cital_pm", "diplotype": "*2/*2", "rsid": "rs4244285 A/A",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Citalopram: 50% dose reduction or alternative SSRI"},
            {"id": "cyp2c19_cital_im", "diplotype": "*1/*2", "rsid": "rs4244285 G/A",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Citalopram: standard dosing; monitor for adverse effects"},
            {"id": "cyp2c19_cital_nm", "diplotype": "*1/*1", "rsid": "All CYP2C19 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Citalopram: standard dosing"},
            {"id": "cyp2c19_cital_rm", "diplotype": "*1/*17", "rsid": "rs12248560 C/T",
             "phenotype": "Rapid Metaboliser", "activity": "1.5",
             "recommendation": "Citalopram: alternative SSRI (subtherapeutic likely)"},
            {"id": "cyp2c19_cital_um", "diplotype": "*17/*17", "rsid": "rs12248560 T/T",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "1.5+",
             "recommendation": "Citalopram: alternative SSRI (treatment failure risk)"},
        ],
        "pop_note": {
            "EUR": "CYP2C19*2/*17 distribution well characterised in psychiatric SSRI dosing.",
            "AMR": "CYP2C19 phenotypes variable in admixed Latin Americans.",
            "AFR": "CYP2C19*2 ~15-18%, *17 ~16-26% — population-specific dosing relevant.",
        },
    },
    {
        "gene": "CYP2C19",
        "drug": "omeprazole",
        "tiers": [
            {"id": "cyp2c19_omep_pm", "diplotype": "*2/*2", "rsid": "rs4244285 A/A",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Omeprazole: standard dose may be excessive; consider 25-50% reduction"},
            {"id": "cyp2c19_omep_im", "diplotype": "*1/*2", "rsid": "rs4244285 G/A",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Omeprazole: standard dosing"},
            {"id": "cyp2c19_omep_nm", "diplotype": "*1/*1", "rsid": "All CYP2C19 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Omeprazole: standard dosing"},
            {"id": "cyp2c19_omep_rm", "diplotype": "*1/*17", "rsid": "rs12248560 C/T",
             "phenotype": "Rapid Metaboliser", "activity": "1.5",
             "recommendation": "Omeprazole: increase dose 50-100% or alternative PPI; reduced efficacy"},
            {"id": "cyp2c19_omep_um", "diplotype": "*17/*17", "rsid": "rs12248560 T/T",
             "phenotype": "Ultra-rapid Metaboliser", "activity": "1.5+",
             "recommendation": "Omeprazole: increase dose or alternative PPI (treatment failure for H. pylori, reflux)"},
        ],
        "pop_note": {
            "EUR": "CYP2C19 phenotyping common for H. pylori eradication regimens.",
            "AMR": "CYP2C19*17 RM/UM ~15%+ — affects PPI choice.",
            "AFR": "CYP2C19*17 ~16-26% — treatment-failure risk for H. pylori with omeprazole.",
        },
    },

    # ============================== CYP2C9 ==============================
    {
        "gene": "CYP2C9",
        "drug": "warfarin",
        "tiers": [
            {"id": "cyp2c9_warf_pm", "diplotype": "*3/*3", "rsid": "rs1057910 C/C (homozygous alt)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm (e.g., IWPC or Gage) incorporating CYP2C9, VKORC1, and (in African ancestry) CYP2C9*5/*6/*8/*11; expected maintenance dose substantially lower than standard; consider alternative anticoagulant if algorithm unavailable"},
            {"id": "cyp2c9_warf_im", "diplotype": "*1/*3", "rsid": "rs1057910 A/C (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm (e.g., IWPC or Gage) incorporating CYP2C9, VKORC1, and (in African ancestry) CYP2C9*5/*6/*8/*11"},
            {"id": "cyp2c9_warf_nm", "diplotype": "*1/*1", "rsid": "All CYP2C9 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm (e.g., IWPC or Gage); standard starting dose acceptable if algorithm unavailable"},
        ],
        "pop_note": {
            "EUR": "CYP2C9*2 ~13%, *3 ~7% in EUR; combined with VKORC1 explains ~30-40% of dose variance.",
            "AMR": "CYP2C9*2 ~5-10%, *3 ~3-5% in admixed Latin Americans.",
            "AFR": "CYP2C9*3 rare (<1%) in African populations; *5, *6, *8, *11 more relevant — often missed.",
        },
    },
    {
        "gene": "CYP2C9",
        "drug": "phenytoin",
        "tiers": [
            {"id": "cyp2c9_phen_pm", "diplotype": "*3/*3", "rsid": "rs1057910 C/C",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Phenytoin: 50% dose reduction; monitor levels closely"},
            {"id": "cyp2c9_phen_im", "diplotype": "*1/*3", "rsid": "rs1057910 A/C",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Phenytoin: 25% dose reduction with TDM"},
            {"id": "cyp2c9_phen_nm", "diplotype": "*1/*1", "rsid": "All CYP2C9 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Phenytoin: standard dosing"},
        ],
        "pop_note": {
            "EUR": "CYP2C9*3 ~7% — phenytoin toxicity risk in PM/IM.",
            "AMR": "CYP2C9*3 ~3-5% in admixed Latin Americans.",
            "AFR": "CYP2C9*3 rare in Africans; *8, *11 more frequent and clinically actionable.",
        },
    },
    {
        "gene": "CYP2C9",
        "drug": "celecoxib",
        "tiers": [
            {"id": "cyp2c9_celec_pm", "diplotype": "*3/*3", "rsid": "rs1057910 C/C",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Celecoxib: 50% dose reduction or alternative NSAID"},
            {"id": "cyp2c9_celec_im", "diplotype": "*1/*3", "rsid": "rs1057910 A/C",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Celecoxib: standard dosing with caution"},
            {"id": "cyp2c9_celec_nm", "diplotype": "*1/*1", "rsid": "All CYP2C9 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Celecoxib: standard dosing"},
        ],
        "pop_note": {
            "EUR": "CYP2C9*3 PM ~0.4% — relevant for chronic NSAID therapy.",
            "AMR": "CYP2C9*3 ~3-5% in Latin Americans — IM common.",
            "AFR": "CYP2C9 African-specific alleles (*5, *8, *11) affect NSAID metabolism.",
        },
    },

    # ============================== CYP3A5 ==============================
    {
        "gene": "CYP3A5",
        "drug": "tacrolimus",
        "tiers": [
            {"id": "cyp3a5_tacr_nm_expressor", "diplotype": "*1/*1", "rsid": "rs776746 A/A (homozygous *1)",
             "phenotype": "Normal Metaboliser (Expressor)", "activity": "2.0",
             "recommendation": "Tacrolimus: increase dose 1.5-2x; rapid clearance"},
            {"id": "cyp3a5_tacr_im", "diplotype": "*1/*3", "rsid": "rs776746 A/G (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Tacrolimus: increase dose 1.5x; partial expression"},
            {"id": "cyp3a5_tacr_pm_nonexpressor", "diplotype": "*3/*3", "rsid": "rs776746 G/G (homozygous)",
             "phenotype": "Poor Metaboliser (Non-expressor)", "activity": "0.0",
             "recommendation": "Tacrolimus: standard dosing (slow clearance, normal target levels)"},
        ],
        "pop_note": {
            "EUR": "CYP3A5*3 ~85-95% in EUR — most are non-expressors with standard tacrolimus dosing.",
            "AMR": "CYP3A5*3 ~70-80% in admixed Latin Americans.",
            "AFR": "CYP3A5*1 (expressor) common in Africans (~50-70%) — higher tacrolimus dose required.",
        },
    },

    # ============================== CYP2B6 ==============================
    {
        "gene": "CYP2B6",
        "drug": "efavirenz",
        "tiers": [
            {"id": "cyp2b6_efa_pm", "diplotype": "*6/*6", "rsid": "rs3745274 T/T (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Efavirenz: 50% dose reduction or alternative ART; CNS toxicity risk"},
            {"id": "cyp2b6_efa_im", "diplotype": "*1/*6", "rsid": "rs3745274 G/T (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Efavirenz: standard dosing with monitoring; CNS adverse effect risk"},
            {"id": "cyp2b6_efa_nm", "diplotype": "*1/*1", "rsid": "All CYP2B6 SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Efavirenz: standard dosing"},
            {"id": "cyp2b6_efa_rm", "diplotype": "*1/*4", "rsid": "rs2279343 A/G (heterozygous)",
             "phenotype": "Rapid Metaboliser", "activity": "1.5",
             "recommendation": "Efavirenz: standard dosing; possible subtherapeutic levels"},
        ],
        "pop_note": {
            "EUR": "CYP2B6*6 ~25% in EUR. Efavirenz largely replaced by integrase inhibitors.",
            "AMR": "CYP2B6*6 ~30-40% in admixed Latin Americans. Efavirenz still used in some HIV programmes.",
            "AFR": "CYP2B6*6 ~40-50% in East Africa. Efavirenz remains first-line in many African HIV programmes — HIGH IMPACT.",
        },
    },

    # ============================== DPYD ==============================
    {
        "gene": "DPYD",
        "drug": "fluorouracil",
        "tiers": [
            {"id": "dpyd_fu_pm", "diplotype": "*2A/*2A", "rsid": "rs3918290 T/T (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Fluorouracil: AVOID (potentially LETHAL — complete DPD deficiency causes severe toxicity)"},
            {"id": "dpyd_fu_im_2a", "diplotype": "*1/*2A", "rsid": "rs3918290 C/T (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Fluorouracil: reduce starting dose by 50%"},
            {"id": "dpyd_fu_im_hapb3", "diplotype": "*1/HapB3", "rsid": "rs56038477 C/T (HapB3 marker)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Fluorouracil: reduce starting dose by 50%; HapB3 reduced function"},
            {"id": "dpyd_fu_nm", "diplotype": "*1/*1", "rsid": "All DPYD risk variants reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Fluorouracil: standard dosing"},
        ],
        "pop_note": {
            "EUR": "DPYD*2A ~1% in EUR. Pre-treatment testing mandated by EMA (2020) in many EU countries.",
            "AMR": "DPYD*2A poorly characterised in Latin American populations. Testing infrastructure limited.",
            "AFR": "DPYD variant frequencies poorly studied in African populations — risk of undetected DPD deficiency.",
        },
    },
    {
        "gene": "DPYD",
        "drug": "capecitabine",
        "tiers": [
            {"id": "dpyd_cape_pm", "diplotype": "*2A/*2A", "rsid": "rs3918290 T/T",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Capecitabine: AVOID (potentially LETHAL — same DPD deficiency mechanism as 5-FU)"},
            {"id": "dpyd_cape_im", "diplotype": "*1/*2A", "rsid": "rs3918290 C/T",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Capecitabine: reduce starting dose by 50%"},
            {"id": "dpyd_cape_nm", "diplotype": "*1/*1", "rsid": "All DPYD risk variants reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Capecitabine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "DPYD pre-treatment testing recommended for capecitabine in EU.",
            "AMR": "[CHECK] DPYD frequencies in admixed Latin Americans — testing infrastructure varies.",
            "AFR": "[CHECK] DPYD genotyping rare in African oncology — risk of unmonitored toxicity.",
        },
    },

    # ============================== TPMT ==============================
    {
        "gene": "TPMT",
        "drug": "azathioprine",
        "tiers": [
            {"id": "tpmt_aza_pm", "diplotype": "*3A/*3A",
             "rsid": "rs1800460 T/T + rs1142345 C/C (forward-genomic; transcript: c.460G>A + c.719A>G; TPMT on minus strand)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Azathioprine: AVOID (potentially LETHAL myelosuppression) or reduce dose by 90%"},
            {"id": "tpmt_aza_im_3b", "diplotype": "*1/*3B",
             "rsid": "rs1800460 C/T (heterozygous, forward-genomic; transcript: c.460G>A; TPMT on minus strand)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Azathioprine: reduce dose by 30-70%; monitor blood counts"},
            {"id": "tpmt_aza_im_3c", "diplotype": "*1/*3C", "rsid": "rs1142345 T/C (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Azathioprine: reduce dose by 30-70%; monitor blood counts"},
            {"id": "tpmt_aza_nm", "diplotype": "*1/*1", "rsid": "All TPMT SNPs homozygous reference",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Azathioprine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "TPMT*3A ~5% in EUR; PM ~0.3%. Pre-treatment testing routine in IBD/transplant.",
            "AMR": "TPMT variant frequencies similar to EUR in admixed populations.",
            "AFR": "TPMT*3C more common than *3A in African populations (~5-8%) — *3A-only panels miss this.",
        },
    },

    # ============================== NUDT15 ==============================
    {
        "gene": "NUDT15",
        "drug": "thiopurines",
        "tiers": [
            {"id": "nudt15_thio_pm", "diplotype": "*3/*3", "rsid": "rs116855232 T/T (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Thiopurines: AVOID (potentially LETHAL myelosuppression) or reduce dose by 90%"},
            {"id": "nudt15_thio_im", "diplotype": "*1/*3", "rsid": "rs116855232 C/T (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Thiopurines: reduce dose by 30-50%"},
            {"id": "nudt15_thio_nm", "diplotype": "*1/*1", "rsid": "rs116855232 C/C (reference)",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Thiopurines: standard dosing"},
        ],
        "pop_note": {
            "EUR": "NUDT15*3 rare in EUR (~0.2%) — TPMT testing prioritised.",
            "AMR": "NUDT15*3 ~3-5% in admixed Latin Americans (Indigenous ancestry).",
            "AFR": "NUDT15*3 rare in Africans — TPMT remains primary marker.",
        },
    },

    # ============================== SLCO1B1 ==============================
    {
        "gene": "SLCO1B1",
        "drug": "simvastatin",
        "tiers": [
            {"id": "slco1b1_simva_poor", "diplotype": "*5/*5", "rsid": "rs4149056 C/C (homozygous)",
             "phenotype": "Poor Function", "activity": "0.0",
             "recommendation": "Simvastatin: avoid high dose (≥40mg); myopathy risk; alternative statin"},
            {"id": "slco1b1_simva_decr", "diplotype": "*1a/*5", "rsid": "rs4149056 T/C (heterozygous)",
             "phenotype": "Decreased Function", "activity": "1.0",
             "recommendation": "Simvastatin: limit to 20mg/day or alternative statin"},
            {"id": "slco1b1_simva_normal", "diplotype": "*1a/*1a", "rsid": "rs4149056 T/T (reference)",
             "phenotype": "Normal Function", "activity": "2.0",
             "recommendation": "Simvastatin: standard dosing"},
        ],
        "pop_note": {
            "EUR": "SLCO1B1*5 ~15% in EUR; homozygous ~2%.",
            "AMR": "SLCO1B1*5 ~10-15% in admixed Latin Americans.",
            "AFR": "SLCO1B1*5 ~2-3% in Africans — lower statin myopathy risk at population level.",
        },
    },
    {
        "gene": "SLCO1B1",
        "drug": "atorvastatin",
        "tiers": [
            {"id": "slco1b1_atorv_poor", "diplotype": "*5/*5", "rsid": "rs4149056 C/C",
             "phenotype": "Poor Function", "activity": "0.0",
             "recommendation": "Atorvastatin: standard dose acceptable; lower myopathy risk than simvastatin"},
            {"id": "slco1b1_atorv_decr", "diplotype": "*1a/*5", "rsid": "rs4149056 T/C",
             "phenotype": "Decreased Function", "activity": "1.0",
             "recommendation": "Atorvastatin: standard dosing"},
            {"id": "slco1b1_atorv_normal", "diplotype": "*1a/*1a", "rsid": "rs4149056 T/T",
             "phenotype": "Normal Function", "activity": "2.0",
             "recommendation": "Atorvastatin: standard dosing"},
        ],
        "pop_note": {
            "EUR": "SLCO1B1*5 carriers can use atorvastatin instead of simvastatin.",
            "AMR": "Atorvastatin preferred over simvastatin in *5 carriers.",
            "AFR": "Atorvastatin generally safe across SLCO1B1 phenotypes.",
        },
    },

    # ============================== UGT1A1 ==============================
    {
        "gene": "UGT1A1",
        "drug": "irinotecan",
        "tiers": [
            {"id": "ugt1a1_iri_pm", "diplotype": "*28/*28", "rsid": "rs8175347 TA7/TA7 (homozygous)",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Irinotecan: reduce dose by 30%; severe neutropenia/diarrhoea risk"},
            {"id": "ugt1a1_iri_im", "diplotype": "*1/*28", "rsid": "rs8175347 TA6/TA7 (heterozygous)",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Irinotecan: standard dosing with monitoring; possible dose reduction at high doses"},
            {"id": "ugt1a1_iri_nm", "diplotype": "*1/*1", "rsid": "rs8175347 TA6/TA6 (reference)",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Irinotecan: standard dosing"},
        ],
        "pop_note": {
            "EUR": "UGT1A1*28 ~35% in EUR; homozygous ~10%.",
            "AMR": "UGT1A1*28 ~30-40% in admixed Latin Americans.",
            "AFR": "UGT1A1*28 ~40-56% in African populations — higher homozygous frequency.",
        },
    },
    {
        "gene": "UGT1A1",
        "drug": "atazanavir",
        "tiers": [
            {"id": "ugt1a1_ataz_pm", "diplotype": "*28/*28", "rsid": "rs8175347 TA7/TA7",
             "phenotype": "Poor Metaboliser", "activity": "0.0",
             "recommendation": "Atazanavir: alternative protease inhibitor; hyperbilirubinaemia / jaundice risk"},
            {"id": "ugt1a1_ataz_im", "diplotype": "*1/*28", "rsid": "rs8175347 TA6/TA7",
             "phenotype": "Intermediate Metaboliser", "activity": "1.0",
             "recommendation": "Atazanavir: standard dosing; monitor bilirubin"},
            {"id": "ugt1a1_ataz_nm", "diplotype": "*1/*1", "rsid": "rs8175347 TA6/TA6",
             "phenotype": "Normal Metaboliser", "activity": "2.0",
             "recommendation": "Atazanavir: standard dosing"},
        ],
        "pop_note": {
            "EUR": "UGT1A1*28 ~35% — hyperbilirubinaemia common in HIV patients on atazanavir.",
            "AMR": "UGT1A1*28 ~30-40% — relevant for HIV protease inhibitor choice.",
            "AFR": "UGT1A1*28 ~40-56% — high prevalence affects atazanavir tolerability.",
        },
    },

    # ============================== HLA-B*57:01 ==============================
    {
        "gene": "HLA-B*57:01",
        "drug": "abacavir",
        "tiers": [
            {"id": "hlab5701_aba_pos", "diplotype": "HLA-B*57:01 positive",
             "rsid": "rs2395029 G/T or G/G (any positive)",
             "phenotype": "Positive", "activity": "carrier",
             "recommendation": "Abacavir: AVOID (potentially LETHAL hypersensitivity reaction with rechallenge)"},
            {"id": "hlab5701_aba_neg", "diplotype": "HLA-B*57:01 negative",
             "rsid": "rs2395029 T/T (no risk allele)",
             "phenotype": "Negative", "activity": "non-carrier",
             "recommendation": "Abacavir: standard dosing"},
        ],
        "pop_note": {
            "EUR": "HLA-B*57:01 ~5-8% in EUR — pre-treatment testing standard of care.",
            "AMR": "HLA-B*57:01 ~3-5% in admixed Latin Americans.",
            "AFR": "HLA-B*57:01 <1% in most African populations — abacavir hypersensitivity rare.",
        },
    },

    # ============================== HLA-B*15:02 ==============================
    {
        "gene": "HLA-B*15:02",
        "drug": "carbamazepine",
        "tiers": [
            {"id": "hlab1502_carba_pos", "diplotype": "HLA-B*15:02 positive",
             "rsid": "HLA-B typing positive for *15:02",
             "phenotype": "Positive", "activity": "carrier",
             "recommendation": "Carbamazepine: AVOID (potentially LETHAL — Stevens-Johnson syndrome / TEN risk)"},
            {"id": "hlab1502_carba_neg", "diplotype": "HLA-B*15:02 negative",
             "rsid": "HLA-B typing negative for *15:02",
             "phenotype": "Negative", "activity": "non-carrier",
             "recommendation": "Carbamazepine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "HLA-B*15:02 <0.1% in EUR — testing not routine.",
            "AMR": "HLA-B*15:02 <0.5% in admixed Latin Americans.",
            "AFR": "HLA-B*15:02 absent in most African populations.",
        },
    },
    {
        "gene": "HLA-B*15:02",
        "drug": "oxcarbazepine",
        "tiers": [
            {"id": "hlab1502_oxcarba_pos", "diplotype": "HLA-B*15:02 positive",
             "rsid": "HLA-B typing positive for *15:02",
             "phenotype": "Positive", "activity": "carrier",
             "recommendation": "Oxcarbazepine: AVOID (potentially LETHAL — SJS/TEN risk, cross-reactivity with carbamazepine)"},
            {"id": "hlab1502_oxcarba_neg", "diplotype": "HLA-B*15:02 negative",
             "rsid": "HLA-B typing negative for *15:02",
             "phenotype": "Negative", "activity": "non-carrier",
             "recommendation": "Oxcarbazepine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "HLA-B*15:02 testing not routine in EUR populations.",
            "AMR": "HLA-B*15:02 <0.5% in admixed Latin Americans.",
            "AFR": "HLA-B*15:02 absent in most African populations.",
        },
    },

    # ============================== HLA-A*31:01 ==============================
    {
        "gene": "HLA-A*31:01",
        "drug": "carbamazepine",
        "tiers": [
            {"id": "hlaa3101_carba_pos", "diplotype": "HLA-A*31:01 positive",
             "rsid": "HLA-A typing positive for *31:01",
             "phenotype": "Positive", "activity": "carrier",
             "recommendation": "Carbamazepine: AVOID (potentially LETHAL — DRESS/SJS/TEN risk)"},
            {"id": "hlaa3101_carba_neg", "diplotype": "HLA-A*31:01 negative",
             "rsid": "HLA-A typing negative for *31:01",
             "phenotype": "Negative", "activity": "non-carrier",
             "recommendation": "Carbamazepine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "HLA-A*31:01 ~2-5% in EUR — testing recommended in some guidelines.",
            "AMR": "HLA-A*31:01 ~5-10% in indigenous Latin Americans (Quechua, Aymara).",
            "AFR": "HLA-A*31:01 ~1-3% in African populations.",
        },
    },

    # ============================== HLA-B*58:01 ==============================
    {
        "gene": "HLA-B*58:01",
        "drug": "allopurinol",
        "tiers": [
            {"id": "hlab5801_allo_pos", "diplotype": "HLA-B*58:01 positive",
             "rsid": "HLA-B typing positive for *58:01",
             "phenotype": "Positive", "activity": "carrier",
             "recommendation": "Allopurinol: AVOID (potentially LETHAL SJS/TEN/DRESS risk)"},
            {"id": "hlab5801_allo_neg", "diplotype": "HLA-B*58:01 negative",
             "rsid": "HLA-B typing negative for *58:01",
             "phenotype": "Negative", "activity": "non-carrier",
             "recommendation": "Allopurinol: standard dosing"},
        ],
        "pop_note": {
            "EUR": "HLA-B*58:01 ~1-2% in EUR.",
            "AMR": "HLA-B*58:01 ~3-7% in admixed Latin Americans.",
            "AFR": "HLA-B*58:01 ~3-5% in African populations.",
        },
    },

    # ============================== G6PD ==============================
    # CPIC G6PD nomenclature uses named variants (B, A, A-, Mediterranean, etc.)
    # not star alleles. G6PD is X-linked: males are hemizygous (single allele),
    # females are homozygous (one variant on each X) or heterozygous.
    {
        "gene": "G6PD",
        "drug": "rasburicase",
        "tiers": [
            {"id": "g6pd_rasb_def_severe",
             "diplotype": "G6PD Mediterranean (hemizygous male or homozygous female)",
             "rsid": "rs5030868 A (hemizygous male) or A/A (homozygous female) [forward-genomic; transcript c.563C>T, p.Ser188Phe; G6PD on minus strand of X]",
             "phenotype": "Deficient", "activity": "0.0",
             "recommendation": "Rasburicase: AVOID (potentially LETHAL acute haemolytic anaemia; Class II severe deficiency)"},
            {"id": "g6pd_rasb_def_mild",
             "diplotype": "G6PD A- haplotype (hemizygous male or homozygous female)",
             "rsid": "rs1050828 T + rs1050829 C in cis (compound A- haplotype; hemizygous male) or rs1050828 T/T + rs1050829 C/C (homozygous female) [forward-genomic; transcript c.202G>A + c.376A>G; G6PD on minus strand of X]",
             "phenotype": "Deficient", "activity": "0.5",
             "recommendation": "Rasburicase: AVOID (haemolysis risk; Class III moderate deficiency, lower severity than Mediterranean)"},
            {"id": "g6pd_rasb_normal",
             "diplotype": "G6PD B (wildtype, hemizygous male or homozygous female)",
             "rsid": "rs1050828 C + rs5030868 G (no deficiency variants; hemizygous male) or C/C + G/G (homozygous female) [forward-genomic; G6PD on minus strand of X]",
             "phenotype": "Normal", "activity": "2.0",
             "recommendation": "Rasburicase: standard dosing"},
        ],
        "pop_note": {
            "EUR": "G6PD Mediterranean ~5-8% in Mediterranean populations; rare in northern EUR.",
            "AMR": "G6PD A- ~5-10% in admixed Latin Americans with African ancestry; testing recommended.",
            "AFR": "G6PD A- ~10-15% in West/Central African populations; A- haplotype (202A in cis with 376G) most common.",
        },
    },
    {
        "gene": "G6PD",
        "drug": "primaquine",
        "tiers": [
            {"id": "g6pd_prim_def",
             "diplotype": "G6PD A- haplotype (hemizygous male or homozygous female)",
             "rsid": "rs1050828 T + rs1050829 C in cis (compound A- haplotype) [forward-genomic; transcript c.202G>A + c.376A>G; G6PD on minus strand of X]",
             "phenotype": "Deficient", "activity": "0.5",
             "recommendation": "Primaquine: AVOID (potentially LETHAL haemolysis in deficient individuals)"},
            {"id": "g6pd_prim_int",
             "diplotype": "G6PD B/A- (heterozygous female)",
             "rsid": "rs1050828 C/T + rs1050829 T/C or T/T (heterozygous female; X-inactivation creates variable enzyme activity) [forward-genomic]",
             "phenotype": "Variable", "activity": "1.0",
             "recommendation": "Primaquine: caution; assess quantitative G6PD enzyme activity before treatment"},
            {"id": "g6pd_prim_normal",
             "diplotype": "G6PD B (wildtype, hemizygous male or homozygous female)",
             "rsid": "rs1050828 C + rs5030868 G (no deficiency variants) [forward-genomic; G6PD on minus strand of X]",
             "phenotype": "Normal", "activity": "2.0",
             "recommendation": "Primaquine: standard dosing"},
        ],
        "pop_note": {
            "EUR": "G6PD deficiency rare in northern EUR; Mediterranean variants higher in southern EUR.",
            "AMR": "G6PD A- ~5-10% in Latin Americans with African ancestry; relevant for malaria treatment.",
            "AFR": "G6PD A- ~10-15% in West/Central Africa; primaquine for vivax malaria requires testing.",
        },
    },

    # ============================== IFNL3 (formerly IL28B) ==============================
    {
        "gene": "IFNL3",
        "drug": "PEG-IFN-α",
        "tiers": [
            {"id": "ifnl3_pegifn_fav", "diplotype": "rs12979860 C/C (favourable)",
             "rsid": "rs12979860 C/C (homozygous reference)",
             "phenotype": "Favourable Response", "activity": "responder",
             "recommendation": "PEG-IFN-α: higher response rate; standard dosing"},
            {"id": "ifnl3_pegifn_unf", "diplotype": "rs12979860 T/T (unfavourable)",
             "rsid": "rs12979860 T/T (homozygous variant)",
             "phenotype": "Unfavourable Response", "activity": "poor responder",
             "recommendation": "PEG-IFN-α: lower response rate; consider direct-acting antivirals if available"},
        ],
        "pop_note": {
            "EUR": "IFNL3 rs12979860 C ~70% in EUR.",
            "AMR": "IFNL3 rs12979860 C ~50-70% in admixed Latin Americans.",
            "AFR": "IFNL3 rs12979860 T-allele ~50-60% — lower PEG-IFN response in African populations.",
        },
    },

    # ============================== VKORC1 ==============================
    {
        "gene": "VKORC1",
        "drug": "warfarin",
        "tiers": [
            {"id": "vkorc1_warf_sens", "diplotype": "rs9923231 A/A (homozygous)",
             "rsid": "rs9923231 A/A (homozygous, -1639 promoter)",
             "phenotype": "Sensitive", "activity": "low expression",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm (e.g., IWPC or Gage) incorporating VKORC1 + CYP2C9; expected lower maintenance dose"},
            {"id": "vkorc1_warf_normal", "diplotype": "rs9923231 G/G (wildtype)",
             "rsid": "rs9923231 G/G (homozygous reference)",
             "phenotype": "Normal", "activity": "normal expression",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm (e.g., IWPC or Gage); standard dosing acceptable if algorithm unavailable"},
        ],
        "pop_note": {
            "EUR": "VKORC1 -1639A allele ~40% in EUR.",
            "AMR": "VKORC1 -1639A ~30-50% in admixed Latin Americans.",
            "AFR": "VKORC1 -1639A ~10-15% in African populations — explains higher warfarin doses.",
        },
    },

    # ============================== CYP4F2 ==============================
    {
        "gene": "CYP4F2",
        "drug": "warfarin",
        "tiers": [
            {"id": "cyp4f2_warf_decr", "diplotype": "*3/*3", "rsid": "rs2108622 T/T (homozygous)",
             "phenotype": "Decreased Function", "activity": "0.5",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm incorporating CYP2C9 + VKORC1 + CYP4F2; CYP4F2*3 contributes ~5-10% higher maintenance dose only when combined with the other markers"},
            {"id": "cyp4f2_warf_normal", "diplotype": "*1/*1", "rsid": "rs2108622 C/C (reference)",
             "phenotype": "Normal Function", "activity": "2.0",
             "recommendation": "Warfarin: use validated pharmacogenetic dosing algorithm; standard dosing acceptable if algorithm unavailable"},
        ],
        "pop_note": {
            "EUR": "CYP4F2*3 ~30% in EUR — minor warfarin dose effect.",
            "AMR": "CYP4F2*3 ~20-30% in admixed Latin Americans.",
            "AFR": "CYP4F2*3 ~10% in African populations.",
        },
    },

    # ============================== CFTR ==============================
    {
        "gene": "CFTR",
        "drug": "ivacaftor",
        "tiers": [
            {"id": "cftr_iva_resp", "diplotype": "G551D mutation present",
             "rsid": "rs75527207 G/A or A/A (G551D positive)",
             "phenotype": "Responsive", "activity": "responder",
             "recommendation": "Ivacaftor: indicated; gating mutation responsive"},
            {"id": "cftr_iva_nonresp", "diplotype": "F508del homozygous (no gating mutation)",
             "rsid": "rs113993960 (F508del homozygous)",
             "phenotype": "Non-responsive", "activity": "non-responder",
             "recommendation": "Ivacaftor: not indicated (use combination therapy with elexacaftor/tezacaftor)"},
        ],
        "pop_note": {
            "EUR": "CFTR F508del most common in EUR (~70% of CF mutations).",
            "AMR": "CFTR mutation spectrum diverse in admixed Latin Americans.",
            "AFR": "CFTR mutations less characterised in African populations; F508del less common.",
        },
    },

    # ============================== MT-RNR1 ==============================
    {
        "gene": "MT-RNR1",
        "drug": "aminoglycosides",
        "tiers": [
            {"id": "mtrnr1_amino_var", "diplotype": "m.1555A>G",
             "rsid": "m.1555A>G (homoplasmic mitochondrial variant)",
             "phenotype": "Variant Carrier", "activity": "carrier",
             "recommendation": "Aminoglycosides: AVOID (potentially LETHAL ototoxicity — irreversible deafness)"},
            {"id": "mtrnr1_amino_wt", "diplotype": "m.1555A (wildtype)",
             "rsid": "m.1555A (no variant)",
             "phenotype": "Normal", "activity": "wildtype",
             "recommendation": "Aminoglycosides: standard dosing"},
        ],
        "pop_note": {
            "EUR": "MT-RNR1 m.1555A>G ~0.2% in EUR populations.",
            "AMR": "MT-RNR1 m.1555A>G frequency unknown in most Latin American populations.",
            "AFR": "MT-RNR1 m.1555A>G poorly characterised in African populations.",
        },
    },

    # ============================== RYR1 ==============================
    {
        "gene": "RYR1",
        "drug": "volatile-anaesthetics",
        "tiers": [
            {"id": "ryr1_anaes_susc", "diplotype": "RYR1 pathogenic variant present",
             "rsid": "RYR1 mutation (e.g. p.Arg614Cys, c.1840C>T)",
             "phenotype": "Malignant Hyperthermia Susceptible", "activity": "carrier",
             "recommendation": "Volatile anaesthetics + succinylcholine: AVOID (potentially LETHAL malignant hyperthermia)"},
            {"id": "ryr1_anaes_norm", "diplotype": "RYR1 wildtype",
             "rsid": "RYR1 no pathogenic variant",
             "phenotype": "Normal", "activity": "non-carrier",
             "recommendation": "Volatile anaesthetics: standard use"},
        ],
        "pop_note": {
            "EUR": "RYR1 pathogenic variants ~1:2000 in EUR — anaesthesia screening relevant.",
            "AMR": "RYR1 pathogenic variant frequency variable in Latin Americans.",
            "AFR": "RYR1 pathogenic variants poorly characterised in African populations.",
        },
    },
]


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

def generate_cases() -> list[dict]:
    """Expand GUIDELINES into flat case-list compatible with 02-run-benchmark."""
    cases = []
    for g in GUIDELINES:
        for tier in g["tiers"]:
            cases.append({
                "id": tier["id"],
                "gene": g["gene"],
                "drug": g["drug"],
                "genotype": tier["rsid"],
                "gt_diplotype": tier["diplotype"],
                "gt_phenotype": tier["phenotype"],
                "gt_drug": f"{g['drug']}: {tier['recommendation'].split(': ', 1)[-1] if ': ' in tier['recommendation'] else tier['recommendation']}",
                "gt_activity": tier["activity"],
                "pop_note": g["pop_note"],
            })
    return cases


def write_spot_check(cases: list[dict], path: Path, n: int = 10) -> None:
    """Pick 10 cases stratified across phenotype categories for diverse manual review."""
    by_id = {c["id"]: c for c in cases}
    target_ids = [
        "cyp2d6_codeine_pm",        # classic PM, high-impact opioid
        "cyp2d6_codeine_um",        # UM (lethal morphine accumulation)
        "cyp2c19_clop_rm",          # Rapid metaboliser (rare tier)
        "cyp2c9_warf_im",           # Heterozygous IM
        "cyp3a5_tacr_nm_expressor", # Expressor (CYP3A5*1/*1 — common in AFR, missed in EUR panels)
        "dpyd_fu_pm",               # Lethal DPD deficiency
        "tpmt_aza_pm",              # Lethal myelosuppression
        "slco1b1_simva_decr",       # Function-tier (not metaboliser)
        "hlab5701_aba_pos",         # HLA carrier (binary)
        "g6pd_rasb_def_severe",     # X-linked deficiency (Mediterranean variant)
    ]
    sample = [by_id[i] for i in target_ids if i in by_id]
    # Backfill if any target id didn't generate
    if len(sample) < n:
        for c in cases:
            if c not in sample:
                sample.append(c)
            if len(sample) >= n:
                break

    lines = [
        "# Spot-validation sample — BiB benchmark v3 test cases",
        "",
        f"Sample of {len(sample)} cases across distinct gene-drug pairs. Cross-check each entry against the corresponding CPIC guideline (https://cpicpgx.org/guidelines/). Flag any mismatch in `gt_phenotype`, `gt_diplotype`, or `gt_drug`.",
        "",
    ]
    for i, c in enumerate(sample, 1):
        lines.extend([
            f"## {i}. {c['id']}",
            f"- **Gene:** {c['gene']}",
            f"- **Drug:** {c['drug']}",
            f"- **Genotype (rsID):** {c['genotype']}",
            f"- **GT diplotype:** {c['gt_diplotype']}",
            f"- **GT phenotype:** {c['gt_phenotype']}",
            f"- **GT activity:** {c['gt_activity']}",
            f"- **GT recommendation:** {c['gt_drug']}",
            f"- **Pop note (AFR):** {c['pop_note']['AFR']}",
            "",
        ])
    path.write_text("\n".join(lines))


def write_log(cases: list[dict], path: Path) -> None:
    by_gene = {}
    for c in cases:
        by_gene.setdefault(c["gene"], {}).setdefault(c["drug"], 0)
        by_gene[c["gene"]][c["drug"]] += 1
    lines = [f"# Case generation log — total: {len(cases)} cases", ""]
    for gene in sorted(by_gene):
        for drug, n in sorted(by_gene[gene].items()):
            lines.append(f"  {gene:<15} × {drug:<25} = {n} tier(s)")
    pending = sum(1 for c in cases if "[CHECK]" in str(c["pop_note"]))
    lines.append("")
    lines.append(f"Cases with [CHECK] population notes (pending Heinner/Segun review): {pending}")
    path.write_text("\n".join(lines))


def main():
    cases = generate_cases()
    out_json = SPECS / "test_cases_v3.json"
    out_md = SPECS / "spot_check_sample.md"
    out_log = SPECS / "generation_log.txt"

    out_json.write_text(json.dumps(cases, indent=2))
    write_spot_check(cases, out_md)
    write_log(cases, out_log)

    print(f"Wrote {out_json}  ({len(cases)} cases)")
    print(f"Wrote {out_md}    (10-case spot-check sample)")
    print(f"Wrote {out_log}")
    print("\n--- Generation log ---")
    print(out_log.read_text())


if __name__ == "__main__":
    main()
