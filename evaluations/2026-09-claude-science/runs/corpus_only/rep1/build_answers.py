"""
build_answers.py -- resolve the 110 blinded pharmacogenomics cases.

Pipeline
  1. load cases_blinded.json and cpic_guideline_corpus.json (local files only)
  2. verify every case gene has a guideline excerpt in the corpus
  3. genotype string -> diplotype + phenotype  (pgx_tables.GENOTYPE_RULES)
  4. (gene, drug, phenotype) -> recommendation (pgx_tables.RECOMMENDATIONS)
  5. abstain where step 3 or step 4 has no confident answer
  6. write claude_science_answers.csv and the supporting tables

No network access is used anywhere in this pipeline.

Run:  python build_answers.py
"""
import csv
import hashlib
import json
import re
import sys
from collections import Counter

import pgx_tables as T

OUT_ANSWERS = "claude_science_answers.csv"
OUT_ALLELES = "table_allele_definitions.csv"
OUT_DIP2PHENO = "table_diplotype_to_phenotype.csv"
OUT_PHENO2REC = "table_phenotype_to_recommendation.csv"
OUT_RATIONALE = "case_rationale.csv"


def load_inputs():
    with open(T.CASES_PATH) as fh:
        cases = json.load(fh)
    with open(T.CORPUS_PATH) as fh:
        corpus = json.load(fh)
    return cases, corpus


def check_corpus_coverage(cases, corpus):
    """Every case gene must have guideline text, and the named drug should
    appear in that gene's drugs_covered list. Report anything missing."""
    genes = corpus["genes"]
    problems = []
    for c in cases:
        g = c["gene"]
        if g not in genes:
            problems.append((c["id"], g, c["drug"], "gene absent from corpus"))
            continue
        if not genes[g]["guideline_excerpt"].strip():
            problems.append((c["id"], g, c["drug"], "empty guideline excerpt"))
        if c["drug"] not in genes[g]["drugs_covered"]:
            problems.append((c["id"], g, c["drug"], "drug not in drugs_covered"))
    return problems


def resolve_genotype(gene, genotype):
    """Apply the gene's ordered regex rules. Returns
    (diplotype, activity_score, phenotype, note) or None."""
    for rule in T.GENOTYPE_RULES.get(gene, []):
        pattern, diplo, act, pheno, note = rule
        if re.search(pattern, genotype, flags=re.IGNORECASE):
            return diplo, act, pheno, note
    return None


def lethal_class(gene, drug, phenotype):
    key = (gene, phenotype)
    if key not in T.LETHAL_CLASS:
        return False, ""
    restrict = T.LETHAL_DRUG_RESTRICTION.get(key)
    if restrict is not None and drug not in restrict:
        return False, ""
    return True, T.LETHAL_CLASS[key]


def main():
    cases, corpus = load_inputs()
    assert len(cases) == 110, "expected 110 cases, got %d" % len(cases)

    problems = check_corpus_coverage(cases, corpus)
    if problems:
        print("corpus coverage problems:")
        for p in problems:
            print("  ", p)
    else:
        print("corpus coverage: all %d cases have gene-level guideline text and an "
              "in-scope drug" % len(cases))

    answers, rationale = [], []
    for c in cases:
        gene, drug, gt, cid = c["gene"], c["drug"], c["genotype"], c["id"]
        res = resolve_genotype(gene, gt)

        if res is None:
            answers.append({
                "case_id": cid, "diplotype": "abstained", "phenotype": "abstained",
                "recommendation": "abstained", "abstained": "true"})
            rationale.append({
                "case_id": cid, "gene": gene, "drug": drug, "genotype": gt,
                "diplotype": "", "activity_score": "", "phenotype": "",
                "recommendation": "", "guideline_strength": "", "source_pa_id": "",
                "abstained": "true", "abstain_reason":
                    "genotype string could not be mapped to a diplotype with confidence",
                "lethal_class": "false", "lethal_basis": "", "mapping_note": ""})
            continue

        diplo, act, pheno, note = res
        rec_entry = T.RECOMMENDATIONS.get((gene, drug, pheno))

        if rec_entry is None:
            # diplotype/phenotype resolvable, but the supplied corpus carries no
            # phenotype-indexed recommendation for this gene-drug pair.
            reason = ("the supplied corpus excerpt for %s contains no phenotype "
                      "assignment or %s-specific dosing action for this gene "
                      "(only the CYP2C9 allele table and the VKORC1 x CYP2C9 dose "
                      "table are present), so no recommendation can be grounded "
                      "in the available text" % (gene, drug))
            answers.append({
                "case_id": cid, "diplotype": diplo, "phenotype": pheno,
                "recommendation": "abstained - no %s recommendation in supplied corpus" % gene,
                "abstained": "true"})
            rationale.append({
                "case_id": cid, "gene": gene, "drug": drug, "genotype": gt,
                "diplotype": diplo, "activity_score": act or "", "phenotype": pheno,
                "recommendation": "", "guideline_strength": "", "source_pa_id": "",
                "abstained": "true", "abstain_reason": reason,
                "lethal_class": "false", "lethal_basis": "", "mapping_note": note})
            continue

        rec, strength, pa = rec_entry
        leth, leth_basis = lethal_class(gene, drug, pheno)
        answers.append({
            "case_id": cid, "diplotype": diplo, "phenotype": pheno,
            "recommendation": rec, "abstained": "false"})
        rationale.append({
            "case_id": cid, "gene": gene, "drug": drug, "genotype": gt,
            "diplotype": diplo, "activity_score": act or "", "phenotype": pheno,
            "recommendation": rec, "guideline_strength": strength, "source_pa_id": pa,
            "abstained": "false", "abstain_reason": "",
            "lethal_class": "true" if leth else "false", "lethal_basis": leth_basis,
            "mapping_note": note})

    # --- write answers -----------------------------------------------------
    with open(OUT_ANSWERS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "case_id", "diplotype", "phenotype", "recommendation", "abstained"])
        w.writeheader()
        w.writerows(answers)

    with open(OUT_RATIONALE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "case_id", "gene", "drug", "genotype", "diplotype", "activity_score",
            "phenotype", "recommendation", "guideline_strength", "source_pa_id",
            "abstained", "abstain_reason", "lethal_class", "lethal_basis",
            "mapping_note"])
        w.writeheader()
        w.writerows(rationale)

    # --- write the knowledge tables ---------------------------------------
    with open(OUT_ALLELES, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "marker", "variant_allele", "defines",
                    "allele_function", "source"])
        w.writerows(T.ALLELE_DEFS)

    with open(OUT_DIP2PHENO, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "genotype_pattern_regex", "diplotype", "activity_score",
                    "phenotype", "note"])
        for gene, rules in T.GENOTYPE_RULES.items():
            for pattern, diplo, act, pheno, note in rules:
                w.writerow([gene, pattern, diplo, act or "", pheno, note])

    with open(OUT_PHENO2REC, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "drug", "phenotype", "recommendation",
                    "guideline_strength", "source_pa_id"])
        for (gene, drug, pheno), (rec, strength, pa) in T.RECOMMENDATIONS.items():
            w.writerow([gene, drug, pheno, rec, strength, pa])

    # --- summary -----------------------------------------------------------
    n_abs = sum(1 for a in answers if a["abstained"] == "true")
    n_leth = sum(1 for r in rationale if r["lethal_class"] == "true")
    print("cases resolved      : %d / %d" % (len(answers) - n_abs, len(answers)))
    print("abstentions         : %d" % n_abs)
    print("lethal-class cases  : %d" % n_leth)
    print("phenotypes assigned :")
    for p, n in sorted(Counter(r["phenotype"] for r in rationale if r["phenotype"]).items()):
        print("   %-70s %d" % (p, n))
    print("input checksums:")
    for path in (T.CASES_PATH, T.CORPUS_PATH):
        with open(path, "rb") as fh:
            print("   %-28s sha256 %s" % (path, hashlib.sha256(fh.read()).hexdigest()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
