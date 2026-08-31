#!/usr/bin/env python3
"""Build Supplementary v45 from v44 and canonical V4 Flash outputs."""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree as ET


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_NS = "http://www.w3.org/XML/1998/namespace"
W = f"{{{W_NS}}}"

CELLS = [
    ("free_generation", "free generation"),
    ("rag_generation", "RAG generation"),
    ("rag_execution", "RAG call + authored execution"),
    ("skill_generation", "authored rules, generated"),
    ("skill_execution", "authored rules, executed"),
]
MODELS = [
    "Claude Opus 4.5",
    "Claude Sonnet 4.5",
    "GPT-5.2",
    "GPT-4.1",
    "o3",
    "o4-mini",
    "DeepSeek V4 Flash",
    "Gemini 2.5 Flash",
]


def _load(path: Path):
    return json.loads(path.read_text())


def _table_s1(repo: Path) -> list[list[str]]:
    lines = (repo / "data" / "v4_flash_supp_tables_s1_s3.txt").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("cell "))
    rows = []
    for line in lines[start + 1 : start + 6]:
        fields = re.split(r"\s{2,}", line.strip())
        if len(fields) != 9:
            raise ValueError(f"cannot parse Table S1 row: {line}")
        fields[1] = fields[1].replace(" ", "")
        rows.append(fields)
    return [[
        "configuration", "n", "coverage", "phenotype, baseline",
        "phenotype, equivalence", "recommendation, baseline",
        "recommendation, equivalence", "lethal errors", "parse failures",
    ], *rows]


def _lethal_tables(repo: Path) -> tuple[list[list[str]], list[list[str]]]:
    rows = _load(repo / "data" / "v4_flash_matched_scored_rows_all5.json")
    lethal = [row for row in rows if row.get("lethal_action") is not None]
    by_cell_gene: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in lethal:
        by_cell_gene[(row["cell"], row["gene"])].append(row)

    genes = [
        "HLA-A*31:01", "HLA-B*15:02", "HLA-B*57:01", "HLA-B*58:01",
        "DPYD", "TPMT", "NUDT15", "G6PD", "CYP2D6", "MT-RNR1", "RYR1",
    ]
    s3 = [["gene", *[label for _, label in CELLS]]]
    totals: dict[str, tuple[int, int]] = {}
    for gene in genes:
        row = [gene]
        for cell, _ in CELLS:
            subset = by_cell_gene[(cell, gene)]
            errors = sum(item["lethal_action"] == 0 for item in subset)
            row.append(f"{errors}/{len(subset)}")
        s3.append(row)
    total_row = ["all loci"]
    for cell, _ in CELLS:
        subset = [row for row in lethal if row["cell"] == cell]
        errors = sum(row["lethal_action"] == 0 for row in subset)
        totals[cell] = (errors, len(subset))
        total_row.append(f"{errors}/{len(subset)}")
    s3.append(total_row)

    s2 = [["configuration", "HLA risk-allele loci", "other lethal-class loci"]]
    for cell, label in CELLS:
        hla = [row for row in lethal if row["cell"] == cell and row["gene"].startswith("HLA-")]
        other = [row for row in lethal if row["cell"] == cell and not row["gene"].startswith("HLA-")]
        hla_errors = sum(row["lethal_action"] == 0 for row in hla)
        other_errors = sum(row["lethal_action"] == 0 for row in other)
        s2.append([
            label,
            f"{hla_errors}/{len(hla)} ({100 * hla_errors / len(hla):.1f}%)",
            f"{other_errors}/{len(other)} ({100 * other_errors / len(other):.1f}%)",
        ])
    return s2, s3


def _table_s4(repo: Path) -> list[list[str]]:
    data = _load(repo / "data" / "v4_flash_hierarchical_stats.json")
    out = [[
        "configuration", "mean", "case-clustered 95% CI", "model-clustered 95% CI",
        "case variance", "model variance", "framing variance", "replicate variance", "residual",
    ]]
    for cell, label in CELLS:
        item = data[cell]
        case = item["intervals"]["case_clustered"]
        model = item["intervals"]["model_clustered"]
        var = item["variance_components"]
        out.append([
            label,
            f"{100 * case['mean']:.1f}%",
            f"{100 * case['ci95'][0]:.1f}% to {100 * case['ci95'][1]:.1f}%",
            f"{100 * model['ci95'][0]:.1f}% to {100 * model['ci95'][1]:.1f}%",
            f"{var['case']:.4f}",
            f"{var['model']:.5f}",
            f"{var['framing']:.5f}",
            f"{var['replicate']:.6f}",
            f"{var['residual']:.5f}",
        ])
    return out


def _table_s5(repo: Path) -> list[list[str]]:
    data = _load(repo / "data" / "v4_flash_crossed_mixed_model.json")
    out = [[
        "configuration", "mean", "fitted var(case)", "fitted var(model)",
        "moments var(case)", "moments var(model)",
    ]]
    for cell, label in CELLS:
        item = data[cell]
        out.append([
            label,
            f"{item['intercept']:.3f}",
            f"{item['var_case']:.4f}",
            f"{item['var_model']:.5f}",
            f"{item['moments_case']:.4f}",
            f"{item['moments_model']:.5f}",
        ])
    return out


def _table_s6(repo: Path) -> list[list[str]]:
    def counts(name: str) -> tuple[int, int, int]:
        rows = _load(repo / "data" / name)
        c = Counter(row.get("classification") for row in rows)
        echo, hedge = c["ECHO_SCRAMBLED"], c["HEDGE"]
        return echo, hedge, len(rows) - echo - hedge

    forward = counts("v4_flash_adversarial_forward_panel.json")
    reverse = counts("v4_flash_adversarial_reverse_panel.json")
    combined = tuple(a + b for a, b in zip(forward, reverse))
    cases = "DPYD, CYP2D6, HLA-B*57:01, TPMT, CYP2C19"
    return [
        ["Direction", "Cases (n=5 each)", "ECHO scrambled", "HEDGE (still unsafe)", "Override to CPIC truth"],
        ["Forward (lethal -> safe)", cases, f"{forward[0]} / 45", f"{forward[1]} / 45", f"{forward[2]} / 45"],
        ["Reverse (safe -> dangerous)", cases, f"{reverse[0]} / 45", f"{reverse[1]} / 45", f"{reverse[2]} / 45"],
        ["Combined", "Both directions", f"{combined[0]} / 90", f"{combined[1]} / 90", f"{combined[2]} / 90"],
    ]


def _table_s7(repo: Path) -> list[list[str]]:
    main = _load(repo / "data" / "v4_flash_hierarchical_stats.json")
    keyed = _load(repo / "data" / "v4_flash_rag_drugkeyed_summary.json")["baseline"]
    return [
        ["retrieval indexing", "lethal-class errors, model generates", "lethal-class errors, code executes"],
        ["gene-keyed", str(main["rag_generation"]["lethal_errors"]), str(main["rag_execution"]["lethal_errors"])],
        ["(gene, drug)-keyed", str(keyed["rag_generation"]["lethal_errors"]), str(keyed["rag_execution"]["lethal_errors"])],
    ]


def _table_s9(repo: Path) -> list[list[str]]:
    freeze = _load(repo / "data" / "v4_flash_input_normalisation_eight_model_freeze.json")
    rows = [[
        "arm and model", "attempted", "call", "abstention", "provider error",
        "truncated", "behavioural coverage",
    ]]
    for arm_key, arm_label in (("without_definitions", "No definitions"), ("with_definitions", "With definitions")):
        for model in MODELS:
            op = freeze["models_analysis"][model][arm_key]["operational"]
            denom = op["call"] + op["abstain"]
            coverage = 100 * op["call"] / denom if denom else 0
            rows.append([
                f"{arm_label}: {model}", str(op["attempted"]), str(op["call"]),
                str(op["abstain"]), str(op["error"]), str(op["truncated_output"]),
                f"{coverage:.1f}%",
            ])
    return rows


def canonical_tables(repo: Path) -> list[list[list[str]] | None]:
    s2, s3 = _lethal_tables(repo)
    return [
        _table_s1(repo), s2, s3, _table_s4(repo), _table_s5(repo),
        _table_s6(repo), _table_s7(repo), None, _table_s9(repo),
    ]


def paragraph_updates() -> dict[str, str]:
    return {
        "Table S1.": "Table S1. Headline metrics for the five matched configurations under both scorers. All rates use 2,640 attempted evaluations per configuration. The n column reports returned records, so authored-rule execution has n = 2,639 because one attempt returned no record; that missing attempt remains a failure in the reported rates. Coverage is the fraction of attempted evaluations producing a usable answer. The clinical-equivalence layer modifies phenotype identification only. Recommendation accuracy and lethal-class error counts are identical under both scorers. Generated by code/89-supp-headline-and-per-gene.py from the canonical V4 Flash eight-model panel.",
        "Table S2.": "Table S2. Lethal-class errors split by locus class. The split is by locus class, not by the action required. The three configurations that apply authored structured rules at the decision point make no lethal-class errors on HLA loci. Counts are regenerated from v4_flash_matched_scored_rows_all5.json.",
        "Table S3.": "Table S3. Per-gene lethal-class errors across the five matched configurations, using 336 attempted lethal-class evaluations per configuration. RAG generation contributes 137 errors, including 45 of 48 HLA-B*15:02 attempts. Both execution configurations apply the same authored mapping; their residual errors arise from the model's diplotype call. Generated by code/89-supp-headline-and-per-gene.py from the canonical V4 Flash panel.",
        "Table S4.": "Table S4. Case-clustered and model-clustered confidence intervals and variance decomposition for phenotype accuracy. Each of the eight models evaluated the same 110 cases. Case clustering resamples cases while retaining the model panel; model clustering resamples models while retaining the cases. Variance components are reported on the proportion scale and are regenerated from v4_flash_matched_scored_rows_all5.json by code/65-hierarchical-stats.py.",
        "Table S5.": "Table S5. Crossed random-effects model against the method-of-moments decomposition. Random intercepts for case and model were fitted by REML, one model per configuration. Near-ceiling components should be interpreted alongside the clustered intervals in Table S4. Generated by code/73-crossed-mixed-model.py from the canonical V4 Flash panel.",
        "Table S6.": "Table S6. Contract-faithfulness is symmetric. Forward and reverse tests each contain five gene-drug pairs, three models (Claude Opus 4, GPT-5.2 and DeepSeek V4 Flash) and three replicates per case. Across 90 responses, 87 echoed the corrupted contract, three hedged while retaining unsafe fields, and none reverted to the canonical CPIC recommendation. Generated from v4_flash_adversarial_forward_panel.json and v4_flash_adversarial_reverse_panel.json.",
        "Table S7.": "Table S7. Gene-keyed and gene-drug-keyed retrieval in the eight-model panel. Drug-keying changed lethal-class errors from 137 to 134 when the model generated the recommendation and from 22 to 29 when code executed the authored mapping after the model's call. The result does not support a consistent improvement from either indexing strategy. Generated from the canonical five-cell and drug-keyed V4 Flash panels.",
        "Table S9.": "Table S9. Attempted requests, calls, explicit abstentions and apparatus attrition for input normalisation. Both arms use the same 527 variant-call inputs for each of eight explicitly versioned models, including DeepSeek V4 Flash. Behavioural coverage is call divided by call plus abstention; provider errors and output-budget truncations are reported but excluded. Generated by code/85-freeze-eight-model-normalisation.py.",
        "Population-stratified phenotype accuracy": "Population-stratified phenotype accuracy in the earlier three-arm protocol showed a small aggregate EUR to AFR difference under free prompting: 78.0%, 76.9% and 77.1% for EUR, AMR and AFR on the eight-model panel. Retrieval produced 86.0%, 85.1% and 85.0%, respectively. Because the genotypes are identical across framings, these differences measure prompt-framing sensitivity rather than performance in a population. Coverage of the validated vocabulary in real cohorts, which is a population result, is reported in main-text Table 3 and Figure S6.",
        "Figure S2.": "Figure S2. Information without action. Each point is one lethal-class gene, plotted at its mean phenotype accuracy and mean rate of correct lethal-class action, one panel per configuration. Points below the diagonal identify loci where the phenotype is named but the required action is not taken. Computed from 336 lethal-class evaluations per configuration in the canonical V4 Flash eight-model panel.",
        "Figure S3.": "Figure S3. Correctness by coincidence. The 336 lethal-class evaluations per configuration are decomposed by phenotype and required-action correctness. Correctness by coincidence falls from 122 of 336 under free generation to zero wherever authored rules are applied at the decision point. Information without action rises from 7 under free generation to 58 under RAG generation. The two action-wrong quadrants sum to 46, 137, 22, 24 and 23 lethal-class errors.",
        "Figure S4.": "Figure S4. Within the two authored-execution configurations, residual benchmark error is localised to the input call. End-to-end phenotype accuracy is plotted against input-call accuracy per gene. Across 5,272 usable calls, call concordance is 95.79%; end-to-end accuracy is 96.5% and 96.7%. This localisation is conditional on the authored mapping and benchmark vocabulary and is not independent clinical validation.",
        "Figure S5.": "Figure S5. Between-model spread across the five matched configurations. Per-model phenotype accuracy for eight models. Spread falls from 21.5 percentage points in free generation to 6.1 percentage points in authored-rule execution. The RAG-assisted execution configuration executes the authored mapping and is labelled accordingly.",
        "Figure S7.": "Figure S7. Models follow a corrupted contract. Bidirectional adversarial testing across five high-severity gene-drug pairs, three models and three replicates per direction produced no reversion to the canonical CPIC answer: 87 of 90 responses echoed the corrupted contract and three hedged while retaining unsafe fields. The Wilson 95% upper bound on reversion is 4.1%. This motivates signed, version-pinned and independently reviewed artefacts with fail-closed integrity checks.",
        "Figure S8.": "Figure S8. Conditional invariance of the executed clinical mapping: a scoped control. Identical curated genotypes were presented under European, Latin American and East African framing across five arms. Bars show phenotype accuracy with Wilson 95% intervals. EUR minus AFR differences are 0.83, 0.98, 1.74, 1.86 and 0.00 percentage points for free prompting, RAG, authored-rule generation, authored-rule execution and control, respectively. The control is invariant by construction conditional on identical valid inputs and is not evidence of end-to-end equity. Sources are v4_flash_population_three_arm_panel.json and v4_flash_population_skill_arms_panel.json; unparsed responses count as failures. Generated by code/90-figure-s8-conditional-invariance.py.",
        "Figure S9.": "Figure S9. One lethal-class case traced through all five configurations. Case hlab1502_carba_pos (HLA-B*15:02 positive, carbamazepine), GPT-5.2, replicate 1. Every configuration receives the same patient text, target drug and output schema; only the knowledge supplied and the decision agent differ. Responses come from v4_flash_eight_model_five_cell.json, and code-executed recommendations are computed by the authored mapping. Under the original endpoint, RAG generation produces 45 errors among 48 evaluations for this locus, while each other configuration produces 0 of 48. The result remains endpoint-sensitive as described in the main text.",
    }


def table_titles() -> dict[str, str]:
    return {
        "Table S1.": "Table S1. Headline metrics for the five matched configurations",
        "Table S2.": "Table S2. Lethal-class errors by locus class",
        "Table S3.": "Table S3. Per-gene lethal-class errors",
        "Table S4.": "Table S4. Case- and model-clustered confidence intervals and variance decomposition",
        "Table S5.": "Table S5. Crossed random-effects model",
        "Table S6.": "Table S6. Bidirectional adversarial scrambled-specification test",
        "Table S7.": "Table S7. Retrieval indexing in the eight-model rerun",
        "Table S8.": "Table S8. Input normalisation by gene",
        "Table S9.": "Table S9. Input-normalisation operational outcomes",
    }


def _text(element: ET.Element) -> str:
    return "".join(node.text or "" for node in element.iter() if node.tag in {W + "t", W + "delText"})


def _new_run(text: str, rpr: ET.Element | None = None, bold: bool = False) -> ET.Element:
    run = ET.Element(W + "r")
    if rpr is not None:
        run.append(copy.deepcopy(rpr))
    elif bold:
        props = ET.SubElement(run, W + "rPr")
        ET.SubElement(props, W + "b")
    text_node = ET.SubElement(run, W + "t")
    if text[:1].isspace() or text[-1:].isspace():
        text_node.set(f"{{{XML_NS}}}space", "preserve")
    text_node.text = text
    return run


def _set_cell_text(cell: ET.Element, value: str) -> None:
    tcpr = cell.find(W + "tcPr")
    para = cell.find(W + "p")
    ppr = para.find(W + "pPr") if para is not None else None
    first_run = para.find(W + "r") if para is not None else None
    rpr = first_run.find(W + "rPr") if first_run is not None else None
    for child in list(cell):
        if child is not tcpr:
            cell.remove(child)
    new_para = ET.SubElement(cell, W + "p")
    if ppr is not None:
        new_para.append(copy.deepcopy(ppr))
    new_para.append(_new_run(value, rpr))


def _replace_paragraph(para: ET.Element, replacement: str) -> None:
    ppr = para.find(W + "pPr")
    for child in list(para):
        if child is not ppr:
            para.remove(child)
    label = None
    match = re.match(r"((?:Table|Figure) S\d+\.)\s*(.*)", replacement)
    if match:
        label, replacement = match.groups()
        para.append(_new_run(label, bold=True))
        if replacement:
            para.append(_new_run(" " + replacement))
    else:
        para.append(_new_run(replacement))


def build(source: Path, output: Path, repo: Path) -> None:
    tables_data = canonical_tables(repo)
    updates = paragraph_updates()
    titles = table_titles()
    media = {
        "word/media/image2.png": repo / "figures" / "FigureS2-information-without-action.png",
        "word/media/image3.png": repo / "figures" / "FigureS3-correctness-by-coincidence.png",
        "word/media/image4.png": repo / "figures" / "FigureR2_error_localisation.png",
        "word/media/image5.png": repo / "figures" / "FigureR3_model_spread.png",
        "word/media/image7.png": repo / "figures" / "Figure6_adversarial_contract.png",
        "word/media/image8.png": repo / "figures" / "FigureS8_conditional_invariance.png",
        "word/media/image9.png": repo / "figures" / "FigureS9.png",
    }
    for path in media.values():
        if not path.exists():
            raise FileNotFoundError(path)

    with ZipFile(source) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
        tables = list(root.iter(W + "tbl"))
        if len(tables) != 9:
            raise ValueError(f"expected 9 supplementary tables, found {len(tables)}")
        for index, matrix in enumerate(tables_data):
            if matrix is None:
                continue
            rows = tables[index].findall(W + "tr")
            if len(rows) != len(matrix):
                raise ValueError(f"Table S{index + 1}: {len(rows)} rows, expected {len(matrix)}")
            for row_node, row_values in zip(rows, matrix):
                cells = row_node.findall(W + "tc")
                if len(cells) != len(row_values):
                    raise ValueError(
                        f"Table S{index + 1}: {len(cells)} columns, expected {len(row_values)}"
                    )
                for cell, value in zip(cells, row_values):
                    _set_cell_text(cell, value)

        found: Counter[str] = Counter()
        for para in root.iter(W + "p"):
            current = _text(para).strip()
            for prefix, replacement in updates.items():
                if current.startswith(prefix):
                    if prefix in titles and found[prefix] == 0:
                        _replace_paragraph(para, titles[prefix])
                    else:
                        _replace_paragraph(para, replacement)
                    found[prefix] += 1
                    break
        missing = set(updates) - set(found)
        if missing:
            raise ValueError(f"paragraph replacements not found: {sorted(missing)}")
        wrong_counts = {
            prefix: count
            for prefix, count in found.items()
            if count != (2 if prefix in titles else 1)
        }
        if wrong_counts:
            raise ValueError(f"unexpected paragraph replacement counts: {wrong_counts}")

        document_xml = ET.tostring(
            root, encoding="UTF-8", xml_declaration=True, standalone=True
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix="supp-v45-", suffix=".docx", dir=output.parent)
        os.close(fd)
        tmp = Path(tmp_name)
        try:
            with ZipFile(tmp, "w", ZIP_DEFLATED) as target:
                for item in archive.infolist():
                    if item.filename == "word/document.xml":
                        payload = document_xml
                    elif item.filename in media:
                        payload = media[item.filename].read_bytes()
                    else:
                        payload = archive.read(item.filename)
                    target.writestr(item, payload)
            os.replace(tmp, output)
        finally:
            if tmp.exists():
                tmp.unlink()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    build(args.source, args.output, args.repo)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
