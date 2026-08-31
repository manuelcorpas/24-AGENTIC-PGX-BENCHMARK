#!/usr/bin/env python3
"""Build tracked manuscript v88 and consistently formatted Supplementary v46."""

from __future__ import annotations

import argparse
import copy
from importlib.util import module_from_spec, spec_from_file_location
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree as ET


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_NS = "http://www.w3.org/XML/1998/namespace"
W = f"{{{W_NS}}}"


def manuscript_replacements() -> list[tuple[str, str]]:
    return [
        ("97.2%", "97.3%"),
        (
            "Stacked bars pool all eight version-matched models.",
            "Stacked bars pool the seven version-matched models with both curated and real-genome results.",
        ),
        (
            "They contribute equally because each receives the same states.",
            "They contribute equally because each receives the same states. DeepSeek V4 Flash is shown in Panel A but excluded from Panel B because no curated V4 benchmark result is available.",
        ),
        (
            "all eight version-matched models pooled in Figure 3B",
            "the seven version-matched models pooled in Figure 3B",
        ),
        (
            "Figure 3B pools all eight version-matched models",
            "Figure 3B pools the seven version-matched models with matched curated benchmark results",
        ),
        (
            "; equal state counts give each model equal weight.",
            "; equal state counts give each model equal weight. DeepSeek V4 Flash is shown in Figure 3A but excluded from Figure 3B because no curated V4 benchmark result is available.",
        ),
        (
            "88 echoed the corrupted contract and two hedged",
            "87 echoed the corrupted contract and three hedged",
        ),
        ("5,271 contained a usable", "5,272 contained a usable"),
    ]


def _zip_replace(source: Path, output: Path, replacements: dict[str, bytes], prefix: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=prefix, suffix=".docx", dir=output.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with ZipFile(source) as archive, ZipFile(tmp, "w", ZIP_DEFLATED) as target:
            for item in archive.infolist():
                target.writestr(item, replacements.get(item.filename, archive.read(item.filename)))
        os.replace(tmp, output)
    finally:
        if tmp.exists():
            tmp.unlink()


def _max_revision_id(root: ET.Element) -> int:
    values = []
    for tag in (W + "ins", W + "del"):
        for node in root.iter(tag):
            try:
                values.append(int(node.get(W + "id", "0")))
            except ValueError:
                pass
    return max(values, default=0)


def _run_with_text(template: ET.Element, text: str, deletion: bool = False) -> ET.Element:
    run = ET.Element(W + "r")
    props = template.find(W + "rPr")
    if props is not None:
        run.append(copy.deepcopy(props))
    node = ET.SubElement(run, W + ("delText" if deletion else "t"))
    if text[:1].isspace() or text[-1:].isspace():
        node.set(f"{{{XML_NS}}}space", "preserve")
    node.text = text
    return run


def _revision(tag: str, revision_id: int, run: ET.Element) -> ET.Element:
    node = ET.Element(W + tag)
    node.set(W + "id", str(revision_id))
    node.set(W + "author", "Codex")
    node.set(W + "date", datetime.now(timezone.utc).replace(microsecond=0).isoformat())
    node.append(run)
    return node


def _replace_in_paragraph(para: ET.Element, old: str, new: str, next_id: list[int]) -> bool:
    for text_node in para.iter(W + "t"):
        ancestor = text_node.getparent()
        skip = False
        while ancestor is not None and ancestor is not para:
            if ancestor.tag == W + "ins" and ancestor.get(W + "author") == "Codex":
                skip = True
                break
            ancestor = ancestor.getparent()
        if skip:
            continue
        if old not in (text_node.text or ""):
            continue
        run = text_node.getparent()
        if run.tag != W + "r":
            continue
        full = text_node.text or ""
        before, after = full.split(old, 1)
        top = run
        while top.getparent() is not para:
            top = top.getparent()
        parent = para
        position = parent.index(top)
        sequence: list[ET.Element] = []

        if top.tag == W + "r":
            if before:
                sequence.append(_run_with_text(run, before))
            next_id[0] += 1
            sequence.append(_revision("del", next_id[0], _run_with_text(run, old, deletion=True)))
            next_id[0] += 1
            sequence.append(_revision("ins", next_id[0], _run_with_text(run, new)))
            if after:
                sequence.append(_run_with_text(run, after))
        elif top.tag == W + "ins":
            children = list(top)
            run_index = children.index(run)
            if run_index or before:
                retained_before = copy.deepcopy(top)
                for child in list(retained_before):
                    retained_before.remove(child)
                for child in children[:run_index]:
                    retained_before.append(copy.deepcopy(child))
                if before:
                    retained_before.append(_run_with_text(run, before))
                sequence.append(retained_before)
            next_id[0] += 1
            sequence.append(_revision("del", next_id[0], _run_with_text(run, old, deletion=True)))
            next_id[0] += 1
            sequence.append(_revision("ins", next_id[0], _run_with_text(run, new)))
            if after or run_index + 1 < len(children):
                retained_after = copy.deepcopy(top)
                for child in list(retained_after):
                    retained_after.remove(child)
                if after:
                    retained_after.append(_run_with_text(run, after))
                for child in children[run_index + 1 :]:
                    retained_after.append(copy.deepcopy(child))
                sequence.append(retained_after)
        else:
            raise ValueError(f"unsupported tracked-change container: {top.tag}")

        parent.remove(top)
        for offset, node in enumerate(sequence):
            parent.insert(position + offset, node)
        return True
    return False


def build_manuscript(source: Path, output: Path) -> None:
    with ZipFile(source) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
        settings = ET.fromstring(archive.read("word/settings.xml"))

    next_id = [_max_revision_id(root)]
    counts: dict[str, int] = {}
    for old, new in manuscript_replacements():
        count = 0
        while True:
            changed = False
            for para in root.iter(W + "p"):
                if _replace_in_paragraph(para, old, new, next_id):
                    changed = True
                    count += 1
                    break
            if not changed:
                break
        if count == 0:
            raise ValueError(f"manuscript replacement not found: {old}")
        counts[old] = count

    expected = {
        "97.2%": 4,
        "Stacked bars pool all eight version-matched models.": 1,
        "They contribute equally because each receives the same states.": 1,
        "all eight version-matched models pooled in Figure 3B": 1,
        "Figure 3B pools all eight version-matched models": 1,
        "; equal state counts give each model equal weight.": 1,
        "88 echoed the corrupted contract and two hedged": 1,
        "5,271 contained a usable": 1,
    }
    if counts != expected:
        raise ValueError(f"unexpected manuscript replacement counts: {counts}, expected {expected}")

    if settings.find(W + "trackRevisions") is None:
        settings.insert(0, ET.Element(W + "trackRevisions"))
    payloads = {
        "word/document.xml": ET.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True),
        "word/settings.xml": ET.tostring(settings, encoding="UTF-8", xml_declaration=True, standalone=True),
    }
    _zip_replace(source, output, payloads, "manuscript-v88-")


def _ensure_props(run: ET.Element) -> ET.Element:
    props = run.find(W + "rPr")
    if props is None:
        props = ET.Element(W + "rPr")
        run.insert(0, props)
    return props


def _set_run_style(run: ET.Element, size: int, bold: bool | None = None) -> None:
    props = _ensure_props(run)
    fonts = props.find(W + "rFonts")
    if fonts is None:
        fonts = ET.Element(W + "rFonts")
        props.insert(0, fonts)
    for key in ("ascii", "hAnsi", "eastAsia", "cs"):
        fonts.set(W + key, "Arial")
    for name in ("sz", "szCs"):
        node = props.find(W + name)
        if node is None:
            node = ET.SubElement(props, W + name)
        node.set(W + "val", str(size))
    if bold is not None:
        for name in ("b", "bCs"):
            node = props.find(W + name)
            if bold:
                if node is None:
                    node = ET.SubElement(props, W + name)
                node.attrib.pop(W + "val", None)
            elif not bold and node is not None:
                props.remove(node)


def _set_paragraph_spacing(para: ET.Element, before: int, after: int, line: int = 240) -> None:
    ppr = para.find(W + "pPr")
    if ppr is None:
        ppr = ET.Element(W + "pPr")
        para.insert(0, ppr)
    spacing = ppr.find(W + "spacing")
    if spacing is None:
        spacing = ET.SubElement(ppr, W + "spacing")
    spacing.set(W + "before", str(before))
    spacing.set(W + "after", str(after))
    spacing.set(W + "line", str(line))
    spacing.set(W + "lineRule", "auto")


def _paragraph_text(para: ET.Element) -> str:
    return "".join(node.text or "" for node in para.iter(W + "t"))


def _next_content_sibling(para: ET.Element) -> ET.Element | None:
    parent = para.getparent()
    index = parent.index(para)
    for sibling in list(parent)[index + 1 :]:
        if sibling.tag == W + "tbl":
            return sibling
        if sibling.tag == W + "p" and _paragraph_text(sibling).strip():
            return sibling
    return None


def _style_caption(para: ET.Element, short_table_title: bool) -> None:
    first_text_run = True
    for run in para.iter(W + "r"):
        text = "".join(node.text or "" for node in run.iter(W + "t"))
        if not text:
            continue
        label_run = first_text_run and bool(re.match(r"(?:Table|Figure) S\d+\.", text))
        _set_run_style(run, 18, bold=(short_table_title or label_run))
        first_text_run = False
    _set_paragraph_spacing(para, 120 if short_table_title else 60, 60 if short_table_title else 120)


def build_supplement(source: Path, output: Path) -> None:
    with ZipFile(source) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))

    canonical_path = Path(__file__).with_name("build_supplement_v45.py")
    spec = spec_from_file_location("build_supplement_v45", canonical_path)
    canonical_module = module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"cannot load canonical supplement builder: {canonical_path}")
    spec.loader.exec_module(canonical_module)
    tables = list(root.iter(W + "tbl"))
    canonical_tables = canonical_module.canonical_tables(Path(__file__).resolve().parents[2])
    if len(tables) != len(canonical_tables):
        raise ValueError(f"expected {len(canonical_tables)} supplementary tables, found {len(tables)}")
    for table_index, matrix in enumerate(canonical_tables):
        if matrix is None:
            continue
        rows = tables[table_index].findall(W + "tr")
        if len(rows) != len(matrix):
            raise ValueError(f"Table S{table_index + 1}: row count mismatch")
        for row, values in zip(rows, matrix):
            cells = row.findall(W + "tc")
            if len(cells) != len(values):
                raise ValueError(f"Table S{table_index + 1}: column count mismatch")
            for cell, value in zip(cells, values):
                canonical_module._set_cell_text(cell, value)

    in_results = False
    for para in root.iter(W + "p"):
        text = _paragraph_text(para).strip()
        if not text:
            continue
        if text.startswith(("Table S", "Figure S")):
            short_title = text.startswith("Table S") and _next_content_sibling(para) is not None and _next_content_sibling(para).tag == W + "tbl"
            _style_caption(para, short_title)
            continue
        if text == "Supplementary Results":
            in_results = True
            for run in para.iter(W + "r"):
                _set_run_style(run, 24, bold=True)
            _set_paragraph_spacing(para, 240, 80)
            continue
        if in_results and text in {"Population context", "Per-gene robustness"}:
            for run in para.iter(W + "r"):
                _set_run_style(run, 20, bold=True)
            _set_paragraph_spacing(para, 120, 60)
            continue
        if in_results:
            for run in para.iter(W + "r"):
                _set_run_style(run, 18, bold=False)
            _set_paragraph_spacing(para, 0, 120)

    for table in root.iter(W + "tbl"):
        for row_index, row in enumerate(table.findall(W + "tr")):
            for run in row.iter(W + "r"):
                if "".join(node.text or "" for node in run.iter(W + "t")):
                    _set_run_style(run, 16, bold=(row_index == 0))

    for text_node in root.iter(W + "t"):
        if text_node.text and "97.2%" in text_node.text:
            text_node.text = text_node.text.replace("97.2%", "97.3%")

    payload = ET.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True)
    _zip_replace(source, output, {"word/document.xml": payload}, "supplement-v46-")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--main-source", type=Path, required=True)
    parser.add_argument("--main-output", type=Path, required=True)
    parser.add_argument("--supp-source", type=Path, required=True)
    parser.add_argument("--supp-output", type=Path, required=True)
    args = parser.parse_args()
    build_manuscript(args.main_source, args.main_output)
    build_supplement(args.supp_source, args.supp_output)
    print(f"wrote {args.main_output}")
    print(f"wrote {args.supp_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
