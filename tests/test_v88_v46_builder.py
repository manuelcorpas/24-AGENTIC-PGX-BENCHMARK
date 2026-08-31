from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from zipfile import ZipFile
import re

from lxml import etree


REPO = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO
    / "docs"
    / "deepseek-v4-flash-rerun-audit-2026-08-30"
    / "build_v88_v46.py"
)
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def load_module():
    spec = spec_from_file_location("build_v88_v46", SCRIPT)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def current_text(path: Path) -> str:
    with ZipFile(path) as archive:
        root = etree.fromstring(archive.read("word/document.xml"))
    for deleted in list(root.iter(W + "del")):
        deleted.getparent().remove(deleted)
    return " ".join(node.text or "" for node in root.iter(W + "t"))


def test_correction_manifest_uses_canonical_values():
    module = load_module()
    replacements = module.manuscript_replacements()
    replacement_text = "\n".join(new for old, new in replacements)

    assert "97.3%" in replacement_text
    assert "seven version-matched models" in replacement_text
    assert "87 echoed" in replacement_text
    assert "three hedged" in replacement_text
    assert "5,272 contained a usable" in replacement_text

    table_script = REPO / "code" / "89-supp-headline-and-per-gene.py"
    spec = spec_from_file_location("supp_tables", table_script)
    table_module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(table_module)
    assert table_module.pct(0.9725) == "97.3%"


def test_builds_real_tracked_v88_and_formatted_v46(tmp_path):
    module = load_module()
    current = (
        Path.home()
        / "Library"
        / "Mobile Documents"
        / "com~apple~CloudDocs"
        / "PUBLICATIONS"
        / "00-CLAWBIO-BENCHMARK"
        / "BENCHMARK"
        / ("D" + "OCS")
        / "CELL-GENOMICS"
        / "00-CURRENT"
    )
    main_source = current / "00-Manuscript-v87.docx"
    supp_source = current / "01-Supplementary-v45.docx"
    if not main_source.exists() or not supp_source.exists():
        return

    main_out = tmp_path / "00-Manuscript-v88.docx"
    supp_out = tmp_path / "01-Supplementary-v46.docx"
    module.build_manuscript(main_source, main_out)
    module.build_supplement(supp_source, supp_out)

    main_text = re.sub(r"\s+", " ", current_text(main_out))
    assert "97.2% under authored-rule execution" not in main_text
    assert "97.3% under authored-rule execution" in main_text
    assert "all eight version-matched models pooled in Figure 3B" not in main_text
    assert "seven version-matched models pooled in Figure 3B" in main_text
    assert "88 echoed the corrupted contract and two hedged" not in main_text
    assert "87 echoed the corrupted contract and three hedged" in main_text
    assert "5,271 contained a usable" not in main_text
    assert "5,272 contained a usable" in main_text

    with ZipFile(main_out) as archive:
        main_root = etree.fromstring(archive.read("word/document.xml"))
        settings = etree.fromstring(archive.read("word/settings.xml"))
    assert main_root.findall(".//" + W + "ins")
    assert main_root.findall(".//" + W + "del")
    assert settings.find(W + "trackRevisions") is not None

    with ZipFile(supp_out) as archive:
        supp_root = etree.fromstring(archive.read("word/document.xml"))
    paragraphs = list(supp_root.iter(W + "p"))
    for para in paragraphs:
        text = "".join(node.text or "" for node in para.iter(W + "t"))
        if text.strip() in {"Population context", "Per-gene robustness"}:
            runs = [run for run in para.iter(W + "r") if "".join(t.text or "" for t in run.iter(W + "t"))]
            assert runs
            for run in runs:
                props = run.find(W + "rPr")
                assert props is not None
                assert props.find(W + "sz").get(W + "val") == "20"
                assert props.find(W + "b") is not None
                assert props.find(W + "b").get(W + "val") not in {"0", "false", "off"}
        if text.startswith(("Table S", "Figure S")):
            runs = [run for run in para.iter(W + "r") if "".join(t.text or "" for t in run.iter(W + "t"))]
            assert runs
            for run in runs:
                props = run.find(W + "rPr")
                assert props is not None
                assert props.find(W + "sz").get(W + "val") == "18"
                assert props.find(W + "rFonts").get(W + "ascii") == "Arial"

    tables = list(supp_root.iter(W + "tbl"))
    assert "97.3%" in "".join(t.text or "" for t in tables[0].iter(W + "t"))
    canonical_script = (
        REPO
        / "docs"
        / "deepseek-v4-flash-rerun-audit-2026-08-30"
        / "build_supplement_v45.py"
    )
    canonical_spec = spec_from_file_location("build_supplement_v45", canonical_script)
    canonical_module = module_from_spec(canonical_spec)
    assert canonical_spec.loader is not None
    canonical_spec.loader.exec_module(canonical_module)
    expected_tables = canonical_module.canonical_tables(REPO)
    for table_index, expected in enumerate(expected_tables):
        if expected is None:
            continue
        actual = []
        for row in tables[table_index].findall(W + "tr"):
            actual.append(
                [
                    re.sub(
                        r"\s+",
                        " ",
                        "".join(t.text or "" for t in cell.iter(W + "t")),
                    ).strip()
                    for cell in row.findall(W + "tc")
                ]
            )
        assert actual == expected
    for table in tables:
        for run in table.iter(W + "r"):
            if "".join(t.text or "" for t in run.iter(W + "t")):
                props = run.find(W + "rPr")
                assert props is not None
                assert props.find(W + "sz").get(W + "val") == "16"
                assert props.find(W + "rFonts").get(W + "ascii") == "Arial"
