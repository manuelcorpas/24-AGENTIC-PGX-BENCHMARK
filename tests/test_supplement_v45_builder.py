from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from zipfile import ZipFile

from lxml import etree


REPO = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO
    / "docs"
    / "deepseek-v4-flash-rerun-audit-2026-08-30"
    / "build_supplement_v45.py"
)


def load_module():
    spec = spec_from_file_location("build_supplement_v45", SCRIPT)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_canonical_tables_bind_to_v4_flash_outputs():
    module = load_module()
    tables = module.canonical_tables(REPO)

    assert len(tables) == 9
    assert tables[0][1][-2:] == ["46", "6"]
    assert tables[0][-1][1:3] == ["2639", "99.96%"]
    assert tables[1][1][1:] == ["12/120 (10.0%)", "34/216 (15.7%)"]
    assert tables[2][-1][1:] == ["46/336", "137/336", "22/336", "24/336", "23/336"]
    assert tables[5][-1][2:] == ["87 / 90", "3 / 90", "0 / 90"]
    assert tables[6][1][1:] == ["137", "22"]
    assert tables[6][2][1:] == ["134", "29"]

    deepseek_rows = [row for row in tables[8] if row[0].endswith("DeepSeek V4 Flash")]
    assert deepseek_rows == [
        ["No definitions: DeepSeek V4 Flash", "527", "267", "259", "0", "1", "50.8%"],
        ["With definitions: DeepSeek V4 Flash", "527", "439", "88", "0", "0", "83.3%"],
    ]
    assert "DeepSeek V3" not in repr(tables)


def test_updated_captions_remove_stale_model_and_counts():
    module = load_module()
    updates = module.paragraph_updates()
    text = "\n".join(updates.values())
    assert "DeepSeek V4 Flash" in text
    assert "DeepSeek V3" not in text
    assert "54, 137, 23, 26 and 25" not in text
    assert "46, 137, 22, 24 and 23" in text
    assert "v4_flash_population_three_arm_panel.json" in text
    assert "pending synchronized evidence release" not in text
    assert module.table_titles()["Table S1."] == "Table S1. Headline metrics for the five matched configurations"


def test_build_preserves_word_compatibility_namespaces(tmp_path):
    module = load_module()
    source = (
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
        / "01-Supplementary-v44.docx"
    )
    if not source.exists():
        return
    output = tmp_path / "supplement-v45.docx"
    module.build(source, output, REPO)
    with ZipFile(output) as archive:
        root = etree.fromstring(archive.read("word/document.xml"))
    ignorable = root.get(
        "{http://schemas.openxmlformats.org/markup-compatibility/2006}Ignorable", ""
    ).split()
    assert ignorable
    assert all(prefix in root.nsmap for prefix in ignorable)
    w = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    paragraphs = [
        "".join(node.text or "" for node in para.iter(w + "t"))
        for para in root.iter(w + "p")
    ]
    s1 = [text for text in paragraphs if text.startswith("Table S1.")]
    assert len(s1) == 2
    assert s1[0] == "Table S1. Headline metrics for the five matched configurations"
    assert "All rates use 2,640 attempted evaluations" in s1[1]


def test_revision_figures_use_readable_compact_labels():
    path = REPO / "code" / "70-figures-revision.py"
    spec = spec_from_file_location("revision_figures", path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    assert all("\n" in label for label in module.FIGURE_R3_LABELS)
    assert module.FIGURE_R2_CLUSTER_CUTOFF == (0.97, 0.98)
    assert module.FIGURE_R2_LABEL_OFFSETS["HLA-B*15:02"] != (5, -3)
    assert module.FIGURE_R2_LABEL_OFFSETS["NUDT15"] != (5, -3)
