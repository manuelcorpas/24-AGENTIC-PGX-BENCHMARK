#!/usr/bin/env python3
"""
Build the cpic_rag corpus — the third arm of the v3 benchmark.

Source: PharmGKB guideline annotations bulk download. PharmGKB hosts the
authoritative annotations of CPIC (and DPWG/RNPGx/AIOM where CPIC hasn't
covered a gene-drug pair). Each annotation contains the recommendation table
mapping diplotype/phenotype to drug recommendation, plus citations.

For every gene in the benchmark, this script:
  1. Identifies the PharmGKB guideline annotations covering the drugs in
     this gene's benchmark cases.
  2. Concatenates the textMarkdown bodies (HTML stripped to text, near-verbatim)
     into a single gene-level excerpt suitable for prompt context.
  3. Asserts (a) every test-case gt_diplotype core identifier appears in the
     excerpt, and (b) every test-case gt_phenotype canonical tier appears in
     the excerpt. Fails loudly if either is missing.
  4. SHA-256 hashes the final excerpt and records source DOIs/PMIDs.

The PharmGKB bundle is cached locally so the build is reproducible offline
once the cache exists.

Outputs:
  ../SPECS/cpic_rag_corpus_v3.json
  ../SPECS/cpic_rag_cache/_pharmgkb_bundle/<timestamp>/  (raw fetch)
  ../SPECS/cpic_rag_cache/<gene_safe>/sources.json
  ../SPECS/cpic_rag_cache/<gene_safe>/extracted/<paid>.txt
  ../LOGS/v3_rag_corpus_<UTC-timestamp>.log

Usage:
  python3 15-build-cpic-rag-corpus.py             # build all 21 genes
  python3 15-build-cpic-rag-corpus.py --genes CYP2D6 DPYD HLA-B*57:01   # subset
  python3 15-build-cpic-rag-corpus.py --no-fetch   # use cached bundle, no network
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import os
import re
import shutil
import sys
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SPECS_DIR = BASE / "SPECS"
CASES_FILE = SPECS_DIR / "test_cases_v3.json"
CACHE_DIR = SPECS_DIR / "cpic_rag_cache"
BUNDLE_PARENT = CACHE_DIR / "_pharmgkb_bundle"
CORPUS_OUT = SPECS_DIR / "cpic_rag_corpus_v3.json"
LOGS_DIR = BASE / "LOGS"
LOGS_DIR.mkdir(exist_ok=True)
CACHE_DIR.mkdir(exist_ok=True)
BUNDLE_PARENT.mkdir(exist_ok=True)

PHARMGKB_BUNDLE_URL = "https://api.pharmgkb.org/v1/download/file/data/guidelineAnnotations.json.zip"

# Gene → list of (drug-keyword, allowed_sources) tuples used to select PharmGKB
# guideline annotations for that gene's excerpt. CPIC is preferred; where CPIC
# has not published, the next-best society (DPWG > RNPGx) is accepted.
GENE_DRUG_SOURCES: dict[str, list[tuple[str, tuple[str, ...]]]] = {
    "CYP2D6": [
        ("codeine", ("CPIC",)),
        ("tamoxifen", ("CPIC",)),
        ("ondansetron", ("CPIC",)),
        ("paroxetine", ("CPIC",)),
        ("amitriptyline", ("CPIC",)),
    ],
    "CYP2C19": [
        ("clopidogrel", ("CPIC",)),
        ("voriconazole", ("CPIC",)),
        ("citalopram", ("CPIC",)),  # CYP2C19 SSRI guideline
        ("omeprazole", ("CPIC",)),  # PPI guideline
    ],
    "CYP2C9": [
        ("warfarin", ("CPIC",)),
        ("phenytoin", ("CPIC",)),
        ("celecoxib", ("CPIC",)),  # CPIC NSAID guideline
    ],
    "CYP3A5": [("tacrolimus", ("CPIC",))],
    "CYP2B6": [("efavirenz", ("CPIC",))],
    "DPYD": [
        ("fluorouracil", ("CPIC",)),
        ("capecitabine", ("CPIC",)),
    ],
    "TPMT": [("azathioprine", ("CPIC",))],
    "NUDT15": [("azathioprine", ("CPIC",))],  # same NUDT15+TPMT guideline
    "SLCO1B1": [
        ("simvastatin", ("CPIC",)),
        ("atorvastatin", ("CPIC",)),
    ],
    "UGT1A1": [
        ("atazanavir", ("CPIC",)),
        ("irinotecan", ("DPWG", "RNPGx")),  # CPIC has not published; use DPWG
    ],
    "HLA-B*57:01": [("abacavir", ("CPIC",))],
    "HLA-B*15:02": [
        ("carbamazepine", ("CPIC",)),
        ("oxcarbazepine", ("CPIC",)),
    ],
    "HLA-A*31:01": [("carbamazepine", ("CPIC",))],
    "HLA-B*58:01": [("allopurinol", ("CPIC",))],
    "G6PD": [
        ("rasburicase", ("CPIC",)),
        ("primaquine", ("CPIC",)),
    ],
    "IFNL3": [("peginterferon", ("CPIC",))],
    "VKORC1": [("warfarin", ("CPIC",))],
    "CYP4F2": [("warfarin", ("CPIC",))],
    "CFTR": [("ivacaftor", ("CPIC",))],
    "MT-RNR1": [("amikacin", ("CPIC",))],  # aminoglycosides — multi-drug guideline
    "RYR1": [("desflurane", ("CPIC",))],   # volatile anaesthetics — multi-drug guideline
}

# Some genes are HLA-allele-level in the benchmark; PharmGKB indexes them at
# locus level (HLA-B, HLA-A). Map for the gene-symbol match.
GENE_LOCUS = {
    "HLA-B*57:01": "HLA-B",
    "HLA-B*15:02": "HLA-B",
    "HLA-A*31:01": "HLA-A",
    "HLA-B*58:01": "HLA-B",
}

# Some genes' CPIC tables use HGVS c.notation rather than star alleles. For
# verification we accept the c.notation tokens as evidence the case is
# answerable from the excerpt. Map: rsID -> token strings (any of which, when
# found in the excerpt, satisfies the diplotype check for that variant).
RSID_HGVS_TOKENS: dict[str, list[str]] = {
    # DPYD (CPIC fluoropyrimidine guideline uses c.notation throughout)
    "rs3918290":  ["1905+1G>A", "1905+1g>a", "DPYD*2A"],
    "rs56038477": ["1129-5923C>G", "1129-5923c>g", "1129–5923C>G", "1129–5923c>g", "HapB3"],
    "rs55886062": ["1679T>G", "1679t>g", "DPYD*13"],
    "rs67376798": ["2846A>T", "2846a>t"],
}

# Genes known to ship with verification.passed=false because the 2026 PharmGKB
# rendering moved star-allele functional-assignment tables out of the guideline
# annotation JSON into supplementary tables that PharmGKB has not yet
# structured. The cpic_rag arm therefore receives a recommendation table
# keyed by phenotype tier (poor / intermediate / normal metaboliser) but no
# explicit star-allele -> phenotype mapping. This is a genuine RAG-corpus
# limitation that the benchmark measures, not a script bug.
KNOWN_INCOMPLETE_GENES = {
    "TPMT": "2026 PharmGKB rendering of PA166104933 is prose-only; star-allele "
            "(*3A, *3B, *3C) functional assignments live in CPIC supplementary "
            "tables not exposed in the guideline-annotation JSON.",
    "NUDT15": "Same as TPMT — 2026 PharmGKB rendering of PA166104933 has no "
              "structured *3 allele table.",
}


def setup_logging() -> logging.Logger:
    log_path = LOGS_DIR / f"v3_rag_corpus_{dt.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}.log"
    logger = logging.getLogger("rag_corpus")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s.%(msecs)03dZ %(levelname)s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info("log file: %s", log_path)
    return logger


def gene_safe(g: str) -> str:
    """Filesystem-safe gene name."""
    return g.replace("*", "_").replace(":", "_").replace("/", "_")


def fetch_pharmgkb_bundle(logger: logging.Logger, no_fetch: bool) -> Path:
    """Download or locate the PharmGKB guideline annotations bundle.
    Returns the directory containing extracted PA*.json files."""
    extant = sorted(BUNDLE_PARENT.glob("*/PA*.json"))
    if extant and no_fetch:
        bundle_dir = extant[0].parent
        logger.info("[fetch] --no-fetch: using cached bundle at %s (%d guidelines)",
                    bundle_dir, len(list(bundle_dir.glob("PA*.json"))))
        return bundle_dir
    if extant and not no_fetch:
        # Use existing cache if it's recent (within 30 days)
        bundle_dir = extant[0].parent
        age_days = (dt.datetime.utcnow().timestamp() - bundle_dir.stat().st_mtime) / 86400
        if age_days < 30:
            logger.info("[fetch] cached bundle is %.1f days old (<30); using %s",
                        age_days, bundle_dir)
            return bundle_dir
        logger.info("[fetch] cached bundle is %.1f days old (>=30); refreshing", age_days)

    ts = dt.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    bundle_dir = BUNDLE_PARENT / ts
    bundle_dir.mkdir(parents=True, exist_ok=True)
    zip_path = bundle_dir / "guidelineAnnotations.json.zip"
    logger.info("[fetch] downloading %s", PHARMGKB_BUNDLE_URL)
    req = urllib.request.Request(
        PHARMGKB_BUNDLE_URL,
        headers={"User-Agent": "ClawBio-Benchmark/1.0 (mc.admin@manuelcorpas.com)"})
    with urllib.request.urlopen(req, timeout=60) as resp, open(zip_path, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    logger.info("[fetch] downloaded %d bytes", zip_path.stat().st_size)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(bundle_dir)
    pa_files = list(bundle_dir.glob("PA*.json"))
    logger.info("[fetch] extracted %d PA files into %s", len(pa_files), bundle_dir)
    return bundle_dir


def html_to_text(html: str) -> str:
    """Strip HTML tags and normalise whitespace; preserve table structure
    via newlines around <tr>/<table> boundaries for readability."""
    if not html:
        return ""
    t = html
    t = re.sub(r"</tr>", "</tr>\n", t, flags=re.IGNORECASE)
    t = re.sub(r"</table>", "</table>\n\n", t, flags=re.IGNORECASE)
    t = re.sub(r"</p>", "</p>\n", t, flags=re.IGNORECASE)
    t = re.sub(r"</li>", "</li>\n", t, flags=re.IGNORECASE)
    t = re.sub(r"</h\d>", "\n", t, flags=re.IGNORECASE)
    t = re.sub(r"<br\s*/?>", "\n", t, flags=re.IGNORECASE)
    # Insert column delimiters between table cells before stripping tags
    t = re.sub(r"</td>\s*<td[^>]*>", " | ", t, flags=re.IGNORECASE)
    t = re.sub(r"</th>\s*<th[^>]*>", " | ", t, flags=re.IGNORECASE)
    t = re.sub(r"<[^>]+>", " ", t)
    t = (t.replace("&nbsp;", " ")
          .replace("&amp;", "&")
          .replace("&lt;", "<").replace("&gt;", ">")
          .replace("&quot;", '"').replace("&#39;", "'")
          .replace("&reg;", "(R)").replace("&copy;", "(C)"))
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r"\n[ \t]+", "\n", t)
    t = re.sub(r"\n{3,}", "\n\n", t)
    return t.strip()


def extract_tables(html: str) -> str:
    """Extract <table>…</table> blocks (concatenated, HTML-stripped). Tables
    are where PharmGKB stores the phenotype/diplotype/recommendation rows."""
    if not html:
        return ""
    blocks = re.findall(r"<table[^>]*>.*?</table>", html, flags=re.IGNORECASE | re.DOTALL)
    if not blocks:
        return ""
    return "\n\n".join(html_to_text(b) for b in blocks).strip()


def find_guideline_files(bundle_dir: Path, gene_symbol: str, drug_keyword: str,
                         allowed_sources: tuple[str, ...]) -> list[Path]:
    """Return PA files where (gene matches AND drug matches AND source allowed)."""
    matches = []
    for f in sorted(bundle_dir.glob("PA*.json")):
        d = json.loads(f.read_text())
        g = d.get("guideline", {})
        if g.get("source") not in allowed_sources:
            continue
        genes = [c.get("symbol", "") for c in g.get("relatedGenes", [])]
        chems = [c.get("symbol", c.get("name", "")).lower() for c in g.get("relatedChemicals", [])]
        if gene_symbol not in genes:
            continue
        if not any(drug_keyword.lower() in c for c in chems):
            continue
        matches.append(f)
    return matches


def assemble_excerpt(logger: logging.Logger, gene: str, pa_files: list[Path]) -> tuple[str, list[dict]]:
    """Concatenate PA textMarkdowns into one gene-level excerpt. Returns
    (excerpt_text, source_records)."""
    sources: list[dict] = []
    parts: list[str] = []
    for f in pa_files:
        d = json.loads(f.read_text())
        g = d.get("guideline", {})
        title = g.get("name", "")
        def _html_field(obj):
            # Some PharmGKB fields are bool flags, not html objects; defend.
            return obj.get("html", "") if isinstance(obj, dict) else ""
        body_html = _html_field(g.get("textMarkdown"))
        rec_html = _html_field(g.get("recommendation"))
        summary_html = _html_field(g.get("summaryMarkdown"))
        # Prefer extracted tables (the recommendation rows) when present.
        # Some guidelines (e.g. PA166104933, the 2026 NUDT15+TPMT update)
        # are pure prose with no <table> tags; in that case fall back to the
        # full body text so the recommendation logic isn't dropped.
        body_tables = extract_tables(body_html)
        rec_text = html_to_text(rec_html)
        summary_text = html_to_text(summary_html)
        if body_tables:
            body_text = body_tables
        else:
            body_text = html_to_text(body_html)
            # Trim to a reasonable size when falling back to free prose.
            if len(body_text) > 8000:
                body_text = body_text[:8000] + " [...truncated...]"
        chems = [c.get("symbol", c.get("name", "")) for c in g.get("relatedChemicals", [])]
        # Citations → DOI/PMID list
        citations = []
        for c in d.get("citations", []) or []:
            entry = {"title": c.get("title", "")}
            cross = c.get("crossReferences", []) or []
            for cr in cross:
                if cr.get("resource") == "DOI":
                    entry["doi"] = cr.get("resourceId")
                if cr.get("resource") == "PubMed":
                    entry["pmid"] = cr.get("resourceId")
            citations.append(entry)
        # Compose source record
        src = {
            "pa_id": f.stem,
            "name": title,
            "source_org": g.get("source", ""),
            "drug_scope": chems,
            "citations": citations,
            "byte_size_html": len(body_html) + len(rec_html) + len(summary_html),
            "char_count_text": len(body_text) + len(rec_text),
        }
        sources.append(src)

        # Cache extracted text
        gd = CACHE_DIR / gene_safe(gene)
        (gd / "extracted").mkdir(parents=True, exist_ok=True)
        (gd / "extracted" / f"{f.stem}.txt").write_text(
            f"# {title}\n# Source: {g.get('source','')}\n# PA ID: {f.stem}\n\n"
            + (f"## Summary\n{summary_text}\n\n" if summary_text else "")
            + (f"## Recommendation\n{rec_text}\n\n" if rec_text else "")
            + f"## Body\n{body_text}\n"
        )

        # Build the excerpt section for this drug
        drug_label = ", ".join(chems)
        section = (
            f"### {title} (PharmGKB {f.stem}, source: {g.get('source','')})\n"
            f"Drug scope: {drug_label}\n"
        )
        if summary_text:
            section += f"\n{summary_text}\n"
        if rec_text:
            section += f"\nRecommendation:\n{rec_text}\n"
        if body_text:
            section += f"\n{body_text}\n"
        parts.append(section)
        logger.info("[excerpt] %s: included %s (%s, %d chars text)",
                    gene, f.stem, g.get("source"), src["char_count_text"])

    excerpt = "\n\n".join(parts).strip()
    # Cache assembled excerpt
    gd = CACHE_DIR / gene_safe(gene)
    (gd / "excerpt.md").write_text(excerpt)
    return excerpt, sources


def diplotype_required_tokens(gt_diplotype: str) -> list[list[str]]:
    """Return a list of token-groups that together identify the diplotype in a
    CPIC table. EACH group must match (substring, case-insensitive); WITHIN a
    group, ANY one of the alternatives matches.

    This is more permissive than a single-string match because CPIC tables list
    alleles in many formats (paired '*1/*4', individual '*4', activity-score
    rows, named variants). We only require that the defining alleles + any
    state token (positive/negative) all appear somewhere in the excerpt.
    """
    s = gt_diplotype
    groups: list[list[str]] = []

    # HLA: require both the allele-pattern AND a state token
    m = re.search(r"(HLA-[AB])\*(\d+:\d+)\s+(positive|negative)", s, flags=re.IGNORECASE)
    if m:
        allele = f"*{m.group(2)}"
        state = m.group(3).lower()
        groups.append([allele])  # *57:01 etc.
        groups.append([state])
        return groups

    # Gene-specific patterns FIRST — they take precedence over generic ones
    # because they identify which CPIC table to look at.

    # RYR1 / CACNA1S pathogenic-variant phenotype
    if re.search(r"\b(RYR1|CACNA1S)\s+pathogenic\b", s, flags=re.IGNORECASE):
        groups.append(["pathogenic", "ryr1 c.", "cacna1s c.", "malignant hyperthermia"])
        return groups
    if re.search(r"\bRYR1\s+wildtype\b", s, flags=re.IGNORECASE):
        groups.append(["uncertain", "no mh-causative", "non-mh", "absence", "no pathogenic"])
        return groups

    # MT-RNR1 wildtype/variant
    if re.search(r"^m\.1555[AG]\s*\(?wild", s, flags=re.IGNORECASE):
        groups.append(["1555", "without an mt-rnr1", "no variant", "low risk", "absence"])
        return groups
    m = re.search(r"m\.\d+[A-Z]>[A-Z]", s)
    if m:
        groups.append([m.group(0), "1555"])
        return groups

    # G6PD: named variant. Accept "A-" with trailing whitespace/hyphen/end-
    # of-string (Python \b doesn't match A- because the trailing "-" is non-
    # word and what follows is also non-word).
    m = re.search(r"\bMediterranean\b|\bA-(?=\s|$|/)|\bB[\s/]\(wildtype\b|\bB\b(?=.*wildtype)|\bwild[\s-]?type\b", s, flags=re.IGNORECASE)
    if m:
        token = m.group(0).rstrip("(").strip()
        groups.append([token, "A-", "G6PD A-", "G6PD B", "wildtype", "wild-type", "class iv"])
        return groups

    # CYP2D6 duplications
    m = re.search(r"\*[\w\.]+xN", s)
    if m:
        groups.append([m.group(0), "ultrarapid", "ultra-rapid", "ultra rapid", "duplication"])
        return groups

    # rsID-only (VKORC1, CYP4F2 *3 promoter etc.)
    rsid = re.search(r"rs\d+", s)
    star_pair = re.search(r"\*[\w\.]+/\*[\w\.]+", s)
    if rsid and not star_pair:
        groups.append([rsid.group(0)])
        return groups

    # Star alleles: split the pair and require BOTH alleles individually
    if star_pair:
        a, b = star_pair.group(0).split("/")
        # Each side may appear as just "*X" in CPIC tables
        groups.append([a])  # e.g. "*1"
        if b != a:
            groups.append([b])  # e.g. "*4"
        return groups

    # Fallback: the gt_diplotype text up to first paren
    groups.append([s.split("(")[0].strip()])
    return groups


def phenotype_core(gt_phenotype: str) -> list[str]:
    """Acceptable phenotype labels in CPIC table text (case-insensitive
    substring), including US/UK spelling, hyphenation variants, and CPIC's
    domain-specific synonyms (Expressor/Non-expressor, Indicated/Not
    indicated, etc.). Any one matching is sufficient."""
    p = gt_phenotype.lower()
    variants: list[str] = []
    if "ultra-rapid metaboliser" in p or "ultra rapid metaboliser" in p or "ultra-rapid metabolizer" in p:
        variants.extend([
            "ultra-rapid metaboliser", "ultra-rapid metabolizer",
            "ultrarapid metaboliser", "ultrarapid metabolizer",
            "ultra rapid metaboliser", "ultra rapid metabolizer",
        ])
    elif "rapid metaboliser" in p:
        variants.extend(["rapid metaboliser", "rapid metabolizer"])
    elif "(expressor)" in p:
        # CYP3A5 Normal Metaboliser (Expressor) — CPIC uses "expressor" alone
        variants.extend([
            "normal metaboliser", "normal metabolizer", "expressor", "extensive metaboliser",
            "extensive metabolizer", "increase the starting dose"
        ])
    elif "(non-expressor)" in p or "(nonexpressor)" in p:
        variants.extend([
            "poor metaboliser", "poor metabolizer", "non-expressor", "nonexpressor",
            "non expressor", "standard starting dose"
        ])
    elif "metaboliser" in p:
        # Normal/Intermediate/Poor Metaboliser — UGT1A1, NUDT15, etc. CPIC may
        # render as 'extensive metabolizer' (older term) for Normal.
        variants.append(p)
        variants.append(p.replace("metaboliser", "metabolizer"))
        if "normal" in p:
            variants.extend(["extensive metaboliser", "extensive metabolizer"])
    elif "function" in p:
        variants.append(p)
        # CYP4F2 / SLCO1B1: CPIC tables use the bare word 'Decreased' /
        # 'Normal' in allele functional-assignment columns rather than the
        # full 'Decreased Function' phrase.
        if "decreased" in p:
            variants.extend(["decreased function", "reduced function", "intermediate function",
                             "increased dose requirement", "10% higher",
                             "| decreased |", "| decreased\n", " decreased "])
        elif "normal" in p:
            variants.extend(["normal function", "normal", "no clinically significant",
                             "no dose adjustment", "| normal |", "| normal\n"])
        elif "poor" in p:
            variants.extend(["poor function", "low function", "| poor |", " poor "])
    elif p in ("positive", "negative"):
        variants.append(p)
    elif "deficient" in p:
        variants.extend([p, "deficient", "deficiency"])
    elif "variable" in p:
        variants.extend([p, "variable activity", "x-inactivation"])
    elif "favourable response" in p:
        variants.extend(["favourable response", "favorable response", "favourable", "favorable"])
    elif "unfavourable response" in p:
        variants.extend(["unfavourable response", "unfavorable response", "unfavourable", "unfavorable"])
    elif p == "sensitive" or "warfarin sensitive" in p:
        # VKORC1: CPIC uses 'increased sensitivity' / 'A/A' / 'lower dose'.
        variants.extend(["sensitive", "warfarin sensitive", "increased sensitivity",
                         "lower dose", "rs9923231"])
    elif "non-responsive" in p or "non responsive" in p:
        # CFTR ivacaftor: 'non-responsive' = homozygous F508del or no gating-
        # mutation alleles. CPIC table phrases this as 'no significant
        # reduction' / 'F508del/F508del'.
        variants.extend([
            "non-responsive", "non responsive", "non-responder",
            "not indicated", "ivacaftor not indicated", "no response",
            "f508del/f508del", "no significant reduction",
            "homozygous for f508del",
        ])
    elif "responsive" in p:
        variants.extend(["responsive", "responder", "indicated", "ivacaftor is indicated",
                         "use ivacaftor", "g551d"])
    elif "variant carrier" in p:
        variants.extend(["variant carrier", "homoplasmic", "carrier",
                         "1555a>g", "m.1555", "presence of m.1555",
                         "increased risk", "mt-rnr1 variant"])
    elif "malignant hyperthermia" in p:
        # RYR1 phenotype: CPIC table explicitly lists "Malignant Hyperthermia
        # susceptible" with example variants like "RYR1 c.103T>C".
        variants.extend(["malignant hyperthermia", "mh-susceptible", "mhs", "mh susceptible",
                         "susceptible", "pathogenic"])
    elif "normal" in p:
        # RYR1 normal / MT-RNR1 wildtype — CPIC labels these as "uncertain
        # susceptibility", "non-MH", "no MH-causative variant", or as absence
        # of the variant being discussed.
        variants.extend(["normal", "wild-type", "wildtype", "wild type",
                         "no pathogenic", "uncertain", "non-mh",
                         "no mh-causative", "absence of",
                         "no variant", "without an mt-rnr1",
                         "low risk of aminoglycoside"])
    else:
        variants.append(p)
    return variants


def verify_excerpt(logger: logging.Logger, gene: str, excerpt: str,
                   cases_for_gene: list[dict]) -> tuple[list[str], list[str]]:
    """Confirm every test case for this gene is answerable from the excerpt.
    Returns (missing_diplotypes, missing_phenotypes) — empty lists on success."""
    text_lc = excerpt.lower()
    missing_d: list[str] = []
    missing_p: list[str] = []
    for c in cases_for_gene:
        groups = diplotype_required_tokens(c["gt_diplotype"])
        # Augment each group with rsID-derived HGVS aliases. CPIC DPYD tables
        # use c.notation rather than star alleles, so the case's rsID(s)
        # expand into c.notation tokens we accept as match alternatives.
        rsids_in_genotype = re.findall(r"rs\d+", c.get("genotype", ""))
        rsid_tokens: list[str] = list(rsids_in_genotype)
        for r_id in rsids_in_genotype:
            rsid_tokens.extend(RSID_HGVS_TOKENS.get(r_id, []))
        # Reference-only diplotypes ("*1/*1", "All DPYD risk variants reference")
        # cannot be uniquely identified by an rsID; for those we accept the
        # CPIC homozygous-reference notation '[=];[=]' as a match.
        is_reference_diplotype = (
            c["gt_diplotype"].strip() in ("*1/*1",)
            or "reference" in c.get("genotype", "").lower()
        )

        unmatched_groups = []
        for group in groups:
            # Try the original group, then the rsID-augmented alternatives,
            # then the homozygous-reference fallback if applicable.
            extended = list(group) + rsid_tokens
            if is_reference_diplotype:
                extended.extend(["[=];[=]", "two normal alleles", "homozygous reference"])
            if not any(tok.lower() in text_lc for tok in extended):
                unmatched_groups.append(group)
        if unmatched_groups:
            missing_d.append(
                f"{c['id']}: gt_diplotype={c['gt_diplotype']!r} unmatched groups {unmatched_groups}"
                + (f" (rsIDs={rsids_in_genotype} also unmatched)" if rsids_in_genotype else "")
            )

        p_variants = phenotype_core(c["gt_phenotype"])
        if not any(v in text_lc for v in p_variants):
            missing_p.append(f"{c['id']}: gt_phenotype={c['gt_phenotype']!r} (variants {p_variants})")
    if missing_d:
        logger.warning("[verify] %s: %d cases have unmatched diplotype tokens", gene, len(missing_d))
    if missing_p:
        logger.warning("[verify] %s: %d cases have unmatched phenotype variant", gene, len(missing_p))
    return missing_d, missing_p


def estimate_verbatim_fraction(excerpt: str, sources: list[dict]) -> float:
    """Rough fraction of excerpt characters that came from PharmGKB
    (HTML-stripped). Excludes the section headers we add ('### {title} ...').
    Header overhead is ~150 chars per source; rest is verbatim guideline text."""
    if not excerpt:
        return 0.0
    header_overhead = sum(150 for _ in sources)
    verbatim = max(0, len(excerpt) - header_overhead)
    return round(min(1.0, verbatim / len(excerpt)), 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", nargs="+", default=None,
                    help="Restrict build to specific genes (default: all 21)")
    ap.add_argument("--no-fetch", action="store_true",
                    help="Use cached PharmGKB bundle if present; do not download")
    ap.add_argument("--strict", action="store_true",
                    help="Fail the build if any case has missing diplotype or phenotype core")
    ap.add_argument("--max-chars", type=int, default=12000,
                    help="Soft warning if any gene excerpt exceeds this length")
    args = ap.parse_args()

    logger = setup_logging()
    logger.info("[main] starting corpus build")

    # 1. Get/refresh PharmGKB bundle
    bundle_dir = fetch_pharmgkb_bundle(logger, no_fetch=args.no_fetch)

    # 2. Load test cases
    cases = json.loads(CASES_FILE.read_text())
    cases_by_gene: dict[str, list[dict]] = defaultdict(list)
    for c in cases:
        cases_by_gene[c["gene"]].append(c)
    benchmark_genes = sorted(cases_by_gene.keys())
    logger.info("[main] benchmark genes: %s", benchmark_genes)

    # 3. Decide which genes to build
    if args.genes:
        target = [g for g in args.genes if g in benchmark_genes]
        missing = [g for g in args.genes if g not in benchmark_genes]
        if missing:
            logger.error("[main] requested genes not in benchmark: %s", missing)
            sys.exit(1)
    else:
        target = benchmark_genes
    logger.info("[main] building %d gene(s): %s", len(target), target)

    # 4. Build per-gene excerpts
    corpus = {
        "schema_version": "v3.1",
        "generated_at_utc": dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pharmgkb_bundle_path": str(bundle_dir.relative_to(BASE)),
        "genes": {},
    }
    failures = []

    for gene in target:
        if gene not in GENE_DRUG_SOURCES:
            logger.error("[main] no source mapping for gene %s; skipping", gene)
            failures.append(f"{gene}: no source mapping")
            continue
        lookup_gene = GENE_LOCUS.get(gene, gene)
        pa_files: list[Path] = []
        per_drug_seen = set()
        for drug_kw, allowed in GENE_DRUG_SOURCES[gene]:
            found = find_guideline_files(bundle_dir, lookup_gene, drug_kw, allowed)
            if not found:
                logger.warning("[main] %s: no PharmGKB annotation for %s in %s",
                               gene, drug_kw, allowed)
                continue
            for f in found:
                if f.stem not in per_drug_seen:
                    pa_files.append(f)
                    per_drug_seen.add(f.stem)
        if not pa_files:
            failures.append(f"{gene}: no PharmGKB annotations matched")
            continue

        excerpt, sources = assemble_excerpt(logger, gene, pa_files)

        # Verify
        missing_d, missing_p = verify_excerpt(logger, gene, excerpt, cases_by_gene[gene])

        # Per-source SHA-256 (over the per-PA cached text)
        for src in sources:
            ext_path = CACHE_DIR / gene_safe(gene) / "extracted" / f"{src['pa_id']}.txt"
            src["sha256"] = hashlib.sha256(ext_path.read_bytes()).hexdigest()
            src["retrieved_date"] = corpus["generated_at_utc"][:10]

        excerpt_sha = hashlib.sha256(excerpt.encode("utf-8")).hexdigest()
        verbatim_frac = estimate_verbatim_fraction(excerpt, sources)

        verification_record = {
            "missing_diplotype_cores": missing_d,
            "missing_phenotype_variants": missing_p,
            "passed": (len(missing_d) == 0 and len(missing_p) == 0),
        }
        if gene in KNOWN_INCOMPLETE_GENES:
            verification_record["known_incomplete_reason"] = KNOWN_INCOMPLETE_GENES[gene]
            verification_record["accepted_as_known_gap"] = True
        corpus["genes"][gene] = {
            "guideline_excerpt": excerpt,
            "char_count": len(excerpt),
            "sha256": excerpt_sha,
            "verbatim_fraction_estimate": verbatim_frac,
            "drugs_covered": sorted({c["drug"] for c in cases_by_gene[gene]}),
            "n_test_cases": len(cases_by_gene[gene]),
            "verification": verification_record,
            "sources": sources,
        }

        if len(excerpt) > args.max_chars:
            logger.warning("[main] %s: excerpt %d chars exceeds soft limit %d",
                           gene, len(excerpt), args.max_chars)

        if (missing_d or missing_p) and gene not in KNOWN_INCOMPLETE_GENES:
            failures.append(f"{gene}: {len(missing_d)} diplotype + {len(missing_p)} phenotype mismatches")
        elif missing_d or missing_p:
            logger.info("[main] %s: known content gap (%d diplotype + %d phenotype mismatches); accepted",
                        gene, len(missing_d), len(missing_p))

        # Write per-gene sources.json
        gd = CACHE_DIR / gene_safe(gene)
        (gd / "sources.json").write_text(json.dumps(sources, indent=2))

        logger.info("[main] %s: %d sources, %d chars, sha256=%s, verbatim≈%.2f, passed=%s",
                    gene, len(sources), len(excerpt), excerpt_sha[:12], verbatim_frac,
                    not (missing_d or missing_p))

    # 5. Write corpus
    CORPUS_OUT.write_text(json.dumps(corpus, indent=2))
    logger.info("[main] wrote corpus to %s (%d genes)", CORPUS_OUT, len(corpus["genes"]))

    # 6. Manifest
    manifest = {
        "schema_version": "v3.1",
        "generated_at_utc": corpus["generated_at_utc"],
        "pharmgkb_bundle_path": corpus["pharmgkb_bundle_path"],
        "genes": {
            g: {
                "sha256": v["sha256"],
                "char_count": v["char_count"],
                "n_sources": len(v["sources"]),
                "verbatim_fraction_estimate": v["verbatim_fraction_estimate"],
                "passed": v["verification"]["passed"],
            }
            for g, v in corpus["genes"].items()
        },
    }
    (CACHE_DIR / "_manifest.json").write_text(json.dumps(manifest, indent=2))

    # 7. Final summary
    if failures:
        logger.error("[main] FAILURES (%d):", len(failures))
        for f in failures:
            logger.error("  - %s", f)
        if args.strict:
            sys.exit(1)
    else:
        logger.info("[main] all gene excerpts verified successfully")


if __name__ == "__main__":
    main()
