#!/usr/bin/env python3
"""
PharmCAT comparator: an independent implementation of CPIC (revision item N2;
Reviewer 2 points 2 and 7).

WHAT THE REVIEWERS SAID
R2.2 and R2.7: the executed mapping is validated only against the authors' own
encoded rules, so "the skill computes the right answer" risks being circular.
An implementation of CPIC that we did not write is needed.

WHY PHARMCAT, AND WHY VIA OUTSIDE CALLS
PharmCAT is the reference open-source CPIC implementation, maintained by
PharmGKB. Its phenotyper accepts "outside calls" (gene plus diplotype), which
lets us compare the MAPPING step directly: the same diplotype goes to our
validated skill and to PharmCAT, and we compare the phenotypes. Using outside
calls rather than VCFs deliberately isolates the mapping from the calling step,
which is the comparison R2.2 asks for. The separate 23andMe-to-VCF preprocessing
dependency for the Corpas family affects only a VCF-based comparison and is
reported as a limitation there, not here.

WHAT IS COMPARED
For every benchmark case, three phenotypes:
    ours      the validated skill's executed mapping (the paper's claim)
    pharmcat  PharmCAT 3.4.0's CPIC phenotype for the same diplotype
    truth     the benchmark ground truth encoded from the CPIC tables
Agreement between `ours` and `pharmcat` is the independent check. Disagreement
is a finding to report, not an error to silence.

SPELLING
PharmCAT emits US spelling ("Metabolizer"); the benchmark uses UK spelling
("Metaboliser"). Normalisation is explicit and tested, because a spelling
mismatch silently scored as disagreement would fabricate a discordance rate.

PINNED VERSION
PharmCAT 3.4.0, tools/pharmcat-3.4.0-all.jar
sha256 9317ef632bf6c9786ff0d9d455d4c9f6d2882ebd66ad7256b4ae958ddf454741

USAGE
    python code/63-pharmcat-comparator.py            # compare on all 110 cases
    python code/63-pharmcat-comparator.py --fetch    # print the download command
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
JAR = BASE / "tools" / "pharmcat-3.4.0-all.jar"
JAR_SHA256 = "9317ef632bf6c9786ff0d9d455d4c9f6d2882ebd66ad7256b4ae958ddf454741"
JAR_URL = ("https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/"
           "pharmcat-3.4.0-all.jar")
OUT = BASE / "data" / "v3_pharmcat_comparison.json"
REPORT = BASE / "data" / "v3_pharmcat_comparison_report.txt"

_spec = spec_from_file_location("_pgx_rules", CODE / "_pgx_rules.py")
rules = module_from_spec(_spec)
_spec.loader.exec_module(rules)

CASES = rules.load_cases()
DIP2PHEN, REC, RULES_DIP, RULES_REC = rules.build_rules(CASES)
execute_skill = rules.make_executor(DIP2PHEN, REC)

# The benchmark names HLA loci by allele; PharmCAT reports them by gene.
GENE_TO_PHARMCAT = {
    "HLA-B*57:01": "HLA-B", "HLA-B*15:02": "HLA-B", "HLA-B*58:01": "HLA-B",
    "HLA-A*31:01": "HLA-A",
}


def pharmcat_gene(gene: str) -> str:
    return GENE_TO_PHARMCAT.get(gene, gene)


# PharmCAT enforces a strict diplotype grammar (allele1/allele2) but accepts an
# allele-status form for HLA. Seventeen of the benchmark's 62 distinct states are
# not diplotypes at all: HLA presence/absence, mitochondrial variants, G6PD
# haplotypes and CFTR/RYR1 variant-presence calls. Each translation is declared
# here rather than inferred, and anything not translatable is reported as
# "not expressible in PharmCAT notation" instead of being quietly dropped.
_HLA_STATUS = re.compile(r"^HLA-[AB](\*\d+:\d+)\s+(positive|negative)$", re.IGNORECASE)

# PharmCAT names some alleles differently from the CPIC guideline tables the
# benchmark encodes. Each mapping was established by probing the pinned jar, is
# declared here rather than guessed, and is pinned by tests. Without this layer
# the comparison reports notation mismatches as scientific disagreements, which
# would be a fabricated discordance rate.
ALLELE_TRANSLATIONS: dict[str, dict[str, str]] = {
    # PharmCAT calls the DPYD reference allele "Reference" and uses HGVS names
    # for the variants; "*1" and "HapB3" return Indeterminate.
    "DPYD": {"*1": "Reference", "HapB3": "c.1129-5923C>G", "*2A": "c.1905+1G>A"},
    # PharmCAT requires an explicit copy number; "*1xN" does not parse at all.
    "CYP2D6": {"*1xN": "*1x2"},
}


def translate_alleles(gene: str, core: str) -> str:
    table = ALLELE_TRANSLATIONS.get(gene, {})
    if not table:
        return core
    return "/".join(table.get(a.strip(), a.strip()) for a in core.split("/"))


def to_outside_call(gene: str, diplotype: str) -> tuple[str | None, str]:
    """Translate a benchmark state to PharmCAT notation.

    Returns (call, route) where route is "diplotype", "hla_status" or
    "not_expressible". The route is reported, so comparator coverage is a
    stated number rather than an implied 100%.
    """
    d = (diplotype or "").strip()
    m = _HLA_STATUS.match(d)
    if m:
        return f"{m.group(1)} {m.group(2).lower()}", "hla_status"
    core = re.sub(r"\(.*?\)", "", d).strip()
    if core.count("/") == 1 and all(part.strip() for part in core.split("/")):
        return translate_alleles(gene, core), "diplotype"
    return None, "not_expressible"


def batch_states(pairs: list[tuple[str, str]]) -> list[list[tuple[str, str]]]:
    """Group states so no gene appears twice in one PharmCAT run.

    A single phenotyper run reports every gene once, so submitting two CYP2D6
    diplotypes together would silently lose one. Batching keeps each run
    unambiguous.
    """
    batches: list[list[tuple[str, str]]] = []
    seen: list[set[str]] = []
    for gene, call in pairs:
        for i, used in enumerate(seen):
            if gene not in used:
                batches[i].append((gene, call))
                used.add(gene)
                break
        else:
            batches.append([(gene, call)])
            seen.append({gene})
    return batches


def normalise_phenotype(text: str) -> str:
    """Compare meaning, not orthography.

    PharmCAT writes "Intermediate Metabolizer"; the benchmark writes
    "Intermediate Metaboliser". Scoring that as a disagreement would invent a
    discordance rate out of a spelling convention, so normalisation is explicit
    here and pinned by tests rather than left to chance.
    """
    t = (text or "").strip().lower()
    # PharmCAT reports uncalled genes as "No Result" and unscored ones as "n/a".
    # Treating either as a phenotype would fabricate a disagreement.
    if t in {"no result", "n/a", "unknown", ""}:
        return ""
    t = t.replace("metabolizer", "metaboliser")
    # HLA: PharmCAT says "*57:01 positive"; the benchmark truth says "Positive".
    m = re.match(r"^\*\d+:\d+\s+(positive|negative)$", t)
    if m:
        return m.group(1)
    t = t.replace("ultrarapid", "ultra-rapid").replace("ultra rapid", "ultra-rapid")
    # "Normal Metaboliser (Expressor)" and "Normal Metabolizer" are the same
    # CPIC phenotype; the qualifier is descriptive, not a different call.
    t = re.sub(r"\s*\(.*?\)", "", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def write_outside_calls(pairs: list[tuple[str, str]], path: Path) -> int:
    """Write PharmCAT's outside-call TSV. Returns the number of rows written."""
    lines = [f"{gene}\t{diplotype}" for gene, diplotype in pairs]
    path.write_text("\n".join(lines) + "\n")
    return len(lines)


def parse_pharmcat_json(payload: dict) -> dict[str, list[str]]:
    """Extract {gene: [phenotype, ...]} from a phenotyper JSON report."""
    out: dict[str, list[str]] = {}
    for gene, report in (payload.get("geneReports") or {}).items():
        if not isinstance(report, dict):
            continue
        phenotypes: list[str] = []
        for dip in report.get("sourceDiplotypes") or []:
            phenotypes.extend(dip.get("phenotypes") or [])
        if phenotypes:
            out[gene] = phenotypes
    return out


def run_pharmcat_batch(pairs: list[tuple[str, str]], workdir: Path,
                       tag: str = "batch") -> dict[str, list[str]]:
    """Run the pinned PharmCAT phenotyper on one gene-unique batch."""
    if not JAR.exists():
        raise FileNotFoundError(
            f"PharmCAT jar not found at {JAR}. Fetch it with:\n  curl -sL -o {JAR} {JAR_URL}")
    workdir.mkdir(parents=True, exist_ok=True)
    tsv = workdir / f"{tag}.tsv"
    write_outside_calls(pairs, tsv)
    subprocess.run(
        ["java", "-jar", str(JAR), "-phenotyper", "-po", str(tsv),
         "-o", str(workdir), "-bf", tag],
        check=True, capture_output=True, text=True, timeout=900)
    return parse_pharmcat_json(json.loads((workdir / f"{tag}.phenotype.json").read_text()))


def hla_default_statuses(workdir: Path) -> dict[tuple[str, str], str]:
    """PharmCAT's status for HLA loci that were NOT asserted present.

    PharmCAT's outside-call semantics for HLA are "list the alleles present":
    submitting "*57:01 negative" is read as asserting *57:01, and the report
    comes back positive. A negative therefore cannot be stated directly. What
    PharmCAT does provide is its default for any locus not asserted, and it
    reports every locus of the gene in each run. So we assert one locus and read
    the others' defaults, which is PharmCAT's own answer for absence rather than
    an assumption of ours.
    """
    out: dict[tuple[str, str], str] = {}
    probes = [("HLA-B", "*57:01 positive"), ("HLA-B", "*15:02 positive"),
              ("HLA-A", "*31:01 positive")]
    for i, (gene, call) in enumerate(probes):
        by_gene = run_pharmcat_batch([(gene, call)], workdir, tag=f"hla{i}")
        asserted = call.split()[0]
        for status in by_gene.get(gene, []):
            allele = status.split()[0]
            if allele != asserted and status.strip().lower().endswith("negative"):
                out[(gene, f"{allele} negative")] = status
    return out


def phenotypes_for_states(pairs: list[tuple[str, str]],
                          workdir: Path | None = None) -> dict[tuple[str, str], str]:
    """{(gene, outside_call): phenotype} for every state, run in gene-safe batches."""
    tmp = Path(workdir or tempfile.mkdtemp(prefix="pharmcat-"))
    out: dict[tuple[str, str], str] = dict(hla_default_statuses(tmp))
    # HLA negatives are never submitted: PharmCAT reads any asserted allele as
    # present, so submitting "*57:01 negative" returns "*57:01 positive". Their
    # only valid source is hla_default_statuses() above.
    pairs = [(g, c) for g, c in pairs
             if not (g.startswith("HLA-") and c.strip().lower().endswith("negative"))]
    for i, batch in enumerate(batch_states(pairs)):
        by_gene = run_pharmcat_batch(batch, tmp, tag=f"batch{i:02d}")
        for gene, call in batch:
            phenos = [p for p in by_gene.get(gene, []) if normalise_phenotype(p)]
            if gene.startswith("HLA-"):
                # PharmCAT returns every locus status; keep the queried one.
                allele = call.split()[0]
                phenos = [p for p in phenos if p.startswith(allele)] or phenos
            if phenos and (gene, call) not in out:
                out[(gene, call)] = phenos[0]
    return out


def compare(cases: list[dict], state_phenotypes: dict[tuple[str, str], str]) -> list[dict]:
    # An HLA negative is only obtainable when the gene carries another locus we
    # can assert; HLA-A*31:01 is the only HLA-A locus in the benchmark, so its
    # negative cannot be expressed at all. That is a coverage limit of the
    # comparator, and it is reported as such rather than as a disagreement.
    """Three-way comparison per case: ours, PharmCAT, ground truth."""
    rows = []
    for c in cases:
        ours, _ = execute_skill(c["gene"], c["drug"], c["gt_diplotype"])
        pc_gene = pharmcat_gene(c["gene"])
        call, route = to_outside_call(c["gene"], c["gt_diplotype"])
        pc = state_phenotypes.get((pc_gene, call)) if call else None
        if route == "hla_status" and call.endswith("negative") and pc is None:
            route = "not_expressible_hla_negative"
        n_ours, n_pc, n_truth = (normalise_phenotype(ours),
                                 normalise_phenotype(pc),
                                 normalise_phenotype(c["gt_phenotype"]))
        # PharmCAT returning "Indeterminate" is a refusal to call, not a rival
        # call. Counting it as disagreement would report an abstention as a
        # contradiction and understate concordance; it is the same distinction
        # the paper draws for its own executed pipeline.
        pc_abstained = n_pc == "indeterminate"
        rows.append({
            "case_id": c["id"], "gene": c["gene"], "pharmcat_gene": pc_gene,
            "diplotype": c["gt_diplotype"],
            "outside_call": call, "translation_route": route,
            "ours": ours, "pharmcat": pc, "truth": c["gt_phenotype"],
            "pharmcat_available": pc is not None,
            "pharmcat_abstained": pc_abstained,
            "ours_vs_pharmcat": (None if (pc is None or pc_abstained) else n_ours == n_pc),
            "ours_vs_truth": n_ours == n_truth,
            "pharmcat_vs_truth": (None if pc is None else n_pc == n_truth),
        })
    return rows


def summarise(rows: list[dict]) -> dict:
    covered = [r for r in rows if r["pharmcat_available"]]
    abstained = [r for r in covered if r["pharmcat_abstained"]]
    both_emitted = [r for r in covered if not r["pharmcat_abstained"]]
    agree = [r for r in both_emitted if r["ours_vs_pharmcat"]]
    return {
        "n_cases": len(rows),
        "pharmcat_covered": len(covered),
        "pharmcat_not_covered": len(rows) - len(covered),
        "pharmcat_abstained": len(abstained),
        "both_emitted": len(both_emitted),
        "ours_vs_pharmcat_agree": len(agree),
        "ours_vs_pharmcat_disagree": len(both_emitted) - len(agree),
        "concordance_among_both_emitted": round(len(agree) / len(both_emitted), 4) if both_emitted else None,
        "ours_vs_truth_correct": sum(1 for r in rows if r["ours_vs_truth"]),
        "pharmcat_vs_truth_correct": sum(1 for r in both_emitted if r["pharmcat_vs_truth"]),
        "disagreements": [
            {k: r[k] for k in ("case_id", "gene", "diplotype", "ours", "pharmcat", "truth")}
            for r in both_emitted if not r["ours_vs_pharmcat"]
        ],
        "abstention_states": sorted({
            (r["gene"], r["diplotype"]) for r in abstained}),
        "genes_not_covered": sorted({r["gene"] for r in rows if not r["pharmcat_available"]}),
        "translation_routes": dict(Counter(r["translation_route"] for r in rows)),
        "not_expressible": sorted({
            (r["gene"], r["diplotype"]) for r in rows
            if r["translation_route"] == "not_expressible"}),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="Print the download command and exit")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    if args.fetch:
        print(f"mkdir -p {JAR.parent}")
        print(f"curl -sL -o {JAR} {JAR_URL}")
        print(f"shasum -a 256 {JAR}   # expect {JAR_SHA256}")
        return 0

    if not shutil.which("java"):
        sys.stderr.write("java not found; PharmCAT needs a JRE (17+)\n")
        return 1

    # One outside call per distinct (gene, diplotype) in the benchmark.
    states = sorted({(pharmcat_gene(c["gene"]),
                      to_outside_call(c["gene"], c["gt_diplotype"])[0] or "",
                      to_outside_call(c["gene"], c["gt_diplotype"])[1])
                     for c in CASES})
    translatable = [(g, call) for g, call, route in states if call]
    skipped = [(g, route) for g, call, route in states if not call]
    print(f"PharmCAT 3.4.0: {len(translatable)} of {len(states)} distinct states are "
          f"expressible in PharmCAT notation ({len(skipped)} are not)")
    state_phenotypes = phenotypes_for_states(translatable)

    rows = compare(CASES, state_phenotypes)
    summary = summarise(rows)
    args.out.write_text(json.dumps({"summary": summary, "rows": rows}, indent=2))

    lines = [
        "PHARMCAT COMPARISON (independent CPIC implementation, v3.4.0)",
        "",
        f"  cases                          {summary['n_cases']}",
        f"  covered by PharmCAT            {summary['pharmcat_covered']}",
        f"  not covered                    {summary['pharmcat_not_covered']}  "
        f"{summary['genes_not_covered']}",
        f"  PharmCAT abstained             {summary['pharmcat_abstained']}  (Indeterminate)",
        f"  both engines emitted           {summary['both_emitted']}",
        f"  ours agrees with PharmCAT      {summary['ours_vs_pharmcat_agree']}",
        f"  ours disagrees with PharmCAT   {summary['ours_vs_pharmcat_disagree']}",
        f"  concordance among both-emitted {summary['concordance_among_both_emitted']}",
        f"  ours correct vs ground truth   {summary['ours_vs_truth_correct']}/{summary['n_cases']}",
        f"  PharmCAT correct vs truth      {summary['pharmcat_vs_truth_correct']}/{summary['pharmcat_covered']}",
        "",
    ]
    if summary["disagreements"]:
        lines.append("DISAGREEMENTS (reported, not silenced):")
        for d in summary["disagreements"]:
            lines.append(f"  {d['gene']:<14} {d['diplotype']:<22} "
                         f"ours={d['ours']!r} pharmcat={d['pharmcat']!r} truth={d['truth']!r}")
        by_gene_counts = Counter(d["gene"] for d in summary["disagreements"])
        lines.append(f"  by gene: {dict(by_gene_counts)}")
    report = "\n".join(lines)
    REPORT.write_text(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
