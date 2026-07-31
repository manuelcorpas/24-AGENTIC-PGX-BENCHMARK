#!/usr/bin/env python3
"""
Input normalisation: the experiment for reviewer point R1.1.

WHY THIS EXISTS
The paper argues that the model belongs at input normalisation and that
validated code belongs at the mapping. Every experiment in the paper measures
the second half of that claim. None measures the first. The Discussion of v30
says so in as many words, disclaiming the architecture rather than showing it,
and R1.1 is the reviewer objecting that the "agentic" framing exceeds what is
tested. A disclaimer concedes the point permanently; this measures it.

WHAT THE TASK IS
A model is shown the variant-level genotypes a real sample carries across a
gene region, rendered in one of three surface forms a clinical service actually
receives, and asked for the canonical star-allele diplotype. That call is then
handed to the validated skill, which computes phenotype and recommendation in
code. This is the configuration the paper argues for, run end to end.

WHAT IT IS SCORED AGAINST
PyPGx, the deterministic caller, on the same genotypes. Where GeT-RM publishes
a consensus genotype for the same (sample, gene) pair, that external truth is
reported alongside, because the caller is not itself infallible and the paper
already reports its 0.761 raw concordance.

TWO ARMS
The default arm gives the model the genotypes and the gene name and nothing
else. It measures whether a model can produce a diplotype from real input, and
it cannot separate a failure to normalise from a failure to recall allele
definitions from memory. It does not claim to.

`--definitions` adds the allele-definition table PyPGx itself uses, in the same
GRCh37 coordinate space as the genotypes. That arm puts model and caller on the
same information footing and isolates normalisation from recall. It is the arm
that decides whether the paper's architectural proposal survives.

WHAT NEITHER ARM DOES
Neither shows the model a candidate diplotype list or the caller's answer. The
definitions block names single alleles only, never a diplotype, so the answer is
never in the prompt: the model must still work out which alleles this individual
carries and how they pair. The tests enforce this, and the rendered patient
genotypes are asserted to contain no "*" in either arm. Supplying the answer is
what made the earlier extraction claim vacuous.

The renderings are VCF-native, HGVS genomic and a prose report fragment. They
are NOT rsID-keyed: the 1000 Genomes phase 3 slices used here carry no rsIDs in
the ID column, and joining dbSNP to manufacture them would add a dependency
whose failures would be scored as model failures. That limit is reported rather
than papered over.

USAGE
    python code/74-input-normalisation.py --build-inputs
    python code/74-input-normalisation.py --estimate
    python code/74-input-normalisation.py --run --max-spend 40
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
CODE = Path(__file__).resolve().parent
RGA = BASE / "real-genome-arm"
VCF = RGA / "work" / "getrm-1000g" / "cohort.vcf.gz"
CALLS = RGA / "work" / "getrm-1000g" / "pypgx_calls.tsv"
CONSENSUS = RGA / "getrm" / "getrm_consensus.tsv"
REGIONS = RGA / "config" / "pgx_regions_grch37.bed"

INPUTS_OUT = BASE / "data" / "v3_normalisation_inputs.json"
RAW_OUT = BASE / "data" / "v3_input_normalisation.json"
REPORT_OUT = BASE / "data" / "v3_input_normalisation.txt"

FORMS = ("vcf", "hgvs", "prose")
DEFAULT_MODELS = ("GPT-5.2", "Claude Opus 4.5", "o3")

# A cap on how many variants are shown, set ABOVE the observed maximum (681) so
# that in practice nothing is truncated. An earlier value of 60 truncated 318 of
# 527 pairs, which would have withheld the allele-defining variant from the model
# and then scored the resulting failure as a model failure. That is the same
# class of error as the HLA regex in CORRECTIONS.md: a harness limit reported as
# a finding. The cap is kept only as a guard against an unbounded prompt, and the
# truncation count is reported so a reader can see it is zero.
MAX_VARIANTS = 800


def _load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- rendering

def render(variants: list[dict], form: str, truncated: bool = False) -> str:
    """Render genotypes in one surface form. Never emits a star allele.

    `truncated` is not cosmetic. The prose form otherwise closes by asserting
    that no other non-reference genotype was observed, which is a claim about
    the input and is false once the list has been cut. Every form discloses
    truncation so that a model's answer is conditioned on what it was actually
    shown rather than on a completeness the harness invented.
    """
    if form not in FORMS:
        raise ValueError(f"unknown rendering {form!r}; expected one of {FORMS}")
    if not variants:
        raise ValueError(
            "refusing to render an empty variant set: a sample with no non-reference "
            "genotype is a decision to make upstream, not a prompt asking a model to "
            "call a diplotype from nothing")

    note = ("This list was truncated and is not the complete set of "
            "non-reference genotypes in the region.")

    if form == "vcf":
        lines = ["#CHROM\tPOS\tREF\tALT\tGT"]
        lines += [f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gt']}"
                  for v in variants]
        if truncated:
            lines.append(f"# {note}")
        return "\n".join(lines)

    if form == "hgvs":
        out = []
        for v in variants:
            zygosity = "homozygous" if v["gt"] in ("1|1", "1/1") else "heterozygous"
            out.append(f"NC_0000{int(v['chrom']):02d}.10:g.{v['pos']}{v['ref']}>{v['alt']} ({zygosity})")
        if truncated:
            out.append(f"({note})")
        return "\n".join(out)

    # prose: how a report fragment actually reads
    out = ["Sequencing report extract. The following non-reference genotypes were "
           "observed across the region:"]
    for v in variants:
        zygosity = "two copies" if v["gt"] in ("1|1", "1/1") else "one copy"
        out.append(f"  - position {v['pos']} on chromosome {v['chrom']}: reference base "
                   f"{v['ref']}, {zygosity} of the {v['alt']} allele.")
    out.append(note if truncated
               else "No other non-reference genotype was observed in the region examined.")
    return "\n".join(out)


ALLELE_TABLE = (BASE / ".venv-pypgx" / "lib" / "python3.10" / "site-packages" /
                "pypgx" / "api" / "data" / "allele-table.csv")


def load_definitions(gene: str, path: Path = None) -> list[dict]:
    """PyPGx's own allele-definition rows for one gene.

    This is the table the deterministic caller uses, in the same GRCh37
    coordinate space as the rendered genotypes, so supplying it puts the model
    and the caller on the same information footing.
    """
    import csv
    src = path or ALLELE_TABLE
    with src.open() as fh:
        return [r for r in csv.DictReader(fh) if r["Gene"] == gene]


def definition_block(definitions: list[dict]) -> str:
    """Allele -> defining variants. Single alleles only, never a diplotype.

    Naming the alleles is fair; the caller knows them. Naming a diplotype would
    hand over the answer and rebuild the circularity that made the extraction
    claim vacuous, so this renders one allele per line and nothing else.
    """
    lines = ["allele\tfunction\tdefining variants (GRCh37, chrom-pos-ref-alt)"]
    for d in definitions:
        core = d.get("GRCh37Core") or "N/A"
        lines.append(f"{d['StarAllele']}\t{d.get('Function', '')}\t{core}")
    return "\n".join(lines)


def build_prompt(gene: str, variants: list[dict], form: str,
                 truncated: bool = False,
                 definitions: list[dict] | None = None) -> str:
    """The model is given genotypes and the gene, and nothing else.

    With `definitions`, it is additionally given the allele-definition table the
    deterministic caller uses. That arm separates a failure to normalise from a
    failure to recall allele definitions from memory, which the first arm could
    not distinguish and did not claim to.
    """
    defs = ""
    if definitions:
        defs = (f"\nThese are the {gene} allele definitions. An allele is present when the "
                f"individual carries its defining variants.\n\n"
                f"{definition_block(definitions)}\n")
    return f"""You are performing input normalisation for a clinical pharmacogenomics pipeline.

Below are the non-reference genotypes observed for one individual across the {gene} region, on the GRCh37 assembly. Determine the {gene} star-allele diplotype this individual carries.

{render(variants, form, truncated)}
{defs}

Report the diplotype using standard star-allele nomenclature. If these data are not sufficient to determine a diplotype, answer ABSTAIN. An abstention is recorded as such and is not penalised as a wrong answer; a guess is.

## Output (1 line only)
DIPLOTYPE: <diplotype, or ABSTAIN>"""


# ------------------------------------------------------------------ parsing

_DIP = re.compile(r"^\s*DIPLOTYPE:\s*(.+?)\s*$", re.IGNORECASE | re.MULTILINE)
_VALID = re.compile(r"^\*[\w.]+(?:x\d+)?/\*[\w.]+(?:x\d+)?$")


# A gene symbol qualifying a diplotype, with or without a separator. The
# separator is optional because "CYP2D6*4/*4" is the standard PharmVar
# rendering; requiring it cost 744 calls (C14). The symbol itself must look
# like a gene symbol, so a lowercase hedge ("unknown*1/*2") is not stripped
# into a call, which would be C1 in reverse.
_GENE_PREFIX = re.compile(r"^[A-Z][A-Z0-9]{2,9}[- ]*(?=\*)")


def parse_call(text: str | None, gene: str | None = None) -> str | None:
    """Return a diplotype, or None for an abstention, refusal or non-answer.

    Deliberately strict. A response that does not carry a well-formed diplotype
    is an abstention, never a salvage attempt: coverage is the headline of this
    experiment, so a lenient parser would move the number it is measuring.

    `gene` is the gene the row asked about. When supplied, that symbol is
    stripped first and case-insensitively; otherwise any string shaped like a
    gene symbol is stripped. Both paths leave a hedge word in place.
    """
    if not text:
        return None
    m = _DIP.search(text)
    if not m:
        return None
    value = m.group(1).strip()
    if value.upper().startswith("ABSTAIN"):
        return None
    # strip a gene prefix ("CYP2D6 *4/*4" or "CYP2D6*4/*4") if the model
    # supplied one
    if gene:
        value = re.sub(rf"^{re.escape(gene)}[- ]*(?=\*)", "", value,
                       flags=re.IGNORECASE)
    value = _GENE_PREFIX.sub("", value.strip())
    value = value.split()[0] if value.split() else ""
    return value if _VALID.match(value) else None


def classify(text: str | None, out_tokens: int, cap: int = None,
             error: str | None = None, gene: str | None = None) -> str:
    """'call', 'abstain', 'truncated_output' or 'error'.

    Three ways of not getting an answer, which must never be pooled.

    `truncated_output` exists because a model that reasons past its output
    budget never emits the DIPLOTYPE line, and a parser cannot tell that from a
    refusal by looking at the text alone.

    `error` is checked FIRST and overrides everything, including text that looks
    like a valid call. A failed request leaves empty text and zero tokens, which
    is indistinguishable from a terse refusal on the text alone. Recording it as
    an abstention is entry C3 in CORRECTIONS.md, and it recurred here: 94 o3
    calls that returned HTTP 429 on an exhausted quota were first written out as
    94 abstentions, at a total cost of $0.00 that should have given it away.
    """
    if error:
        return "error"
    cap = cap if cap is not None else MAX_OUT_TOKENS
    if parse_call(text, gene) is not None:
        return "call"
    if _DIP.search(text or "") is None and out_tokens >= cap:
        return "truncated_output"
    return "abstain"


# ------------------------------------------------------------------ scoring

def _norm_dip(d: str) -> tuple[str, str]:
    parts = [p.strip() for p in d.split("/")]
    return tuple(sorted(parts))


def score(rows: list[dict]) -> dict:
    """Coverage, abstention and accuracy AMONG EMITTED.

    Accuracy is divided by what was emitted, not by the total. Dividing by the
    total rewards abstention, which is the failure this paper accuses other
    evaluations of.
    """
    for r in rows:
        if "reference" not in r:
            raise ValueError(f"row without a reference cannot be scored: {r}")
    n = len(rows)
    emitted = [r for r in rows if r.get("call")]
    correct = [r for r in emitted
               if r["reference"] and _norm_dip(r["call"]) == _norm_dip(r["reference"])]
    return {
        "n": n,
        "emitted": len(emitted),
        "coverage": len(emitted) / n if n else 0.0,
        "abstention": 1 - (len(emitted) / n) if n else 0.0,
        "correct": len(correct),
        "accuracy_among_emitted": (len(correct) / len(emitted)) if emitted else None,
    }


def stratified_subsample(inputs: dict, n: int, seed: int = 20260730) -> list[str]:
    """Deterministic gene-stratified draw of n keys.

    o3 costs more than the other two models combined because of its reasoning
    tokens, so it runs on a subsample rather than the full set. The draw is
    stratified because the gene counts are very uneven (CYP2D6 72 pairs against
    NUDT15 11), and a proportional draw would drop the small genes; it is seeded
    because a subsample that changes between runs cannot be republished.
    """
    import random
    if n >= len(inputs):
        return sorted(inputs)
    by_gene: dict[str, list[str]] = {}
    for k in sorted(inputs):
        by_gene.setdefault(inputs[k]["gene"], []).append(k)
    rng = random.Random(seed)
    # one from every gene first, so no gene is lost, then proportional
    picked = []
    for g in sorted(by_gene):
        picked.append(rng.choice(by_gene[g]))
    remaining = [k for k in sorted(inputs) if k not in set(picked)]
    rng.shuffle(remaining)
    picked += remaining[: max(0, n - len(picked))]
    return sorted(picked[:n])


def freeze(inputs: dict) -> str:
    """SHA-256 over the exact genotypes shown, so the scored run is provably
    the run whose inputs are published."""
    return hashlib.sha256(
        json.dumps(inputs, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


# ------------------------------------------------------------- input build

def _regions() -> dict[str, tuple[str, int, int]]:
    out = {}
    for line in REGIONS.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, name = line.split()[:4]
        for gene in name.split("_"):
            out[gene] = (chrom, int(start), int(end))
    return out


def build_inputs() -> dict:
    """Extract per (sample, gene) non-reference genotypes for the evaluable pairs.

    Evaluable means both PyPGx and the GeT-RM consensus table report the pair,
    which is the 527-pair set already reported in Results. Using that set means
    every row carries both a deterministic reference and an external truth.
    """
    import csv
    regions = _regions()
    calls = {(r["sample"], r["gene"]): r["diplotype"]
             for r in csv.DictReader(CALLS.open(), delimiter="\t")}
    cons = {(r["sample"], r["gene"]): r["diplotype"]
            for r in csv.DictReader(CONSENSUS.open(), delimiter="\t")}
    pairs = sorted(set(calls) & set(cons))

    inputs = {}
    skipped = {"no_region": 0, "no_variants": 0}
    by_gene: dict[str, list[str]] = {}
    for sample, gene in pairs:
        by_gene.setdefault(gene, []).append(sample)

    for gene, samples in sorted(by_gene.items()):
        if gene not in regions:
            skipped["no_region"] += len(samples)
            continue
        chrom, start, end = regions[gene]
        cmd = ["bcftools", "query", "-r", f"{chrom}:{start}-{end}",
               "-s", ",".join(samples),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(VCF)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"bcftools failed for {gene}: {res.stderr[:400]}")
        per_sample: dict[str, list[dict]] = {s: [] for s in samples}
        for line in res.stdout.splitlines():
            f = line.split("\t")
            chrom_, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
            for s, gt in zip(samples, f[4:]):
                if gt in ("0|0", "0/0", ".|.", "./.", "."):
                    continue
                per_sample[s].append({"chrom": chrom_, "pos": pos, "ref": ref,
                                      "alt": alt, "gt": gt})
        for s in samples:
            vs = per_sample[s][:MAX_VARIANTS]
            if not vs:
                skipped["no_variants"] += 1
                continue
            inputs[f"{s}|{gene}"] = {
                "sample": s, "gene": gene, "variants": vs,
                "n_variants_total": len(per_sample[s]),
                "truncated": len(per_sample[s]) > MAX_VARIANTS,
                "reference": calls[(s, gene)],
                "getrm": cons[(s, gene)],
            }
    payload = {
        "inputs": inputs,
        "n_pairs_evaluable": len(pairs),
        "n_pairs_with_variants": len(inputs),
        "skipped": skipped,
        "max_variants_shown": MAX_VARIANTS,
        "sha256": None,
    }
    payload["sha256"] = freeze(inputs)
    return payload


# ------------------------------------------------------------------ clients

# The matched factorial caps Anthropic output at 320 tokens, which is right for
# its one-line task. Here it is not: normalisation makes the model reason about
# which variants define which haplotype, and a 320-token cap truncated Claude
# mid-reasoning BEFORE it reached the DIPLOTYPE line. The parser then read a
# missing line as an abstention, so a harness cap would have been published as a
# model declining to answer. The cap is raised for every model in this
# experiment, and 60-matched-factorial.py is left untouched so the published
# factorial runs are not disturbed.
MAX_OUT_TOKENS = 6000


def load_models(names: list[str]) -> dict:
    """Clients for this experiment only, with an output budget that fits the task."""
    import anthropic
    import openai

    env_path = BASE / ".env"
    keys = dict(os.environ)
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, _, v = line.partition("=")
                keys.setdefault(k.strip(), v.strip().strip('"').strip("'"))

    def key(*nm):
        for n in nm:
            if keys.get(n):
                return keys[n]
        raise KeyError(f"Missing API key: set one of {nm}")

    built = {}
    if any(n.startswith("Claude") for n in names):
        ant = anthropic.Anthropic(api_key=key("ANTHROPIC_API_KEY"))

        def _ant(model, p):
            r = ant.messages.create(model=model, max_tokens=MAX_OUT_TOKENS,
                                    messages=[{"role": "user", "content": p}])
            text = "".join(b.text for b in r.content if getattr(b, "type", "") == "text")
            return text, r.usage.input_tokens, r.usage.output_tokens

        built["Claude Opus 4.5"] = lambda p: _ant("claude-opus-4-5-20251101", p)
        built["Claude Sonnet 4.5"] = lambda p: _ant("claude-sonnet-4-5-20250929", p)

    if any(n in ("GPT-5.2", "GPT-4.1", "o3", "o4-mini") for n in names):
        oai = openai.OpenAI(api_key=key("OPENAI_API_KEY"))

        def _oai(model, p, reasoning=False):
            kw = ({"max_completion_tokens": MAX_OUT_TOKENS} if reasoning
                  else {"max_tokens": MAX_OUT_TOKENS})
            r = oai.chat.completions.create(
                model=model, messages=[{"role": "user", "content": p}], **kw)
            return (r.choices[0].message.content,
                    r.usage.prompt_tokens, r.usage.completion_tokens)

        built["GPT-5.2"] = lambda p: _oai("gpt-5.2", p, reasoning=True)
        built["GPT-4.1"] = lambda p: _oai("gpt-4.1", p)
        built["o3"] = lambda p: _oai("o3", p, reasoning=True)
        built["o4-mini"] = lambda p: _oai("o4-mini", p, reasoning=True)

    missing = [n for n in names if n not in built]
    if missing:
        raise KeyError(f"no client built for {missing}")
    return {n: built[n] for n in names}


# ----------------------------------------------------------------- running

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build-inputs", action="store_true")
    ap.add_argument("--estimate", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    ap.add_argument("--forms", nargs="+", default=list(FORMS))
    ap.add_argument("--max-spend", type=float, default=40.0)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap pairs, for a smoke run")
    ap.add_argument("--subsample", type=int, default=None,
                    help="gene-stratified deterministic subsample of pairs")
    ap.add_argument("--out", type=Path, default=RAW_OUT)
    ap.add_argument("--max-out-tokens", type=int, default=None,
                    help="override the output budget for this run")
    ap.add_argument("--definitions", action="store_true",
                    help="supply the PyPGx allele-definition table for the gene")
    ap.add_argument("--skip-completed", type=Path, default=None,
                    help="exclude (sample,gene,form,model) units already in a prior output")
    ap.add_argument("--only-truncated", type=Path, default=None,
                    help="rerun only the units a prior run left truncated_output")
    args = ap.parse_args(argv)

    if args.build_inputs:
        payload = build_inputs()
        INPUTS_OUT.write_text(json.dumps(payload, indent=2))
        print(f"evaluable pairs        {payload['n_pairs_evaluable']}")
        print(f"pairs with variants    {payload['n_pairs_with_variants']}")
        print(f"skipped                {payload['skipped']}")
        print(f"sha256                 {payload['sha256']}")
        print(f"written                {INPUTS_OUT}")
        return 0

    if not INPUTS_OUT.exists():
        print("build inputs first: --build-inputs", file=sys.stderr)
        return 2
    payload = json.loads(INPUTS_OUT.read_text())
    inputs = payload["inputs"]
    if args.subsample:
        keys = stratified_subsample(inputs, args.subsample)
    else:
        keys = sorted(inputs)[: args.limit] if args.limit else sorted(inputs)

    runner = _load("matched_factorial", CODE / "60-matched-factorial.py")

    units = [(k, f, m) for k in keys for f in args.forms for m in args.models]
    if args.skip_completed and args.skip_completed.exists():
        prior = json.loads(args.skip_completed.read_text())["rows"]
        done = {(f"{r['sample']}|{r['gene']}", r["form"], r["model"])
                for r in prior if r.get("status") != "error"}
        units = [u for u in units if u not in done]
    if args.only_truncated:
        prior = json.loads(args.only_truncated.read_text())["rows"]
        want = {(f"{r['sample']}|{r['gene']}", r["form"], r["model"])
                for r in prior if r.get("status") == "truncated_output"}
        units = [u for u in units if u in want]

    if args.estimate:
        # price from the real prompts, not from a guess
        import statistics
        _dcache = {}
        def _defs(g):
            if args.definitions and g not in _dcache:
                _dcache[g] = load_definitions(g)
            return _dcache.get(g)
        sample_prompts = [build_prompt(inputs[k]["gene"], inputs[k]["variants"], f, inputs[k]["truncated"], _defs(inputs[k]["gene"]))
                          for k in keys[:40] for f in args.forms]
        mean_chars = statistics.mean(len(p) for p in sample_prompts)
        in_tok = mean_chars / 3.7          # conservative chars-per-token
        total = 0.0
        for m in args.models:
            out_tok = runner.DEFAULT_OUT_TOKENS.get(m, 320)
            per = runner.Spend.cost(m, in_tok, out_tok)
            n = len(keys) * len(args.forms)
            total += per * n
            print(f"  {m:20s} {n:6d} calls  ~{in_tok:6.0f} in  {out_tok:5d} out  ${per*n:8.2f}")
        print(f"\n  pairs {len(keys)}  forms {len(args.forms)}  models {len(args.models)}")
        print(f"  total calls {len(units)}")
        print(f"  PROJECTED   ${total:.2f}")
        print(f"  inputs sha256 {payload['sha256']}")
        return 0

    if not args.run:
        ap.print_help()
        return 1

    global MAX_OUT_TOKENS
    if args.max_out_tokens:
        MAX_OUT_TOKENS = args.max_out_tokens
    models = load_models(args.models)
    spend = runner.Spend(args.max_spend)
    rows = []

    def one(unit):
        key, form, model_name = unit
        rec = inputs[key]
        defs = load_definitions(rec["gene"]) if args.definitions else None
        prompt = build_prompt(rec["gene"], rec["variants"], form, rec["truncated"], defs)
        text, error, in_tok, out_tok = "", None, 0, 0
        for attempt in range(3):
            try:
                spend.check()
                r = models[model_name](prompt)
                text, in_tok, out_tok = r if isinstance(r, tuple) else (r, 0, 0)
                spend.add(model_name, in_tok, out_tok)
                break
            except runner.BudgetExceeded:
                raise
            except Exception as exc:                    # noqa: BLE001
                error = f"{type(exc).__name__}: {exc}"
                if "429" in str(exc) or "rate" in str(exc).lower():
                    import time
                    time.sleep(5 * (attempt + 1))
        return {
            "sample": rec["sample"], "gene": rec["gene"], "form": form,
            "model": model_name, "raw": text, "error": error,
            "call": parse_call(text, rec["gene"]),
            "status": classify(text, out_tok, error=error, gene=rec["gene"]),
            "reference": rec["reference"], "getrm": rec["getrm"],
            "definitions_supplied": bool(args.definitions),
            "n_variants_shown": len(rec["variants"]),
            "truncated": rec["truncated"],
            "in_tokens": in_tok, "out_tokens": out_tok,
            "cost_usd": round(runner.Spend.cost(model_name, in_tok, out_tok), 6),
        }

    try:
        with ThreadPoolExecutor(max_workers=8) as ex:
            for i, row in enumerate(ex.map(one, units), 1):
                rows.append(row)
                if i % 100 == 0:
                    print(f"  {i}/{len(units)}  ${spend.total:.2f}", file=sys.stderr)
    except runner.BudgetExceeded as exc:
        print(f"budget stop: {exc}", file=sys.stderr)

    out = {"rows": rows, "inputs_sha256": payload["sha256"],
           "models": args.models, "forms": args.forms,
           "total_cost_usd": round(sum(r["cost_usd"] for r in rows), 4)}
    out["subsample"] = args.subsample
    out["n_pairs"] = len(keys)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}  ({len(rows)} rows, ${out['total_cost_usd']:.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
