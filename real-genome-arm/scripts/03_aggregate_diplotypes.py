#!/usr/bin/env python3
"""
Step 3 — aggregate PyPGx per-gene calls into a tidy table of REAL diplotypes
observed in a cohort. Reads each <gene>/results.zip produced by step 2 and emits
one row per (gene, diplotype) with the PyPGx phenotype and the number of carriers.

The (gene, diplotype) -> phenotype pairs are the real-genome test cases: the agent
must reproduce the phenotype for each diplotype actually observed in real people.

Usage: 03_aggregate_diplotypes.py <pypgx_out_dir> <cohort_label> <out.tsv>
"""
import sys
import io
import zipfile
import glob
import os
import csv
from collections import defaultdict

out_dir, cohort, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []  # (cohort, gene, diplotype, phenotype, n_carriers)
for gene_dir in sorted(glob.glob(os.path.join(out_dir, "*"))):
    gene = os.path.basename(gene_dir)
    zpath = os.path.join(gene_dir, "results.zip")
    if not os.path.exists(zpath):
        continue
    zf = zipfile.ZipFile(zpath)
    tsv = [n for n in zf.namelist() if n.endswith("data.tsv")]
    if not tsv:
        continue
    counts = defaultdict(int); phen = {}
    reader = csv.DictReader(io.TextIOWrapper(io.BytesIO(zf.read(tsv[0])), encoding="utf-8"), delimiter="\t")
    for r in reader:
        dip = (r.get("Genotype") or "").strip()
        ph = (r.get("Phenotype") or "").strip()
        if not dip:
            continue
        counts[dip] += 1
        phen[dip] = ph
    for dip, n in counts.items():
        rows.append((cohort, gene, dip, phen.get(dip, ""), n))

with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["cohort", "gene", "diplotype", "phenotype", "n_carriers"])
    w.writerows(sorted(rows))
print(f"{cohort}: {len(rows)} unique (gene, diplotype) observed -> {out_tsv}")
