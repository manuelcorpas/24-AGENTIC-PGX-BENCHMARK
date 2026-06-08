#!/usr/bin/env bash
# Step 1 — extract the PGx-gene regions (GRCh37) from a PLINK1 binary fileset
# (.bed/.bim/.fam) to a bgzipped, tabix-indexed VCF. Reduces a genome-wide set
# (e.g. ~20M variants) to the ~12k variants in the PGx genes.
#
# Used for: Uganda Genome Resource (UGR) and the Peruvian chip set.
# Requires: plink2 (>= v2.0.0-a.7), python3 with pysam.
#
# Usage: 01_extract_pgx_regions.sh <plink_prefix> <out_prefix>
#   <plink_prefix>  path prefix of the .bed/.bim/.fam set (GRCh37, numeric chroms)
#   <out_prefix>    output prefix; writes <out_prefix>.vcf.gz (+ .tbi)
set -euo pipefail
PREFIX="${1:?plink prefix required}"
OUT="${2:?out prefix required}"
BED="$(cd "$(dirname "$0")/.." && pwd)/config/pgx_regions_grch37.bed"

plink2 --bfile "$PREFIX" --extract bed1 "$BED" --export vcf bgz --out "$OUT"
python3 -c "import pysam; pysam.tabix_index('${OUT}.vcf.gz', preset='vcf', force=True)"
echo "wrote ${OUT}.vcf.gz (+ .tbi)"
