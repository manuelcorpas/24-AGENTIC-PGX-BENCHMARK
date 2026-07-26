#!/usr/bin/env bash
# Call PGx diplotypes for the 1000 Genomes EUR superpopulation (revision item N5).
#
# WHY: the European point in the cross-cohort comparison is currently a
# five-member family, which cannot carry the weight R2.4 puts on it. 1000
# Genomes EUR (503 samples, permissively licensed, GRCh37, callable by the
# pinned PyPGx) replaces it with an adequately sized cohort, and its Coriell
# identifiers join to the GeT-RM consensus truth set for the calling-step
# validation.
#
# This is a long compute job (tens of GB of downloads), so it is a script you
# run deliberately, not something imported or triggered by a test.
#
# Requires: the pinned PyPGx environment (see SOFTWARE.md), bcftools, samtools.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARM="$(dirname "$HERE")"
WORK="${1:-$ARM/work/1000g-eur}"
PANEL_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
VCF_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

mkdir -p "$WORK"
cd "$WORK"

echo "[1/4] sample panel -> EUR sample list"
[ -f panel.txt ] || curl -sL "$PANEL_URL" -o panel.txt
awk '$3=="EUR" {print $1}' panel.txt > eur_samples.txt
echo "      EUR samples: $(wc -l < eur_samples.txt)"

echo "[2/4] per-chromosome PGx region subsets (BED from the repo config)"
BED="$ARM/config/pgx_regions_grch37.bed"
[ -f "$BED" ] || { echo "missing $BED" >&2; exit 1; }
cut -f1 "$BED" | sort -u | while read -r chr; do
  out="chr${chr}.pgx.vcf.gz"
  if [ ! -f "$out" ]; then
    echo "      chr${chr}"
    bcftools view \
      -R "$BED" \
      -S eur_samples.txt --force-samples \
      -Oz -o "$out" \
      "${VCF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    bcftools index -t "$out"
  fi
done

echo "[3/4] merge"
[ -f eur.pgx.vcf.gz ] || {
  bcftools concat -Oz -o eur.pgx.vcf.gz chr*.pgx.vcf.gz
  bcftools index -t eur.pgx.vcf.gz
}

echo "[4/4] PyPGx calls (same genes and pinned version as the other cohorts)"
while read -r gene; do
  [ -z "$gene" ] && continue
  [ -d "pypgx-$gene" ] || pypgx run-ngs-pipeline "$gene" "pypgx-$gene" \
      --variants eur.pgx.vcf.gz --assembly GRCh37
done < "$ARM/config/genes.txt"

echo
echo "Done. Aggregate with:"
echo "  python $ARM/scripts/03_aggregate_diplotypes.py --input $WORK --cohort 1000G_EUR"
echo "Then validate the calling step against GeT-RM with:"
echo "  python $ARM/scripts/08_caller_truth_eval.py --truth $ARM/getrm/getrm_consensus.tsv --called <aggregated tsv>"
