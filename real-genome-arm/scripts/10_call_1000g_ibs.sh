#!/usr/bin/env bash
# Call PGx diplotypes for the 1000 Genomes IBS cohort (revision item N5).
#
# WHY THIS COHORT
# R2.4 objects that the European point in the cross-cohort comparison is a
# five-member family, which cannot carry the weight placed on it. This adds 93
# 1000 Genomes IBS (Iberian, Spain) samples: per-sample GATK HaplotypeCaller
# VCFs on GRCh37 with AD/DP retained, already held locally. That is 19x the
# family and, unlike a genotype panel, carries the read-depth information
# copy-number calling needs.
#
# LABEL IT HONESTLY. This is IBS, one Iberian population, not the 503-sample EUR
# superpopulation. The manuscript must say IBS. Anything broader would be the
# same over-claim R2.4 is objecting to.
#
# ASSEMBLY. GRCh37 throughout, matching the family, Peru and Uganda arms and the
# single-caller requirement in SOFTWARE.md. No liftover: a build difference in
# one arm would make every cross-cohort coverage number uninterpretable, because
# an allele-table gap could not be told apart from a coordinate mismatch.
#
# CHROMOSOME NAMING. The source VCFs use chr-prefixed contigs; PyPGx requires
# numeric names. The rename happens here, in the pipeline, not by hand.
#
# Requires: bcftools, and the pinned PyPGx env (see SOFTWARE.md); the repo venv
# at .venv-pypgx satisfies it. Also ~/pypgx-bundle (GRCh37 1KGP panels).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARM="$(dirname "$HERE")"
REPO="$(dirname "$ARM")"

SRC="${1:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/ARCHIVE/CoV-MadrID/DATA/1000G/VCF}"
WORK="${2:-$ARM/work/1000g-ibs}"
PYPGX="${PYPGX:-$REPO/.venv-pypgx/bin/pypgx}"
BED="$ARM/config/pgx_regions_grch37.bed"
GENES="$ARM/config/genes.txt"

[ -x "$PYPGX" ] || { echo "pypgx not found at $PYPGX" >&2; exit 1; }
[ -d "${PYPGX_BUNDLE:-$HOME/pypgx-bundle}" ] || { echo "pypgx-bundle missing" >&2; exit 1; }
[ -f "$BED" ] || { echo "missing $BED" >&2; exit 1; }

mkdir -p "$WORK/samples" "$WORK/calls"
cd "$WORK"

echo "[1/4] PGx regions with chr-prefixed names, and a contig rename map"
awk '{print "chr"$1"\t"$2"\t"$3}' "$BED" > pgx_chr.bed
for c in $(seq 1 22); do printf "chr%s\t%s\n" "$c" "$c"; done > rename.txt
printf "chrX\tX\nchrY\tY\nchrM\tMT\n" >> rename.txt

echo "[2/4] subset each sample to the PGx regions and renumber contigs"
n=0
for f in "$SRC"/*.raw.vcf; do
  s="$(basename "$f" .raw.vcf)"
  out="samples/${s}.vcf.gz"
  if [ ! -f "$out" ]; then
    bcftools view -T pgx_chr.bed -Oz -o "samples/.tmp.${s}.vcf.gz" "$f"
    bcftools annotate --rename-chrs rename.txt -Oz -o "$out" "samples/.tmp.${s}.vcf.gz"
    rm -f "samples/.tmp.${s}.vcf.gz"
    bcftools index -t "$out"
  fi
  n=$((n+1))
done
echo "      prepared $n samples"

echo "[3/4] merge into one cohort VCF"
# -0 so a site absent in a sample is reference, not missing: these are per-sample
# call sets, and treating absence as no-call would inflate the abstention rate
# for reasons that have nothing to do with the allele tables.
if [ ! -f cohort.vcf.gz ]; then
  ls samples/*.vcf.gz > merge_list.txt
  bcftools merge -0 -l merge_list.txt -Oz -o merged_raw.vcf.gz
  # Strip the per-sample depth fields. Merging leaves AD as "." wherever a
  # sample had no call at a site, and PyPGx's VCF reader does int(".") on it and
  # dies. Stripping is the honest fix rather than filling in a fabricated depth:
  # it also makes this cohort's input format identical to the Peru and Uganda
  # chip data, which carry no depth fields at all, so every cohort reaches the
  # caller in the same shape. Genotypes, the only thing the chip pipeline uses,
  # are untouched.
  bcftools annotate -x FORMAT/AD,FORMAT/PL,FORMAT/GQ,FORMAT/DP \
      -Oz -o cohort.vcf.gz merged_raw.vcf.gz
  bcftools index -t cohort.vcf.gz
  rm -f merged_raw.vcf.gz
fi
echo "      cohort: $(bcftools query -l cohort.vcf.gz | wc -l | tr -d ' ') samples, \
$(bcftools view -H cohort.vcf.gz | wc -l | tr -d ' ') sites"

# run-chip-pipeline, NOT run-ngs-pipeline. Two reasons, and the second is the
# important one. First, bcftools merge -0 leaves AD as "." at sites a sample
# lacks, and the NGS pipeline's dependency does int(".") on it and dies; that is
# a merge artefact, not a gene that cannot be called. Second and decisive: the
# Peru and Uganda cohorts were called with the chip pipeline, so using it here
# keeps a single calling mode across every cohort. A per-cohort difference in
# calling mode would make cross-cohort coverage differences uninterpretable,
# which is exactly the confound R2.4 objects to.
echo "[4/4] PyPGx per gene (GRCh37, chip pipeline as for Peru and Uganda)"
grep -v '^#' "$GENES" | tr ' ' '\n' | grep -v '^$' | while read -r gene; do
  [ -d "calls/$gene" ] && { echo "      $gene already called, skipping"; continue; }
  echo "      $gene"
  if "$PYPGX" run-chip-pipeline "$gene" "calls/$gene" cohort.vcf.gz \
        --assembly GRCh37 > "calls/${gene}.log" 2>&1; then
    echo "        ok"
  else
    # A gene that cannot be called is a result, not a crash: record and continue.
    echo "        FAILED (see calls/${gene}.log)"
  fi
done

echo
echo "Done. Aggregate with:"
echo "  python $ARM/scripts/03_aggregate_diplotypes.py --input $WORK/calls --cohort 1000G_IBS"
