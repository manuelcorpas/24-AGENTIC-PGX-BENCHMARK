#!/usr/bin/env bash
# Step 16 — call PGx star-allele diplotypes for the Corpas family from WGS.
#
# Replaces the 23andMe SNP-chip arm (step 00 + step 02) for the family cohort.
# Four members have a GRCh37 Sentieon build; the aunt (PT00010A) is GRCh38-only
# at 4.4x and is excluded, so star-allele calling is not attempted on her.
#
# TWO ARMS ARE CALLED ON THE SAME INPUT, deliberately:
#   ngs   — run-ngs-pipeline, the mode a production WGS deployment would use.
#   chip  — run-chip-pipeline, the mode used for the Iberian, Peruvian and
#           Ugandan cohorts. Called here only so the family stays like-for-like
#           with them, making the assay effect a measured quantity rather than
#           a caveat in prose.
#
# Input is each member's OWN single-sample VCF, not a merged cohort VCF.
# In a single-sample VCF the absence of a site already means reference, whereas
# `bcftools merge` writes those sites as missing and a `-0` fill then produces
# AD='.' which PyPGx cannot parse. Per-sample input avoids inventing genotypes.
#
# Requires:
#   - the pinned caller env: conda create -n pypgx python=3.10 && pip install -r
#     real-genome-arm/requirements.txt  (PyPGx 0.26.0, pandas 2.0.3, numpy 1.26.4)
#   - the pypgx-bundle at $PYPGX_BUNDLE or ~/pypgx-bundle
#   - Java 17 on PATH (Beagle phasing)
#   - per-member PGx-region VCFs produced by 15_extract_family_wgs_regions.sh
#
# Usage: 16_call_family_wgs.sh <vcf_dir> <out_dir>
set -euo pipefail
VCFDIR="${1:?directory holding <SAMPLE>.pgx.vcf.gz required}"
OUTDIR="${2:?output dir required}"
GENES="$(cd "$(dirname "$0")/.." && pwd)/config/genes.txt"
CONDA="${CONDA_BIN:-conda}"
ENVNAME="${PYPGX_ENV:-pypgx}"

# keg-only JDKs are not on the default PATH, and launchd/non-interactive shells
# inherit an even smaller one. Probe rather than assume.
command -v java >/dev/null 2>&1 || {
  for J in /opt/homebrew/opt/openjdk@17/bin \
           /opt/homebrew/opt/openjdk/bin \
           /usr/lib/jvm/*/bin /Library/Java/JavaVirtualMachines/*/Contents/Home/bin; do
    [ -x "$J/java" ] && export PATH="$J:$PATH" && break
  done
}
java -version >/dev/null 2>&1 || { echo "ERROR: java not found on PATH"; exit 1; }

SAMPLES="$(cd "$VCFDIR" && ls *.pgx.vcf.gz 2>/dev/null | sed 's/\.pgx\.vcf\.gz$//')"
[ -n "$SAMPLES" ] || { echo "ERROR: no <SAMPLE>.pgx.vcf.gz in $VCFDIR"; exit 1; }

for ARM in ngs chip; do
  for S in $SAMPLES; do
    VCF="$VCFDIR/$S.pgx.vcf.gz"
    grep -v '^#' "$GENES" | sed '/^[[:space:]]*$/d' | while read -r GENE; do
      DEST="$OUTDIR/$ARM/$S/$GENE"
      rm -rf "$DEST"; mkdir -p "$(dirname "$DEST")"
      if [ "$ARM" = "ngs" ]; then
        "$CONDA" run -n "$ENVNAME" pypgx run-ngs-pipeline "$GENE" "$DEST" \
          --variants "$VCF" --assembly GRCh37 --platform WGS \
          >"$DEST.log" 2>&1 || echo "  SKIP $ARM/$S/$GENE (see $DEST.log)"
      else
        "$CONDA" run -n "$ENVNAME" pypgx run-chip-pipeline "$GENE" "$DEST" "$VCF" \
          --assembly GRCh37 \
          >"$DEST.log" 2>&1 || echo "  SKIP $ARM/$S/$GENE (see $DEST.log)"
      fi
    done
    echo "done: $ARM/$S"
  done
done
echo "calls under $OUTDIR/<arm>/<sample>/<gene>/results.zip"
