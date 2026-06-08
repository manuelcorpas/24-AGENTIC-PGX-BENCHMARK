#!/usr/bin/env bash
# Step 2 — call PGx star-allele diplotypes from a chip/array VCF using PyPGx,
# one gene at a time, on GRCh37 (no liftover). PyPGx is the SINGLE caller used
# for every cohort (UGR, Peru, Corpas) so cross-population comparison is fair.
#
# Requires:
#   - PyPGx (pip install pypgx==0.26.0)
#   - the PyPGx reference bundle: git clone https://github.com/sbslee/pypgx-bundle
#     and set PYPGX_BUNDLE or place at ~/pypgx-bundle (GRCh37 1KGP phasing panels)
#   - Java (>= 17) on PATH (Beagle phasing); Beagle jar ships with PyPGx
#   - the input VCF bgzipped + tabix-indexed (step 1 does this)
#
# Usage: 02_call_diplotypes_pypgx.sh <variants.vcf.gz> <out_dir> [assembly]
set -euo pipefail
VCF="${1:?variants.vcf.gz required}"
OUTDIR="${2:?output dir required}"
ASM="${3:-GRCh37}"
GENES="$(cd "$(dirname "$0")/.." && pwd)/config/genes.txt"

# ensure Java is reachable (keg-only JDKs are not on the default PATH)
command -v java >/dev/null 2>&1 || {
  for J in /opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home/bin \
           /usr/lib/jvm/*/bin /Library/Java/JavaVirtualMachines/*/Contents/Home/bin; do
    [ -x "$J/java" ] && export PATH="$J:$PATH" && break
  done
}
java -version >/dev/null 2>&1 || { echo "ERROR: java not found on PATH"; exit 1; }

mkdir -p "$OUTDIR"
grep -v '^#' "$GENES" | sed '/^[[:space:]]*$/d' | while read -r GENE; do
  echo ">>> $GENE"
  rm -rf "${OUTDIR:?}/$GENE"
  pypgx run-chip-pipeline "$GENE" "$OUTDIR/$GENE" "$VCF" --assembly "$ASM" \
    || echo "  (PyPGx could not call $GENE on this dataset; skipped)"
done
echo "diplotype calls in $OUTDIR/<GENE>/results.zip"
