#!/usr/bin/env bash
# Call PGx diplotypes for the GeT-RM / 1000 Genomes overlap (revision item N5b).
#
# WHY THIS EXISTS
# N0 validated the executed skill's MAPPING step against PyPGx. It did not
# validate the CALLING step, because PyPGx both called the diplotype and supplied
# the phenotype it was scored against, so input-call error was unobservable
# rather than absent. The plan of revision says so to the editor in those terms.
#
# The GeT-RM Consolidated PGx and HLA Table supplies external consensus
# genotypes. 138 of its samples are in 1000 Genomes phase 3 and are called here;
# 113 of those 138 yield at least one (sample, gene) pair the consensus table and
# the caller both report, giving 527 such pairs across the ten genes this arm also
# calls. The remaining 25 samples contribute no evaluable pair: GeT-RM
# characterises every one of them for CYP3A4 alone, and this arm calls CYP3A5,
# not CYP3A4. The exclusion is therefore a vocabulary boundary, not a call
# failure, and it is not a filter applied after seeing concordance. That is enough to
# make "100% correctness among emitted answers" genuinely falsifiable: a call
# that disagrees with the GeT-RM consensus is an end-to-end error attributable
# to input interpretation.
#
# WHY REMOTE AND NOT A LOCAL ARCHIVE
# The local 1000G archive used by the IBS arm holds 93 Iberian samples, of which
# exactly 2 are in GeT-RM. The phase 3 release is public and tabix-indexed, so
# the region slices for these samples are fetched over HTTPS. Nothing here needs
# a human in the loop, which is the point: the truth set is an independent
# software oracle, not an expert opinion.
#
# WHAT THIS DOES AND DOES NOT MEASURE. Read before quoting a number.
# 1000G phase 3 genotypes come from low-coverage sequencing with imputation and
# carry no per-sample read depth. This therefore evaluates PyPGx on phase 3
# genotype data, which is the same calling mode (run-chip-pipeline) already used
# for the Peru, Uganda and IBS arms, so cross-arm comparisons stay meaningful.
# It does NOT evaluate PyPGx on high-depth sequencing, and CYP2D6 structural
# variation cannot be resolved without depth. Both limits are reported, not
# worked around.
#
# ASSEMBLY. GRCh37 throughout, matching every other arm. No liftover.
#
# Requires: bcftools (1.21 used), and the pinned PyPGx env (SOFTWARE.md); the
# repo venv at .venv-pypgx satisfies it. Also ~/pypgx-bundle (GRCh37 panels).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARM="$(dirname "$HERE")"
REPO="$(dirname "$ARM")"

WORK="${1:-$ARM/work/getrm-1000g}"
PYPGX="${PYPGX:-$REPO/.venv-pypgx/bin/pypgx}"
BED="$ARM/config/pgx_regions_grch37.bed"
GENES="$ARM/config/genes.txt"
TRUTH="$ARM/getrm/getrm_consensus.tsv"

KGP_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
KGP_TPL="ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

[ -x "$PYPGX" ] || { echo "pypgx not found at $PYPGX" >&2; exit 1; }
[ -d "${PYPGX_BUNDLE:-$HOME/pypgx-bundle}" ] || { echo "pypgx-bundle missing" >&2; exit 1; }
[ -f "$BED" ] || { echo "missing $BED" >&2; exit 1; }
[ -f "$TRUTH" ] || { echo "missing truth set $TRUTH; see $ARM/getrm/README.md" >&2; exit 1; }

mkdir -p "$WORK/regions" "$WORK/calls"
cd "$WORK"

echo "[1/5] sample list: GeT-RM samples present in 1000G phase 3"
if [ ! -f samples.txt ]; then
  # The intersection is computed here rather than hard-coded, so the list moves
  # if the truth set is updated.
  cut -f1 "$TRUTH" | tail -n +2 | sort -u > getrm_all.txt
  bcftools query -l "${KGP_BASE}/$(printf "$KGP_TPL" 22)" | sort -u > kgp_all.txt
  comm -12 getrm_all.txt kgp_all.txt > samples.txt
fi
echo "      $(wc -l < samples.txt | tr -d ' ') samples"

echo "[2/5] fetch PGx region slices per chromosome, subset to those samples"
# One remote query per chromosome carrying all of that chromosome's PGx regions.
# tabix fetches only the indexed byte ranges, so this is a few tens of MB, not
# the whole release.
for chrom in $(cut -f1 "$BED" | sort -u -n); do
  out="regions/chr${chrom}.vcf.gz"
  [ -f "$out" ] && { echo "      chr${chrom} present, skipping"; continue; }
  regions="$(awk -v c="$chrom" '$1==c {printf "%s%s:%s-%s", (n++?",":""), $1, $2+1, $3}' "$BED")"
  echo "      chr${chrom}: $regions"
  bcftools view -r "$regions" -S samples.txt --force-samples \
      -Oz -o "regions/.tmp.chr${chrom}.vcf.gz" \
      "${KGP_BASE}/$(printf "$KGP_TPL" "$chrom")"
  mv "regions/.tmp.chr${chrom}.vcf.gz" "$out"
  bcftools index -t "$out"
done

echo "[3/5] concatenate into one cohort VCF"
if [ ! -f cohort.vcf.gz ]; then
  ls regions/chr*.vcf.gz | sort -V > concat_list.txt
  bcftools concat -a -Oz -o concat_raw.vcf.gz -f concat_list.txt
  # Phase 3 carries no per-sample AD/DP, but strip defensively so this cohort
  # reaches the caller in exactly the shape the Peru, Uganda and IBS cohorts do.
  # A per-cohort difference in input format would make a caller-accuracy
  # difference uninterpretable.
  bcftools annotate -x FORMAT/AD,FORMAT/PL,FORMAT/GQ,FORMAT/DP \
      -Oz -o cohort.vcf.gz concat_raw.vcf.gz
  bcftools index -t cohort.vcf.gz
  rm -f concat_raw.vcf.gz
fi
echo "      cohort: $(bcftools query -l cohort.vcf.gz | wc -l | tr -d ' ') samples, \
$(bcftools view -H cohort.vcf.gz | wc -l | tr -d ' ') sites"

echo "[4/5] PyPGx per gene (GRCh37, chip pipeline as for every other cohort)"
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

echo "[5/5] next steps"
echo "  python $HERE/13_calls_to_sample_table.py --calls $WORK/calls --out $WORK/pypgx_calls.tsv"
echo "  python $HERE/08_caller_truth_eval.py --truth $TRUTH --called $WORK/pypgx_calls.tsv"
