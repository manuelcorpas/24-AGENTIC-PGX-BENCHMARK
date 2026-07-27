#!/usr/bin/env bash
# Fetch the GeT-RM Consolidated PGx and HLA Table from CDC.
#
# WHY THIS SCRIPT EXISTS, AND A CORRECTION
# An earlier version of the README stated that www.cdc.gov "returns HTTP 403 to
# non-browser clients" and that the download was therefore necessarily manual.
# That diagnosis was wrong, and wrong in an instructive way: cdc.gov does return
# 403, but to requests missing a full browser header set, not to automated
# clients as such. Supplying Accept, Accept-Language, Sec-Fetch-* and
# Upgrade-Insecure-Requests returns 200 for both the page and the asset.
#
# The original conclusion was reached from two failed fetches with a bare
# user-agent. Two negative results with one method were treated as a property of
# the server. That is the same error class this paper is about.
#
# The URL below is a dated CDC asset path and may move when the table is
# revised. If it 404s, open the landing page and take the current link:
#   https://www.cdc.gov/lab-quality/php/get-rm/reference-materials.html
#
# Source: Scheinfeldt L et al. J Mol Diagn 2025;27(6):457-464.
#         doi:10.1016/j.jmoldx.2025.02.008
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../getrm/GeT-RM_consolidated.xlsx"
URL="https://www.cdc.gov/lab-quality/media/files/2025/08/Consolidated_PGx-HLA_table_1-22-25-V4.xlsx"

# SHA-256 of the copy retrieved 2026-07-27. A changed digest means CDC revised
# the table: that is information, not an error. Re-record it deliberately rather
# than letting the truth set change underneath a reported number.
KNOWN_SHA256="c66174befa7a7bbb6cfaf174cc500a269565407807467204bd7bc6eaae3be109"

UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.4 Safari/605.1.15'

echo "fetching $URL"
curl -sSL --compressed --fail \
  -H "User-Agent: $UA" \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' \
  -H 'Accept-Language: en-GB,en;q=0.9' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Referer: https://www.cdc.gov/lab-quality/php/get-rm/reference-materials.html' \
  "$URL" -o "$OUT"

if ! file "$OUT" | grep -qi 'excel'; then
  echo "ERROR: downloaded file is not a spreadsheet (CDC may have served an" >&2
  echo "       error page). Refusing to continue." >&2
  head -c 200 "$OUT" >&2; echo >&2
  exit 1
fi

GOT="$(shasum -a 256 "$OUT" | awk '{print $1}')"
echo "  wrote $OUT"
echo "  sha256 $GOT"
if [ "$GOT" != "$KNOWN_SHA256" ]; then
  echo "  NOTE: digest differs from the recorded 2026-07-27 copy." >&2
  echo "        CDC has revised the table. Re-run the ingester, re-run the" >&2
  echo "        evaluation, and update KNOWN_SHA256 in this script." >&2
fi

echo
echo "next:"
echo "  python real-genome-arm/scripts/11_ingest_getrm_consolidated.py \\"
echo "      --input $OUT \\"
echo "      --out   $HERE/../getrm/getrm_consensus.tsv"
