#!/usr/bin/env bash
# Make the public repository consistent with manuscript v30, so that every file
# the response letter cites resolves at the tag the manuscript points to.
#
# WHY THIS EXISTS
# The response cites 11 files as evidence for the referee's points. At the tag
# cg-revision-2026-07-27 only three of them exist, because the extraction
# experiment, the provenance tooling, the case-level lethal analysis, the
# deterministic baseline and data/MANIFEST.md were all committed after the tag
# was cut. Sending in that state reproduces the referee's "the package does not
# match the manuscript" complaint for a third consecutive round, on the exact
# files offered as the remedy.
#
# DRY RUN BY DEFAULT. Nothing is pushed and no tag moves unless you pass
# --confirm. Read the plan it prints first.
#
#   bash prepare-send.sh              # show what would happen
#   bash prepare-send.sh --confirm    # do it
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

TAG=cg-revision-2026-07-27
LETTER="$HOME/Library/Mobile Documents/com~apple~CloudDocs/PUBLICATIONS/00-CLAWBIO-BENCHMARK/BENCHMARK/DOCS/CELL-GENOMICS/Response-to-Reviewers-v30.md"
CONFIRM=${1:-}

cited() {
  /opt/homebrew/bin/python3 - "$LETTER" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
for f in sorted(set(re.findall(r'`([A-Za-z0-9_./-]+\.(?:py|md|sh|txt|json|tsv))`', t))):
    print(f)
PY
}

echo "== current state =="
echo "  origin/main : $(git rev-parse --short origin/main 2>/dev/null || echo unknown)"
echo "  local HEAD  : $(git rev-parse --short HEAD)"
echo "  $TAG -> $(git rev-parse --short "$TAG^{commit}" 2>/dev/null || echo 'absent')"
echo "  unpushed    : $(git rev-list --count origin/main..HEAD 2>/dev/null || echo '?') commits"
echo

echo "== cited files that do NOT resolve at $TAG =="
miss_tag=0
while read -r f; do
  [ -z "$f" ] && continue
  git cat-file -e "$TAG:$f" 2>/dev/null || { echo "  $f"; miss_tag=$((miss_tag+1)); }
done < <(cited)
[ "$miss_tag" -eq 0 ] && echo "  none"
echo

echo "== of those, still missing at local HEAD (must be committed first) =="
miss_head=0
while read -r f; do
  [ -z "$f" ] && continue
  git cat-file -e "HEAD:$f" 2>/dev/null || { echo "  $f   <-- commit this"; miss_head=$((miss_head+1)); }
done < <(cited)
[ "$miss_head" -eq 0 ] && echo "  none: everything cited exists at HEAD"
echo

if [ "$miss_head" -gt 0 ]; then
  echo "STOP. Commit the file(s) above first; this script will not commit files it"
  echo "does not own. CORRECTIONS.md belongs to the other session and is cited by"
  echo "both the response letter and the manuscript's Limitations."
  exit 1
fi

echo "== plan =="
echo "  1. git push origin main            ($(git rev-list --count origin/main..HEAD) commits)"
echo "  2. move $TAG to $(git rev-parse --short HEAD) and force-push the tag"
echo "  3. re-verify every cited file resolves at the tag"
echo

if [ "$CONFIRM" != "--confirm" ]; then
  echo "DRY RUN. Re-run with --confirm to execute."
  exit 0
fi

echo "== executing =="
git push origin main || { echo "push failed"; exit 1; }
git tag -f -a "$TAG" -m "Manuscript v30 revision snapshot

Retagged so that every file cited by the response to the referee resolves here:
the extraction experiment, per-cell provenance, case-level lethal analysis, the
deterministic baseline on real cohorts, the fabrication firewall and
data/MANIFEST.md were all committed after the tag was first cut." || exit 1
git push --force origin "$TAG" || { echo "tag push failed"; exit 1; }

echo
echo "== verification =="
git fetch origin --tags --force >/dev/null 2>&1
bad=0
while read -r f; do
  [ -z "$f" ] && continue
  git cat-file -e "$TAG:$f" 2>/dev/null || { echo "  STILL MISSING  $f"; bad=$((bad+1)); }
done < <(cited)
if [ "$bad" -eq 0 ]; then
  echo "  all cited files resolve at $TAG -> $(git rev-parse --short "$TAG^{commit}")"
  echo
  echo "Safe to send. The manuscript's citation now matches what a referee will find."
else
  echo "  $bad file(s) still missing; do not send."
  exit 1
fi
