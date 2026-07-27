# Superseded result files

These are earlier scorings of the matched factorial, kept for audit and clearly
separated so that a reader auditing a cell cannot open one by accident and
arrive at a different number from the manuscript.

A referee reported exactly that: three similarly named five-cell reports sat in
`data/` carrying materially different numbers for the retrieval-executed cell,
with nothing identifying which was authoritative. That was our fault, and this
directory is the fix.

**The authoritative file is `data/v3_five_cell_live_report.txt`.** It is the one
that reproduces Table 1 of the manuscript. See `data/MANIFEST.md` for the full
map from every manuscript table and figure to the file and script that generate
it.

| file | why it is superseded | retrieval-executed A1 |
|---|---|---|
| `v3_five_cell_report.txt` | pre-parser-fix (89186f1) and pre-scorer-fix: uses the retired `drug_match` metric, which required the answer to name the drug and so penalised correct-but-terse prose | 0.8998 |
| `v3_five_cell_common6_report.txt` | six-model subset from an intermediate run, not the eight-model common set the manuscript reports | 0.8864 |

`data/v3_matched_factorial_report.txt` is left in place because
`code/61-rescore-matched.py` writes it as its default report path; it is the
report of whatever was last scored, not a published result, and it must not be
cited.

Both defects behind these numbers are documented in `CORRECTIONS.md`.
