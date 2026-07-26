# GeT-RM caller-truth set (revision item N5)

## Why this exists

N0 validated the executed skill's **mapping** step against PyPGx. It did not
validate the **calling** step, because PyPGx both called the diplotype and
supplied the phenotype it was scored against. Input-call error is therefore
unobservable in N0 rather than eliminated, and the plan of revision says so
explicitly to the editor.

GeT-RM supplies external truth for the calling step. With it, the claim "100%
correctness among emitted answers" becomes genuinely falsifiable: a call that
disagrees with the GeT-RM consensus is an end-to-end error attributable to
input interpretation.

## Acquisition is manual, and this is not an oversight

The CDC GeT-RM pages are JavaScript shells with no tabular payload, the NCBI
FTP path returns 404, and the consensus diplotypes are published as journal
supplementary tables. There is no stable machine-readable endpoint to pin, so
an automated fetch would be a scraper that breaks silently and produces a
truth set nobody can reconstruct. The download is therefore a documented human
step and the parser validates whatever arrives.

Sources (cite whichever is used):

- Pratt VM et al. *Characterization of 137 Genomic DNA Reference Materials for
  28 Pharmacogenetic Genes: A GeT-RM Collaborative Project.* J Mol Diagn 2016.
- Pratt VM et al. *Characterization of 137 Genomic DNA Reference Materials...
  (expanded).* J Mol Diagn 2010.
- Gaedigk A et al. GeT-RM CYP2D6 characterisation updates.

## Expected schema

Place a TSV here named `getrm_consensus.tsv` with at least these columns:

    sample	gene	diplotype

`sample` must be the Coriell identifier (NA12878 style) so it joins to the
1000 Genomes calls produced by `scripts/09_call_1000g_eur.sh`.

## Validation

    python real-genome-arm/scripts/08_caller_truth_eval.py \
        --truth real-genome-arm/getrm/getrm_consensus.tsv \
        --called <pypgx calls tsv>

The evaluator refuses to run on an empty or malformed truth set rather than
reporting a vacuous 100%.
