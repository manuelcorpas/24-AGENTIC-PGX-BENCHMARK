# Corpas family, whole-genome arm

Supersedes the 23andMe SNP-chip family arm for every cohort number in the
manuscript. Four members with a GRCh37 Sentieon build (PT00002A son, PT00007A
father, PT00008A mother, PT00009A daughter), sequenced on BGISEQ-500, observed
control-gene depth 23.8x, 25.2x, 30.0x and 43.9x. A fifth member (PT00010A) has
only a 4.4x GRCh38 build and is excluded: star-allele calling is unreliable at
that depth and CYP2D6 copy number is uncallable.

Regenerate:
  15/16  extract PGx regions per member, then call both arms
  19     re-run the model panel over the resulting states
  17     aggregate, Mendelian check, NGS-vs-chip concordance

Files
  family_wgs_ngs.tsv    headline arm: NGS pipeline, copy number for CYP2B6,
                        CYP2D6 and CYP4F2
  family_wgs_chip.tsv   control arm: chip pipeline, the calling mode used for
                        the Iberian, Peruvian and Ugandan cohorts
  family_wgs_qc.txt     Mendelian consistency and cross-arm concordance

Key results
  68 diplotype calls, 37 distinct (gene, diplotype) states
  NGS vs chip: 66/68 identical (97.1%); both differences are CYP2D6, the gene
    where copy number is decisive
  Mendelian: 2 inconsistent of 32 checked, 2 not checkable, under the NGS arm;
    3 of 34 under the chip arm. Copy-number analysis found a Tandem1A
    arrangement in mother and daughter, which explains the CYP2D6 inconsistency
    and converts the mother's call from a confident *10/*41 to an abstention.
