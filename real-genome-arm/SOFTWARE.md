# Software and versions (real-genome arm)

Exact tool versions used to produce the real-genome results. The pipeline runs on
macOS (Apple Silicon) and Linux; commands assume a POSIX shell.

## External tools

| Tool | Version | Purpose | Install |
|---|---|---|---|
| PLINK2 | v2.0.0-a.7.1 (4 May 2026) | extract PGx gene regions from PLINK sets; export VCF | https://www.cog-genomics.org/plink/2.0/ (static binary, no sudo) |
| PyPGx | 0.26.0 | star-allele diplotype calling (single consistent caller, GRCh37) | `pip install pypgx==0.26.0` |
| pypgx-bundle | git HEAD (1KGP GRCh37 phasing panels) | reference panels PyPGx requires | `git clone https://github.com/sbslee/pypgx-bundle` → set `PYPGX_BUNDLE` or place at `~/pypgx-bundle` |
| Beagle | 22Jul22.46e | haplotype phasing (invoked by PyPGx) | ships inside the PyPGx package |
| Java (OpenJDK) | 17 | runs Beagle | any JDK 17+; ensure `java` is on PATH (keg-only JDKs are not by default) |
| pysam | >= 0.22 | bgzip/tabix indexing of VCFs | `pip install pysam` |
| rclone | v1.74.3 | authenticated download of controlled-access data (no public sharing) | https://rclone.org/downloads/ |
| bcftools | >= 1.19 (optional) | VCF inspection | https://samtools.github.io/bcftools/ |

CrossMap 0.7.0 is **not required** for this pipeline: PyPGx supports GRCh37 natively,
so the data is called on its original build with no liftover (avoiding coordinate-
conversion errors and keeping one caller for every cohort).

## Python

Python >= 3.11 (results produced under 3.13.12). Install Python deps:

```bash
pip install -r requirements.txt
```

## Environment / configuration

- API keys for the model panel are read from the environment, never stored:
  `ANTHROPIC_API_KEY`, `OPENAI_API_KEY`, `GOOGLE_API_KEY`, `DEEPSEEK_API_KEY`.
- PyPGx bundle location: `PYPGX_BUNDLE=/path/to/pypgx-bundle` (or `~/pypgx-bundle`).
- All calling is on **GRCh37**; inputs must use numeric chromosome names (`1`..`22`).

## Determinism notes

- PLINK2 extraction and PyPGx allele matching are deterministic.
- Beagle phasing uses a fixed reference panel (pypgx-bundle); PyPGx invokes it with
  `em=true impute=false` for chip data.
- The agent step is stochastic (LLM generation); the manuscript reports replicate
  behaviour. Reasoning models are queried with provider defaults (see
  `scripts/04_run_agent_realgenome.py`).
