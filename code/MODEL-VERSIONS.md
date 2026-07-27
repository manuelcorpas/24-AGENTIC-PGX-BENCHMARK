# Model versions and API snapshots

Corrected 2026-07-27. The previous version of this file listed the panel used at
original submission and was not updated when three models were substituted for
the revision. A referee who audited it, as requested in the first round, was
therefore given a panel that no longer matched the manuscript. That is the file's
one job, so the correction is recorded here rather than made silently.

## The panel that produced the matched-factorial numbers in v30

Regenerate this table with `python code/69-cell-provenance.py`, which reads the
evaluation rows rather than trusting this document.

| Manuscript name  | API identifier issued              | Pinned? |
|------------------|------------------------------------|---------|
| Claude Opus 4.5  | `claude-opus-4-5-20251101`         | yes (dated) |
| Claude Sonnet 4.5| `claude-sonnet-4-5-20250929`       | yes (dated) |
| GPT-5.2          | `gpt-5.2`                          | yes (fixed release id) |
| GPT-4.1          | `gpt-4.1`                          | yes (fixed release id) |
| o3               | `o3`                               | yes (fixed release id) |
| o4-mini          | `o4-mini`                          | yes (fixed release id) |
| Gemini 2.5 Flash | `gemini-2.5-flash`                 | yes (fixed release id) |
| DeepSeek V3      | `deepseek-chat`                    | yes (fixed release id) |

All five cells were issued against this same eight-model panel. Every row in
every cell carries its own token counts and cost; `code/69-cell-provenance.py`
reports per-cell totals and confirms the panel is identical across cells.

Mistral is reported separately as a paced sensitivity check and is not in the
headline eight-model set. The version issued for the revision is
`mistral-large-2512`, a dated snapshot, replacing the moving pointer used at
submission.

## Substitutions, and why

| Submitted        | Status at revision                  | Replaced by |
|------------------|-------------------------------------|-------------|
| Claude Opus 4 (`claude-opus-4-20250514`)   | retired by vendor, returns HTTP 404 | Claude Opus 4.5 |
| Claude Sonnet 4 (`claude-sonnet-4-20250514`)| retired by vendor, returns HTTP 404 | Claude Sonnet 4.5 |
| Mistral Large 2 (`mistral-large-latest`)    | moving pointer; `mistral-large-2411` withdrawn | `mistral-large-2512`, reported separately |

The submitted run called Mistral through `mistral-large-latest`, which resolves
to whatever the provider serves at request time. The served version was not
recorded in the response metadata and cannot be recovered from the logs. That is
why Mistral is reported separately rather than substituted into the headline: the
original comparison is not reconstructible, and pretending otherwise would be a
provenance claim we cannot support.

## The wider point

A benchmark whose result is a property of specific model versions decays as
vendors retire them, and this one partially did inside a single review cycle.
That is an argument for correctness executed in versioned code, and it is made
in the manuscript's Discussion.
