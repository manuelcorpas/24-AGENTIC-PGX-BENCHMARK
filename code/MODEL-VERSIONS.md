# Model versions and API snapshots

The nine-model panel and the exact API identifier each condition issues. All models
except Mistral are pinned to a dated or fixed identifier. This file exists to answer
the reviewer request to audit model identifiers.

| Manuscript name  | API identifier issued by the code | Pinned? |
|------------------|------------------------------------|---------|
| Claude Opus 4    | `claude-opus-4-20250514`           | yes (dated) |
| Claude Sonnet 4  | `claude-sonnet-4-20250514`         | yes (dated) |
| GPT-5.2          | `gpt-5.2`                          | yes (fixed release id) |
| GPT-4.1          | `gpt-4.1`                          | yes (fixed release id) |
| o3               | `o3`                               | yes (fixed release id) |
| o4-mini          | `o4-mini`                          | yes (fixed release id) |
| Gemini 2.5 Flash | `gemini-2.5-flash`                 | yes (fixed release id) |
| DeepSeek V3      | `deepseek-chat`                    | yes (fixed release id) |
| Mistral Large 2  | `mistral-large-latest`             | **NO: moving pointer** |

## The Mistral moving-pointer caveat

The code calls Mistral through `mistral-large-latest`, which resolves to whatever
version Mistral serves at request time, rather than a dated snapshot such as
`mistral-large-2411`. The runs used the API key's default served version at the time
of testing (early 2026); the served version was not recorded in the response
metadata, so it cannot be recovered from the raw logs after the fact.

This does not affect any headline number. Mistral returned an elevated non-response
rate under this pointer (about 18% empty cells under the paced protocol; 3.7% usable
in the primary run), so it is **excluded from the headline skill-arm aggregate and
reported separately** in the manuscript. The primary comparison is the common
eight-model set.

Reconciliation performed for the revision:
1. All headline figures are reported on the common eight-model set as the primary
   comparison (per reviewer request).
2. Mistral is re-run under a dated, pinned identifier (`mistral-large-2411`) as a
   sensitivity check; both the pointer value and the dated identifier are recorded in
   the response metadata of the re-run so the served version is auditable.
3. The manuscript model table is reconciled to the exact identifiers in this file.
