# DeepSeek V4 Flash paper-wide rerun audit

Date: 2026-08-30

## Decision

The paper cannot support the label `DeepSeek V3`. The relevant runners issued
the rolling alias `deepseek-chat`, and most historical rows did not retain the
model identifier returned by the provider. DeepSeek documents that this alias
routed to V4 Flash non-thinking mode from 24 April 2026. The corrected
real-genome run was the first paper experiment to retain enough provider
metadata to expose the returned identifier, `deepseek-v4-flash`.

To preserve the eight-model design and use one DeepSeek model throughout the
paper, every currently reported DeepSeek result must be rerun with:

```text
model = deepseek-v4-flash
thinking.type = disabled
```

Disabling thinking is essential. Explicit V4 Flash now defaults to thinking
enabled, whereas DeepSeek documents the retired `deepseek-chat` alias as the
non-thinking route.

## Scope

The complete scope is 8,901 calls across ten jobs:

| Experiment | Calls | Current reported outputs |
|---|---:|---|
| Five matched cells | 1,650 | Table 1, Figure 2, Figures S2 to S5 and S9, Tables S1 to S5 |
| Extracted rule tables | 21 | Table 2 |
| Adversarial contract, forward | 15 | Figure S7, Table S6 |
| Adversarial contract, reverse | 15 | Figure S7, Table S6 |
| Drug-keyed retrieval | 660 | Table S7 |
| Population framing, three arms | 2,970 | Figure S8 |
| Population framing, two skill arms | 1,980 | Figure S8 |
| Input normalisation without definitions | 527 | Figure 4, Tables S8 and S9 |
| Input normalisation with definitions | 527 | Figure 4, Tables S8 and S9 |
| Four-cohort real-genome interpretation | 536 | Figure 3 |
| **Total** | **8,901** | |

The 4,950 population-framing calls and 660 drug-keyed calls are easy to miss
because the manuscript does not print a DeepSeek-specific number for them.
DeepSeek nevertheless contributes to the pooled results, so leaving those rows
unchanged would still mix model versions inside the paper.

## Not current and therefore not rerun

The following DeepSeek-containing files are historical or superseded and are
not authoritative inputs to a current submitted display item:

- `data/v3_family_wgs_preds.json`
- `data/v3_rag_genedrug_chunking.json`
- `data/v3_armA9_armBv2_POP.json`
- `data/v3_five_cell_common6.json`
- `data/v3_matched_factorial_clean.json`
- sync-conflict files and files under `data/SUPERSEDED/`

They must remain available for provenance and must not be silently relabelled.

## Execution order

1. Patch and test every affected runner before issuing calls. Each runner must
   use the explicit identifier, disable thinking, retain the full provider
   response and refuse a returned-model mismatch.
2. Run one pilot for each of the 18 distinct prompt families against the exact
   frozen inputs. Confirm the request payload, returned model, output budget,
   parser compatibility and projected spend for each prompt family.
3. Execute jobs DSV4-01 through DSV4-10 into new `data/v4_flash_*` files. Never
   overwrite historical data.
4. Validate 8,901 unique completed units, with no missing or duplicate unit and
   no request through `deepseek-chat`.
5. Replace only the DeepSeek rows in each canonical analysis dataset.
6. Regenerate all statistics, figures, supplementary tables and model audits.
7. Update the manuscript, supplement and response letter. Remove every
   `DeepSeek V3` label and the temporary Figure 3 V3/V4 mismatch explanation.
8. Run the registered-number validator and submission-package validator.
9. Cut a new immutable repository tag and Zenodo version. Do not modify the
   cited v2.8 release or any existing public URL.

## Non-averageable gates

One wrong returned model identifier fails the run. One missing expected unit
fails the run. Historical rows cannot be relabelled as V4. A mixture of old and
new DeepSeek rows cannot be pooled. The paper is synchronised only when the
main manuscript, supplement, response letter, figures, data, code, tag and
Zenodo deposit all identify DeepSeek V4 Flash consistently.

## Files

- `RERUN-MANIFEST.json`: machine-checkable job and dependency manifest.
- `SOURCE-HASHES.tsv`: frozen input and historical-source hashes recorded at
  audit time.
- `pilots/deepseek-v4-flash-pilots.jsonl`: 18 live pilot requests with prompts,
  hashes, complete request and response bodies, usage, timing and attempt logs.
- `PILOT-VALIDATION.json`: deterministic validation summary for the pilot file.
- `tests/test_deepseek_v4_rerun_manifest.py`: manifest completeness checks.

## Pilot result

All 18 prompt-family pilots passed on 30 August 2026. Every request explicitly
used `deepseek-v4-flash` with thinking disabled. Every provider response
returned `deepseek-v4-flash`, completed with stop reason `stop`, and succeeded
on its first attempt. The pilots consumed 29,019 input tokens and 1,768 output
tokens. They validate routing, output budgets and provenance capture. They do
not establish the scientific results of the full rerun.

## Full rerun result

Manuel authorised the exact 8,901-call rerun on 30 August 2026. All 8,901
canonical units completed and passed the model, thinking-mode, provenance,
error and uniqueness gates. The canonical set contains 10,095,507 input tokens
and 1,021,874 output tokens, with an estimated recorded cost of $3.85 under the
runner price table.

The first rule-extraction pass exposed an output-budget artefact in G6PD,
MT-RNR1 and CYP2C19. Their incomplete responses reached the configured token
ceiling and were not treated as scientific results. The affected genes were
rerun with a larger ceiling. All 21 canonical rule tables then parsed with
`stop` completion. The failed attempts remain preserved outside the canonical
file.

`FULL-RERUN-VALIDATION.json` records the ten canonical files, exact unit counts
and SHA-256 hashes. Statistical analysis, figures, tables and manuscript text
have not yet been regenerated from these new rows.

## Official provider evidence

- DeepSeek change log: <https://api-docs.deepseek.com/updates/>
- Thinking mode: <https://api-docs.deepseek.com/guides/thinking_mode/>
- Chat Completions schema: <https://api-docs.deepseek.com/api/create-chat-completion/>
