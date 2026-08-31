# Model versions used in the V4 Flash revision

The primary matched five-configuration comparison uses the same eight-model
panel in every configuration. The DeepSeek arm was rerun across every reported
prompt family on 30 August 2026. All 8,901 canonical requests were explicitly
issued to `deepseek-v4-flash`, returned `deepseek-v4-flash`, and used disabled
thinking. No historical DeepSeek response was relabelled.

| Manuscript name | API identifier issued | Identifier returned |
|---|---|---|
| Claude Opus 4.5 | `claude-opus-4-5-20251101` | retained per provider row |
| Claude Sonnet 4.5 | `claude-sonnet-4-5-20250929` | retained per provider row |
| DeepSeek V4 Flash | `deepseek-v4-flash` | `deepseek-v4-flash` |
| GPT-4.1 | `gpt-4.1` | `gpt-4.1-2025-04-14` where returned |
| GPT-5.2 | `gpt-5.2` | `gpt-5.2-2025-12-11` where returned |
| Gemini 2.5 Flash | `gemini-2.5-flash` | `gemini-2.5-flash` |
| o3 | `o3` | `o3-2025-04-16` where returned |
| o4-mini | `o4-mini` | `o4-mini-2025-04-16` where returned |

The matched analysis contains 13,200 attempted evaluations. Exactly 13,199
records returned because one prespecified Gemini 2.5 Flash request in authored
rule execution returned no record. It remains a failure in the attempted
denominator.

Mistral Large 2512 (`mistral-large-2512`) is reported only as a paced
sensitivity analysis and is not part of the common eight-model comparison.

## Figure 3 version boundary

Figure 3A reports all eight requested models on the four real-genome cohorts.
Figure 3B pools the seven models with version-matched curated and real-genome
results. No matched curated V4 Flash result exists, so the DeepSeek real-genome
row is shown in panel A and is not included in the pooled panel B comparison.

## Input normalisation

The reported input-normalisation analysis uses the common eight-model panel.
The definition-supplied arm contains all 527 variant-call pairs for every
model. In the paired no-definition comparison, o3 contains its prespecified
150-pair gene-stratified arm and the other models contain all 527 pairs.

## Canonical provenance

The ten canonical V4 files and their SHA-256 hashes are fixed in
`docs/deepseek-v4-flash-rerun-audit-2026-08-30/FULL-RERUN-VALIDATION.json`.
Provider response identifiers, request and returned model identifiers, thinking
mode, stop reasons, usage, timestamps, latency, retry history, raw text and
error bodies are retained where supplied by the provider and runner.
