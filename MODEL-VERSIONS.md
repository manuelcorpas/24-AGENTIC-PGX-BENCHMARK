# Model versions used in the corrected revision

The matched five-cell comparison used the same eight-model panel in every cell.
The mapping below is regenerated from `data/v3_cell_provenance.json` by
`code/69-cell-provenance.py`.

| Manuscript name | API identifier issued |
|---|---|
| Claude Opus 4.5 | `claude-opus-4-5-20251101` |
| Claude Sonnet 4.5 | `claude-sonnet-4-5-20250929` |
| DeepSeek V3 | `deepseek-chat` |
| GPT-4.1 | `gpt-4.1` |
| GPT-5.2 | `gpt-5.2` |
| Gemini 2.5 Flash | `gemini-2.5-flash` |
| o3 | `o3` |
| o4-mini | `o4-mini` |

Each model contributed 330 returned records to each cell except Gemini 2.5
Flash in authored-rule execution, where one prespecified request returned no
record. That request remains part of the 13,200 attempted denominator and is
counted as a failure.

Mistral Large 2512 (`mistral-large-2512`) is reported only as a paced
sensitivity analysis and is not part of the common eight-model comparison.

The input-normalisation experiment used seven models in July 2026: Claude Opus
4.5, Claude Sonnet 4.5, GPT-5.2, GPT-4.1, o3, o4-mini and DeepSeek V3. The
definition-supplied arm used all 527 variant-call pairs for every model. The
paired no-definition comparison used all 527 variant-call pairs except for o3,
whose prespecified gene-stratified arm contains 150 pairs. Claude Opus 4.5,
GPT-5.2 and o3 also have the original three-rendering no-definition experiment.

Gemini 2.5 Flash receives no input-normalisation performance estimate because a
complete comparable run was not collected. Its pilot returned an explicitly
marked diplotype inside a LaTeX box. The frozen C15 parser accepts that wrapper
under the same policy used for every model; Gemini is not excluded because of
formatting, and no result is inferred from the incomplete pilot.
