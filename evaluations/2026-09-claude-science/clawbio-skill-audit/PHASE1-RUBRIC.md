# Phase 1 classification rubric

The criteria every reader applies, recorded here so the classification is reproducible and so
disagreement between readers can be audited later. Readers follow this file; they do not
invent categories.

## The question per skill

> If a frontier model with tool and database access were handed this skill's task tomorrow,
> **what would be lost by deleting the skill?**

Classify by what would be lost, not by what the skill contains.

## Primary class, one per skill, by centre of gravity

**D1 knowledge-encoding.** What would be lost is a mapping, threshold or criteria set that
already exists publicly and authoritatively elsewhere: CPIC tables, ACMG criteria, ClinGen
dosage ratings, published formulas, guideline cut-offs. Test: could a capable agent with the
source open reconstruct the decision logic in a session? If yes, D1.

**D2 tool-orchestration.** What would be lost is correct invocation and sequencing of external
software or services: nextflow and nf-core pipelines, VEP, PLINK, samtools, scanpy, REST APIs.
Includes handling their formats, failure modes and version quirks. Test: is the hard part
driving something else correctly? If yes, D2.

**D3 contract and policy.** What would be lost is what the skill refuses, validates, records or
guarantees: input validation, abstention when evidence is insufficient, provenance, determinism,
reporting standards, safety gates, governance. Test: does the skill's value survive the model
getting smarter? If yes, D3.

**D0 unclear.** The SKILL.md is too thin to support a judgement, typically spec-only stubs. Use
this rather than guessing. D0 is a finding, not a failure.

Many skills contain all three. Pick the one whose loss would hurt most, and say why in the
justification.

## Modifiers, one value each

- `verifiability`: `ground-truth` (a correct answer exists and is checkable), `partial`
  (checkable in places), `none` (no reference answer exists).
- `consequence`: `clinical` (output could reach a patient decision), `research` (feeds a
  scientific claim), `exploratory` (browsing, discovery, convenience).
- `locality`: `local-required` (needs private data, a named cluster, or a governance regime a
  consumer product cannot reach) or `portable`.

## Evidence discipline

- Read `SKILL.md`. If the class is still unclear, read the main script's first 60 lines.
- The justification must cite something specific from the file: a named table, a command, a
  guideline, a validation step. "Encodes clinical logic" is not a justification.
- Classify what the file actually says, not what the name suggests. `article-data-fetcher`
  sounds like orchestration and downloads deposited files instead.
- If the inventory's signal label contradicts what you read, set `signal_misleading: true`.
  That flag is itself a measurement of the heuristic's error rate.
- Low confidence is an acceptable answer. An invented class is not.

## Output schema, strict

One JSON object per skill, no prose outside the JSON:

```json
{
  "name": "clinpgx",
  "primary_class": "D1",
  "confidence": "high",
  "verifiability": "ground-truth",
  "consequence": "clinical",
  "locality": "portable",
  "justification": "SKILL.md encodes the CPIC diplotype-to-phenotype table and cites CPIC v1.9 directly.",
  "signal_misleading": false
}
```

Controlled vocabularies: `primary_class` D1|D2|D3|D0, `confidence` high|medium|low,
`verifiability` ground-truth|partial|none, `consequence` clinical|research|exploratory,
`locality` portable|local-required. Anything else is rejected by the validator.

## Why these classes

D1 is the wasting asset: a general agent with CPIC's API and no skill scored 0.942 (sd 0.004,
n=3) against 0.966 for the authored skill on the 110-case pharmacogenomics benchmark. D3 is
what did not transfer: the same agent abstained 0, 5 and 0 times across identical runs, filled
missing genotypes with reference alleles, and returned Ultrarapid Metaboliser for an input
stating "no duplication detected". D2 sits between and is being absorbed by harnesses that ship
their own connectors.
