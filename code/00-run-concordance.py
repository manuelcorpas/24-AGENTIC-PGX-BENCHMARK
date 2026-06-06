import json
import os
import time

SPEC = open("/tmp/clawbio-analysis/concordance_spec.md").read()

PROMPT = f"""You are executing a ClawBio pharmacogenomics skill. Follow the specification below EXACTLY. Do not add interpretation, caveats, or additional text. Return ONLY the four required output lines.

{SPEC}
"""

results = {}

# 1. Anthropic Claude
def run_anthropic():
    import anthropic
    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=200,
        messages=[{"role": "user", "content": PROMPT}]
    )
    return response.content[0].text

# 2. OpenAI GPT-4o
def run_openai():
    import openai
    client = openai.OpenAI(api_key=os.environ["OPENAI_API_KEY"])
    response = client.chat.completions.create(
        model="gpt-4o",
        max_tokens=200,
        messages=[{"role": "user", "content": PROMPT}]
    )
    return response.choices[0].message.content

# 3. Google Gemini
def run_gemini():
    import google.generativeai as genai
    genai.configure(api_key=os.environ["GEMINI_API_KEY"])
    model = genai.GenerativeModel("gemini-2.0-flash")
    response = model.generate_content(PROMPT)
    return response.text

# 4. Mistral
def run_mistral():
    from mistralai import Mistral
    client = Mistral(api_key=os.environ["MISTRAL_API_KEY"])
    response = client.chat.complete(
        model="mistral-large-latest",
        messages=[{"role": "user", "content": PROMPT}]
    )
    return response.choices[0].message.content

# 5. DeepSeek
def run_deepseek():
    import openai
    client = openai.OpenAI(
        api_key=os.environ["DEEPSEEK_API_KEY"],
        base_url="https://api.deepseek.com"
    )
    response = client.chat.completions.create(
        model="deepseek-chat",
        max_tokens=200,
        messages=[{"role": "user", "content": PROMPT}]
    )
    return response.choices[0].message.content

models = [
    ("Claude_Sonnet", run_anthropic),
    ("GPT-4o", run_openai),
    ("Gemini_2.0_Flash", run_gemini),
    ("Mistral_Large", run_mistral),
    ("DeepSeek_Chat", run_deepseek),
]

for name, func in models:
    print(f"\n=== {name} ===")
    try:
        result = func()
        results[name] = result
        print(result)
    except Exception as e:
        results[name] = f"ERROR: {e}"
        print(f"ERROR: {e}")

# Save results
with open("/tmp/clawbio-analysis/results/concordance_results.json", "w") as f:
    json.dump(results, f, indent=2)

# Parse and check concordance
print("\n" + "=" * 60)
print("CONCORDANCE ANALYSIS")
print("=" * 60)

def parse_result(text):
    fields = {}
    for line in text.strip().split("\n"):
        line = line.strip()
        for key in ["DIPLOTYPE:", "ACTIVITY_SCORE:", "PHENOTYPE:", "CODEINE_RECOMMENDATION:"]:
            if key in line.upper():
                val = line[line.upper().index(key) + len(key):].strip()
                fields[key.rstrip(":")] = val
    return fields

parsed = {}
for name, result in results.items():
    if not result.startswith("ERROR"):
        parsed[name] = parse_result(result)

print("\nParsed outputs:")
for name, fields in parsed.items():
    print(f"  {name}:")
    for k, v in fields.items():
        print(f"    {k}: {v}")

# Check agreement
if len(parsed) >= 2:
    all_diplotypes = [p.get("DIPLOTYPE", "?") for p in parsed.values()]
    all_phenotypes = [p.get("PHENOTYPE", "?") for p in parsed.values()]
    all_scores = [p.get("ACTIVITY_SCORE", "?") for p in parsed.values()]

    print(f"\nDiplotype agreement: {len(set(all_diplotypes)) == 1} ({set(all_diplotypes)})")
    print(f"Phenotype agreement: {len(set(all_phenotypes)) == 1} ({set(all_phenotypes)})")
    print(f"Activity score agreement: {len(set(all_scores)) == 1} ({set(all_scores)})")
    print(f"\nModels producing output: {len(parsed)}/{len(models)}")

    agree_count = sum(1 for s in [all_diplotypes, all_phenotypes, all_scores] if len(set(s)) == 1)
    print(f"Fields with full agreement: {agree_count}/3")
