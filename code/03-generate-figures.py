"""
Generate publication-quality figures for the ClawBio PGx benchmark paper.
Corpas, Fatumo & Guio (2026) PLOS Computational Biology.

Figures:
  2. Clinical correctness (Tier A) with vs without specification
  3. Multidimensional radar (6 dimensions, 2 tiers)
  4. Consistency heatmap (3-run agreement per model x test case)
  5. Population-specific performance (requires v2 data)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import json
import os
import warnings
warnings.filterwarnings('ignore')

# ===== PATHS (self-relative) =====
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(BENCHMARK_DIR, 'RESULTS')
FIGURES_DIR = os.path.join(BENCHMARK_DIR, 'FIGURES')

# ===== STYLE =====
plt.style.use('default')
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
})

# ===== LOAD DATA =====
def load_data():
    v2_path = os.path.join(RESULTS_DIR, 'v2_raw.json')
    v1_path = os.path.join(RESULTS_DIR, 'raw_results.json')
    data_path = v2_path if os.path.exists(v2_path) else v1_path
    with open(data_path) as f:
        raw = json.load(f)
    print(f"Loaded {len(raw)} results from {os.path.basename(data_path)}")
    return raw, data_path

# ===== DIMENSION MAPPING =====
# v1 data uses different dimension names than our two-tier framework.
# Map them for consistent access.
V1_DIM_MAP = {
    "A1": "factual_correctness",
    "A2": "specification_compliance",
    "A3": "clinical_safety",
    "B1": "semantic_accuracy",
    "B2": "reasoning_quality",
    "B3": "domain_knowledge",
}

def get_score(result, dim):
    """Get score for a dimension, handling v1 vs v2 naming."""
    scores = result.get("scores", {})
    if dim in scores:
        return scores[dim]
    mapped = V1_DIM_MAP.get(dim)
    if mapped and mapped in scores:
        return scores[mapped]
    return 0

def get_condition(result):
    return result.get("cond", result.get("condition", ""))

def get_test_case(result):
    return result.get("tc", result.get("test_case", ""))

def is_valid(result):
    """Check if result parsed correctly (handles both dict and bool parsed)."""
    p = result.get("parsed")
    if isinstance(p, dict):
        return bool(p)
    return bool(p)

# ===== MODEL CONFIG =====
MODEL_ORDER = [
    "Claude Opus 4", "Claude Sonnet 4", "GPT-5.2", "GPT-4.1",
    "o3", "o4-mini", "Gemini 2.5 Flash", "DeepSeek V3", "Mistral Large 2"
]

COLORS = {
    "Claude Opus 4": "#8B5CF6",
    "Claude Sonnet 4": "#A78BFA",
    "GPT-5.2": "#10B981",
    "GPT-4.1": "#34D399",
    "o3": "#F59E0B",
    "o4-mini": "#FBBF24",
    "Gemini 2.5 Flash": "#3B82F6",
    "DeepSeek V3": "#EF4444",
    "Mistral Large 2": "#EC4899",
    "Mistral Large": "#EC4899",
}

SHORT_NAMES = {
    "Claude Opus 4": "Claude\nOpus 4",
    "Claude Sonnet 4": "Claude\nSonnet 4",
    "GPT-5.2": "GPT-5.2",
    "GPT-4.1": "GPT-4.1",
    "o3": "o3",
    "o4-mini": "o4-mini",
    "Gemini 2.5 Flash": "Gemini\n2.5 Flash",
    "DeepSeek V3": "DeepSeek\nV3",
    "Mistral Large 2": "Mistral\nLarge 2",
    "Mistral Large": "Mistral\nLarge",
}

def get_models_present(raw):
    fallback = list(MODEL_ORDER)
    if not any(r.get("model") == "Mistral Large 2" for r in raw):
        fallback = [("Mistral Large" if m == "Mistral Large 2" else m) for m in fallback]
    return [m for m in fallback if any(r.get("model") == m for r in raw)]


# ===== FIGURE 2: TIER A CLINICAL CORRECTNESS =====
def figure2_clinical_correctness(raw, models):
    """Grouped bar chart: Tier A scores by model, no_spec vs with_spec.

    Uses error bars showing stdev across test cases to reveal variance
    that the ceiling-effect means hide.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), sharey=True)

    tier_a = [("A1", "Phenotype Accuracy"), ("A2", "Drug Recommendation"), ("A3", "Clinical Safety")]

    for ax_idx, (dim, label) in enumerate(tier_a):
        ax = axes[ax_idx]
        no_spec_means, no_spec_stds = [], []
        wi_spec_means, wi_spec_stds = [], []

        for m in models:
            for cond, means_list, stds_list in [
                ("no_spec", no_spec_means, no_spec_stds),
                ("with_spec", wi_spec_means, wi_spec_stds)
            ]:
                vals = [get_score(r, dim) for r in raw
                        if r.get("model") == m and get_condition(r) == cond and is_valid(r)]
                means_list.append(np.mean(vals) if vals else 0)
                stds_list.append(np.std(vals) if vals else 0)

        x = np.arange(len(models))
        width = 0.35

        ax.bar(x - width/2, no_spec_means, width, yerr=no_spec_stds,
               label='Without specification', color='#94A3B8', edgecolor='white',
               linewidth=0.5, capsize=2, error_kw={'linewidth': 0.8})
        ax.bar(x + width/2, wi_spec_means, width, yerr=wi_spec_stds,
               label='With ClawBio specification', color='#22C55E', edgecolor='white',
               linewidth=0.5, capsize=2, error_kw={'linewidth': 0.8})

        ax.set_title(label, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([SHORT_NAMES.get(m, m) for m in models],
                           rotation=45, ha='right', fontsize=7)
        ax.set_ylim(0, 1.15)
        ax.set_ylabel('Mean score' if ax_idx == 0 else '')
        ax.axhline(y=1.0, color='#CBD5E1', linestyle='--', linewidth=0.5)

        if ax_idx == 0:
            ax.legend(loc='lower left', framealpha=0.9, fontsize=8)

    fig.suptitle('Figure 2. Tier A Clinical Correctness: With vs Without ClawBio Specification',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['png', 'tiff']:
        plt.savefig(os.path.join(FIGURES_DIR, f'Figure2_clinical_correctness.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Figure 2 saved")
    plt.close()


# ===== FIGURE 3: MULTIDIMENSIONAL RADAR =====
def figure3_radar(raw, models):
    """Radar chart: 6 dimensions per model, side-by-side no_spec vs with_spec."""
    dims = ["A1", "A2", "A3", "B1", "B2", "B3"]
    dim_names = ["Phenotype\nAccuracy", "Drug\nRecommendation", "Clinical\nSafety",
                 "Dataset\nSpecificity", "Reasoning\nChain", "Domain\nGrounding"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), subplot_kw=dict(polar=True))

    for ax_idx, cond in enumerate(["no_spec", "with_spec"]):
        ax = axes[ax_idx]
        cond_label = "Without Specification" if cond == "no_spec" else "With ClawBio Specification"

        angles = np.linspace(0, 2 * np.pi, len(dims), endpoint=False).tolist()
        angles += angles[:1]

        for m in models:
            valid = [r for r in raw if r.get("model") == m and get_condition(r) == cond and is_valid(r)]
            if not valid:
                continue
            scores = [np.mean([get_score(r, d) for r in valid]) for d in dims]
            scores += scores[:1]

            ax.plot(angles, scores, 'o-', linewidth=1.5, markersize=4,
                    label=m, color=COLORS.get(m, '#666'), alpha=0.8)
            ax.fill(angles, scores, alpha=0.05, color=COLORS.get(m, '#666'))

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(dim_names, size=8)
        ax.set_ylim(0, 1.1)
        ax.set_title(cond_label, fontsize=11, fontweight='bold', pad=25)

        if ax_idx == 1:
            ax.legend(loc='upper right', bbox_to_anchor=(1.5, 1.05), fontsize=8)

    fig.suptitle('Figure 3. Multidimensional Performance Across 6 Evaluation Dimensions',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    for ext in ['png', 'tiff']:
        plt.savefig(os.path.join(FIGURES_DIR, f'Figure3_radar.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Figure 3 saved")
    plt.close()


# ===== FIGURE 4: CONSISTENCY HEATMAP =====
def figure4_consistency(raw, models):
    """Heatmap: 3-run consistency per model x test case, side-by-side."""
    test_cases = sorted(set(get_test_case(r) for r in raw))

    # Prettier test case labels
    tc_labels = {
        'cyp2b6_pm': 'CYP2B6\nPM',
        'cyp2c19_pm': 'CYP2C19\nPM',
        'cyp2c19_rapid': 'CYP2C19\nRapid',
        'cyp2c9_pm': 'CYP2C9\nPM',
        'cyp2d6_normal': 'CYP2D6\nNormal',
        'cyp2d6_pm': 'CYP2D6\nPM',
        'dpyd_het': 'DPYD\nhet',
        'dpyd_hom': 'DPYD\nhom *',
        'dpyd_normal': 'DPYD\nnormal',
        'slco1b1': 'SLCO1B1',
        'tpmt_het': 'TPMT\nhet',
        'ugt1a1_hom': 'UGT1A1\nhom',
    }

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    for ax_idx, cond in enumerate(["no_spec", "with_spec"]):
        ax = axes[ax_idx]
        cond_label = "Without Specification" if cond == "no_spec" else "With ClawBio Specification"

        matrix = np.full((len(models), len(test_cases)), np.nan)

        for i, m in enumerate(models):
            for j, tc in enumerate(test_cases):
                runs = [r for r in raw
                        if r.get("model") == m
                        and get_test_case(r) == tc
                        and get_condition(r) == cond]
                valid_runs = [r for r in runs if is_valid(r)]
                if valid_runs:
                    correct = sum(1 for r in valid_runs if get_score(r, "A1") >= 1.0)
                    matrix[i, j] = correct / len(valid_runs)
                elif runs:
                    matrix[i, j] = 0  # all failed to parse

        im = ax.imshow(matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')

        ax.set_xticks(range(len(test_cases)))
        ax.set_xticklabels([tc_labels.get(tc, tc) for tc in test_cases],
                           fontsize=7, ha='center')
        ax.set_yticks(range(len(models)))
        ax.set_yticklabels(models, fontsize=9)
        ax.set_title(cond_label, fontweight='bold', fontsize=11)

        for i in range(len(models)):
            for j in range(len(test_cases)):
                val = matrix[i, j]
                if not np.isnan(val):
                    text = f"{val:.0%}"
                    color = 'white' if val < 0.5 else 'black'
                    fontw = 'bold' if val < 1.0 else 'normal'
                    ax.text(j, i, text, ha='center', va='center',
                            fontsize=7, color=color, fontweight=fontw)

    fig.colorbar(im, ax=axes, shrink=0.6, label='Fraction correct (3 runs)',
                 orientation='vertical', pad=0.02)
    fig.suptitle('Figure 4. Consistency: Fraction of Correct Runs per Model and Test Case',
                 fontsize=13, fontweight='bold')
    fig.text(0.5, -0.02, '* DPYD hom (rs3918290 T/T): standard fluorouracil dosing is potentially lethal',
             ha='center', fontsize=8, fontstyle='italic', color='#DC2626')
    plt.tight_layout()
    for ext in ['png', 'tiff']:
        plt.savefig(os.path.join(FIGURES_DIR, f'Figure4_consistency.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Figure 4 saved")
    plt.close()


# ===== FIGURE 5: SUMMARY COMPARISON =====
def figure5_summary(raw, models):
    """
    Summary bar chart: aggregate metrics with vs without specification.
    Shows the key numbers from the paper in a single visual.
    If v2 data with population dimension is available, shows population breakdown.
    """
    has_pop = any("pop" in r for r in raw)

    if has_pop:
        figure5_population(raw, models)
        return

    # Aggregate summary figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Mean Tier A score by model
    ax = axes[0]
    for cond_idx, (cond, color, label) in enumerate([
        ("no_spec", "#94A3B8", "Without spec"),
        ("with_spec", "#22C55E", "With spec")
    ]):
        tier_a_means = []
        for m in models:
            valid = [r for r in raw if r.get("model") == m and get_condition(r) == cond and is_valid(r)]
            if valid:
                a1 = np.mean([get_score(r, "A1") for r in valid])
                a2 = np.mean([get_score(r, "A2") for r in valid])
                a3 = np.mean([get_score(r, "A3") for r in valid])
                tier_a_means.append(np.mean([a1, a2, a3]))
            else:
                tier_a_means.append(0)

        x = np.arange(len(models))
        width = 0.35
        offset = -width/2 if cond_idx == 0 else width/2
        ax.bar(x + offset, tier_a_means, width, label=label, color=color,
               edgecolor='white', linewidth=0.5)

    ax.set_title('Mean Tier A Score (Clinical Correctness)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_NAMES.get(m, m) for m in models], rotation=45, ha='right', fontsize=7)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel('Mean Tier A score')
    ax.axhline(y=1.0, color='#CBD5E1', linestyle='--', linewidth=0.5)
    ax.legend(loc='lower left', fontsize=8)

    # Right: Consistency (perfect 3/3 rate)
    ax = axes[1]
    test_cases = sorted(set(get_test_case(r) for r in raw))

    for cond_idx, (cond, color, label) in enumerate([
        ("no_spec", "#94A3B8", "Without spec"),
        ("with_spec", "#22C55E", "With spec")
    ]):
        consistency_rates = []
        for m in models:
            perfect = 0
            total = 0
            for tc in test_cases:
                runs = [r for r in raw if r.get("model") == m and get_test_case(r) == tc
                        and get_condition(r) == cond]
                valid = [r for r in runs if is_valid(r)]
                if valid:
                    total += 1
                    correct = sum(1 for r in valid if get_score(r, "A1") >= 1.0)
                    if correct == len(valid):
                        perfect += 1
            consistency_rates.append(perfect / total if total > 0 else 0)

        x = np.arange(len(models))
        width = 0.35
        offset = -width/2 if cond_idx == 0 else width/2
        ax.bar(x + offset, consistency_rates, width, label=label, color=color,
               edgecolor='white', linewidth=0.5)

    ax.set_title('Perfect Consistency (3/3 correct runs)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT_NAMES.get(m, m) for m in models], rotation=45, ha='right', fontsize=7)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel('Fraction of test cases with 3/3 correct')
    ax.axhline(y=1.0, color='#CBD5E1', linestyle='--', linewidth=0.5)
    ax.legend(loc='lower left', fontsize=8)

    fig.suptitle('Figure 5. Aggregate Performance: Clinical Correctness and Consistency',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    for ext in ['png', 'tiff']:
        plt.savefig(os.path.join(FIGURES_DIR, f'Figure5_summary.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Figure 5 saved")
    plt.close()


def figure5_population(raw, models):
    """Population-specific performance (v2 data only)."""
    pops = ["EUR", "AMR", "AFR"]
    pop_labels = ["European\n(Corpasome)", "Latin American\n(Peru)", "East African\n(Uganda)"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax_idx, cond in enumerate(["no_spec", "with_spec"]):
        ax = axes[ax_idx]
        cond_label = "Without Specification" if cond == "no_spec" else "With ClawBio Specification"

        for m in models:
            scores = []
            for pop in pops:
                vals = [get_score(r, "A1") for r in raw
                        if r.get("model") == m and r.get("pop") == pop
                        and get_condition(r) == cond and is_valid(r)]
                scores.append(np.mean(vals) if vals else 0)

            ax.plot(pop_labels, scores, 'o-', label=m,
                    color=COLORS.get(m, '#666'), alpha=0.7, linewidth=1.5, markersize=6)

        ax.set_ylabel('Phenotype Accuracy (A1)')
        ax.set_ylim(0, 1.15)
        ax.set_title(cond_label, fontweight='bold')
        ax.axhline(y=1.0, color='#CBD5E1', linestyle='--', linewidth=0.5)

        if ax_idx == 1:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)

    fig.suptitle('Figure 5. Phenotype Accuracy by Population Context',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    for ext in ['png', 'tiff']:
        plt.savefig(os.path.join(FIGURES_DIR, f'Figure5_population.{ext}'),
                    dpi=300, bbox_inches='tight')
    print("  Figure 5 saved")
    plt.close()


# ===== MAIN =====
if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    raw, data_path = load_data()
    models = get_models_present(raw)

    print(f"Generating figures from {len(raw)} results...")
    print(f"Models: {models}")
    print()

    figure2_clinical_correctness(raw, models)
    figure3_radar(raw, models)
    figure4_consistency(raw, models)
    figure5_summary(raw, models)

    print(f"\nAll figures saved to {FIGURES_DIR}")
