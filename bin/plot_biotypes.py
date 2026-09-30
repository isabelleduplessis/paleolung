#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

# ── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="featurecounts_summary.csv")
parser.add_argument("--output-prefix", default="biotype_composition")
parser.add_argument("--level", choices=["gene", "exon"], required=True,
                     help="Whether this summary is gene-level or exon-level counts")
args = parser.parse_args()

TITLES = {
    "exon": "FeatureCounts Biotype Composition - Exon Level and Repetitive Elements",
    "gene": "FeatureCounts Biotype Composition - Gene Level",
}

df = pd.read_csv(args.input)

# ── CONFIG ────────────────────────────────────────────────────────────────────
meta_cols = ["sample", "condition", "total_alignments", "assigned_alignments", "pct_assigned"]
label_cols = [c for c in df.columns if c not in meta_cols]

GROUP_MAP = {
    "protein_coding":                     "protein_coding",
    "Mt_protein_coding":                  "Mt_protein_coding",
    "Mt_rRNA":                            "Mt_rRNA",
    "Mt_tRNA":                            "Mt_tRNA",
    "lncRNA":                             "noncoding_RNA",
    "miRNA":                              "noncoding_RNA",
    "misc_RNA":                           "noncoding_RNA",
    "rRNA":                               "rRNA",
    "rRNA_pseudogene":                    "noncoding_RNA",
    "scaRNA":                             "noncoding_RNA",
    "scRNA":                              "noncoding_RNA",
    "snoRNA":                             "noncoding_RNA",
    "snRNA":                              "noncoding_RNA",
    "sRNA":                               "noncoding_RNA",
    "srpRNA":                             "noncoding_RNA",
    "vault_RNA":                          "noncoding_RNA",
    "ribozyme":                           "noncoding_RNA",
    "tRNA":                               "tRNA",
    "SINE":                               "SINEs",
    "SINE?":                              "SINEs",
    "LINE":                               "LINEs",
    "LTR":                                "LTRs",
    "LTR?":                               "LTRs",
    "DNA":                                "Other_repetitive",
    "DNA?":                               "Other_repetitive",
    "RC":                                 "Other_repetitive",
    "RC?":                                "Other_repetitive",
    "Retroposon":                         "Other_repetitive",
    "RNA":                                "Other_repetitive",
    "Satellite":                          "Other_repetitive",
    "Unknown":                            "Other_repetitive",
    "Low_complexity":                     "Other_repetitive",
    "Simple_repeat":                      "Other_repetitive",
    "artifact":                           "Other_repetitive",
    "processed_pseudogene":               "Other_assigned",
    "unprocessed_pseudogene":             "Other_assigned",
    "transcribed_processed_pseudogene":   "Other_assigned",
    "transcribed_unprocessed_pseudogene": "Other_assigned",
    "transcribed_unitary_pseudogene":     "Other_assigned",
    "translated_processed_pseudogene":    "Other_assigned",
    "translated_unprocessed_pseudogene":  "Other_assigned",
    "unitary_pseudogene":                 "Other_assigned",
    "polymorphic_pseudogene":             "Other_assigned",
    "pseudogene":                         "Other_assigned",
    "TEC":                                "Other_assigned",
    "IG_C_gene":                          "Other_assigned",
    "IG_C_pseudogene":                    "Other_assigned",
    "IG_D_gene":                          "Other_assigned",
    "IG_J_gene":                          "Other_assigned",
    "IG_J_pseudogene":                    "Other_assigned",
    "IG_pseudogene":                      "Other_assigned",
    "IG_V_gene":                          "Other_assigned",
    "IG_V_pseudogene":                    "Other_assigned",
    "TR_C_gene":                          "Other_assigned",
    "TR_D_gene":                          "Other_assigned",
    "TR_J_gene":                          "Other_assigned",
    "TR_J_pseudogene":                    "Other_assigned",
    "TR_V_gene":                          "Other_assigned",
    "TR_V_pseudogene":                    "Other_assigned",
}

COLORS = {
    "protein_coding":    "#6BB5D6",
    "Mt_protein_coding": "#1A7A2E",
    "noncoding_RNA":     "#ADD8E6",
    "tRNA":              "#3E7CB1",
    "Mt_tRNA":           "#E8C61E",
    "rRNA":              "#F28E2B", 
    "Mt_rRNA":           "#8C5A2B",
    "SINEs":             "#8B9A2D",
    "LINEs":             "#D4C97A",
    "LTRs":              "#E8897A",
    "Other_repetitive":  "#7B2D5E",
    "Other_assigned":    "#D46BB5",
    "unassigned":        "#2D2D8B",
}

STACK_ORDER = [
    "protein_coding", "Mt_protein_coding", "noncoding_RNA",
    "tRNA", "Mt_tRNA", "rRNA", "Mt_rRNA",
    "SINEs", "LINEs", "LTRs", "Other_repetitive", "Other_assigned",
    "unassigned",
]

legend_labels = {
    "protein_coding":    "Protein Coding",
    "Mt_protein_coding": "Mitochondrial Protein Coding",
    "noncoding_RNA":     "Noncoding RNA",
    "tRNA":              "tRNA",
    "Mt_tRNA":           "Mitochondrial tRNA",
    "rRNA":              "rRNA",
    "Mt_rRNA":           "Mitochondrial rRNA",
    "SINEs":             "SINEs",
    "LINEs":             "LINEs",
    "LTRs":              "LTRs",
    "Other_repetitive":  "Other Repetitive",
    "Other_assigned":    "Other Assigned",
    "unassigned":        "Unassigned",
}

# ── PREPARE DATA ──────────────────────────────────────────────────────────────
plot_df = df.copy()

group_data = {}
for col in label_cols:
    group = GROUP_MAP.get(col, col)
    group_data[group] = group_data.get(group, 0) + plot_df[col].values

group_df = pd.DataFrame(group_data, index=plot_df.index)
group_df["unassigned"] = plot_df["total_alignments"] - plot_df["assigned_alignments"]
group_df["sample"] = plot_df["sample"]
group_df["total_alignments"] = plot_df["total_alignments"]

biotype_cols = [c for c in group_df.columns if c not in ["sample", "total_alignments"]]
for col in biotype_cols:
    group_df[col] = group_df[col] / group_df["total_alignments"] * 100

STACK_ORDER = [c for c in STACK_ORDER if c in group_df.columns]

gdf = group_df.set_index("sample")

# ── FORMAT HELPER ─────────────────────────────────────────────────────────────
def fmt_reads(n):
    return f"{int(n):,}"

# ── PLOT ──────────────────────────────────────────────────────────────────────
samples = list(gdf.index)
x_positions = np.arange(len(samples))

fig_width = max(6, len(samples) * 0.9)
fig, ax = plt.subplots(figsize=(fig_width, 5))

for i, sample in enumerate(samples):
    x = x_positions[i]
    bottom = 0
    width = 0.7

    for biotype in STACK_ORDER:
        val = gdf.loc[sample, biotype]
        ax.bar(x, val, bottom=bottom,
               color=COLORS.get(biotype, "#AAAAAA"),
               width=width,
               edgecolor="none")
        bottom += val

    total_reads = gdf.loc[sample, "total_alignments"]
    ax.text(x, 101.5, fmt_reads(total_reads), ha="center", va="bottom",
            fontsize=8, fontweight="bold", color="black")

ax.set_xticks(x_positions)
ax.set_xticklabels(samples, fontsize=9, rotation=45, ha="right")

ax.set_ylim(0, 108)
ax.set_ylabel("Percentage (%) of Total Aligned Reads", fontsize=11)
ax.set_title(TITLES[args.level], fontsize=12, pad=25)
ax.set_xlim(x_positions[0] - 0.6, x_positions[-1] + 0.6)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

handles = [
    mpatches.Patch(color=COLORS[b], label=legend_labels[b])
    for b in STACK_ORDER if b in COLORS
]
ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc="upper left",
          fontsize=9, frameon=False, title_fontsize=10)

plt.tight_layout()
plt.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight")
plt.savefig(f"{args.output_prefix}.png", dpi=150, bbox_inches="tight")
print(f"Saved: {args.output_prefix}.pdf / .png")