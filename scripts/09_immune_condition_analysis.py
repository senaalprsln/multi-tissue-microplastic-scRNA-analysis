# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:32:13 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata_immune = sc.read_h5ad(
    "processed/immune_cells_condition.h5ad"
)

cond_counts = pd.crosstab(adata_immune.obs["immune_celltype"], adata_immune.obs["condition_simple"])
print("\nCell type x condition_simple Distribution (Counts):")
print(cond_counts)

cond_pct = pd.crosstab(
    adata_immune.obs["immune_celltype"], adata_immune.obs["condition_simple"], normalize="index"
) * 100
print("\nCell type x condition_simple Distribution (Percentage):")
print(cond_pct.round(1))

fig, ax = plt.subplots(figsize=(10, 6))
cond_pct.plot(kind="bar", stacked=True, ax=ax, colormap="Set2")
plt.title("Condition Ratios Across Immune Cell Types")
plt.ylabel("Percentage (%)")
plt.xlabel("Immune Cell Type")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("immune_celltype_condition_barplot.png", dpi=300)
plt.show()

adata_immune.write(
    "processed/immune_cells_condition_analysis.h5ad"
)