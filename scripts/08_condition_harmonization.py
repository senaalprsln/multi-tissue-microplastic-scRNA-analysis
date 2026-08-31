# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:05:17 2026

@author: senaa
"""
import scanpy as sc
import pandas as pd

adata_immune = sc.read_h5ad(
    "processed/immune_cells_annotated.h5ad"
)

condition_map = {
    "Ctrl": "control", "S": "nano", "M": "nano", "L": "micro",
    "CR": "control", "PE01": "nano", "PET01": "nano", "PET1": "micro",
    "Control": "control", "Microplastic": "micro",
}
adata_immune.obs["condition_simple"] = adata_immune.obs["condition"].map(condition_map)

n_missing = adata_immune.obs["condition_simple"].isna().sum()
if n_missing > 0:
    print(f"WARNING: {n_missing} cells did not match!")
    print(adata_immune.obs.loc[adata_immune.obs["condition_simple"].isna(), "condition"].unique())

print(pd.crosstab(adata_immune.obs["dataset"], adata_immune.obs["condition_simple"]))

adata_immune.write(
    "processed/immune_cells_condition.h5ad"
)