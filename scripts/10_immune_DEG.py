# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:32:44 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd

adata_immune = sc.read_h5ad(
    "processed/immune_cells_condition.h5ad"
)

MIN_CELLS_PER_GROUP = 40
PADJ_THRESHOLD = 0.05

condition_col = "condition_simple"
immune_celltypes = adata_immune.obs["immune_celltype"].unique().tolist()
datasets = adata_immune.obs["dataset"].unique().tolist()

all_deg_results = []

for tissue in datasets:
    for celltype in immune_celltypes:

        subset = adata_immune[
            (adata_immune.obs["dataset"] == tissue) &
            (adata_immune.obs["immune_celltype"] == celltype)
        ].copy()

        if subset.shape[0] == 0:
            continue

        available_conditions = subset.obs[condition_col].value_counts()
        valid_groups = available_conditions[available_conditions >= MIN_CELLS_PER_GROUP].index.tolist()

        if "control" not in valid_groups or len(valid_groups) < 2:
            print(f"[ATLANDI] {tissue} - {celltype}: yetersiz grup/hücre "
                  f"({dict(available_conditions)})")
            continue

        test_groups = [g for g in valid_groups if g != "control"]

        for test_group in test_groups:
            sub_pair = subset[subset.obs[condition_col].isin(["control", test_group])].copy()

            print(f"\n>>> DEG: {tissue} | {celltype} | control vs {test_group} "
                  f"(n_control={sum(sub_pair.obs[condition_col]=='control')}, "
                  f"n_{test_group}={sum(sub_pair.obs[condition_col]==test_group)})")

            try:
                sc.tl.rank_genes_groups(
                    sub_pair,
                    groupby=condition_col,
                    groups=[test_group],
                    reference="control",
                    method="wilcoxon",
                    pts=True
                )
                de_df = sc.get.rank_genes_groups_df(sub_pair, group=test_group)
                de_df["tissue"] = tissue
                de_df["celltype"] = celltype
                de_df["comparison"] = f"{test_group}_vs_control"
                all_deg_results.append(de_df)

            except Exception as e:
                print(f"HATA ({tissue}, {celltype}, {test_group}): {e}")

if all_deg_results:
    all_deg_df = pd.concat(all_deg_results, ignore_index=True)
    all_deg_df.to_csv("processed/all_deg_results.csv", index=False)
    print("Tüm DEG sonuçları processed/all_deg_results.csv olarak kaydedildi.")
    
adata_immune.write(
    "processed/immune_cells_deg.h5ad"
)