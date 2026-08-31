# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:35:04 2026

@author: senaa
"""
import scanpy as sc
import pandas as pd


adata_immune = sc.read_h5ad(
    "processed/immune_subtype_annotation.h5ad"
)

PADJ_THRESHOLD = 0.05
all_deg_df = pd.read_csv("immune_DEG_all_comparisons.csv")

for df in all_deg_df:
    if df.empty:
        continue

    tissue = df["tissue"].iloc[0]
    celltype = df["celltype"].iloc[0]
    comparison = df["comparison"].iloc[0]
    key_base = f"{tissue}_{celltype}_{comparison}"

    print(f"\n{'='*60}")
    print(f"{key_base}")
    print(f"{'='*60}")

    significant = df[df["pvals_adj"] < PADJ_THRESHOLD].sort_values(by="pvals_adj", ascending=True)
    print(f"Total significant gene count: {significant.shape[0]}")

    down_table = df[
        (df["logfoldchanges"] < 0) & (df["pvals_adj"] < PADJ_THRESHOLD)
    ].sort_values(by="logfoldchanges", ascending=True)
    down_genes = down_table["names"].tolist()
    print(f"\nDownregulated gene count: {len(down_genes)}")
    print(" ".join(down_genes))

    up_table = df[
        (df["logfoldchanges"] > 0) & (df["pvals_adj"] < PADJ_THRESHOLD)
    ].sort_values(by="logfoldchanges", ascending=False)
    up_genes = up_table["names"].tolist()
    print(f"\nUpregulated gene count: {len(up_genes)}")
    print(" ".join(up_genes))
    
    
adata_immune.write(
    "processed/deg_summary.h5ad"
)