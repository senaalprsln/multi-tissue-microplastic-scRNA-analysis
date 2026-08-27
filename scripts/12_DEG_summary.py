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

    anlamli = df[df["pvals_adj"] < PADJ_THRESHOLD].sort_values(by="pvals_adj", ascending=True)
    print(f"Toplam anlamlı gen sayısı: {anlamli.shape[0]}")

    down_tablosu = df[
        (df["logfoldchanges"] < 0) & (df["pvals_adj"] < PADJ_THRESHOLD)
    ].sort_values(by="logfoldchanges", ascending=True)
    down_genler = down_tablosu["names"].tolist()
    print(f"\nDownregüle gen sayısı: {len(down_genler)}")
    print(" ".join(down_genler))

    up_tablosu = df[
        (df["logfoldchanges"] > 0) & (df["pvals_adj"] < PADJ_THRESHOLD)
    ].sort_values(by="logfoldchanges", ascending=False)
    up_genler = up_tablosu["names"].tolist()
    print(f"\nUpregüle gen sayısı: {len(up_genler)}")
    print(" ".join(up_genler))
    
    
adata_immune.write(
    "processed/deg_summary.h5ad"
)