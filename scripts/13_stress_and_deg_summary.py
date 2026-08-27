# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 13:32:42 2026

@author: senaa
"""
import scanpy as sc
import pandas as pd

adata_immune = sc.read_h5ad(
    "processed/deg_summary.h5ad"
)

all_deg_df = pd.read_csv("immune_DEG_all_comparisons.csv")

print(f"{'='*70}")
print("ETKİ BÜYÜKLÜĞÜ SIRALAMASI (DEG sayısına göre)")
print(f"{'='*70}")

for df in all_deg_df:
    if df.empty:
        continue
    tissue = df["tissue"].iloc[0]
    celltype = df["celltype"].iloc[0]
    comparison = df["comparison"].iloc[0]
    n_sig = (df["pvals_adj"] < 0.05).sum()
    print(f"{tissue:8s} | {celltype:22s} | {comparison:20s} | {n_sig} anlamlı gen")
    

stress_markers = {
    "Apoptosis": ["Bax", "Casp3", "Casp8", "Bak1", "Trp53"],
    "Inflammation": ["Tnf", "Il1b", "Il6", "Nfkb1", "Ccl2"],
    "Oxidative_stress": ["Nfe2l2", "Hmox1", "Sod1", "Gpx1", "Cat"],
}

for pathway, genes in stress_markers.items():
    valid = [g for g in genes if g in adata_immune.var_names]
    if valid:
        sc.tl.score_genes(adata_immune, valid, score_name=f"module_{pathway}", use_raw=False)

for tissue in adata_immune.obs["dataset"].unique():
    subset = adata_immune[adata_immune.obs["dataset"] == tissue]
    print(f"\n{tissue}:")
    print(subset.obs.groupby("condition_simple", observed=True)[
        [f"module_{p}" for p in stress_markers]
    ].mean())
    

