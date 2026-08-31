# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 21:37:24 2026

@author: senaa
"""
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

adata_immune = sc.read_h5ad(
    "processed/immune_cells_clustered.h5ad"
)

immune_classes = [
    "T_Cells", "NK_Cytotoxic_T", "B_Cells",  
    "Plasma_Cells", "Myeloid_Macrophages", "Dendritic_APC", "Neutrophil"
]

marker_sets = {
   
    "T_Cells":             ["Cd3d", "Cd3e", "Trac", "Il7r", "Lck", "Cd2", "Cd247"],
    "NK_Cytotoxic_T":      ["Nkg7", "Gnly", "Gzmb", "Gzma", "Ccl5", "Klrd1", "Klrk1", "Prf1"],
    "B_Cells":             ["Ms4a1", "Cd79a", "Cd79b"],
    "Plasma_Cells":        ["Mzb1", "Jchain", "Igkc"],
    "Myeloid_Macrophages": ["Lyz2", "Cd68", "Aif1", "C1qa", "C1qb", "Adgre1", "Csf1r", "Tyrobp", "Fcgr3"],
    "Dendritic_APC":       ["H2-Ab1", "H2-Aa", "Cd74", "Fcer1a"],
    "Neutrophil":          ["S100a8", "S100a9", "Cxcr2", "Retnlg"],

    
    "Endothelial":         ["Pecam1", "Vwf", "Kdr", "Cldn5", "Egfl7"],
    "Fibroblasts":         ["Col1a1", "Col1a2", "Dcn", "Lum", "Col3a1"],
    "Schwann_Neuronal":    ["Kcna1", "Sox10", "Erbb3", "Mpz", "Plp1"],

    
    "PT":                  ["Slc27a2", "Lrp2", "Gatm", "Slc5a12", "Slc22a6"],
    "TAL":                 ["Umod", "Slc12a1", "Cldn10"],
    "CNT":                 ["Calb1", "Slc8a1"],
    "DCT":                 ["Slc12a3", "Pvalb"],
    "CD_PC":               ["Aqp2", "Aqp3", "Hsd11b2"],
    "CD_IC":               ["Atp6v1b1", "Atp6v0d2", "Foxi1"],
    "Podocyte":            ["Nphs1", "Nphs2", "Podxl", "Synpo"],
    "Thin_Limb":           ["Bst1"],

    
    "Lung_AT1":            ["Ager", "Hopx", "Aqp5", "Cldn18"],
    "Lung_AT2":            ["Sftpb", "Sftpc", "Sftpd", "Slc34a2"],
    "Lung_Club":           ["Scgb1a1", "Cyp2f2"],
    "Lung_Ciliated":       ["Foxj1", "Dnai1", "Tppp3"],

    
    "SMC_Contractile":     ["Acta2", "Myh11", "Tagln", "Cnn1"],
    "Adventitial_Fib":     ["Pdgfra", "Dpt", "Pi16", "Cd34"],
    "Pericyte":            ["Rgs5", "Kcnj8", "Ndufa4l2"],
    "Mesothelial":         ["Msln", "Upk3b"],   
    "PVAT_Adipocyte":      ["Adipoq", "Plin1", "Scd1"],
}

sc.pp.neighbors(adata_immune, use_rep="X_scVI_selective", n_neighbors=15)
sc.tl.umap(adata_immune)
sc.tl.leiden(adata_immune, resolution=0.4, key_added="immune_leiden")

print("\n=== Step 1 (immune subset): Raw cluster numbers ===")
sc.pl.umap(
    adata_immune,
    color=["immune_leiden", "condition", "dataset"],
    frameon=False
)

plt.show()

cluster_dataset_comp = pd.crosstab(adata_immune.obs["immune_leiden"], adata_immune.obs["dataset"])
print(cluster_dataset_comp)

cluster_dataset_pct = pd.crosstab(
    adata_immune.obs["immune_leiden"], adata_immune.obs["dataset"], normalize="index"
) * 100
print(cluster_dataset_pct.round(1))

if adata_immune.X.max() > 50:
    sc.pp.normalize_total(adata_immune, target_sum=1e4)
    sc.pp.log1p(adata_immune)

immune_only_marker_sets = {
    k: v for k, v in marker_sets.items()
    if k in immune_classes
}

for celltype, genes in immune_only_marker_sets.items():
    missing = [g for g in genes if g not in adata_immune.var_names]
    if missing:
        print(f"WARNING - {celltype}: these genes were not found -> {missing}")

immune_scores = {}
for celltype, genes in immune_only_marker_sets.items():
    valid = [g for g in genes if g in adata_immune.var_names]
    if valid:
        sc.tl.score_genes(adata_immune, valid, score_name=f"score_{celltype}", use_raw=False)
        immune_scores[celltype] = f"score_{celltype}"

immune_score_cols = list(immune_scores.values())
immune_cluster_mean = adata_immune.obs.groupby("immune_leiden", observed=True)[immune_score_cols].mean()
print("\nMarker score means for immune_leiden clusters:")
print(immune_cluster_mean)

immune_cluster_mean_z = (immune_cluster_mean - immune_cluster_mean.mean(axis=0)) / immune_cluster_mean.std(axis=0)

immune_cluster_to_celltype = {}
for cluster in immune_cluster_mean_z.index:
    best_col = immune_cluster_mean_z.loc[cluster].idxmax()
    immune_cluster_to_celltype[cluster] = best_col.replace("score_", "")

print("\n=== immune_leiden -> Cell Type mapping ===")
for cluster, celltype in immune_cluster_to_celltype.items():
    print(f"Küme {cluster} ==> {celltype}")

adata_immune.obs["immune_celltype"] = (
    adata_immune.obs["immune_leiden"].map(immune_cluster_to_celltype).astype("category")
)

print("\nImmune cell type distribution:")
print(adata_immune.obs["immune_celltype"].value_counts())

print("\n=== Step 2 (immune subset): Named cell types ===")
sc.pl.umap(
    adata_immune,
    color=["immune_leiden", "immune_celltype"],
    legend_loc="on data",
    frameon=False,
    size=25
)

plt.show()

adata_immune.write(
    "processed/immune_cells_annotated.h5ad"
)