# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 18:25:04 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad(
    "processed/annotated_all_tissues.h5ad"
)

print(adata)

sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon'
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False
)

marker_table = pd.DataFrame(
    adata.uns["rank_genes_groups"]["names"]
)

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

for celltype, genes in marker_sets.items():

    genes_present = [
        g for g in genes if g in adata.var_names
    ]

    if len(genes_present) == 0:
        print(f"{celltype}: marker genler bulunamadı.")
        continue

    print(f"{celltype} için bulunan genler:",
          genes_present)

    sc.pl.umap(
        adata,
        color=genes_present,
        title=[
            f"{celltype}: {g}"
            for g in genes_present
        ],
        frameon=False,
        size=20,
        ncols=3
    )

scores = {}
for celltype, genes in marker_sets.items():
    valid = [g for g in genes if g in adata.var_names]
    if valid:
        sc.tl.score_genes(adata, valid, score_name=f"score_{celltype}", use_raw=False)
        scores[celltype] = f"score_{celltype}"

score_cols = list(scores.values())
cluster_mean = adata.obs.groupby("leiden", observed=True)[score_cols].mean()

cluster_to_celltype = {}
for cluster in cluster_mean.index:
    best_col = cluster_mean.loc[cluster].idxmax()
    cluster_to_celltype[cluster] = best_col.replace("score_", "")

adata.obs["celltype"] = adata.obs["leiden"].map(cluster_to_celltype).astype("category")

adata.write(
    "processed/immune_cells_pre_scvi.h5ad"
)

