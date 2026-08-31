# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:34:50 2026

@author: senaa
"""
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata_immune = sc.read_h5ad(
    "processed/immune_cells_condition.h5ad"
)

subtype_marker_dict = {
    "T_Cells": {
        "CD4_Helper": ["Cd4", "Cd5", "Cxcr3", "Lef1"],
        "Naive_T": ["Ccr7", "Klf2", "Sell", "Tcf7"],
        "CD8_Cytotoxic": ["Cd8a", "Cd8b1", "Gzmb", "Nkg7"],
        "Effector_Memory_T": ["Ccl5", "Ccl4", "Gzmk"],
        "Regulatory_T": ["Foxp3", "Ctla4", "Il2ra"]
    },
    "Myeloid_Macrophages": {
        "Resident_Mac": ["Cx3cr1", "Mrc1", "Folr2"],
        "Inflammatory_Mac": ["Ccl2", "Il1b", "Tnf", "S100a8"],
        "Alveolar_Mac": ["Pparg", "Trem2", "Gpnmb"],
        "Monocyte_like": ["Plac8", "Ly6c2", "Cd14"]
    },
    "NK_Cytotoxic_T": {
        "NK_Cells": ["Nkg7", "Gzmb", "Klra8", "Ncam1"],
        "Cytotoxic_Effector": ["Prf1", "Ifng", "Gzmb"],
        "NKT_Cells": ["Cd3d", "Klrb1c"]
    },
    "Dendritic_APC": {
        "Conventional_DC1": ["H2-Ab1", "Cd74", "Xcr1", "Clec9a"],
        "Conventional_DC2": ["Cd1c", "Clec10a"],
        "Migratory_DC": ["Ccr7", "Fscn1", "Lamp3"]
    },
    "Plasma_Cells": {
        "Secretory_Plasma": ["Mzb1", "Jchain", "Xbp1"],
        "Mature_Plasma": ["Sdc1", "Cd38", "Cercam"],
        "Proliferating_Plasmablast": ["Mki67", "Top2a"]
    },
    "B_Cells": {
        "Naive_B": ["Cd19", "Ms4a1", "Cd79a", "Ighd"],
        "Memory_B": ["Cd27", "Cd19", "Ms4a1", "Tnfrsf13b"],
        "Germinal_Center_B": ["Aicda", "Bcl6"],
        "Age_Associated_B": ["Itgax", "Tbx21", "Fcrl2"]
    },
    "Neutrophil": {
        "Mature_Neutrophil": ["S100a8", "S100a9", "Cxcr2", "Retnlg"],
        "Activated_Neutrophil": ["Il1b", "Csf3r", "S100a12"],
        "Immature_Neutrophil": ["Mpo", "Elane"]
    }
}

target_cell_types = list(subtype_marker_dict.keys())

for ct in target_cell_types:
    print(f"\n==========================================")
    print(f"Processing & Subtyping: {ct}")
    print(f"==========================================")
    
    sub_adata = adata_immune[adata_immune.obs["immune_celltype"] == ct].copy()
    if sub_adata.shape[0] < 30: 
        print(f"Insufficient cells ({sub_adata.shape[0]}), atlanıyor.")
        continue
        
    if "log1p" in sub_adata.uns:
        del sub_adata.uns["log1p"]
        
    sc.pp.neighbors(sub_adata, use_rep="X_scVI_selective", n_neighbors=15)
    sc.tl.umap(sub_adata)
    leiden_key = f"{ct}_leiden"
    sc.tl.leiden(sub_adata, resolution=0.4, key_added=leiden_key)
    
    markers = subtype_marker_dict[ct]
    score_cols = []
    for subtype, gene_list in markers.items():
        valid_genes = [g for g in gene_list if g in sub_adata.var_names]
        if valid_genes:
            sc.tl.score_genes(sub_adata, valid_genes, score_name=f"score_{subtype}", use_raw=False)
            score_cols.append(f"score_{subtype}")
            
    if score_cols:
        cluster_mean = sub_adata.obs.groupby(leiden_key, observed=True)[score_cols].mean()
        cluster_to_sub = {}
        for cluster in cluster_mean.index:
            best_col = cluster_mean.loc[cluster].idxmax()
            best_sub = best_col.replace("score_", "")
            cluster_to_sub[cluster] = best_sub
            
        sub_adata.obs[f"{ct}_subtype"] = sub_adata.obs[leiden_key].map(cluster_to_sub).astype("category")
        
        
        plt.figure(figsize=(18, 5))
        sc.pl.umap(
            sub_adata,
            color=[leiden_key, f"{ct}_subtype", f"{ct}_subtype"],
            frameon=False,
            ncols=3,
            size=35,
            show=True,
            title=[f"{ct} Leiden", f"{ct} Automatic Subtypes", f"{ct} Subtypes (Distribution)"]
        )
    else:
        print(f"No valid marker gene match found for {ct}.")
        
adata_immune.write(
    "processed/immune_subtype_annotation.h5ad"
)
