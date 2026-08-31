# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 19:19:02 2026

@author: senaa
"""
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import scvi

adata = sc.read_h5ad(
    "processed/immune_cells_scvi.h5ad"
)

print(adata)
print(adata.obsm.keys())

immune_classes = [
    "T_Cells", "NK_Cytotoxic_T", "B_Cells",  
    "Plasma_Cells", "Myeloid_Macrophages", "Dendritic_APC", "Neutrophil"
]

model_sel = scvi.model.SCVI.load("scvi_model_selective", adata=adata)
adata.obsm["X_scVI_selective"] = model_sel.get_latent_representation()
 
sc.pp.neighbors(adata, use_rep="X_scVI_selective", key_added="neighbors_selective")
sc.tl.umap(adata, neighbors_key="neighbors_selective")
adata.obsm["X_umap_selective"] = adata.obsm["X_umap"].copy()
 
print("\n=== Selective integration UMAP: immune mixed, others tissue-specific ===")
sc.pl.embedding(
    adata, basis="X_umap_selective",
    color=["celltype", "dataset"],
    frameon=False, wspace=0.4
)
plt.show()
 
sc.pp.neighbors(adata, use_rep="X_scVI_selective", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.3, key_added="leiden_selective")

print("\n=== Plotting Single Global UMAP (No separate file/plot) ===")
sc.pl.umap(
    adata,
    color=["celltype", "dataset"],  
    frameon=False,
    title="All Cells in One Plot (Only Immune Cells Integrated)"
)

adata_immune = adata[adata.obs["celltype"].isin(immune_classes)].copy()

print(f"\nTotal cells: {adata.shape[0]}")
print(f"Immune cells: {adata_immune.shape[0]}")
print(f"Ratio: {adata_immune.shape[0] / adata.shape[0] * 100:.1f}%")

print("\nImmune cell type distribution:")
print(adata_immune.obs["celltype"].value_counts())

print("\nImmune cell distribution by tissue:")
print(pd.crosstab(adata_immune.obs["celltype"], adata_immune.obs["dataset"]))

if "log1p" in adata_immune.uns:
    del adata_immune.uns["log1p"]
    
adata_immune.write(
    "processed/immune_cells_clustered.h5ad"
)