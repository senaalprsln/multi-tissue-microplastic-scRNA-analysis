# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 19:17:34 2026

@author: senaa
"""

import scanpy as sc
import scvi

adata = sc.read_h5ad(
    "processed/immune_cells_pre_scvi.h5ad"
)

print(adata)


sc.pp.highly_variable_genes(
    adata,
    layer="counts",
    flavor="seurat_v3",
    batch_key="dataset",
    n_top_genes=3000,
    subset=False
)
adata_hvg = adata[:, adata.var.highly_variable].copy()

immune_classes = [
    "T_Cells", "NK_Cytotoxic_T", "B_Cells",  
    "Plasma_Cells", "Myeloid_Macrophages", "Dendritic_APC", "Neutrophil"
]

adata_hvg.obs["celltype"] = adata.obs["celltype"].values
dataset_col = adata_hvg.obs["dataset"].astype(str)


adata_hvg.obs["batch_for_correction"] = dataset_col.values


non_immune_mask = ~adata_hvg.obs["celltype"].isin(immune_classes)
adata_hvg.obs.loc[non_immune_mask, "batch_for_correction"] = (
    dataset_col[non_immune_mask] + "_" + adata_hvg.obs.loc[non_immune_mask, "celltype"].astype(str)
)

print("\n'batch_for_correction' distribution:")
print(adata_hvg.obs["batch_for_correction"].value_counts())

model_path_selective = "scvi_model_selective"

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="counts",
    batch_key="batch_for_correction"
)

model_sel = scvi.model.SCVI(
    adata_hvg,
    n_layers=2,
    n_latent=30
)

model_sel.train()

model_sel.save(
    "scvi_model_selective",
    overwrite=True
)

adata_hvg.obsm["X_scVI"] = model_sel.get_latent_representation()

print("Selective scVI trained and saved.")

adata.write(
    "processed/immune_cells_scvi.h5ad"
)