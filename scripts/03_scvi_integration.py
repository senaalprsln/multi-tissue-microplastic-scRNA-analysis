# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:11:58 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("processed/all_tissues_combined.h5ad")

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

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=21)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=['leiden'])

adata.write("processed/annotated_all_tissues.h5ad")