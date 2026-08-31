# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:10:59 2026

@author: senaa
"""

import scanpy as sc
import scvi
import anndata as ad
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

lung   = sc.read_h5ad("processed/lung_qc.h5ad")
kidney = sc.read_h5ad("processed/kidney_qc.h5ad")
vessel = sc.read_h5ad("processed/vessel_qc.h5ad")

lung.obs["dataset"]   = "lung"
kidney.obs["dataset"] = "kidney"
vessel.obs["dataset"] = "vessel"

lung.obs["tissue"]   = "lung"
kidney.obs["tissue"] = "kidney"
vessel.obs["tissue"] = "vessel"

adata = ad.concat(
    [lung, kidney, vessel],
    join="inner",
    merge="same",
    index_unique="-"
)

print(adata)
print(adata.obs["dataset"].value_counts())
print(adata.obs["condition"].value_counts())

adata.raw = adata.copy()

os.makedirs("processed", exist_ok=True)

adata.write("processed/all_tissues_combined.h5ad")

print("\nSaved:")
print("processed/all_tissues_combined.h5ad")