# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 20:05:44 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

data_paths = {
    "Ctrl_1": "C:/Users/senaa/Downloads/cr/omix007490/KC1/",
    "Ctrl_2": "C:/Users/senaa/Downloads/cr/omix007490/KC2/",
    "S_50nm_1": "C:/Users/senaa/Downloads/mp/omix007490/KS1/", 
    "S_50nm_2": "C:/Users/senaa/Downloads/mp/omix007490/KS2/",
    "M_500nm_1": "C:/Users/senaa/Downloads/mp/omix007490/KM1/", 
    "M_500nm_2": "C:/Users/senaa/Downloads/mp/omix007490/KM2/",
    "L_5um_1": "C:/Users/senaa/Downloads/mp/omix007490/KLARGE1/", 
    "L_5um_2": "C:/Users/senaa/Downloads/mp/omix007490/KLARGE2/"
}

adatas = {}
for name, path in data_paths.items():
    adatas[name] = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    adatas[name].obs['condition'] = name.split('_')[0]  
    adatas[name].obs['sample_id'] = name

adata = sc.concat(adatas.values(), label="sample_id", keys=adatas.keys())
adata.obs_names_make_unique()
adata.var_names_make_unique()

print(adata)


condition_map = {
    "Ctrl_1": "Control",
    "Ctrl_2": "Control",
    "S_50nm_1": "Nanoplastic",
    "S_50nm_2": "Nanoplastic",
    "M_500nm_1": "Nanoplastic",
    "M_500nm_2": "Nanoplastic",
    "L_5um_1": "Microplastic",
    "L_5um_2": "Microplastic"
}

for name, path in data_paths.items():
    adatas[name] = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        cache=True
    )

    adatas[name].obs["condition"] = condition_map[name]
    adatas[name].obs["sample_id"] = name



adata.var["mt"] = adata.var_names.str.startswith("mt-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.total_counts < 50000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :].copy()

print("Filtreleme sonrası veri boyutu:", adata.shape)

sc.pp.scrublet(adata, batch_key="sample_id")
adata = adata[~adata.obs['predicted_doublet']].copy()
print("Doublet filtreleme sonrası:", adata.shape)


adata.layers["counts"] = adata.X.copy()  

os.makedirs("processed", exist_ok=True)
adata.write_h5ad("processed/kidney_qc.h5ad")


sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, facecolor="white")


adata_ctrl = sc.read_10x_mtx("C:/Users/senaa/Downloads/cr/gse306246/")
adata_mp = sc.read_10x_mtx("C:/Users/senaa/Downloads/mp/gse306246/")


adata_ctrl.obs['condition'] = 'Control'
adata_mp.obs['condition'] = 'Microplastic'

adata = adata_ctrl.concatenate(adata_mp)
adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("mt-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.n_genes_by_counts > 2000, :].copy()
adata = adata[adata.obs.n_genes_by_counts < 10000, :].copy()
adata = adata[adata.obs.pct_counts_mt < 10, :].copy()
adata = adata[adata.obs.total_counts < 250000, :].copy()
print("Filtreleme sonrası veri boyutu:", adata.shape)


sc.pp.scrublet(adata, batch_key="condition")
adata = adata[~adata.obs['predicted_doublet']].copy()
print("Doublet filtreleme sonrası:", adata.shape)


adata.layers["counts"] = adata.X.copy()   

os.makedirs("processed", exist_ok=True)
adata.write_h5ad("processed/vessel_qc.h5ad")


data_paths = {
    "CR": "C:/Users/senaa/Downloads/cr/gse32/CR",
    "PE01": "C:/Users/senaa/Downloads/mp/gse32/PE0_1",
    "PET01": "C:/Users/senaa/Downloads/mp/gse32/PET0_1", 
    "PET1": "C:/Users/senaa/Downloads/mp/gse32/PET1_",
   
}

adatas = {}
for name, path in data_paths.items():
    adatas[name] = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    adatas[name].obs['condition'] = name.split('_')[0]  
    adatas[name].obs['sample_id'] = name

adata = sc.concat(adatas.values(), label="sample_id", keys=adatas.keys())
adata.obs_names_make_unique()
adata.var_names_make_unique()

print(adata)

condition_map = {
    "CR": "Control",
    "PE01": "Nanoplastic",
    "PET01": "Nanoplastic",
    "PET1": "Microplastic"
}

for name, path in data_paths.items():
    adatas[name] = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        cache=True
    )

    adatas[name].obs["condition"] = condition_map[name]
    adatas[name].obs["sample_id"] = name

adata.var["mt"] = adata.var_names.str.startswith("mt-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.n_genes_by_counts < 4000, :]
adata = adata[adata.obs.total_counts < 20000, :]
adata = adata[adata.obs.pct_counts_mt < 10, :].copy()

print("Filtreleme sonrası veri boyutu:", adata.shape)

sc.pp.scrublet(adata, batch_key="condition")
adata = adata[~adata.obs['predicted_doublet']].copy()
print("Doublet filtreleme sonrası:", adata.shape)


adata.layers["counts"] = adata.X.copy()   


os.makedirs("processed", exist_ok=True)
adata.write_h5ad("processed/lung_qc.h5ad")
