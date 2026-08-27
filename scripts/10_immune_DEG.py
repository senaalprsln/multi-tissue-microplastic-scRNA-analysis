# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:32:44 2026

@author: senaa
"""

import scanpy as sc
import pandas as pd

adata_immune = sc.read_h5ad(
    "processed/immune_cells_condition.h5ad"
)

MIN_CELLS_PER_GROUP = 40
PADJ_THRESHOLD = 0.05

condition_col = "condition_simple"

immune_celltypes = (
    adata_immune.obs["immune_celltype"]
    .dropna()
    .unique()
    .tolist()
)

datasets = (
    adata_immune.obs["dataset"]
    .dropna()
    .unique()
    .tolist()
)


comparisons = [
    ("nano", "control"),
    ("micro", "control"),
    ("micro", "nano")
]

all_deg_results = []


for tissue in datasets:

    for celltype in immune_celltypes:

        
        subset = adata_immune[
            (adata_immune.obs["dataset"] == tissue) &
            (adata_immune.obs["immune_celltype"] == celltype)
        ].copy()

        if subset.shape[0] == 0:
            continue

        print("\n" + "=" * 70)
        print(f"DOKU: {tissue} | HÜCRE: {celltype}")
        print("=" * 70)

        
        condition_counts = (
            subset.obs[condition_col]
            .value_counts()
        )

        print("Koşul dağılımı:")
        print(condition_counts)

      
        for test_group, reference_group in comparisons:

            n_test = condition_counts.get(test_group, 0)
            n_reference = condition_counts.get(reference_group, 0)

            if (
                n_test < MIN_CELLS_PER_GROUP
                or n_reference < MIN_CELLS_PER_GROUP
            ):

                print(
                    f"[ATLANDI] "
                    f"{test_group} vs {reference_group} | "
                    f"test={n_test}, reference={n_reference}"
                )

                continue

            print(
                f"\n>>> DEG: {tissue} | {celltype} | "
                f"{test_group} vs {reference_group}"
            )

            print(
                f"n_{test_group} = {n_test} | "
                f"n_{reference_group} = {n_reference}"
            )

            sub_pair = subset[
                subset.obs[condition_col].isin(
                    [test_group, reference_group]
                )
            ].copy()

            
            try:

                sc.tl.rank_genes_groups(
                    sub_pair,
                    groupby=condition_col,
                    groups=[test_group],
                    reference=reference_group,
                    method="wilcoxon",
                    pts=True
                )

                de_df = sc.get.rank_genes_groups_df(
                    sub_pair,
                    group=test_group
                )

                
                de_df["tissue"] = tissue
                de_df["celltype"] = celltype
                de_df["comparison"] = (
                    f"{test_group}_vs_{reference_group}"
                )

                de_df["n_test"] = n_test
                de_df["n_reference"] = n_reference


                all_deg_results.append(de_df)

                
                n_sig = (
                    de_df["pvals_adj"] < PADJ_THRESHOLD
                ).sum()

                n_up = (
                    (de_df["pvals_adj"] < PADJ_THRESHOLD) &
                    (de_df["logfoldchanges"] > 0)
                ).sum()

                n_down = (
                    (de_df["pvals_adj"] < PADJ_THRESHOLD) &
                    (de_df["logfoldchanges"] < 0)
                ).sum()

                print(
                    f"Anlamlı DEG: {n_sig} | "
                    f"Up: {n_up} | "
                    f"Down: {n_down}"
                )

            except Exception as e:

                print(
                    f"HATA: "
                    f"{tissue} | {celltype} | "
                    f"{test_group} vs {reference_group}"
                )

                print(e)


if all_deg_results:

    all_deg_df = pd.concat(
        all_deg_results,
        ignore_index=True
    )

    print("\n" + "=" * 70)
    print("DEG ANALİZİ TAMAMLANDI")
    print("=" * 70)

    print(
        f"Toplam DEG satırı: "
        f"{len(all_deg_df)}"
    )

    print("\nKarşılaştırma dağılımı:")

    print(
        all_deg_df["comparison"]
        .value_counts()
    )

    
    all_deg_df.to_csv(
        "immune_DEG_all_comparisons.csv",
        index=False
    )

    print(
        "\nSonuçlar kaydedildi:"
        " immune_DEG_all_comparisons.csv"
    )

else:

    print("\nHiç DEG sonucu oluşturulamadı.")
    
adata_immune.write(
    "processed/immune_cells_deg.h5ad"
)