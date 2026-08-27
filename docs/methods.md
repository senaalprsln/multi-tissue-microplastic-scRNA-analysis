**Materials and Methods**

Data Acquisition and Preprocessing

Publicly available single-cell RNA sequencing (scRNA-seq) datasets investigating the effects of nano- and micro-scale plastic exposure on mouse tissues were retrieved from public repositories (NCBI GEO accessions GSE306246, GSE324602, and NGDC OMIX accession OMIX007490). Raw 10x Genomics count matrices for kidney, lung, and vessel tissues were processed using Scanpy (v1.10+). Per-tissue quality control (QC) was performed by filtering cells based on the number of expressed genes, total UMI counts, and mitochondrial gene percentage thresholds tailored to each tissue type. Potential cell doublets were identified and removed using Scrublet. Count matrices were subsequently normalized and log-transformed for exploratory analyses while raw counts were preserved in data layers for downstream integration and differential expression testing.



**Data Integration and Batch Correction**

To integrate multi-tissue datasets while mitigating batch effects without over-correcting biological variance, a selective batch-aware integration strategy was implemented using scVI (scvi-tools). A composite batch covariate was constructed where immune cells were grouped primarily by their tissue of origin, whereas non-tissue-immune populations were labeled using a combination of tissue type and coarse cell identity. This selective harmonization allowed immune compartments across different organs to align properly in a shared latent space while preserving transcriptionally distinct, tissue-specific structural and parenchymal cell states.



**Cell-Type Annotation and Subtyping**

Following latent space dimensionality reduction and nearest-neighbor graph construction, Leiden clustering was applied. Initial coarse cell-type identities and subsequent refined immune subsets (including T cells, NK/cytotoxic T cells, B cells, plasma cells, myeloid macrophages, dendritic antigen-presenting cells, and neutrophils) were annotated via marker-gene module scoring. Marker panels and cell-type definitions were standardized, and z-score-standardized scores were used to assign biological labels to clusters.



**Differential Expression Analysis**

Differential expression analysis (DEG) was conducted across harmonized cell populations within each tissue using the Wilcoxon rank-sum test. A comprehensive 3-way pairwise comparison framework was established to evaluate transcriptional shifts between exposure conditions: nanoplastic versus control, microplastic versus control, and microplastic versus nanoplastic. Genes with an adjusted p-value of less than 0.05 were considered statistically significant. Up- and down-regulated gene sets from these comparisons were structured to serve as direct inputs for downstream functional enrichment tools (e.g., g:Profiler).



**Biological Stress Module Scoring**

To quantify cellular stress responses, composite module scores were calculated across cells using predefined gene sets representing specific biological pathways, namely Apoptosis (e.g., Bax, Casp3, Trp53), Inflammation (e.g., Tnf, Il1b, Ccl2), and Oxidative Stress (e.g., Nfe2l2, Hmox1, Sod1). Average module scores were aggregated and compared across exposure groups (control, nano, micro) within each tissue type to assess localized pathway activation.

