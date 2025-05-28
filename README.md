# CMML3_ICA2
# 🧬 CMML3_ICA2: Benchmarking scAtlasVAE for Single-cell Integration of Bone Marrow Datasets

This project benchmarks the performance of **scAtlasVAE**, **scVI**, and **Harmony** for single-cell RNA-seq data integration. Using bone marrow single-cell datasets with batch effects, we evaluate how well each method removes batch artifacts while preserving cell type structure.

---

## 📁 Repository Structure
```
CMML3_ICA2 /
├── Data \ # files are too large, which can be found in https://figshare.com/articles/MCA_DGE_Data/5435866
├── Scripts / # Preprocessing, integration, and evaluation scripts
├── Figure / # UMAP visualizations and metric summary plots
└── README.md # Project documentation
```
---

## 🎯 Project Goal

To compare three popular single-cell integration tools—**scAtlasVAE**, **scVI**, and **Harmony**—in terms of their ability to:
- Correct batch effects across multiple bone marrow datasets
- Retain meaningful biological differences
- Perform well across standard benchmarking metrics

---

## 🧪 Data & Preprocessing

**Input data**: Publicly available bone marrow DGE matrices (e.g. `BoneMarrow1_dge.txt.gz`, etc.)

**Steps:**
- Read and transpose DGE matrices
- Annotate batch metadata from cell barcodes
- Normalize and log-transform counts
- Identify highly variable genes (HVGs)
- Apply PCA for downstream integration (Harmony)

Implemented using **Scanpy**, wrapped into reproducible `Scripts/preprocess/` functions.

---

## 🧬 Integration Methods

Three integration approaches were applied to the same datasets:
```
| Method        | Description                                                                 |
|---------------|-----------------------------------------------------------------------------|
| `Harmony`     | PCA-based correction via `scanpy.external.pp.harmony_integrate`             |
| `scVI`        | Probabilistic deep generative model using `scvi-tools`                      |
| `scAtlasVAE`  | VAE-based method trained on raw counts; extracted latent space from encoder |

Scripts and settings for each method are available in `Scripts/integration/`.
```
---

## 📊 Evaluation Metrics

Integration outputs were evaluated using [scIB](https://github.com/theislab/scib):
```
| Metric                  | Meaning                                        |
|-------------------------|------------------------------------------------|
| **Silhouette (celltype)** | Cell type conservation (higher is better)     |
| **ARI / NMI**            | Clustering agreement with known labels         |
| **graph iLISI**          | Batch mixing assessment (ideal value = 1.0)    |
```
Results are stored as `.csv` tables and visualized as bar plots in `Figure/`.

---

## 📈 Example Results

- **Figure 1**: UMAP visualizations after integration by each method (colored by batch and cell type)
- **Figure 2**: Bar plots showing benchmark metrics (Silhouette, ARI, NMI, graph iLISI)

All images are available in the `Figure/` folder.

---

