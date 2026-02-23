---
name: scanpy
description: Single-cell RNA-seq analysis with Scanpy. QC, normalization, clustering, visualization, marker gene identification.
metadata:
    skill-author: Albert Ying
---

# Scanpy

## When to use

- Single-cell RNA-seq data preprocessing and QC
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Cell clustering (Leiden, Louvain)
- Marker gene identification
- Cell type annotation

## Standard workflow

```python
import scanpy as sc
import anndata as ad

# Read 10X data
adata = sc.read_10x_mtx('path/to/data/', var_names='gene_symbols')

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20].copy()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

# Scale, PCA, neighbors, clustering
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

# Marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
```

## Key decisions

- Filter cells with < 200 genes and > 20% mitochondrial reads
- Use batch-aware HVG selection when integrating multiple samples
- Leiden resolution controls cluster granularity (0.3-1.0 typical)
- Use Wilcoxon rank-sum test for marker gene identification
