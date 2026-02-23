---
name: single-cell-rna-seq
description: End-to-end single-cell RNA-seq analysis including data loading, QC, integration, clustering, cell type annotation, and differential expression.
metadata:
    skill-author: Albert Ying
---

# Single-cell RNA-seq Analysis

## When to use

- Processing 10X Genomics scRNA-seq data
- Multi-sample integration and batch correction
- Cell type identification using marker databases
- Condition-specific differential expression within cell types

## Standard workflow

1. **Data loading**: Read 10X MTX format (barcodes, features, matrix)
2. **QC filtering**: Remove low-quality cells (low gene count, high mito%)
3. **Normalization**: Library size normalization + log transformation
4. **Integration**: Batch correction across samples if needed
5. **Clustering**: PCA → neighbors → Leiden clustering
6. **Cell type annotation**: Use known marker genes or reference databases
7. **Differential expression**: Compare conditions within each cell type

## Cell type annotation strategy

```python
import pandas as pd

# Load marker database (e.g., CellMarker)
markers_db = pd.read_excel('Cell_marker_Seq.xlsx')

# Get cluster markers
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# For each cluster, compare top markers against database
# Assign cell type based on best marker overlap
```

## Differential expression between conditions

```python
# For each cell type cluster
for cluster in adata.obs['leiden'].unique():
    mask = adata.obs['leiden'] == cluster
    sub = adata[mask].copy()
    # Compare pre vs post exercise
    sc.tl.rank_genes_groups(sub, 'condition', method='wilcoxon',
                            groups=['post'], reference='pre')
    result = sc.get.rank_genes_groups_df(sub, group='post')
```

## Key decisions

- Use raw counts for DE analysis (not scaled data)
- Apply Benjamini-Hochberg correction for multiple testing
- Report both significant and non-significant results per cluster
