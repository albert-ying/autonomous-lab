---
name: differential-expression
description: Differential gene expression analysis with PyDESeq2. Design matrices, contrasts, multiple testing correction, volcano plots.
metadata:
    skill-author: Albert Ying
---

# Differential expression

## When to use

- Bulk RNA-seq differential expression
- Multi-factor experimental designs
- Batch effect correction in DE analysis
- Generating volcano and MA plots

## Standard PyDESeq2 workflow

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd

# counts: genes x samples DataFrame; metadata: samples DataFrame
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata_df,
    design="~batch + condition",  # batch correction built into design
)
dds.deseq2()

# Extract results for a contrast
stat_res = DeseqStats(dds, contrast=["condition", "treated", "control"])
stat_res.summary()
results_df = stat_res.results_df

# Filter significant genes
sig = results_df[
    (results_df["padj"] < 0.05) & (results_df["log2FoldChange"].abs() > 1)
]
```

## Volcano plot

```python
import matplotlib.pyplot as plt
import numpy as np

df = results_df.copy()
df["-log10p"] = -np.log10(df["padj"].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(8, 6))
colors = np.where(
    (df["padj"] < 0.05) & (df["log2FoldChange"] > 1), "red",
    np.where((df["padj"] < 0.05) & (df["log2FoldChange"] < -1), "blue", "gray")
)
ax.scatter(df["log2FoldChange"], df["-log10p"], c=colors, s=8, alpha=0.5)
ax.axhline(-np.log10(0.05), ls="--", color="gray", lw=0.8)
ax.axvline(-1, ls="--", color="gray", lw=0.8)
ax.axvline(1, ls="--", color="gray", lw=0.8)
ax.set_xlabel("log2 Fold Change")
ax.set_ylabel("-log10 adjusted p-value")
```

## Design matrix tips

- Always include batch as a covariate if samples were processed in batches.
- Put the variable of interest last in the formula.
- For paired designs: `~patient + treatment`.
- For interaction: `~genotype + treatment + genotype:treatment`.
