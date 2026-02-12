#!/usr/bin/env python3
"""Simple mock project setup — no heavy dependencies."""
import json
import os
import struct
import zlib
from pathlib import Path

DEMO = Path("/tmp/autolab_demo_full")

def tiny_png(w=4, h=4, r=80, g=120, b=180):
    """Create a minimal valid PNG."""
    raw = b""
    for y in range(h):
        raw += b"\x00"
        for x in range(w):
            raw += bytes([r, g, b])
    def chunk(t, d):
        c = t + d
        return struct.pack(">I", len(d)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw))
            + chunk(b"IEND", b""))

def main():
    import shutil
    if DEMO.exists():
        shutil.rmtree(DEMO)

    # Dirs
    for d in [".autolab/profiles", "data", "scripts", "figures", "results",
              "paper/sections", "paper/figures"]:
        (DEMO / d).mkdir(parents=True, exist_ok=True)

    # .autolab/idea.md
    (DEMO / ".autolab/idea.md").write_text("""# scRNA-seq Analysis of NSCLC Tumor Microenvironment

10x Genomics scRNA-seq from 12 NSCLC patients (6 responders, 6 non-responders to anti-PD-1).
Goal: identify immune populations and gene programs predicting immunotherapy response.

Data: data/counts_matrix.h5ad (45,000 cells), data/clinical_metadata.csv
Target journal: Nature Medicine
""")

    # Profiles (plain text, yaml-ish)
    (DEMO / ".autolab/profiles/pi.yaml").write_text(
        "title: Principal Investigator\nexpertise: Cancer immunology\n"
        "goal: Publish in Nature Medicine\npersonality:\n- Visionary\n- Critical\n- Publication-focused\n")
    (DEMO / ".autolab/profiles/trainee.yaml").write_text(
        "title: Postdoctoral Researcher\nexpertise: Single-cell genomics, scanpy\n"
        "goal: Execute with technical excellence\npersonality:\n- Dedicated\n- Technical\n- Proactive\n")
    (DEMO / ".autolab/config.yaml").write_text("max_iterations: 50\ntarget_journal: Nature Medicine\n")

    # Data files
    (DEMO / "data/counts_matrix.h5ad").write_bytes(b"mock")
    (DEMO / "data/clinical_metadata.csv").write_text(
        "patient_id,response,age,sex\nP001,responder,62,M\nP002,non-responder,55,F\n")

    # Scripts
    for s in ["01_qc_preprocessing.py", "02_clustering_umap.py",
              "03_differential_expression.py", "04_gene_signatures.py"]:
        (DEMO / f"scripts/{s}").write_text(f'"""Mock script: {s}"""\nimport scanpy as sc\nprint("Running {s}")\n')

    # Results
    (DEMO / "results/cell_type_proportions.csv").write_text(
        "cell_type,responder_mean,non_responder_mean,p_value\n"
        "CD8+ T cells,0.28,0.12,0.003\nTregs,0.08,0.15,0.012\n"
        "Macrophages M1,0.18,0.09,0.007\nMacrophages M2,0.06,0.14,0.009\n")
    (DEMO / "results/de_genes_responders.csv").write_text(
        "gene,log2FC,pvalue,padj\nCXCL13,2.8,1.2e-12,5.6e-10\nGZMB,2.1,3.4e-9,8.1e-7\n")
    (DEMO / "results/signature_scores.json").write_text(
        json.dumps({"response_signature_auc": 0.87, "n_genes": 42, "ci_95": [0.79, 0.94]}))

    # Figures (tiny PNGs with different colors)
    figs = [
        ("umap_cell_types.png", 60, 120, 200),
        ("cell_type_proportions.png", 100, 80, 160),
        ("volcano_de_genes.png", 140, 60, 100),
        ("gene_signature_heatmap.png", 80, 140, 80),
        ("survival_curve.png", 50, 100, 150),
    ]
    for name, r, g, b in figs:
        (DEMO / f"figures/{name}").write_bytes(tiny_png(120, 90, r, g, b))
        (DEMO / f"figures/{name.replace('.png','.pdf')}").write_bytes(b"%PDF-1.4 mock")

    # Paper sections
    sections = {
        "abstract": r"""\begin{abstract}
Immune checkpoint inhibitors have transformed NSCLC treatment, yet only 20--30\% of patients respond. We performed scRNA-seq on 45,000 cells from 12 NSCLC patients and identified 15 distinct cell populations. CD8+ T cells were enriched in responders ($p=0.003$), while Tregs and M2 macrophages dominated non-responders. A 42-gene signature achieved AUC=0.87 (95\% CI: 0.79--0.94), outperforming PD-L1 and TMB biomarkers.
\end{abstract}""",
        "introduction": r"""\section{Introduction}
Immune checkpoint inhibitors targeting PD-1/PD-L1 have revolutionized NSCLC treatment. Despite remarkable responses in a subset of patients, overall response rates remain at 20--30\%. The tumor microenvironment plays a central role in determining response. We present a comprehensive single-cell atlas from 12 patients stratified by immunotherapy response, revealing distinct immune compositions and a novel predictive gene signature.""",
        "methods": r"""\section{Methods}
\subsection{Patient Cohort} Twelve patients with advanced NSCLC enrolled prior to anti-PD-1 therapy. Response assessed by RECIST v1.1 at 12 weeks.
\subsection{scRNA-seq} 10x Genomics Chromium, Illumina NovaSeq 6000, Cell Ranger v7.0.
\subsection{Analysis} QC, normalization (scran), batch correction (Harmony) via Scanpy v1.9. Leiden clustering, Wilcoxon rank-sum test, BH correction.
\subsection{Signature} LASSO logistic regression, 5-fold CV, bootstrap ROC-AUC.""",
        "results": r"""\section{Results}
After QC, 42,318 cells retained. 15 clusters identified. CD8+ T cells: 28\% vs 12\% ($p=0.003$). Tregs: 8\% vs 15\% ($p=0.012$). 847 DE genes identified. 42-gene signature AUC=0.87 (95\% CI: 0.79--0.94), outperforming PD-L1 TPS (0.68) and TMB (0.72).""",
        "discussion": r"""\section{Discussion}
Our analysis reveals distinct immune profiles in responders vs non-responders. The 42-gene signature captures complex immune interplay beyond single-marker approaches. Limitations include small cohort size (n=12). Validation in larger cohorts is essential. We provide a clinically actionable predictive signature for immunotherapy response.""",
    }
    for name, content in sections.items():
        (DEMO / f"paper/sections/{name}.tex").write_text(content)

    (DEMO / "paper/main.tex").write_text(r"""\documentclass[11pt]{article}
\usepackage{amsmath,graphicx,booktabs,hyperref,natbib}
\title{Single-Cell Dissection of the NSCLC Tumor Immune Microenvironment}
\author{Autonomous Lab}
\begin{document}
\maketitle
\input{sections/abstract}
\input{sections/introduction}
\input{sections/methods}
\input{sections/results}
\input{sections/discussion}
\bibliographystyle{unsrtnat}
\bibliography{references}
\end{document}""")
    (DEMO / "paper/references.bib").write_text(r"""@article{gandhi2018,
  author={Gandhi, L. and others},
  title={Pembrolizumab plus Chemotherapy in Metastatic NSCLC},
  journal={NEJM}, year={2018}, volume={378}, pages={2078--2092}}""")

    # Meeting log
    (DEMO / ".autolab/meeting_log.md").write_text("""# Autonomous Lab Meeting Log

---

## Iteration 0 — PI Turn
**Date:** 2026-02-10 09:00 UTC

### Summary
Reviewed the project idea and planned initial scRNA-seq analysis approach. Targeting three core analyses: cell type composition, differential expression, and predictive gene signature.

### Details
### 1. REVIEW
Excellent research question with clear clinical relevance. Dataset of 12 patients with matched response data is solid.

### 2. FIGURE QA
No figures to evaluate.

### 3. AGENDA
1. QC and preprocessing of scRNA-seq data
2. PCA + UMAP + Leiden clustering
3. Cell type annotation using marker genes
4. Generate first UMAP visualization

PROGRESS: 5
STATUS: continue

---

## Iteration 0 — TRAINEE Turn
**Date:** 2026-02-10 11:30 UTC

### Summary
Completed QC, preprocessing, clustering. 42,318 cells retained, 15 distinct populations identified. Generated UMAP visualization.

### Details
### 1. APPROACH
QC → normalization → PCA → UMAP → Leiden clustering → marker annotation.

### 3. RESULTS
- Median genes/cell: 2,847 (IQR: 1,920–3,650)
- 15 Leiden clusters mapped to 8 major cell types
- CD8+ T cells show interesting subclustering

### 4. FIGURES
- `figures/umap_cell_types.png` — UMAP colored by 15 annotated populations

---

## Iteration 1 — PI Turn
**Date:** 2026-02-10 14:00 UTC

### Summary
Reviewed UMAP results — excellent quality. Now prioritizing differential composition and DE analysis between responders and non-responders.

### Details
### 1. REVIEW
Clean UMAP with well-separated clusters. Annotations are reasonable.

### 2. FIGURE QA
- `umap_cell_types.png`: Score **7/10** — needs larger legend font, point outlines, colorblind-safe palette

### 3. AGENDA
1. Cell type proportions per patient, compare responders vs non-responders
2. Differential expression within each cell type
3. Generate bar chart + volcano plot
4. Invite a **Statistician** and **Immunologist** for consultation

PROGRESS: 25
STATUS: continue

---

## Iteration 1 — TRAINEE Turn
**Date:** 2026-02-10 17:00 UTC

### Summary
CD8+ T cells significantly enriched in responders (28% vs 12%, p=0.003). Tregs and M2 macrophages enriched in non-responders. Generated volcano plot and proportion figures.

### Details
### 3. RESULTS
- **CD8+ T cells**: 28% vs 12% (p = 0.003)
- **Tregs**: 8% vs 15% (p = 0.012)
- **M1 Macrophages**: 18% vs 9% (p = 0.007)
- **M2 Macrophages**: 6% vs 14% (p = 0.009)
- 847 DE genes (|log2FC| > 1, FDR < 0.05)

### 4. FIGURES
- `figures/cell_type_proportions.png` — Bar chart with significance annotations
- `figures/volcano_de_genes.png` — Volcano highlighting top DE genes

---

## Iteration 2 — PI Turn
**Date:** 2026-02-11 09:00 UTC

### Summary
Strong results. Clear story of "hot" vs "cold" tumors. Prioritizing gene signature derivation and paper writing.

### Details
### 2. FIGURE QA
- `cell_type_proportions.png`: **8/10** — good significance annotations
- `volcano_de_genes.png`: **7/10** — needs point outlines, label top 10 genes

### 3. AGENDA
1. Derive 42-gene immunotherapy response signature (LASSO)
2. Cross-validation and ROC-AUC computation
3. Generate heatmap + survival curves
4. **Begin writing all paper sections**

PROGRESS: 50
STATUS: continue

---

## Iteration 2 — TRAINEE Turn
**Date:** 2026-02-11 14:00 UTC

### Summary
42-gene signature derived (AUC=0.87, 95% CI: 0.79–0.94). Outperforms PD-L1 and TMB. All paper sections written.

### Details
### 3. RESULTS
- **42-gene signature AUC: 0.87** (95% CI: 0.79–0.94)
- Outperforms PD-L1 TPS (AUC=0.68) and TMB (AUC=0.72)

### 4. FIGURES
- `figures/gene_signature_heatmap.png` — Clustered heatmap of 42 genes
- `figures/survival_curve.png` — Response curves by signature score

---

## Iteration 3 — PI Turn
**Date:** 2026-02-12 10:00 UTC

### Summary
Paper draft is strong. Final figure polish needed. Nearly ready for submission.

### Details
### 2. FIGURE QA
- `gene_signature_heatmap.png`: **8/10** — add clearer dendrogram
- `survival_curve.png`: **9/10** — add at-risk table
- All figures must pass at 9+/10

PROGRESS: 85
STATUS: continue

---

## Iteration 3 — TRAINEE Turn
**Date:** 2026-02-12 12:00 UTC

### Summary
All figures polished to publication quality. Paper sections refined. Ready for PI's final review.

### Details
All 5 figures updated with point outlines, proper fonts, colorblind-safe palettes. At-risk table added to survival curve. All sections refined for clarity.
""")

    (DEMO / ".autolab/meeting_summaries.md").write_text("# Meeting Summaries\n\n")

    # State
    state = {
        "iteration": 4,
        "next_role": "pi",
        "status": "active",
        "user_feedback": "",
        "progress": 70,
        "experts": [
            {"name": "Dr. Yamamoto", "role": "Statistician", "avatar": "statistician",
             "thought": "BH correction is appropriate for discovery with 15 cell types."},
            {"name": "Dr. Alvarez", "role": "Immunologist", "avatar": "immunologist",
             "thought": "Check for TOX and TCF7 in CD8+ subclusters — exhausted vs effector."},
        ],
        "editorial": {"phase": "none", "cover_letter": "", "reviewers": [],
                       "reviews": {}, "decision": "", "decision_feedback": "", "round": 0},
        "created_at": "2026-02-10T09:00:00+00:00",
        "last_updated": "2026-02-12T12:00:00+00:00",
    }
    with open(DEMO / ".autolab/state.json", "w") as f:
        json.dump(state, f, indent=2)

    print("Mock project created at:", DEMO)
    print("Files:", len(list(DEMO.rglob("*"))))

if __name__ == "__main__":
    main()
