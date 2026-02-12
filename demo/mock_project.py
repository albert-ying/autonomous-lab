#!/usr/bin/env python3
"""
Mock the entire Autonomous Lab workflow for UI testing.

Creates a demo project at /tmp/autolab_demo/ with:
- Full .autolab/ state
- Realistic meeting log entries (5 PI/Trainee iterations)
- Paper sections with LaTeX content
- Mock figures (simple PNGs)
- Scripts, results files
- Expert consultants in state
- Progress at ~70%

Then starts the web UI server for visual inspection.
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

DEMO_DIR = "/tmp/autolab_demo_full"


def create_mock_figure(path: Path, label: str):
    """Create a simple colored PNG for figure previews."""
    try:
        import struct
        import zlib

        width, height = 400, 300
        # Create a simple gradient image
        raw_data = b""
        for y in range(height):
            raw_data += b"\x00"  # filter byte
            for x in range(width):
                r = int(60 + 100 * (x / width))
                g = int(80 + 80 * (y / height))
                b = int(140 + 60 * ((x + y) / (width + height)))
                raw_data += bytes([r, g, b])

        def make_chunk(chunk_type, data):
            chunk = chunk_type + data
            return (
                struct.pack(">I", len(data))
                + chunk
                + struct.pack(">I", zlib.crc32(chunk) & 0xFFFFFFFF)
            )

        png = b"\x89PNG\r\n\x1a\n"
        png += make_chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        png += make_chunk(b"IDAT", zlib.compress(raw_data))
        png += make_chunk(b"IEND", b"")

        path.write_bytes(png)
        print(f"  Created figure: {path.name}")
    except Exception as e:
        # Fallback: create a tiny valid PNG
        path.write_bytes(
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01"
            b"\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00"
            b"\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00\x00\x01\x01\x00"
            b"\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
        )
        print(f"  Created placeholder figure: {path.name}")


def setup_project():
    """Create the full mock project."""
    base = Path(DEMO_DIR)

    # Clean slate
    import shutil
    if base.exists():
        shutil.rmtree(base)

    # Directories
    autolab = base / ".autolab" / "profiles"
    autolab.mkdir(parents=True)
    for d in ("data", "scripts", "figures", "results"):
        (base / d).mkdir()
    sections = base / "paper" / "sections"
    sections.mkdir(parents=True)
    (base / "paper" / "figures").mkdir()

    print("=== Creating mock project ===")

    # --- Idea ---
    (base / ".autolab" / "idea.md").write_text(
        """# Single-Cell RNA-seq Analysis of Tumor Microenvironment

We have 10x Genomics scRNA-seq data from 12 patients with non-small cell
lung cancer (NSCLC). 6 patients responded to anti-PD-1 immunotherapy and
6 did not. We aim to identify the immune cell populations and gene programs
that distinguish responders from non-responders.

## Data
- data/counts_matrix.h5ad — AnnData with 45,000 cells × 22,000 genes
- data/clinical_metadata.csv — Patient demographics, treatment, response

## Goals
1. Characterize the tumor microenvironment at single-cell resolution
2. Identify cell types enriched in responders vs non-responders
3. Discover gene signatures predictive of immunotherapy response
4. Produce a publication-ready manuscript for Nature Medicine
""",
        encoding="utf-8",
    )

    # --- Profiles ---
    import yaml

    pi_profile = {
        "title": "Principal Investigator",
        "expertise": "Cancer immunology and computational biology",
        "goal": "Publish a high-impact paper in Nature Medicine on immunotherapy response biomarkers",
        "personality": [
            "Visionary: sees the big picture and focuses on clinical relevance",
            "Critical: demands rigorous statistical validation",
            "Publication-focused: every figure must be Nature-quality",
            "Efficient: prioritizes the 2-3 analyses that tell the story",
        ],
    }
    trainee_profile = {
        "title": "Postdoctoral Researcher",
        "expertise": "Single-cell genomics, Python, scanpy, statistical modeling",
        "goal": "Execute analyses with technical excellence and write clean code",
        "personality": [
            "Dedicated: completes tasks thoroughly",
            "Technical: writes reproducible, well-documented code",
            "Proactive: identifies additional analyses aligned with PI vision",
            "Clear communicator: reports results with proper statistics",
        ],
    }
    with open(base / ".autolab" / "profiles" / "pi.yaml", "w") as f:
        yaml.dump(pi_profile, f)
    with open(base / ".autolab" / "profiles" / "trainee.yaml", "w") as f:
        yaml.dump(trainee_profile, f)

    # --- Config ---
    with open(base / ".autolab" / "config.yaml", "w") as f:
        yaml.dump({"max_iterations": 50, "target_journal": "Nature Medicine"}, f)

    # --- Mock data files ---
    (base / "data" / "counts_matrix.h5ad").write_text("(mock h5ad)")
    (base / "data" / "clinical_metadata.csv").write_text(
        "patient_id,response,age,sex\nP001,responder,62,M\nP002,non-responder,55,F\n"
    )
    print("  Created data files")

    # --- Mock scripts ---
    (base / "scripts" / "01_qc_preprocessing.py").write_text(
        '"""QC and preprocessing pipeline for scRNA-seq data."""\nimport scanpy as sc\n# ... preprocessing code ...\n'
    )
    (base / "scripts" / "02_clustering_umap.py").write_text(
        '"""Dimensionality reduction and clustering."""\nimport scanpy as sc\nimport matplotlib.pyplot as plt\n# ... clustering ...\n'
    )
    (base / "scripts" / "03_differential_expression.py").write_text(
        '"""Differential expression between responders and non-responders."""\nimport scanpy as sc\n# ... DE analysis ...\n'
    )
    (base / "scripts" / "04_gene_signatures.py").write_text(
        '"""Gene signature scoring and survival analysis."""\nimport pandas as pd\n# ... signature analysis ...\n'
    )
    print("  Created scripts")

    # --- Mock results ---
    (base / "results" / "cell_type_proportions.csv").write_text(
        "cell_type,responder_mean,non_responder_mean,p_value\n"
        "CD8+ T cells,0.28,0.12,0.003\n"
        "Tregs,0.08,0.15,0.012\n"
        "Macrophages M1,0.18,0.09,0.007\n"
        "Macrophages M2,0.06,0.14,0.009\n"
        "B cells,0.11,0.10,0.82\n"
        "NK cells,0.09,0.05,0.04\n"
    )
    (base / "results" / "de_genes_responders.csv").write_text(
        "gene,log2FC,pvalue,padj\nCXCL13,2.8,1.2e-12,5.6e-10\nGZMB,2.1,3.4e-9,8.1e-7\n"
        "IFNG,1.9,1.1e-8,2.0e-6\nPDCD1,1.5,4.5e-7,6.2e-5\n"
    )
    (base / "results" / "signature_scores.json").write_text(
        json.dumps({"response_signature_auc": 0.87, "n_genes": 42, "ci_95": [0.79, 0.94]})
    )
    print("  Created results")

    # --- Mock figures ---
    for fig in [
        "umap_cell_types.png",
        "cell_type_proportions.png",
        "volcano_de_genes.png",
        "gene_signature_heatmap.png",
        "survival_curve.png",
    ]:
        create_mock_figure(base / "figures" / fig, fig)
        # Also create PDF placeholder
        pdf_name = fig.replace(".png", ".pdf")
        (base / "figures" / pdf_name).write_bytes(b"%PDF-1.4 mock")

    # --- Paper sections with real content ---
    (sections / "abstract.tex").write_text(
        r"""\begin{abstract}
Immune checkpoint inhibitors have transformed non-small cell lung cancer (NSCLC)
treatment, yet only 20--30\% of patients respond. Understanding the tumor
microenvironment (TME) composition that predicts response remains a critical
unmet need. Here, we performed single-cell RNA sequencing on 45,000 cells from
12 NSCLC patients (6 responders, 6 non-responders to anti-PD-1 therapy).
We identified 15 distinct cell populations and found that CD8+ T cell abundance
was significantly higher in responders ($p = 0.003$), while regulatory T cells
and M2 macrophages were enriched in non-responders. We derived a 42-gene
immunotherapy response signature (AUC = 0.87, 95\% CI: 0.79--0.94) that
outperforms existing biomarkers. Our findings provide a comprehensive atlas
of the NSCLC immune landscape and a clinically actionable predictive signature.
\end{abstract}
""",
        encoding="utf-8",
    )

    (sections / "introduction.tex").write_text(
        r"""\section{Introduction}
Immune checkpoint inhibitors (ICIs) targeting the PD-1/PD-L1 axis have
revolutionized the treatment of non-small cell lung cancer (NSCLC).
Despite remarkable responses in a subset of patients, overall response
rates remain at 20--30\%, highlighting the urgent need for predictive
biomarkers \cite{gandhi2018pembrolizumab}.

The tumor microenvironment (TME) plays a central role in determining
immunotherapy response. Single-cell RNA sequencing (scRNA-seq) has emerged
as a powerful tool to dissect the cellular heterogeneity of the TME at
unprecedented resolution \cite{zheng2017landscape}.

In this study, we present a comprehensive single-cell atlas of the NSCLC
tumor microenvironment from 12 patients stratified by immunotherapy response.
Our analysis reveals distinct immune cell compositions that distinguish
responders from non-responders and identifies a novel gene signature with
strong predictive performance.
""",
        encoding="utf-8",
    )

    (sections / "methods.tex").write_text(
        r"""\section{Methods}
\subsection{Patient Cohort and Sample Collection}
Twelve patients with advanced NSCLC were enrolled prior to first-line
anti-PD-1 immunotherapy. Response was assessed by RECIST v1.1 criteria
at 12 weeks. Fresh tumor biopsies were collected and processed within
2 hours.

\subsection{Single-Cell RNA Sequencing}
Single-cell suspensions were prepared using the 10x Genomics Chromium
platform. Libraries were sequenced on an Illumina NovaSeq 6000 targeting
50,000 reads per cell. Raw reads were processed with Cell Ranger v7.0.

\subsection{Computational Analysis}
Quality control, normalization (scran), and batch correction (Harmony)
were performed using Scanpy v1.9. Leiden clustering identified 15 cell
populations annotated using known marker genes. Differential expression
used the Wilcoxon rank-sum test with Benjamini-Hochberg correction
($\alpha = 0.05$).

\subsection{Gene Signature Derivation}
A 42-gene response signature was derived using LASSO logistic regression
with 5-fold cross-validation. Performance was assessed by ROC-AUC with
bootstrap confidence intervals (1000 iterations).
""",
        encoding="utf-8",
    )

    (sections / "results.tex").write_text(
        r"""\section{Results}
\subsection{Single-Cell Atlas of the NSCLC Tumor Microenvironment}
After quality control, 42,318 cells from 12 patients were retained
(median: 3,526 cells/patient). Unsupervised clustering identified 15
transcriptionally distinct populations (Figure~\ref{fig:umap}).

\subsection{Differential Immune Composition in Responders}
CD8+ T cells were significantly enriched in responders (28\% vs 12\%,
$p = 0.003$, Wilcoxon). Conversely, regulatory T cells (8\% vs 15\%,
$p = 0.012$) and M2 macrophages (6\% vs 14\%, $p = 0.009$) were
enriched in non-responders (Figure~\ref{fig:proportions}).

\subsection{Gene Expression Signatures of Response}
Differential expression analysis identified 847 genes significantly
altered between responders and non-responders (|log2FC| > 1, FDR < 0.05).
The top upregulated genes in responders included \textit{CXCL13},
\textit{GZMB}, and \textit{IFNG} (Figure~\ref{fig:volcano}).

\subsection{Predictive Immunotherapy Response Signature}
Our 42-gene signature achieved an AUC of 0.87 (95\% CI: 0.79--0.94),
outperforming PD-L1 TPS (AUC = 0.68) and TMB (AUC = 0.72).
""",
        encoding="utf-8",
    )

    (sections / "discussion.tex").write_text(
        r"""\section{Discussion}
Our single-cell analysis reveals that the composition and functional state
of the tumor immune microenvironment are strongly associated with
immunotherapy response in NSCLC. The enrichment of CD8+ T cells and M1
macrophages in responders, coupled with regulatory T cell and M2 macrophage
enrichment in non-responders, paints a clear picture of immunologically
``hot'' versus ``cold'' tumors.

The 42-gene immunotherapy response signature we derived substantially
outperforms existing biomarkers including PD-L1 expression and tumor
mutational burden. This signature captures the complex interplay between
immune cell types and their functional states, which single-marker
approaches cannot.

\subsection{Limitations}
Our study is limited by the relatively small cohort size (n=12).
Validation in larger, independent cohorts is essential. Additionally,
spatial information is lost in dissociated scRNA-seq.

\subsection{Conclusions}
We present a comprehensive single-cell atlas of the NSCLC immune
microenvironment and a clinically actionable 42-gene predictive signature.
These findings advance our understanding of immunotherapy resistance
mechanisms and provide tools for patient stratification.
""",
        encoding="utf-8",
    )

    # References
    (base / "paper" / "references.bib").write_text(
        r"""@article{gandhi2018pembrolizumab,
  author  = {Gandhi, Leena and others},
  title   = {Pembrolizumab plus Chemotherapy in Metastatic NSCLC},
  journal = {New England Journal of Medicine},
  year    = {2018},
  volume  = {378},
  pages   = {2078--2092},
}

@article{zheng2017landscape,
  author  = {Zheng, Chunhong and others},
  title   = {Landscape of Infiltrating T Cells in Liver Cancer},
  journal = {Cell},
  year    = {2017},
  volume  = {169},
  pages   = {1342--1356},
}
""",
        encoding="utf-8",
    )

    # Main tex
    (base / "paper" / "main.tex").write_text(
        r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage[margin=1in]{geometry}
\usepackage{natbib}
\usepackage{xcolor}
\usepackage{float}

\title{Single-Cell Dissection of the Tumor Immune Microenvironment
Reveals Predictive Signatures for Immunotherapy Response in NSCLC}
\author{Autonomous Lab}
\date{\today}

\begin{document}
\maketitle
\input{sections/abstract}
\input{sections/introduction}
\input{sections/methods}
\input{sections/results}
\input{sections/discussion}
\bibliographystyle{unsrtnat}
\bibliography{references}
\end{document}
""",
        encoding="utf-8",
    )
    print("  Created paper sections")

    return base


def create_meeting_log(base: Path):
    """Create a realistic meeting log with 5 iterations."""
    log = """# Autonomous Lab Meeting Log

---

## Iteration 0 — PI Turn
**Date:** 2026-02-10 09:00 UTC

### Summary
Reviewed the project idea and planned the initial approach for scRNA-seq analysis of the NSCLC tumor microenvironment.

### Details
### 1. REVIEW
Excellent research question with clear clinical relevance. The dataset of 12 patients with matched response data provides a solid foundation. We need to focus on three core analyses: (1) cell type composition, (2) differential gene expression, and (3) a predictive signature.

### 2. FIGURE QA
No figures to evaluate.

### 3. AGENDA
1. Perform QC and preprocessing of the scRNA-seq data
2. Run dimensionality reduction (PCA + UMAP) and clustering
3. Annotate cell types using known marker genes
4. Generate the first UMAP visualization colored by cell type

### 4. AGENDA QUESTIONS
1. What is the quality of the raw data (% mitochondrial, genes per cell)?
2. How many distinct clusters do we observe?
3. Can we clearly identify the major immune cell populations?

### 5. AGENDA RULES
1. Use scanpy for all single-cell analysis
2. Filter cells with >20% mitochondrial genes
3. Use Leiden clustering with resolution 1.0 as starting point

### 6. MANUSCRIPT PLANNING
Not yet — focus on analyses first.

PROGRESS: 5
STATUS: continue

---

## Iteration 0 — TRAINEE Turn
**Date:** 2026-02-10 11:30 UTC

### Summary
Completed QC, preprocessing, clustering, and cell type annotation. Generated UMAP visualization with 15 distinct cell populations identified.

### Details
### 1. APPROACH
Followed the PI's agenda: QC → normalization → PCA → UMAP → Leiden clustering → marker-based annotation.

### 2. EXECUTION
- Ran `scripts/01_qc_preprocessing.py`: Filtered 2,682 low-quality cells (5.6%), retained 42,318 cells
- Ran `scripts/02_clustering_umap.py`: Identified 15 clusters at resolution 1.0
- Annotated clusters: CD8+ T cells (8 clusters), CD4+ T cells (3), Macrophages (2), B cells (1), NK cells (1)

### 3. RESULTS
- Median genes/cell: 2,847 (IQR: 1,920–3,650)
- Median UMIs/cell: 8,420
- Mean mitochondrial fraction: 4.2% (post-filter)
- 15 Leiden clusters mapped to 8 major cell types

### 4. FIGURES
- `figures/umap_cell_types.png` — UMAP colored by annotated cell types (15 populations)

### 5. INTERPRETATION
Clear separation of immune populations. CD8+ T cells show interesting subclustering that may reveal exhaustion states.

### 6. ADDITIONAL OBSERVATIONS
Noticed a rare population (~2% of cells) that may be plasmacytoid dendritic cells — worth investigating.

### 7. AGENDA ANSWERS
1. Data quality is excellent: low mitochondrial content, high gene detection
2. 15 distinct clusters at resolution 1.0
3. Yes — all major immune populations clearly identified

---

## Iteration 1 — PI Turn
**Date:** 2026-02-10 14:00 UTC

### Summary
Reviewed UMAP and clustering results. Excellent quality. Now prioritizing differential composition analysis and DE between responders/non-responders.

### Details
### 1. REVIEW
The UMAP looks clean with well-separated clusters. Cell type annotations are reasonable. The 15-cluster resolution captures important subtypes.

### 2. FIGURE QA
- `umap_cell_types.png`: Score 7/10. Good layout but needs: (1) larger font for legend text, (2) point outlines for better visibility, (3) color palette should be colorblind-safe. Required fixes: apply the CLAUDE.md visualization rules.

### 3. AGENDA
1. Compute cell type proportions per patient and compare responders vs non-responders
2. Run differential expression between responders and non-responders within each cell type
3. Generate: (a) stacked bar chart of cell proportions, (b) volcano plot of DE genes
4. Start the gene signature analysis

### 7. EXPERT CONSULTANTS
Inviting a statistician to advise on the multiple testing correction strategy and a clinical immunologist to validate cell type annotations.

PROGRESS: 25
STATUS: continue

---

## Iteration 1 — TRAINEE Turn
**Date:** 2026-02-10 17:00 UTC

### Summary
Completed differential composition analysis and DE. Found CD8+ T cells significantly enriched in responders (p=0.003). Generated volcano plot and proportion figures.

### Details
### 1. APPROACH
Ran cell type proportion analysis per patient, then Wilcoxon tests for composition differences. Ran DE within CD8+ T cells using Wilcoxon rank-sum.

### 2. EXECUTION
- Ran `scripts/03_differential_expression.py`
- Generated `figures/cell_type_proportions.png` and `figures/volcano_de_genes.png`
- Saved results to `results/cell_type_proportions.csv` and `results/de_genes_responders.csv`

### 3. RESULTS
Key findings:
- **CD8+ T cells**: 28% in responders vs 12% in non-responders (p = 0.003)
- **Tregs**: 8% vs 15% (p = 0.012)
- **M1 Macrophages**: 18% vs 9% (p = 0.007)
- **M2 Macrophages**: 6% vs 14% (p = 0.009)
- 847 DE genes (|log2FC| > 1, FDR < 0.05)
- Top genes: CXCL13 (log2FC=2.8), GZMB (2.1), IFNG (1.9)

### 4. FIGURES
- `figures/cell_type_proportions.png` — Bar chart with significance annotations
- `figures/volcano_de_genes.png` — Volcano plot highlighting top DE genes

### 5. INTERPRETATION
Clear immunological distinction: responders have "hot" tumors (high CD8+, M1) while non-responders have immunosuppressive profiles (high Tregs, M2).

---

## Iteration 2 — PI Turn
**Date:** 2026-02-11 09:00 UTC

### Summary
Excellent results. The composition and DE analyses tell a compelling story. Now focusing on the gene signature and paper writing.

### Details
### 1. REVIEW
Strong statistical evidence for differential immune composition. The p-values are convincing given the sample size. The volcano plot clearly shows the key genes.

### 2. FIGURE QA
- `cell_type_proportions.png`: Score 8/10. Good use of significance annotations. Minor: increase bar width slightly.
- `volcano_de_genes.png`: Score 7/10. Needs point outlines per CLAUDE.md rules. Label the top 10 genes.

### 3. AGENDA
1. Derive the 42-gene immunotherapy response signature using LASSO
2. Validate with cross-validation and compute ROC-AUC
3. Generate heatmap of signature genes across patients
4. Generate survival/response curve
5. **Begin writing paper sections**: Abstract, Introduction, Methods, Results
6. Update all figures to meet Nature Medicine standards

### 8. PROGRESS
PROGRESS: 50
STATUS: continue

---

## Iteration 2 — TRAINEE Turn
**Date:** 2026-02-11 14:00 UTC

### Summary
Derived 42-gene signature (AUC=0.87). Generated heatmap and survival curves. Wrote all paper sections.

### Details
### 1. APPROACH
Used LASSO logistic regression with 5-fold CV for signature derivation. Wrote all LaTeX sections.

### 2. EXECUTION
- Ran `scripts/04_gene_signatures.py`
- Generated `figures/gene_signature_heatmap.png` and `figures/survival_curve.png`
- Wrote all 5 paper sections in `paper/sections/`

### 3. RESULTS
- **42-gene signature AUC: 0.87** (95% CI: 0.79–0.94)
- Outperforms PD-L1 TPS (AUC=0.68) and TMB (AUC=0.72)
- Cross-validation stable across all 5 folds

### 4. FIGURES
- `figures/gene_signature_heatmap.png` — Clustered heatmap of 42 signature genes
- `figures/survival_curve.png` — Response curves stratified by signature score

### 5. INTERPRETATION
The signature captures multi-dimensional immune states better than any single biomarker. This is the key finding for the paper.

---

## Iteration 3 — PI Turn
**Date:** 2026-02-11 18:00 UTC

### Summary
Paper draft is strong. Making final revisions and figure polish. Nearly ready for submission.

### Details
### 1. REVIEW
Impressive signature performance. The paper sections are well-written with appropriate statistics. A few refinements needed.

### 2. FIGURE QA
- `gene_signature_heatmap.png`: Score 8/10. Add clearer dendrogram and colorbar label.
- `survival_curve.png`: Score 9/10. Excellent. Minor: add at-risk table below.
- All previous figures: re-confirm they follow CLAUDE.md rules.

### 3. AGENDA
1. Final figure polish — all figures must pass QA at 9+/10
2. Review all paper sections for consistency and flow
3. Add missing references
4. Prepare for submission

PROGRESS: 85
STATUS: continue

---

## Iteration 3 — TRAINEE Turn
**Date:** 2026-02-12 10:00 UTC

### Summary
Polished all figures to publication quality. Updated paper sections with reviewer-ready prose. All 5 figures score 9+/10.

### Details
### 1. APPROACH
Applied CLAUDE.md visualization rules to all figures. Refined paper text for clarity and consistency.

### 2. EXECUTION
- Updated all 5 figures with point outlines, proper fonts, colorblind-safe palettes
- Added at-risk table to survival curve
- Polished all LaTeX sections
- Added 2 new references

### 3. RESULTS
All figures now meet Nature Medicine publication standards.

### 4. FIGURES
- All 5 figures updated: umap_cell_types, cell_type_proportions, volcano_de_genes, gene_signature_heatmap, survival_curve

### 5. INTERPRETATION
The manuscript is ready for PI's final review and submission decision.

---

## Iteration 4 — PI Turn
**Date:** 2026-02-12 12:00 UTC

### Summary
Manuscript is complete and ready for submission. All analyses are rigorous, figures meet top-journal standards, and the narrative is compelling.

### Details
### 1. REVIEW
The manuscript is excellent. The story flows logically from atlas → differential composition → gene signature → clinical utility. All figures are publication-quality.

### 2. FIGURE QA
All figures: 9/10 or above. Ready for Nature Medicine.

PROGRESS: 95
STATUS: ready_for_review
"""

    (base / ".autolab" / "meeting_log.md").write_text(log, encoding="utf-8")

    # Meeting summaries (compressed earlier meetings)
    (base / ".autolab" / "meeting_summaries.md").write_text(
        "# Meeting Summaries\n\n", encoding="utf-8"
    )
    print("  Created meeting log (5 iterations, 10 turns)")


def create_state(base: Path, phase: str = "active"):
    """Create state.json at various phases."""
    now = datetime.now(timezone.utc).isoformat()

    state = {
        "iteration": 4,
        "next_role": "pi",
        "status": phase,
        "user_feedback": "",
        "progress": 70,
        "experts": [
            {
                "name": "Dr. Yamamoto",
                "role": "Statistician",
                "avatar": "statistician",
                "thought": "Multiple testing correction with 15 cell types needs Bonferroni or BH — I'd recommend BH for discovery.",
            },
            {
                "name": "Dr. Alvarez",
                "role": "Immunologist",
                "avatar": "immunologist",
                "thought": "The CD8+ T cell subclusters may include exhausted vs effector phenotypes — check for TOX and TCF7.",
            },
        ],
        "editorial": {
            "phase": "none",
            "cover_letter": "",
            "reviewers": [],
            "reviews": {},
            "decision": "",
            "decision_feedback": "",
            "round": 0,
        },
        "created_at": "2026-02-10T09:00:00+00:00",
        "last_updated": now,
    }

    with open(base / ".autolab" / "state.json", "w") as f:
        json.dump(state, f, indent=2)
    print(f"  State: iteration=4, progress=70%, status={phase}")


def main():
    print("\n" + "=" * 60)
    print("  AUTONOMOUS LAB — Full Mock Setup")
    print("=" * 60 + "\n")

    base = setup_project()
    create_meeting_log(base)
    create_state(base, "active")

    print(f"\n  Demo project created at: {DEMO_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
