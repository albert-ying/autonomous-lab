#!/usr/bin/env python3
"""Advance the mock project through editorial workflow phases."""
import json
import sys
from pathlib import Path

DEMO = Path("/tmp/autolab_demo_full")
STATE_PATH = DEMO / ".autolab" / "state.json"

def read_state():
    return json.loads(STATE_PATH.read_text())

def write_state(state):
    STATE_PATH.write_text(json.dumps(state, indent=2))
    print(f"  State updated: status={state['status']}, editorial.phase={state['editorial']['phase']}")

COVER_LETTER = """Dear Editor,

We are pleased to submit our manuscript entitled **"Single-Cell Dissection of the Tumor Immune Microenvironment Reveals Predictive Signatures for Immunotherapy Response in NSCLC"** for consideration in *Nature Medicine*.

## Key Findings

1. **Comprehensive single-cell atlas**: We profiled 42,318 cells from 12 NSCLC patients, identifying 15 distinct immune cell populations in the tumor microenvironment.

2. **Differential immune composition**: CD8+ T cells were significantly enriched in immunotherapy responders (28% vs 12%, p=0.003), while regulatory T cells and M2 macrophages dominated non-responders.

3. **Novel predictive signature**: We derived a 42-gene immunotherapy response signature achieving AUC=0.87 (95% CI: 0.79–0.94), substantially outperforming existing biomarkers including PD-L1 TPS (AUC=0.68) and TMB (AUC=0.72).

## Significance

These findings advance our understanding of immune-mediated mechanisms of immunotherapy response and provide a clinically actionable tool for patient stratification. The 42-gene signature has the potential to guide treatment decisions for the approximately 70% of NSCLC patients who do not respond to immune checkpoint inhibitors.

We believe this work is well-suited for *Nature Medicine* given its direct clinical implications for precision oncology.

Sincerely,
The Autonomous Lab Team
"""

REVIEWS = {
    "reviewer_0": {
        "reviewer": {"name": "Dr. Sarah Chen", "role": "Computational Biologist"},
        "report": """## Overall Assessment

This manuscript presents a well-executed single-cell RNA-seq analysis of the NSCLC tumor microenvironment with clinically relevant findings. The 42-gene signature is a notable contribution.

## Strengths
- Comprehensive profiling of 42,318 cells with clear cell type annotations
- Strong statistical framework with appropriate multiple testing corrections
- Gene signature outperforms existing clinical biomarkers (PD-L1, TMB)
- Well-written manuscript with publication-quality figures

## Weaknesses
- Small cohort size (n=12) limits generalizability — validation in independent cohorts is essential
- Batch effects between patients could confound results despite Harmony correction
- The LASSO signature may be overfit given the small sample-to-feature ratio

## Specific Comments
1. **Methods**: Please clarify the Harmony integration parameters used
2. **Figure 3**: The volcano plot should include fold-change thresholds as dashed lines
3. **Discussion**: Please discuss spatial context loss from dissociated scRNA-seq

## Recommendation
**Minor Revision** — Address the points above and discuss validation plans

## Confidence
4/5 — High confidence in this area""",
    },
    "reviewer_1": {
        "reviewer": {"name": "Dr. James Morrison", "role": "Clinical Oncologist"},
        "report": """## Overall Assessment

An interesting study that bridges molecular profiling with clinical outcomes in NSCLC immunotherapy. The predictive signature could be clinically useful if validated.

## Strengths
- Clear clinical relevance and well-defined patient cohort
- Response assessment by RECIST v1.1 is appropriate
- The gene signature's superior performance over PD-L1 and TMB is compelling
- Statistical methods are generally appropriate

## Weaknesses
- Sample size of 12 patients is a significant limitation for clinical claims
- No external validation cohort is included
- The 42-gene signature may be impractical for clinical implementation without a companion diagnostic

## Specific Comments
1. Include hazard ratios for the signature in addition to AUC
2. Discuss how the 42-gene signature would be implemented clinically (e.g., NanoString, qPCR panel)
3. Compare with the Ayers et al. interferon-gamma signature and other published scores
4. Table 1 with patient demographics would strengthen the manuscript

## Recommendation
**Major Revision** — External validation and clinical implementation discussion are needed

## Confidence
3/5 — Moderate (bioinformatics not my primary expertise)""",
    },
    "reviewer_2": {
        "reviewer": {"name": "Dr. Mei-Lin Park", "role": "Immunologist"},
        "report": """## Overall Assessment

This study provides valuable insights into the immune landscape of NSCLC in the context of immunotherapy response. The cellular characterization is thorough, and the findings are consistent with known immunology.

## Strengths
- Excellent cell type resolution with 15 distinct populations
- CD8+ T cell / Treg / M2 macrophage findings align with established immunology
- CXCL13 as a top gene is consistent with tertiary lymphoid structures literature
- Robust statistical methodology with BH correction

## Weaknesses
- Functional states of CD8+ T cells (exhausted vs effector) are not resolved
- TCR clonotype analysis is missing — would strengthen the T cell narrative
- Spatial relationships between cell types cannot be inferred

## Specific Comments
1. **Critical**: Analyze TOX/TCF7 expression to distinguish exhausted from progenitor CD8+ T cells
2. Include CD8+ T cell subclustering analysis
3. Discuss the role of tertiary lymphoid structures given CXCL13 prominence
4. Consider SCENIC or similar for TF regulon analysis

## Recommendation
**Minor Revision** — The T cell functional characterization needs expansion

## Confidence
5/5 — Very high confidence; this is my core expertise""",
    },
}


def phase_submitted():
    """PI has submitted manuscript — editor sees cover letter."""
    state = read_state()
    state["status"] = "submitted_to_editor"
    state["progress"] = 95
    state["editorial"] = {
        "phase": "submitted",
        "cover_letter": COVER_LETTER,
        "reviewers": [],
        "reviews": {},
        "decision": "",
        "decision_feedback": "",
        "round": 1,
    }
    write_state(state)

def phase_reviewers_invited():
    """Editor invited 3 reviewers."""
    state = read_state()
    state["status"] = "under_review"
    state["editorial"]["phase"] = "reviewers_invited"
    state["editorial"]["reviewers"] = [
        {"name": "Dr. Sarah Chen", "role": "Computational Biologist", "id": "reviewer_0"},
        {"name": "Dr. James Morrison", "role": "Clinical Oncologist", "id": "reviewer_1"},
        {"name": "Dr. Mei-Lin Park", "role": "Immunologist", "id": "reviewer_2"},
    ]
    write_state(state)

def phase_reviewing():
    """Reviews in progress — first reviewer done."""
    state = read_state()
    state["editorial"]["phase"] = "under_review"
    state["editorial"]["reviews"] = {
        "reviewer_0": REVIEWS["reviewer_0"],
    }
    state["next_role"] = "reviewer_1"
    write_state(state)

def phase_reviews_complete():
    """All 3 reviews complete."""
    state = read_state()
    state["editorial"]["phase"] = "reviews_complete"
    state["editorial"]["reviews"] = REVIEWS
    state["next_role"] = "editor"
    write_state(state)

def phase_decision_minor():
    """Editor decides: minor revision."""
    state = read_state()
    state["editorial"]["phase"] = "decision_made"
    state["editorial"]["decision"] = "minor_revision"
    state["editorial"]["decision_feedback"] = (
        "The reviewers are largely positive about this work. "
        "Please address the specific comments from all three reviewers, "
        "particularly the T cell functional characterization (Reviewer 3) "
        "and the clinical implementation discussion (Reviewer 2). "
        "A revised manuscript addressing these points should be suitable for publication."
    )
    state["status"] = "revision_requested"
    state["next_role"] = "pi"
    write_state(state)


PHASES = {
    "submitted": phase_submitted,
    "invited": phase_reviewers_invited,
    "reviewing": phase_reviewing,
    "complete": phase_reviews_complete,
    "decision": phase_decision_minor,
}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in PHASES:
        print(f"Usage: {sys.argv[0]} <{'|'.join(PHASES.keys())}>")
        sys.exit(1)
    PHASES[sys.argv[1]]()
