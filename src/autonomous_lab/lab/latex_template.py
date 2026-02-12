"""
Autonomous Lab - LaTeX Paper Template

Generates the paper/ directory structure with main.tex,
section files, and references.bib.
"""

from pathlib import Path


MAIN_TEX = r"""\documentclass[11pt]{article}

% --- Packages ---
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

% --- Metadata ---
\title{TITLE}
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
"""

SECTION_TEMPLATES = {
    "abstract": r"""% Abstract
\begin{abstract}
% TODO: Write abstract summarizing the key findings.
\end{abstract}
""",
    "introduction": r"""% Introduction
\section{Introduction}
% TODO: Write introduction establishing context, gap, and contribution.
""",
    "methods": r"""% Methods
\section{Methods}
% TODO: Describe data, analysis methods, and statistical approaches.
""",
    "results": r"""% Results
\section{Results}
% TODO: Present findings with figures, tables, and statistics.
""",
    "discussion": r"""% Discussion
\section{Discussion}
% TODO: Interpret results, compare with literature, discuss limitations.
""",
}

REFERENCES_BIB = r"""% References
% Add BibTeX entries here. Example:
%
% @article{key2025,
%   author  = {Author, A. and Author, B.},
%   title   = {Title of the paper},
%   journal = {Journal Name},
%   year    = {2025},
%   volume  = {1},
%   pages   = {1--10},
%   doi     = {10.1234/example},
% }
"""


def create_paper_structure(project_dir: str, title: str = "TITLE") -> None:
    """
    Create the paper/ directory with main.tex, sections/*.tex,
    figures/, and references.bib.
    """
    base = Path(project_dir) / "paper"
    sections = base / "sections"
    figures = base / "figures"

    sections.mkdir(parents=True, exist_ok=True)
    figures.mkdir(exist_ok=True)

    # Write main.tex
    main_content = MAIN_TEX.replace("TITLE", title)
    (base / "main.tex").write_text(main_content, encoding="utf-8")

    # Write section files
    for name, template in SECTION_TEMPLATES.items():
        (sections / f"{name}.tex").write_text(template, encoding="utf-8")

    # Write references.bib
    (base / "references.bib").write_text(REFERENCES_BIB, encoding="utf-8")
