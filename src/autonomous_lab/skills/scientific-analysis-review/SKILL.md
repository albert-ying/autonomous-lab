---
name: scientific-analysis-review
description: >
  Critically review AI-agent-conducted scientific analyses for correctness, rigor,
  and completeness. Use this skill whenever an analysis session has completed and
  needs validation, when a user asks to "review," "validate," "check," or "audit"
  a computational analysis, or when an agent pipeline produces scientific results
  that require quality control before reporting. Also trigger when the user references
  an execution trace, notebook, or conversation history from a prior analysis session.
  This skill should run as the final step of any autonomous scientific analysis pipeline.
---

# Scientific Analysis Review

## Purpose

This skill validates AI-agent-conducted scientific analyses by systematically checking
for bugs, hallucinations, logical errors, and scientific rigor issues. It produces a
structured review report with severity-graded findings and actionable recommendations.

The review operates on the execution trace of a completed analysis: conversation history,
code cells, intermediate outputs, and final claims. It does NOT re-run the analysis.
It audits what was done.


## When to Run

1. After any autonomous analysis session completes (automatic trigger)
2. When a user uploads or references an execution trace, notebook, or session log
3. When a user asks to validate, audit, or review computational results
4. Before results are incorporated into a manuscript or presentation


## Input Requirements

The reviewer needs access to one or more of the following (in order of preference):

1. Full conversation/execution trace (best: contains code, outputs, and reasoning)
2. Jupyter notebook or script with outputs
3. Summary of analysis steps with key code blocks and results
4. Final report or claims to be checked against available data

If the input is empty or contains no analysis, state this in one sentence and stop.
Do not fill out N/A tables for empty sessions.


## Review Categories

Run every category in order. Each category contains specific checks.
Skip a category only if it is entirely inapplicable (e.g., no code was executed).

### 1. Execution Integrity

Check whether the code ran as intended and produced the expected outputs.

- **Silent failures**: Did any operation return empty results, null values, or default
  outputs without raising an error? (Example: an API call that returns `[]` but the
  analysis continues as if data was retrieved.)
- **Unexecuted branches**: Were any code blocks skipped, commented out, or short-circuited?
  Did conditional logic cause important steps to be bypassed?
- **Data loading**: Were all referenced files, databases, or APIs actually accessed?
  Do row/column counts match expectations after loading?
- **Dependency chain**: Did later steps depend on outputs from failed or incomplete
  earlier steps? Trace the data flow from input to final claim.
- **Duplicate entries**: Were deduplication steps applied where needed? Check for
  repeated identifiers in filtered or merged datasets.

### 2. Numerical and Statistical Verification

Check whether numbers, calculations, and statistical claims are correct.

- **Denominator checking**: For every percentage, fraction, or ratio in the final output,
  verify the base (denominator) is correct and complete. Flag any percentage where the
  denominator was derived from an incomplete or earlier subset of the data.
- **Threshold justification**: For every numeric threshold used in filtering (e.g.,
  expression cutoffs, p-value thresholds, score cutoffs), check whether a justification
  was provided. Flag arbitrary thresholds that materially affect results.
- **Sample size propagation**: Track how many data points survive each filtering step.
  Flag any step that removes >50% of data without explicit justification.
- **Statistical test assumptions**: Were parametric tests applied to non-normal data?
  Were multiple comparison corrections applied when needed?
- **Rounding and precision**: Do reported numbers match computed values? Check for
  truncation errors or rounding inconsistencies between text and code output.

### 3. Logical Chain Validation

Trace whether each conclusion follows from the preceding analysis step.

- **Inferential leaps**: Does the analysis jump from observation A to conclusion B
  without establishing the logical connection? (Example: "A is associated with B
  in database X" -> "A directly regulates B." Association does not imply a specific
  mechanistic relationship.)
- **Proxy validity**: When a proxy measure is used (e.g., aggregate tissue-level
  measurements as proxy for cell-type-specific activity), is the proxy acknowledged
  and its limitations stated? Flag proxies that are treated as ground truth.
- **Causal language**: Does the analysis use causal language ("causes," "drives,"
  "mediates") when only correlational or associational evidence was provided?
- **Negative results**: Were null or negative results reported, or were they silently
  dropped? Check whether the analysis only presents findings that support the narrative.

### 4. Scope Drift Detection

Check whether the agent answered the question that was actually asked.

- **Original question**: Restate the user's original request in one sentence.
- **What was delivered**: Summarize what the analysis actually produced.
- **Gap analysis**: Identify any mismatch. Common failure modes:
  - Agent answered an easier question than the one asked
  - Agent produced a general survey instead of the specific analysis requested
  - Agent introduced scope that was not requested and used it to avoid the hard part
  - Agent conflated related but distinct concepts (e.g., "differentially expressed"
    vs. "causally involved")

### 5. Hallucination and Fabrication Check

Check whether claims in the output are supported by the actual computed results.

- **Uncomputed statistics**: Does the final report cite numbers, p-values, or effect
  sizes that do not appear in any code output? Cross-reference every number in the
  final summary against the execution trace.
- **Invented identifiers**: Are all gene names, protein IDs, database accessions, and
  pathway names real? Flag any identifier that was not retrieved from a database or
  computed from data.
- **Literature claims**: Does the analysis assert novelty ("has not been previously
  described," "first to show") without performing a literature search? Flag unverified
  novelty claims as a distinct category from incorrect claims.
- **Manually injected data**: Was any data point added manually (hardcoded) rather than
  derived from the analysis pipeline? If so, is this disclosed and justified?
- **Source attribution**: For every factual claim, can you trace it to either (a) a code
  output in the trace, or (b) a cited reference? Flag unsourced claims.

### 6. Scientific Rigor

Domain-specific checks for biological and computational analyses.

- **Data source appropriateness**: Is the data source appropriate for the biological
  question? (Example: bulk tissue measurements cannot support cell-type-specific
  claims. A functional association database does not establish physical binding.)
- **Confounders**: Are obvious confounding variables acknowledged? (Example: cell
  composition differences in bulk samples; batch effects in multi-cohort analyses;
  population stratification in genetic studies.)
- **Biological plausibility**: Do the top candidates make biological sense, or do they
  suggest a systematic bias? (Example: if all top hits belong to one functional
  category unrelated to the target biology, the method may be detecting a confound
  rather than the signal of interest.)
- **Validation**: Were computational predictions validated by any orthogonal evidence
  (literature, independent dataset, functional annotation)?
- **Effect size vs. significance**: Were results evaluated for practical significance,
  not just statistical significance?
- **Model performance context**: For ML models, is performance (AUC, accuracy, etc.)
  contextualized against baselines and the difficulty of the task?

### 7. Reproducibility

Check whether the analysis could be independently reproduced.

- **Parameter documentation**: Are all parameters, thresholds, and random seeds recorded?
- **Data availability**: Are input data sources specified with versions or accession dates?
- **Software versions**: Are package versions recorded for non-standard dependencies?
- **Intermediate outputs**: Are key intermediate results (not just final outputs) preserved?


## Severity Classification

Every finding must be assigned exactly one severity level:

| Severity | Definition | Criteria |
|----------|-----------|----------|
| Critical | Invalidates a primary conclusion | A main claim of the analysis is unsupported, fabricated, or logically unfounded |
| Major | Materially affects interpretation | A significant methodological flaw, unjustified threshold, or logical gap that changes what can be concluded, but does not fully invalidate the work |
| Minor | Should be fixed but does not change conclusions | Cosmetic issues, incomplete documentation, missing edge cases, or limitations that are real but acknowledged |
| Note | Observation for improvement | Suggestions for better practice that do not affect current results |


## Capability Boundary Declarations

For every category, the reviewer must state what it CANNOT check.
Do not silently skip a check. Instead, declare:

```
CANNOT VERIFY: [specific claim]. Reason: [would require literature search /
wet-lab validation / access to raw data / domain expertise in X].
```

This is a distinct output from "no issues found." "No issues found" means the check
was performed and passed. "CANNOT VERIFY" means the check could not be performed.


## Output Format

The review produces a single structured report. Use this exact structure:

```
# Analysis Review: [Brief title]

## Session Summary
- **Original question**: [one sentence]
- **Analysis performed**: [two to three sentences]
- **Primary conclusions**: [numbered list of main claims]

## Findings

### [SEVERITY]: [Short title]
**Category**: [which review category]
**Evidence**: [specific code block, output, or claim that triggered this finding]
**Impact**: [what this means for the conclusions]
**Recommendation**: [specific action to resolve]

[Repeat for each finding, ordered by severity: Critical -> Major -> Minor -> Note]

## Capability Boundaries
[List of CANNOT VERIFY declarations]

## Verdict

**Overall assessment**: [One of: VALID, VALID WITH CAVEATS, MAJOR REVISION NEEDED, UNRELIABLE]

- VALID: No critical or major issues. Conclusions are supported.
- VALID WITH CAVEATS: No critical issues. Major issues exist but are acknowledged
  limitations, not errors. Conclusions hold with stated caveats.
- MAJOR REVISION NEEDED: Major issues that change interpretation. Analysis should
  be revised before conclusions are reported.
- UNRELIABLE: Critical issues found. Primary conclusions are not supported by the
  analysis as conducted.

**Summary**: [Two to three sentences stating what can and cannot be concluded from
this analysis.]
```


## Behavioral Rules

1. Be specific. "The threshold is arbitrary" is not useful. "The filtering step
   uses a cutoff of X, which excludes N data points (including borderline cases
   at X+0.1); no justification is provided for this value versus alternative
   cutoffs" is useful.

2. Distinguish between errors and limitations. An error is something the analysis got
   wrong. A limitation is something the analysis cannot address given available data.
   Both matter, but they require different responses (fix vs. acknowledge).

3. Do not grade on effort. An analysis that tried many things but reached unsupported
   conclusions is worse than a simple analysis with well-supported conclusions.

4. Do not invent issues. If a category has no findings, say "No issues identified"
   for that category. Do not manufacture minor concerns to appear thorough.

5. Check the obvious first. Before looking for subtle issues, verify: Did the code
   actually run? Did it load the right data? Do the numbers in the summary match the
   numbers in the output? Most serious errors are mundane, not exotic.

6. Trace claims backward. Start from the final conclusions and trace each one back
   through the analysis to its data source. This catches more issues than reading
   the analysis forward.

7. When the session is empty or contains no analysis, state this in one sentence.
   Do not produce the full report template with N/A entries.
