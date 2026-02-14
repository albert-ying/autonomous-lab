/* ==============================================================
   AUTONOMOUS LAB — Website Main Script
   Tutorial, Marketplace, Navigation, Character System
   ============================================================== */

(function () {
  "use strict";

  // ============================================================
  // DOMAIN DEFINITIONS — Six configurable domains
  // Each domain defines role labels, tutorial content, and UI text
  // ============================================================
  let activeDomain = "research";

  const DOMAINS = {
    research: {
      label: "Research Lab", icon: "&#x1F52C;",
      senior: { label: "Professor", short: "PI", icon: "&#x1F52C;" },
      junior: { label: "Trainee", short: "Trainee", icon: "&#x1F9EA;" },
      overseer: { label: "Editor", short: "Editor", icon: "&#x1F9D1;&#x200D;&#x1F52C;" },
      artifact: "Paper", workspace: "Lab", meetingLog: "Meeting Log",
      inventory: ["data", "scripts", "figures", "results", "paper"],
      artifactSections: ["abstract", "introduction", "methods", "results", "discussion"],
      reviewProcess: "Peer Review", consultants: "Consultants",
      decisions: { accept: "Accept", minor: "Minor Revision", major: "Major Revision", reject: "Reject" },
      heroSubtitle: "You did the research. Now you judge it.",
      heroDescription: 'AI takes the <strong class="dyn-senior">PI</strong> and <strong class="dyn-junior">Trainee</strong> roles\u200a\u2014\u200adesigning experiments, running analysis, writing the paper. You become the <strong class="dyn-overseer">Editor</strong>: the sole decision maker who curates, judges, and steers the science.',
      expert: { name: "Dr. Tanaka", role: "Statistician", avatar: "statistician" },
      reviewers: [
        { name: "Dr. Kim", role: "Comp. Biologist", avatar: "comp_biologist" },
        { name: "Dr. Patel", role: "Immunologist", avatar: "immunologist" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Professor", meta: "Iteration 1",
          text: '<strong>Setting the research agenda.</strong> Welcome to the lab! I\'ve reviewed the literature and identified a promising direction. Let\'s investigate the relationship between <strong>single-cell chromatin accessibility</strong> and <strong>gene expression</strong> in tumor microenvironments. Our goal: identify cell-type-specific regulatory programs that drive treatment resistance.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Trainee" },
        { role: "trainee", label: "Trainee", meta: "Iteration 1",
          text: '<strong>Executing the analysis pipeline.</strong> Project structure set up. Downloaded the scRNA-seq dataset from GEO, ran QC with <strong>scanpy</strong>, performed normalization and batch correction using <strong>scvi-tools</strong>. The initial UMAP shows 12 distinct clusters across 45,000 cells.',
          inv: { data: 2, scripts: 3, figures: 1, results: 1, paper: 0 }, paper: 10, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "PI" },
        { role: "pi", label: "Professor", meta: "Iteration 2",
          text: '<strong>Reviewing results and consulting an expert.</strong> Good first pass. The clustering looks clean, but I want to validate the cell type annotations. Let me consult <strong>Dr. Tanaka</strong> (Statistician) on the differential expression approach. Also, let\'s add a survival analysis component.',
          inv: { data: 2, scripts: 3, figures: 1, results: 1, paper: 0 }, paper: 15, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Trainee" },
        { role: "expert", label: "Dr. Tanaka (Statistician)", meta: "Consultation",
          text: '<strong>Expert consultation.</strong> For DE analysis, use a <strong>pseudo-bulk approach</strong> with pyDESeq2 rather than per-cell tests&mdash;this avoids inflated p-values from pseudoreplication. For survival, use a <strong>Cox proportional hazards model</strong> with C-index. Adjust for multiple comparisons with Benjamini-Hochberg.',
          inv: { data: 2, scripts: 3, figures: 1, results: 1, paper: 0 }, paper: 15, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Trainee" },
        { role: "trainee", label: "Trainee", meta: "Iteration 3",
          text: '<strong>Completing analysis and drafting the paper.</strong> All analyses done. Pseudo-bulk DE identified <strong>347 significant genes</strong>. Survival model shows C-index = 0.78. All 4 figures complete. Full LaTeX draft written: abstract through discussion.',
          inv: { data: 3, scripts: 6, figures: 4, results: 5, paper: 6 }, paper: 85,
          sections: { abstract: "250w", intro: "1,200w", methods: "900w", results: "1,500w", discussion: "800w" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "PI" },
        { role: "pi", label: "Professor", meta: "Submission",
          text: '<strong>Submitting manuscript for editorial review.</strong> The paper is ready. I\'ve written a cover letter highlighting the key contributions and novelty. Submitting to the Editor (you) for review. The editorial office is now open.',
          inv: { data: 3, scripts: 6, figures: 4, results: 5, paper: 6 }, paper: 90,
          sections: { abstract: "250w", intro: "1,200w", methods: "900w", results: "1,500w", discussion: "800w" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "Editor" },
        { role: "editor", label: "You (Editor)", meta: "Editorial Office",
          text: '<strong>Your turn as Editor!</strong> You open the <strong>Editorial Office</strong> overlay. You can <strong>Desk Reject</strong> the manuscript or <strong>Send to Reviewers</strong>. You choose to invite two AI reviewers. The reviews are generated automatically.',
          inv: { data: 3, scripts: 6, figures: 4, results: 5, paper: 6 }, paper: 90,
          sections: { abstract: "250w", intro: "1,200w", methods: "900w", results: "1,500w", discussion: "800w" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "reviewer", label: "Reviewer #1 (Comp. Biologist)", meta: "Peer Review",
          text: '<strong>Strengths:</strong> Novel integration of chromatin accessibility with expression data. Solid statistical framework using pseudo-bulk DE. <strong>Weaknesses:</strong> The survival model lacks external validation. Figure 3 needs higher resolution. <strong>Recommendation: Minor Revision</strong>&mdash;address the validation cohort and figure quality.',
          inv: { data: 3, scripts: 6, figures: 4, results: 5, paper: 6 }, paper: 90,
          sections: { abstract: "250w", intro: "1,200w", methods: "900w", results: "1,500w", discussion: "800w" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Editor" },
        { role: "editor", label: "You (Editor)", meta: "Decision",
          text: '<strong>Editorial decision: Minor Revision.</strong> Based on the reviewer reports, you write editorial feedback requesting the authors add an external validation cohort and improve figure resolution. The decision is sent back to the PI, and the loop continues with revisions. <em>The cycle repeats until you Accept.</em>',
          inv: { data: 3, scripts: 6, figures: 4, results: 5, paper: 6 }, paper: 92,
          sections: { abstract: "250w", intro: "1,200w", methods: "900w", results: "1,500w", discussion: "800w" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "PI" },
      ],
    },
    software: {
      label: "Software Engineering", icon: "&#x1F4BB;",
      senior: { label: "Tech Lead", short: "Lead", icon: "&#x1F468;&#x200D;&#x1F4BB;" },
      junior: { label: "Developer", short: "Dev", icon: "&#x1F469;&#x200D;&#x1F4BB;" },
      overseer: { label: "Code Reviewer", short: "Reviewer", icon: "&#x1F50D;" },
      artifact: "Pull Request", workspace: "Sprint", meetingLog: "Standup Log",
      inventory: ["src", "tests", "docs", "configs", "PR"],
      artifactSections: ["summary", "architecture", "implementation", "test results", "review notes"],
      reviewProcess: "Code Review", consultants: "Specialists",
      decisions: { accept: "Approve", minor: "Request Changes", major: "Needs Redesign", reject: "Reject" },
      heroSubtitle: "You wrote the code. Now you review it.",
      heroDescription: 'AI takes the <strong class="dyn-senior">Tech Lead</strong> and <strong class="dyn-junior">Developer</strong> roles\u200a\u2014\u200aarchitecting systems, writing code, shipping features. You become the <strong class="dyn-overseer">Code Reviewer</strong>: the sole decision maker who approves, shapes, and directs.',
      expert: { name: "Alex Chen", role: "DevOps Engineer", avatar: "engineer" },
      reviewers: [
        { name: "Sarah K.", role: "Security Engineer", avatar: "engineer" },
        { name: "Mike R.", role: "Backend Lead", avatar: "data_scientist" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Tech Lead", meta: "Sprint 1",
          text: '<strong>Setting the sprint agenda.</strong> We need to build a <strong>real-time notification service</strong> for our platform. Requirements: WebSocket connections, message queue integration with Redis, rate limiting, and a REST API for notification preferences. Let\'s target 99.9% delivery rate.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Dev" },
        { role: "trainee", label: "Developer", meta: "Sprint 1",
          text: '<strong>Implementing the core service.</strong> Project scaffolded with <strong>FastAPI</strong> and <strong>Redis</strong>. WebSocket handler supports connection pooling for 10K concurrent users. REST endpoints for CRUD on notification preferences. Initial test suite with 23 unit tests passing.',
          inv: { data: 3, scripts: 5, figures: 0, results: 2, paper: 0 }, paper: 15, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "Lead" },
        { role: "pi", label: "Tech Lead", meta: "Sprint 2",
          text: '<strong>Reviewing implementation.</strong> Solid foundation. The WebSocket handler needs <strong>graceful reconnection</strong> logic. Let me consult <strong>Alex Chen</strong> (DevOps) on the deployment architecture. Also, add load testing with Locust before we ship.',
          inv: { data: 3, scripts: 5, figures: 0, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Dev" },
        { role: "expert", label: "Alex Chen (DevOps Engineer)", meta: "Consultation",
          text: '<strong>Infrastructure review.</strong> Use <strong>Kubernetes HPA</strong> for auto-scaling the WebSocket pods. Redis Cluster with 3 nodes for high availability. Add a <strong>circuit breaker</strong> pattern for downstream service failures. Prometheus metrics for connection count and delivery latency.',
          inv: { data: 3, scripts: 5, figures: 0, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Dev" },
        { role: "trainee", label: "Developer", meta: "Sprint 3",
          text: '<strong>Feature complete.</strong> Added reconnection logic, circuit breaker, and Prometheus metrics. Load test shows <strong>12K concurrent connections</strong> at p99 latency of 45ms. 47 tests passing, 92% coverage. PR description and architecture doc written.',
          inv: { data: 4, scripts: 8, figures: 2, results: 5, paper: 5 }, paper: 85,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "Lead" },
        { role: "pi", label: "Tech Lead", meta: "PR Submission",
          text: '<strong>Submitting PR for code review.</strong> The implementation looks solid. I\'ve written the PR summary highlighting the architecture decisions and performance benchmarks. Submitting to the Code Reviewer (you) for approval.',
          inv: { data: 4, scripts: 8, figures: 2, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "editor", label: "You (Code Reviewer)", meta: "Code Review",
          text: '<strong>Your turn as Code Reviewer!</strong> You open the <strong>Review Board</strong>. You can <strong>Reject</strong> the PR or <strong>Send to Reviewers</strong>. You choose to invite a Security Engineer and a Backend Lead for thorough review.',
          inv: { data: 4, scripts: 8, figures: 2, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "reviewer", label: "Reviewer #1 (Security Engineer)", meta: "Code Review",
          text: '<strong>Security review.</strong> WebSocket auth looks solid. Rate limiting is properly implemented. <strong>Issue:</strong> The Redis connection string is logged in debug mode&mdash;needs to be masked. Add input validation on the notification payload. <strong>Recommendation: Request Changes.</strong>',
          inv: { data: 4, scripts: 8, figures: 2, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "editor", label: "You (Code Reviewer)", meta: "Decision",
          text: '<strong>Decision: Request Changes.</strong> Based on the security review, you request the team mask sensitive log data and add payload validation. The feedback is sent to the Tech Lead, and the sprint continues with fixes. <em>The cycle repeats until you Approve.</em>',
          inv: { data: 4, scripts: 8, figures: 2, results: 5, paper: 5 }, paper: 92,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "Lead" },
      ],
    },
    consulting: {
      label: "Consulting", icon: "&#x1F4BC;",
      senior: { label: "Partner", short: "Partner", icon: "&#x1F454;" },
      junior: { label: "Associate", short: "Associate", icon: "&#x1F4CA;" },
      overseer: { label: "Client", short: "Client", icon: "&#x1F465;" },
      artifact: "Strategy Report", workspace: "Engagement", meetingLog: "Engagement Log",
      inventory: ["data", "analysis", "slides", "exhibits", "report"],
      artifactSections: ["executive summary", "market analysis", "recommendations", "financials", "appendix"],
      reviewProcess: "Client Review", consultants: "Subject Experts",
      decisions: { accept: "Approve", minor: "Revise", major: "Rework", reject: "Reject" },
      heroSubtitle: "You ran the analysis. Now you direct it.",
      heroDescription: 'AI takes the <strong class="dyn-senior">Partner</strong> and <strong class="dyn-junior">Associate</strong> roles\u200a\u2014\u200abuilding strategy, crunching data, crafting deliverables. You become the <strong class="dyn-overseer">Client</strong>: the sole decision maker who steers vision and approves.',
      expert: { name: "Rachel Torres", role: "Industry Analyst", avatar: "data_scientist" },
      reviewers: [
        { name: "David L.", role: "Finance Director", avatar: "statistician" },
        { name: "Nina S.", role: "Strategy VP", avatar: "data_scientist" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Partner", meta: "Week 1",
          text: '<strong>Scoping the engagement.</strong> A fintech startup needs a <strong>market entry strategy</strong> for Southeast Asia. We need competitive landscape analysis, regulatory mapping across 5 markets, and a go-to-market recommendation with financial projections. Target: $50M revenue in 3 years.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Associate" },
        { role: "trainee", label: "Associate", meta: "Week 1",
          text: '<strong>Market research completed.</strong> Analyzed 23 competitors across Singapore, Indonesia, Thailand, Vietnam, and Philippines. Built a <strong>regulatory complexity matrix</strong> scoring each market. Identified Singapore as the beachhead with lowest barriers. Initial financial model projects $12M year-one revenue.',
          inv: { data: 3, scripts: 2, figures: 3, results: 2, paper: 0 }, paper: 15, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "Partner" },
        { role: "pi", label: "Partner", meta: "Week 2",
          text: '<strong>Reviewing analysis.</strong> Strong competitive mapping. The financial model needs <strong>sensitivity analysis</strong> on customer acquisition cost. Let me bring in <strong>Rachel Torres</strong> (Industry Analyst) for deeper market sizing. Also, the Indonesia opportunity may be undervalued.',
          inv: { data: 3, scripts: 2, figures: 3, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Associate" },
        { role: "expert", label: "Rachel Torres (Industry Analyst)", meta: "Consultation",
          text: '<strong>Market insight.</strong> Southeast Asian digital payments grew <strong>42% YoY</strong>. Indonesia\'s new OJK regulations actually favor foreign fintechs with banking partnerships. Recommend a <strong>dual-market launch</strong>: Singapore for credibility, Indonesia for volume. Use a <strong>TAM/SAM/SOM framework</strong> for each market.',
          inv: { data: 3, scripts: 2, figures: 3, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Associate" },
        { role: "trainee", label: "Associate", meta: "Week 3",
          text: '<strong>Report drafted.</strong> Complete strategy report with dual-market approach. Sensitivity analysis shows <strong>break-even at 18 months</strong> under conservative scenario. 8 exhibit slides, executive summary, and 3-year financial model complete. All recommendations backed by market data.',
          inv: { data: 5, scripts: 4, figures: 8, results: 4, paper: 5 }, paper: 85,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "Partner" },
        { role: "pi", label: "Partner", meta: "Delivery",
          text: '<strong>Submitting to client.</strong> The strategy report is compelling. I\'ve prepared an executive briefing highlighting the dual-market approach and key financials. Presenting to the Client (you) for approval.',
          inv: { data: 5, scripts: 4, figures: 8, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "Client" },
        { role: "editor", label: "You (Client)", meta: "Client Review",
          text: '<strong>Your turn as Client!</strong> You review the deliverable package. You can <strong>Reject</strong> the proposal or <strong>Send to Advisors</strong> for additional review. You choose to get feedback from your Finance Director and Strategy VP.',
          inv: { data: 5, scripts: 4, figures: 8, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Advisor" },
        { role: "reviewer", label: "Advisor #1 (Finance Director)", meta: "Client Review",
          text: '<strong>Financial review.</strong> Revenue projections are realistic. <strong>Concern:</strong> The CAC assumptions in Indonesia need a 20% buffer for regulatory compliance costs. Recommend adding a <strong>scenario planning section</strong> with bull/bear cases. <strong>Recommendation: Revise.</strong>',
          inv: { data: 5, scripts: 4, figures: 8, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Client" },
        { role: "editor", label: "You (Client)", meta: "Decision",
          text: '<strong>Decision: Revise.</strong> You ask the team to add scenario planning and adjust the Indonesia CAC projections. The feedback is sent to the Partner, and the engagement continues. <em>The cycle repeats until you Approve.</em>',
          inv: { data: 5, scripts: 4, figures: 8, results: 4, paper: 5 }, paper: 92,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "Partner" },
      ],
    },
    legal: {
      label: "Legal", icon: "&#x2696;&#xFE0F;",
      senior: { label: "Senior Partner", short: "Partner", icon: "&#x2696;&#xFE0F;" },
      junior: { label: "Associate", short: "Associate", icon: "&#x1F4D1;" },
      overseer: { label: "Reviewing Partner", short: "Reviewer", icon: "&#x1F9D1;&#x200D;&#x2696;&#xFE0F;" },
      artifact: "Legal Brief", workspace: "Case", meetingLog: "Case Log",
      inventory: ["evidence", "research", "drafts", "exhibits", "brief"],
      artifactSections: ["summary", "facts", "legal analysis", "arguments", "conclusion"],
      reviewProcess: "Partner Review", consultants: "Specialists",
      decisions: { accept: "File", minor: "Minor Edits", major: "Major Rewrite", reject: "Withdraw" },
      heroSubtitle: "You drafted the briefs. Now you rule on them.",
      heroDescription: 'AI takes the <strong class="dyn-senior">Senior Partner</strong> and <strong class="dyn-junior">Associate</strong> roles\u200a\u2014\u200aresearching precedent, drafting motions, building arguments. You become the <strong class="dyn-overseer">Reviewing Partner</strong>: the sole decision maker who rules on strategy.',
      expert: { name: "Prof. Hayes", role: "Regulatory Specialist", avatar: "generic" },
      reviewers: [
        { name: "M. Walsh", role: "IP Counsel", avatar: "generic" },
        { name: "J. Reeves", role: "Litigation Expert", avatar: "generic" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Senior Partner", meta: "Phase 1",
          text: '<strong>Setting case strategy.</strong> Our client faces a <strong>patent infringement claim</strong> in the Eastern District. We need to build a defense around prior art, conduct claim construction analysis, and prepare an invalidity argument. Discovery deadline is in 8 weeks.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Associate" },
        { role: "trainee", label: "Associate", meta: "Phase 1",
          text: '<strong>Legal research completed.</strong> Identified <strong>12 prior art references</strong> predating the patent filing. Claim construction memo drafted for the 3 contested claims. Found a key <strong>prosecution history estoppel</strong> argument from the patent\'s file wrapper.',
          inv: { data: 4, scripts: 3, figures: 0, results: 2, paper: 0 }, paper: 15, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "Partner" },
        { role: "pi", label: "Senior Partner", meta: "Phase 2",
          text: '<strong>Reviewing research.</strong> Strong prior art collection. The estoppel argument is promising. Let me consult <strong>Prof. Hayes</strong> (Regulatory Specialist) on the technical standard implications. Also, we need expert declarations for the invalidity motion.',
          inv: { data: 4, scripts: 3, figures: 0, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Associate" },
        { role: "expert", label: "Prof. Hayes (Regulatory Specialist)", meta: "Consultation",
          text: '<strong>Technical analysis.</strong> The patent\'s key claims overlap with the <strong>IEEE 802.11 standard</strong> published 2 years earlier. This is strong prior art. Recommend focusing the invalidity argument on <strong>anticipation under 35 U.S.C. 102</strong> rather than obviousness. The standard\'s technical committee minutes are publicly available.',
          inv: { data: 4, scripts: 3, figures: 0, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Associate" },
        { role: "trainee", label: "Associate", meta: "Phase 3",
          text: '<strong>Brief drafted.</strong> Complete invalidity brief with <strong>anticipation and estoppel arguments</strong>. Expert declaration from Dr. Williams attached. Claim charts mapping prior art to all contested claims. 42-page brief with 28 exhibits.',
          inv: { data: 6, scripts: 5, figures: 3, results: 4, paper: 6 }, paper: 85,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "Partner" },
        { role: "pi", label: "Senior Partner", meta: "Filing",
          text: '<strong>Submitting for partner review.</strong> The brief is thorough and well-argued. Submitting to the Reviewing Partner (you) for final approval before filing with the court.',
          inv: { data: 6, scripts: 5, figures: 3, results: 4, paper: 6 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "editor", label: "You (Reviewing Partner)", meta: "Partner Review",
          text: '<strong>Your turn as Reviewing Partner!</strong> You review the filing package. You can <strong>Withdraw</strong> the motion or <strong>Send to Counsel</strong> for additional review. You invite IP Counsel and a Litigation Expert.',
          inv: { data: 6, scripts: 5, figures: 3, results: 4, paper: 6 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "reviewer", label: "Reviewer #1 (IP Counsel)", meta: "Partner Review",
          text: '<strong>IP review.</strong> Anticipation argument is strong. <strong>Concern:</strong> Claim 7 has a dependent limitation not fully addressed by the IEEE standard. Need to add an <strong>obviousness argument</strong> as a fallback for that claim. <strong>Recommendation: Minor Edits.</strong>',
          inv: { data: 6, scripts: 5, figures: 3, results: 4, paper: 6 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "editor", label: "You (Reviewing Partner)", meta: "Decision",
          text: '<strong>Decision: Minor Edits.</strong> You request the team add an obviousness fallback for Claim 7 and strengthen the expert declaration. The feedback is sent to the Senior Partner. <em>The cycle repeats until you approve for Filing.</em>',
          inv: { data: 6, scripts: 5, figures: 3, results: 4, paper: 6 }, paper: 92,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "Partner" },
      ],
    },
    medical: {
      label: "Medical", icon: "&#x1FA7A;",
      senior: { label: "Attending", short: "Attending", icon: "&#x1F468;&#x200D;&#x2695;&#xFE0F;" },
      junior: { label: "Resident", short: "Resident", icon: "&#x1FA7A;" },
      overseer: { label: "Chief", short: "Chief", icon: "&#x1F3E5;" },
      artifact: "Treatment Plan", workspace: "Rounds", meetingLog: "Rounds Log",
      inventory: ["records", "labs", "imaging", "notes", "plan"],
      artifactSections: ["assessment", "history", "workup", "treatment", "follow-up"],
      reviewProcess: "Peer Review", consultants: "Specialists",
      decisions: { accept: "Approve", minor: "Modify Plan", major: "Reassess", reject: "Override" },
      heroSubtitle: "You ran the rounds. Now you call the shots.",
      heroDescription: 'AI takes the <strong class="dyn-senior">Attending</strong> and <strong class="dyn-junior">Resident</strong> roles\u200a\u2014\u200adiagnosing, treating, and documenting care. You become the <strong class="dyn-overseer">Chief</strong>: the sole decision maker who oversees and approves.',
      expert: { name: "Dr. Nakamura", role: "Cardiologist", avatar: "clinician" },
      reviewers: [
        { name: "Dr. Park", role: "Pharmacologist", avatar: "pharmacologist" },
        { name: "Dr. Singh", role: "Radiologist", avatar: "radiologist" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Attending", meta: "Day 1",
          text: '<strong>Setting the treatment plan.</strong> 58-year-old patient presenting with <strong>acute chest pain</strong>, elevated troponin (0.8 ng/mL), and ST-segment changes on ECG. We need a full cardiac workup: serial troponins, echocardiogram, coronary CT angiography. Initiate antiplatelet therapy and continuous monitoring.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Resident" },
        { role: "trainee", label: "Resident", meta: "Day 1",
          text: '<strong>Workup completed.</strong> Serial troponins trending up: 0.8 &rarr; 1.2 &rarr; 2.1 ng/mL. Echo shows <strong>EF 45%</strong> with anterior wall hypokinesis. CT angiography reveals <strong>80% LAD stenosis</strong>. Started dual antiplatelet, heparin drip, statin. Vitals stable on telemetry.',
          inv: { data: 3, scripts: 2, figures: 2, results: 3, paper: 0 }, paper: 15, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "Attending" },
        { role: "pi", label: "Attending", meta: "Day 2",
          text: '<strong>Reviewing results.</strong> Troponin trend and CT findings confirm <strong>NSTEMI with significant LAD disease</strong>. Let me consult <strong>Dr. Nakamura</strong> (Cardiology) on intervention timing. We need to risk-stratify with GRACE score and decide between PCI and medical management.',
          inv: { data: 3, scripts: 2, figures: 2, results: 3, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Resident" },
        { role: "expert", label: "Dr. Nakamura (Cardiologist)", meta: "Consultation",
          text: '<strong>Cardiology consult.</strong> GRACE score 142 = <strong>high-risk category</strong>. Recommend <strong>early invasive strategy</strong> with cardiac catheterization within 24 hours. Given the EF of 45%, start an ACE inhibitor. If PCI is successful, plan discharge on day 3-4 with cardiac rehab referral.',
          inv: { data: 3, scripts: 2, figures: 2, results: 3, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Resident" },
        { role: "trainee", label: "Resident", meta: "Day 3",
          text: '<strong>Treatment plan complete.</strong> Cath lab: successful <strong>PCI with drug-eluting stent</strong> to LAD. Post-procedure EF improved to 50%. Discharge plan: DAPT for 12 months, high-intensity statin, ACE inhibitor, beta-blocker. Cardiac rehab scheduled. Full documentation completed.',
          inv: { data: 5, scripts: 4, figures: 4, results: 5, paper: 5 }, paper: 85,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "Attending" },
        { role: "pi", label: "Attending", meta: "Discharge",
          text: '<strong>Submitting for peer review.</strong> The treatment plan and documentation are thorough. Submitting to the Chief (you) for final approval before patient discharge.',
          inv: { data: 5, scripts: 4, figures: 4, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "Chief" },
        { role: "editor", label: "You (Chief)", meta: "Peer Review",
          text: '<strong>Your turn as Chief!</strong> You review the case and treatment plan. You can <strong>Override</strong> the plan or <strong>Send to Committee</strong> for peer review. You request review from a Pharmacologist and Radiologist.',
          inv: { data: 5, scripts: 4, figures: 4, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "reviewer", label: "Reviewer #1 (Pharmacologist)", meta: "Peer Review",
          text: '<strong>Medication review.</strong> Appropriate DAPT and statin regimen. <strong>Concern:</strong> Check renal function before ACE inhibitor dosing&mdash;creatinine was borderline at admission. Consider adjusting to a lower starting dose. <strong>Recommendation: Modify Plan.</strong>',
          inv: { data: 5, scripts: 4, figures: 4, results: 5, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Chief" },
        { role: "editor", label: "You (Chief)", meta: "Decision",
          text: '<strong>Decision: Modify Plan.</strong> You request a renal function recheck and ACE inhibitor dose adjustment. The feedback is sent to the Attending. <em>The cycle repeats until you Approve the discharge.</em>',
          inv: { data: 5, scripts: 4, figures: 4, results: 5, paper: 5 }, paper: 92,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "Attending" },
      ],
    },
    creative: {
      label: "Creative", icon: "&#x1F3A8;",
      senior: { label: "Art Director", short: "AD", icon: "&#x1F3A8;" },
      junior: { label: "Designer", short: "Designer", icon: "&#x270F;&#xFE0F;" },
      overseer: { label: "Creative Director", short: "CD", icon: "&#x1F451;" },
      artifact: "Campaign", workspace: "Studio", meetingLog: "Creative Brief",
      inventory: ["assets", "sketches", "mockups", "finals", "deck"],
      artifactSections: ["concept", "moodboard", "visual system", "deliverables", "guidelines"],
      reviewProcess: "Creative Review", consultants: "Specialists",
      decisions: { accept: "Approve", minor: "Polish", major: "Rework", reject: "Scrap" },
      heroSubtitle: "You made the art. Now you curate it.",
      heroDescription: 'AI takes the <strong class="dyn-senior">Art Director</strong> and <strong class="dyn-junior">Designer</strong> roles\u200a\u2014\u200aconcepting, creating, and iterating on work. You become the <strong class="dyn-overseer">Creative Director</strong>: the sole decision maker who shapes the vision.',
      expert: { name: "Lena Novak", role: "Brand Strategist", avatar: "psychologist" },
      reviewers: [
        { name: "Tom H.", role: "Copywriter", avatar: "science_writer" },
        { name: "Yuki M.", role: "UX Lead", avatar: "generic" },
      ],
      tutorialSteps: [
        { role: "pi", label: "Art Director", meta: "Week 1",
          text: '<strong>Setting the creative brief.</strong> A sustainable fashion brand needs a <strong>complete visual rebrand</strong>. They\'re pivoting from luxury minimalism to bold eco-activism. Deliverables: new logo system, color palette, typography, social media templates, and a hero campaign for Earth Day launch.',
          inv: { data: 0, scripts: 0, figures: 0, results: 0, paper: 0 }, paper: 0, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 1, nextRole: "Designer" },
        { role: "trainee", label: "Designer", meta: "Week 1",
          text: '<strong>Initial concepts delivered.</strong> Created 3 logo directions exploring <strong>organic forms, geometric boldness, and typographic play</strong>. Moodboard with 40 reference images establishing the eco-activist aesthetic. Color palette: deep forest, sunset orange, recycled paper beige. Two font pairings tested.',
          inv: { data: 5, scripts: 0, figures: 3, results: 2, paper: 0 }, paper: 15, sections: {},
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 1, nextRole: "AD" },
        { role: "pi", label: "Art Director", meta: "Week 2",
          text: '<strong>Reviewing concepts.</strong> Direction B (geometric boldness) has the most potential. The logo needs more <strong>weight and confidence</strong>. Let me bring in <strong>Lena Novak</strong> (Brand Strategist) to validate the positioning against competitors. Also, the color palette needs a vibrant accent.',
          inv: { data: 5, scripts: 0, figures: 3, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 2, nextRole: "Designer" },
        { role: "expert", label: "Lena Novak (Brand Strategist)", meta: "Consultation",
          text: '<strong>Brand analysis.</strong> The geometric direction differentiates from competitors who are all doing <strong>soft, pastel eco-aesthetics</strong>. Add an <strong>electric lime accent</strong> for Gen-Z appeal. The logo should work as both a mark and a pattern system. Recommend a <strong>"visual activism" manifesto</strong> as the campaign\'s anchor.',
          inv: { data: 5, scripts: 0, figures: 3, results: 2, paper: 0 }, paper: 25, sections: {},
          piActive: true, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 2, nextRole: "Designer" },
        { role: "trainee", label: "Designer", meta: "Week 3",
          text: '<strong>Campaign complete.</strong> Final logo system with <strong>12 responsive variants</strong>. Complete visual identity: colors, typography, iconography, and pattern library. Earth Day hero campaign with 6 social templates, 2 billboard concepts, and a brand guidelines PDF. All files production-ready.',
          inv: { data: 8, scripts: 2, figures: 6, results: 4, paper: 5 }, paper: 85,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: true, piStatus: "Idle", traineeStatus: "speaking", iteration: 3, nextRole: "AD" },
        { role: "pi", label: "Art Director", meta: "Presentation",
          text: '<strong>Submitting for creative review.</strong> The campaign is strong and cohesive. I\'ve prepared the presentation deck with the brand story, visual system, and campaign concepts. Presenting to the Creative Director (you) for approval.',
          inv: { data: 8, scripts: 2, figures: 6, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: true, traineeActive: false, piStatus: "speaking", traineeStatus: "Idle", iteration: 4, nextRole: "CD" },
        { role: "editor", label: "You (Creative Director)", meta: "Creative Review",
          text: '<strong>Your turn as Creative Director!</strong> You review the brand package and campaign. You can <strong>Scrap</strong> the direction or <strong>Send to Review</strong>. You invite a Copywriter and UX Lead for feedback.',
          inv: { data: 8, scripts: 2, figures: 6, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "Reviewer" },
        { role: "reviewer", label: "Reviewer #1 (Copywriter)", meta: "Creative Review",
          text: '<strong>Copy review.</strong> Visual identity is powerful. <strong>Concern:</strong> The manifesto tone feels too aggressive for the target demo&mdash;soften the call-to-action language while keeping the bold aesthetic. Social templates need <strong>more whitespace</strong> for copy integration. <strong>Recommendation: Polish.</strong>',
          inv: { data: 8, scripts: 2, figures: 6, results: 4, paper: 5 }, paper: 90,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 4, nextRole: "CD" },
        { role: "editor", label: "You (Creative Director)", meta: "Decision",
          text: '<strong>Decision: Polish.</strong> You request softer manifesto language and more whitespace in social templates. The feedback is sent to the Art Director. <em>The cycle repeats until you Approve the campaign.</em>',
          inv: { data: 8, scripts: 2, figures: 6, results: 4, paper: 5 }, paper: 92,
          sections: { abstract: "done", intro: "done", methods: "done", results: "done", discussion: "done" },
          piActive: false, traineeActive: false, piStatus: "Idle", traineeStatus: "Idle", iteration: 5, nextRole: "AD" },
      ],
    },
  };

  // ============================================================
  // DOMAIN EXTRAS — badges, subtitles, character YAML, steps, code
  // Merged into DOMAINS at startup to keep the main object readable
  // ============================================================
  (function () {
    var extras = {
      research: {
        badges: ["&#x1F4DD; LaTeX Papers", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch a research session unfold, step by step",
        characterYaml: {
          filename: "pi.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">Computational Biology PI</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">single-cell genomics, ML</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">discover cell-type-specific\n  gene regulatory programs</span>\n<span class="yaml-key">skills:</span>  <span class="code-comment"># \u2192 skills/ folder</span>\n  <span class="yaml-val">- single-cell-analysis</span>\n  <span class="yaml-val">- multi-modal-integration</span>\n  <span class="yaml-val">- scientific-writing</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Visionary: novel questions"</span>\n  <span class="yaml-val">- "Rigorous: reproducible"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, expertise, and personality traits. These shape how the AI agent approaches problems and communicates." },
          { title: "Attach Cursor Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills used by Cursor and Claude. Includes <em>scanpy</em>, <em>scientific-writing</em>, and 200+ more.' },
          { title: "Deploy to Your Lab", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your research team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The research cycle:</span>\nautolab_next  <span class="code-comment"># PI sets agenda</span>\nautolab_next  <span class="code-comment"># Trainee executes</span>\nautolab_next  <span class="code-comment"># PI reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When paper is ready:</span>\nautolab_editorial <span class="code-comment"># You decide!</span>\n<span class="code-comment"># Accept / Minor / Major / Reject</span>'
      },
      software: {
        badges: ["&#x1F500; Pull Requests", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch a sprint unfold, step by step",
        characterYaml: {
          filename: "lead.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">Full-Stack Tech Lead</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">distributed systems, APIs</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">ship reliable, well-tested\n  software on schedule</span>\n<span class="yaml-key">skills:</span>\n  <span class="yaml-val">- system-design</span>\n  <span class="yaml-val">- code-review</span>\n  <span class="yaml-val">- testing-strategy</span>\n  <span class="yaml-val">- ci-cd-pipelines</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Pragmatic: ships working code"</span>\n  <span class="yaml-val">- "Mentoring: grows the team"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, expertise, and personality traits. These shape how the AI writes and reviews code." },
          { title: "Attach Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills for <em>system-design</em>, <em>testing-strategy</em>, <em>ci-cd</em>, and more.' },
          { title: "Deploy to Your Sprint", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your engineering team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The sprint cycle:</span>\nautolab_next  <span class="code-comment"># Lead plans</span>\nautolab_next  <span class="code-comment"># Dev implements</span>\nautolab_next  <span class="code-comment"># Lead reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When PR is ready:</span>\nautolab_editorial <span class="code-comment"># You review!</span>\n<span class="code-comment"># Approve / Request Changes / Reject</span>'
      },
      consulting: {
        badges: ["&#x1F4CA; Strategy Reports", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch a consulting engagement unfold, step by step",
        characterYaml: {
          filename: "partner.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">Strategy Partner</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">market entry, financial modeling</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">deliver actionable client\n  recommendations</span>\n<span class="yaml-key">skills:</span>\n  <span class="yaml-val">- market-analysis</span>\n  <span class="yaml-val">- financial-modeling</span>\n  <span class="yaml-val">- presentation-design</span>\n  <span class="yaml-val">- stakeholder-mgmt</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Analytical: data-driven"</span>\n  <span class="yaml-val">- "Client-focused: actionable"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, expertise, and personality traits. These shape how the AI approaches client engagements." },
          { title: "Attach Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills for <em>market-analysis</em>, <em>financial-modeling</em>, <em>presentation-design</em>, and more.' },
          { title: "Deploy to Your Engagement", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your consulting team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The engagement cycle:</span>\nautolab_next  <span class="code-comment"># Partner scopes</span>\nautolab_next  <span class="code-comment"># Associate analyzes</span>\nautolab_next  <span class="code-comment"># Partner reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When report is ready:</span>\nautolab_editorial <span class="code-comment"># You decide!</span>\n<span class="code-comment"># Approve / Revise / Rework / Reject</span>'
      },
      legal: {
        badges: ["&#x2696;&#xFE0F; Legal Briefs", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch a legal case unfold, step by step",
        characterYaml: {
          filename: "senior-partner.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">IP Litigation Partner</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">patent litigation, prior art</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">build winning legal\n  strategies for clients</span>\n<span class="yaml-key">skills:</span>\n  <span class="yaml-val">- legal-research</span>\n  <span class="yaml-val">- claim-construction</span>\n  <span class="yaml-val">- brief-writing</span>\n  <span class="yaml-val">- case-analysis</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Meticulous: every detail"</span>\n  <span class="yaml-val">- "Strategic: long-game thinking"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, expertise, and personality traits. These shape how the AI approaches legal strategy and drafting." },
          { title: "Attach Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills for <em>legal-research</em>, <em>brief-writing</em>, <em>case-analysis</em>, and more.' },
          { title: "Deploy to Your Case", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your legal team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The case cycle:</span>\nautolab_next  <span class="code-comment"># Partner strategizes</span>\nautolab_next  <span class="code-comment"># Associate researches</span>\nautolab_next  <span class="code-comment"># Partner reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When brief is ready:</span>\nautolab_editorial <span class="code-comment"># You review!</span>\n<span class="code-comment"># File / Minor Edits / Rewrite</span>'
      },
      medical: {
        badges: ["&#x1FA7A; Treatment Plans", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch clinical rounds unfold, step by step",
        characterYaml: {
          filename: "attending.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">Cardiology Attending</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">interventional cardiology</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">optimal evidence-based\n  patient outcomes</span>\n<span class="yaml-key">skills:</span>\n  <span class="yaml-val">- clinical-reasoning</span>\n  <span class="yaml-val">- evidence-based-medicine</span>\n  <span class="yaml-val">- treatment-planning</span>\n  <span class="yaml-val">- risk-stratification</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Thorough: systematic workup"</span>\n  <span class="yaml-val">- "Compassionate: patient-first"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, specialty, and personality traits. These shape how the AI approaches clinical reasoning and documentation." },
          { title: "Attach Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills for <em>clinical-reasoning</em>, <em>treatment-planning</em>, <em>risk-stratification</em>, and more.' },
          { title: "Deploy to Your Rounds", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your clinical team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The rounds cycle:</span>\nautolab_next  <span class="code-comment"># Attending assesses</span>\nautolab_next  <span class="code-comment"># Resident executes</span>\nautolab_next  <span class="code-comment"># Attending reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When plan is ready:</span>\nautolab_editorial <span class="code-comment"># You decide!</span>\n<span class="code-comment"># Approve / Modify / Reassess</span>'
      },
      creative: {
        badges: ["&#x1F3A8; Campaigns", "&#x1F9EC; Skill Containers", "&#x23F1;&#xFE0F; 24-Hour Sessions", "&#x2699;&#xFE0F; Fully Configurable"],
        tutorialSubtitle: "Watch a creative project unfold, step by step",
        characterYaml: {
          filename: "art-director.yaml",
          code: '<span class="yaml-key">title:</span> <span class="yaml-val">Brand Art Director</span>\n<span class="yaml-key">expertise:</span> <span class="yaml-val">visual identity, campaigns</span>\n<span class="yaml-key">goal:</span> <span class="yaml-val">create impactful brand\n  experiences</span>\n<span class="yaml-key">skills:</span>\n  <span class="yaml-val">- brand-strategy</span>\n  <span class="yaml-val">- visual-design</span>\n  <span class="yaml-val">- typography</span>\n  <span class="yaml-val">- motion-graphics</span>\n<span class="yaml-key">personality:</span>\n  <span class="yaml-val">- "Bold: pushes boundaries"</span>\n  <span class="yaml-val">- "Detail-oriented: pixel-perfect"</span>'
        },
        characterSteps: [
          { title: "Define a Persona", desc: "Give your character a title, expertise, and personality traits. These shape how the AI approaches creative direction and design." },
          { title: "Attach Skills", desc: 'Each skill maps to a <code>SKILL.md</code>&mdash;agent skills for <em>brand-strategy</em>, <em>visual-design</em>, <em>typography</em>, and more.' },
          { title: "Deploy to Your Studio", desc: 'Drop the YAML into <code>.autolab/profiles/</code>. The character joins your creative team with all skills active.' }
        ],
        gsStep4Code: '<span class="code-comment"># The studio cycle:</span>\nautolab_next  <span class="code-comment"># AD directs</span>\nautolab_next  <span class="code-comment"># Designer creates</span>\nautolab_next  <span class="code-comment"># AD reviews</span>\n  ...          <span class="code-comment"># iterate</span>\n\n<span class="code-comment"># When campaign is ready:</span>\nautolab_editorial <span class="code-comment"># You review!</span>\n<span class="code-comment"># Approve / Polish / Rework / Scrap</span>'
      }
    };
    Object.keys(extras).forEach(function (k) { Object.assign(DOMAINS[k], extras[k]); });
  })();

  // ============================================================
  // INITIALIZATION — Render all sprites on load
  // ============================================================
  document.addEventListener("DOMContentLoaded", () => {
    const S = window.LabSprites;
    if (!S) return console.error("LabSprites not loaded");

    // Hero sprites
    S.renderPI(document.getElementById("hero-pi-sprite"));
    S.renderTrainee(document.getElementById("hero-trainee-sprite"));

    // Tutorial sprites
    S.renderPI(document.getElementById("tut-pi-sprite"));
    S.renderTrainee(document.getElementById("tut-trainee-sprite"));
    S.renderExpert(document.getElementById("tut-expert-sprite"), "statistician");

    // Character explainer sprite
    S.renderExpert(document.getElementById("explainer-sprite"), "neuroscientist");

    // Init modules
    initNavigation();
    initDomainSwitcher();
    initTutorial();
    initMarketplace();
    initLinkRepoModal();
  });

  // ============================================================
  // SLOT-MACHINE SWAP — vertical spin animation helper
  // Uses inline transitions + setTimeout for maximum reliability
  // ============================================================
  function slotSwap(el, newContent, delay, isHTML) {
    if (!el) return;
    delay = delay || 0;
    var dur = 180; // ms per phase
    setTimeout(function () {
      // Phase 1: slide up + fade out
      el.style.transition = "opacity " + dur + "ms ease, transform " + dur + "ms ease";
      el.style.opacity = "0";
      el.style.transform = "translateY(-60%)";
      setTimeout(function () {
        // Phase 2: swap content, position below
        if (isHTML) el.innerHTML = newContent;
        else el.textContent = newContent;
        el.style.transition = "none";
        el.style.transform = "translateY(60%)";
        // Force reflow so the browser registers the new position
        void el.offsetWidth;
        // Phase 3: slide up into place + fade in
        el.style.transition = "opacity " + dur + "ms ease, transform " + dur + "ms ease";
        el.style.opacity = "1";
        el.style.transform = "translateY(0)";
        setTimeout(function () {
          el.style.transition = "";
          el.style.opacity = "";
          el.style.transform = "";
        }, dur + 50);
      }, dur + 20);
    }, delay);
  }

  // ============================================================
  // DOMAIN SWITCHER — Tag-machine style domain selector
  // ============================================================
  const _domainKeys = Object.keys(DOMAINS);
  let _autoLoopTimer = null;
  let _autoLoopIndex = 0;
  let _pauseResumeTimer = null;  // setTimeout ID for auto-resume after manual interaction
  let _tutorialActive = false;   // true while tutorial is in progress — blocks auto-loop

  function initDomainSwitcher() {
    const pills = document.querySelectorAll(".domain-pill");
    if (!pills.length) return;

    pills.forEach(pill => {
      pill.addEventListener("click", () => {
        _tutorialActive = false;   // clear tutorial lock so resume timer can be set
        _pauseAutoLoop(30000);     // pause 30s, then resume
        pills.forEach(p => p.classList.remove("active"));
        pill.classList.add("active");
        switchDomain(pill.dataset.domain);
      });
    });

    // Pause auto-loop on hero hover, resume on leave
    const hero = document.getElementById("hero");
    if (hero) {
      hero.addEventListener("mouseenter", function () { _stopAutoLoop(); });
      hero.addEventListener("mouseleave", function () {
        if (!_tutorialActive) setTimeout(_startAutoLoop, 2000);
      });
    }

    // Start auto-loop after initial delay
    setTimeout(_startAutoLoop, 2500);
  }

  /** Pause auto-loop and optionally auto-resume after `ms` milliseconds.
   *  Pass 0 or omit to pause without auto-resume. */
  function _pauseAutoLoop(ms) {
    _stopAutoLoop();
    if (_pauseResumeTimer) { clearTimeout(_pauseResumeTimer); _pauseResumeTimer = null; }
    if (ms > 0 && !_tutorialActive) {
      _pauseResumeTimer = setTimeout(_startAutoLoop, ms);
    }
  }

  function _startAutoLoop() {
    if (_autoLoopTimer || _tutorialActive) return;
    if (_pauseResumeTimer) { clearTimeout(_pauseResumeTimer); _pauseResumeTimer = null; }
    const pills = document.querySelectorAll(".domain-pill");
    _autoLoopTimer = setInterval(function () {
      _autoLoopIndex = (_autoLoopIndex + 1) % _domainKeys.length;
      const next = _domainKeys[_autoLoopIndex];
      pills.forEach(p => p.classList.toggle("active", p.dataset.domain === next));
      switchDomain(next);
    }, 4000);
  }

  function _stopAutoLoop() {
    if (_autoLoopTimer) { clearInterval(_autoLoopTimer); _autoLoopTimer = null; }
  }

  function switchDomain(domainId) {
    if (!DOMAINS[domainId]) return;
    activeDomain = domainId;
    const d = DOMAINS[domainId];
    const S = window.LabSprites;

    // --- Hero section (slot-machine animation for prominent elements) ---
    slotSwap(document.querySelector(".pi-label"), d.senior.label, 0);
    slotSwap(document.querySelector(".trainee-label"), d.junior.label, 50);
    slotSwap(document.querySelector(".editor-label"), "You (" + d.overseer.label + ")", 100);
    slotSwap(document.querySelector(".hero-subtitle"), d.heroSubtitle, 80, true);
    slotSwap(document.querySelector(".hero-description"), d.heroDescription, 140, true);

    // Hero you-icon (instant, not text)
    const youIcon = document.querySelector(".hero-you-icon");
    if (youIcon) youIcon.innerHTML = d.overseer.icon;

    // --- Lab frame (lighter fade for secondary elements) ---
    const piName = document.querySelector("#tut-pi-card .lf-char-name");
    const traineeName = document.querySelector("#tut-trainee-card .lf-char-name");
    if (piName) piName.textContent = d.senior.label;
    if (traineeName) traineeName.textContent = d.junior.label;

    const convTitle = document.querySelector(".lf-conv-title");
    if (convTitle) convTitle.textContent = d.meetingLog;

    const expertDiv = document.querySelector(".lf-experts-divider span");
    if (expertDiv) expertDiv.textContent = d.consultants;

    const expertName = document.querySelector(".lf-expert-name");
    const expertRole = document.querySelector(".lf-expert-role");
    if (expertName) expertName.textContent = d.expert.name;
    if (expertRole) expertRole.textContent = d.expert.role;
    const expertSprite = document.getElementById("tut-expert-sprite");
    if (expertSprite && S) S.renderExpert(expertSprite, d.expert.avatar);

    const invRows = document.querySelectorAll(".lf-inv-content .lf-inv-row");
    d.inventory.forEach((label, i) => {
      if (invRows[i]) {
        const labelSpan = invRows[i].querySelector("span:nth-child(2)");
        if (labelSpan) labelSpan.textContent = label;
      }
    });

    const paperTitle = document.querySelector(".lf-paper-section .lf-inv-title");
    if (paperTitle) paperTitle.textContent = d.artifact;

    const secRows = document.querySelectorAll(".lf-paper-sections .lf-paper-row");
    d.artifactSections.forEach((label, i) => {
      if (secRows[i]) {
        const labelSpan = secRows[i].querySelector("span:first-child");
        if (labelSpan) labelSpan.textContent = label;
      }
    });

    const roleEl = document.getElementById("tut-role");
    if (roleEl) roleEl.textContent = d.senior.short;

    // --- Editorial overlay ---
    const edTitle = document.querySelector(".tut-editor-title");
    if (edTitle) edTitle.textContent = d.overseer.label + "'s Desk";

    const decAccept = document.querySelector(".tut-dec-accept");
    const decMinor = document.querySelector("#tut-btn-minor");
    const decMajor = document.querySelector(".tut-dec-major");
    const decReject = document.querySelector(".tut-dec-reject");
    if (decAccept) decAccept.innerHTML = "&#x2705; " + d.decisions.accept;
    if (decMinor) decMinor.innerHTML = "&#x1F4DD; " + d.decisions.minor + " &#x2192;";
    if (decMajor) decMajor.innerHTML = "&#x1F6A7; " + d.decisions.major;
    if (decReject) decReject.innerHTML = "&#x274C; " + d.decisions.reject;

    const doneH4 = document.querySelector(".tut-editor-done-msg h4");
    const doneHint = document.querySelector(".tut-done-hint");
    if (doneH4) doneH4.textContent = "Decision sent: " + d.decisions.minor;
    if (doneHint) doneHint.innerHTML = "&#x2714; Tutorial complete! The " + d.senior.short + "&ndash;" + d.junior.short + " loop resumes with revisions.";

    const emptyHint = document.querySelector(".lf-hint");
    if (emptyHint) emptyHint.textContent = "Watch how " + d.senior.label + " and " + d.junior.label + " collaborate on a " + d.workspace.toLowerCase();

    // --- Hero badges ---
    if (d.badges) {
      const badgeContainer = document.querySelector(".hero-badges");
      if (badgeContainer) {
        badgeContainer.innerHTML = d.badges.map(function (b) { return '<span class="badge">' + b + '</span>'; }).join("");
      }
    }

    // --- Section subtitles ---
    if (d.tutorialSubtitle) {
      const tutSub = document.querySelector("#tutorial .section-subtitle");
      if (tutSub) tutSub.textContent = d.tutorialSubtitle;
    }

    // --- Characters section ---
    if (d.characterYaml) {
      const yamlHeader = document.querySelector(".yaml-header span:last-child");
      if (yamlHeader) yamlHeader.textContent = d.characterYaml.filename;
      const yamlCode = document.querySelector(".yaml-code code");
      if (yamlCode) yamlCode.innerHTML = d.characterYaml.code;
    }
    if (d.characterSteps) {
      const stepEls = document.querySelectorAll(".explainer-step");
      d.characterSteps.forEach(function (step, i) {
        if (stepEls[i]) {
          var inner = stepEls[i].querySelector("div");
          if (inner) inner.innerHTML = "<strong>" + step.title + "</strong><p>" + step.desc + "</p>";
        }
      });
    }

    // --- Get Started step 4 ---
    const gsLoop = document.getElementById("gs-loop-desc");
    if (gsLoop) gsLoop.innerHTML = "The " + d.senior.short + "&ndash;" + d.junior.short + " loop runs, then you act as " + d.overseer.label + ":";
    if (d.gsStep4Code) {
      const gsCode = document.querySelector("#gs-step4-code code");
      if (gsCode) gsCode.innerHTML = d.gsStep4Code;
    }

    // --- Marketplace filters ---
    const mpFilters = document.querySelectorAll(".mp-filter");
    mpFilters.forEach(function (btn) {
      var f = btn.dataset.filter;
      if (f === "pi") btn.textContent = d.senior.short;
      else if (f === "trainee") btn.textContent = d.junior.short;
    });

    // --- Reset tutorial to use new domain steps ---
    if (window._resetTutorial) window._resetTutorial();
  }

  // ============================================================
  // NAVIGATION — Scroll spy + smooth scrolling
  // ============================================================
  function initNavigation() {
    const links = document.querySelectorAll(".nav-link");
    const sections = ["hero", "tutorial", "characters", "marketplace", "get-started"];

    function updateActive() {
      const scrollY = window.scrollY + 100;
      let current = "hero";
      for (const id of sections) {
        const el = document.getElementById(id);
        if (el && el.offsetTop <= scrollY) current = id;
      }
      links.forEach(l => {
        l.classList.toggle("active", l.getAttribute("href") === "#" + current);
      });
    }

    window.addEventListener("scroll", updateActive, { passive: true });
    updateActive();
  }

  // ============================================================
  // TUTORIAL — Interactive meeting-log walkthrough (9 steps)
  // Includes: PI/Trainee loop, expert consultation,
  //           editorial submission, reviewer reports, decision
  // ============================================================
  function initTutorial() {
    let steps = DOMAINS[activeDomain].tutorialSteps;

    let currentStep = -1;
    const container = document.getElementById("tut-conversation");
    const emptyState = document.getElementById("tut-empty");
    const prevBtn = document.getElementById("tut-prev");
    const nextBtn = document.getElementById("tut-next");
    const restartBtn = document.getElementById("tut-restart");
    const stepLabel = document.getElementById("tut-step-label");
    const progressFill = document.getElementById("tut-progress-fill");

    const MAX_VISIBLE_BUBBLES = 3;

    function updateStep(dir) {
      const newStep = currentStep + dir;
      if (newStep < 0 || newStep >= steps.length) return;
      currentStep = newStep;
      const step = steps[currentStep];

      // Hide empty state
      if (emptyState) emptyState.style.display = "none";

      // When going backward, remove bubbles after current step
      if (dir < 0) {
        const existing = container.querySelectorAll(".tut-bubble");
        existing.forEach(b => {
          const idx = parseInt(b.id.replace("tut-step-", ""), 10);
          if (idx > currentStep) b.remove();
        });
      }

      // Add bubble if it doesn't already exist
      if (!document.getElementById("tut-step-" + currentStep)) {
        const dom = DOMAINS[activeDomain];
        // Resolve icon from domain config based on role
        const iconMap = { pi: dom.senior.icon, trainee: dom.junior.icon, editor: dom.overseer.icon, expert: "&#x1F4CA;", reviewer: "&#x1F50D;" };
        const icon = iconMap[step.role] || "&#x1F4AC;";
        const bubble = document.createElement("div");
        bubble.className = "tut-bubble " + step.role;
        bubble.id = "tut-step-" + currentStep;
        bubble.innerHTML = `
          <div class="tut-bubble-header">
            <span class="tut-bubble-role ${step.role}">${icon} ${step.label}</span>
            <span class="tut-bubble-meta">${step.meta}</span>
          </div>
          <div class="tut-bubble-text">${step.text}</div>
        `;
        container.appendChild(bubble);
      }

      // Keep only the last MAX_VISIBLE_BUBBLES visible — fade-remove oldest
      const allBubbles = container.querySelectorAll(".tut-bubble");
      if (allBubbles.length > MAX_VISIBLE_BUBBLES) {
        const toRemove = allBubbles.length - MAX_VISIBLE_BUBBLES;
        for (let i = 0; i < toRemove; i++) {
          allBubbles[i].style.transition = "opacity 0.3s, max-height 0.3s";
          allBubbles[i].style.opacity = "0";
          allBubbles[i].style.maxHeight = "0";
          allBubbles[i].style.overflow = "hidden";
          allBubbles[i].style.marginBottom = "0";
          allBubbles[i].style.padding = "0";
          const el = allBubbles[i];
          setTimeout(() => el.remove(), 350);
        }
      }

      // Scroll within the conversation container
      requestAnimationFrame(() => {
        container.scrollTop = container.scrollHeight;
      });

      // Update inventory
      const inv = step.inv;
      for (const [key, val] of Object.entries(inv)) {
        const el = document.getElementById("tut-inv-" + key);
        if (el && el.textContent !== String(val)) {
          el.textContent = val;
          el.classList.add("updated");
          setTimeout(() => el.classList.remove("updated"), 600);
        }
      }

      // Update paper progress
      const paperFill = document.getElementById("tut-paper-fill");
      const paperLabel = document.getElementById("tut-paper-label");
      if (paperFill) paperFill.style.width = step.paper + "%";
      if (paperLabel) paperLabel.textContent = step.paper + "%";

      // Update paper sections
      const secMap = { abstract: "tut-sec-abstract", intro: "tut-sec-intro", methods: "tut-sec-methods", results: "tut-sec-results", discussion: "tut-sec-discussion" };
      for (const [key, id] of Object.entries(secMap)) {
        const el = document.getElementById(id);
        if (el) {
          const val = step.sections[key] || "--";
          el.textContent = val;
          el.classList.toggle("has-content", val !== "--");
        }
      }

      // Update character cards
      const piCard = document.getElementById("tut-pi-card");
      const traineeCard = document.getElementById("tut-trainee-card");
      piCard.className = "lf-char-card" + (step.piActive ? " active-pi" : "");
      traineeCard.className = "lf-char-card" + (step.traineeActive ? " active-trainee" : "");

      // Update status dots
      const piDot = document.getElementById("tut-pi-dot");
      const traineeDot = document.getElementById("tut-trainee-dot");
      piDot.className = "lf-status-dot " + (step.piActive ? "speaking" : "idle");
      traineeDot.className = "lf-status-dot " + (step.traineeActive ? "speaking" : "idle");
      document.getElementById("tut-pi-status").textContent = step.piStatus;
      document.getElementById("tut-trainee-status").textContent = step.traineeStatus;

      // Update top bar
      document.getElementById("tut-iteration").textContent = step.iteration;
      document.getElementById("tut-role").textContent = step.nextRole;

      // Show/hide editorial overlay on editor steps
      const editorOverlay = document.getElementById("tut-editor-overlay");
      const edPhases = {
        submitted: document.getElementById("tut-editor-submitted"),
        selectReviewers: document.getElementById("tut-editor-select-reviewers"),
        decision: document.getElementById("tut-editor-decision"),
        done: document.getElementById("tut-editor-done"),
      };
      if (editorOverlay) {
        if (step.role === "editor") {
          editorOverlay.classList.remove("hidden");
          const isDecision = step.meta === "Decision";
          // Show only the correct phase
          Object.values(edPhases).forEach(el => { if (el) el.style.display = "none"; });
          if (isDecision && edPhases.decision) edPhases.decision.style.display = "";
          else if (edPhases.submitted) edPhases.submitted.style.display = "";
          const roundEl = document.getElementById("tut-editor-round");
          if (roundEl) roundEl.textContent = isDecision ? "Round 1 — Decision" : "Round 1";
          // Render reviewer sprites on first show
          document.querySelectorAll(".tut-reviewer-sprite").forEach(c => {
            const av = c.getAttribute("data-avatar");
            if (av && window.LabSprites) window.LabSprites.renderExpert(c, av);
          });
        } else {
          editorOverlay.classList.add("hidden");
          // Reset phases for next visit
          Object.values(edPhases).forEach(el => { if (el) el.style.display = "none"; });
          if (edPhases.submitted) edPhases.submitted.style.display = "";
        }
      }

      // Update controls
      prevBtn.disabled = currentStep <= 0;
      nextBtn.textContent = currentStep >= steps.length - 1 ? "Done \u2714" : "Next Step \u2192";
      nextBtn.disabled = currentStep >= steps.length - 1;
      stepLabel.textContent = `Step ${currentStep + 1} / ${steps.length}`;
      progressFill.style.width = ((currentStep + 1) / steps.length * 100) + "%";
      // Show restart button once tutorial has started
      if (restartBtn) restartBtn.classList.toggle("hidden", currentStep < 0);
    }

    function resetTutorial() {
      // Allow auto-loop to resume when tutorial is reset
      _tutorialActive = false;
      // Re-read steps from current domain (supports domain switching)
      steps = DOMAINS[activeDomain].tutorialSteps;
      currentStep = -1;
      // Clear all bubbles
      container.querySelectorAll(".tut-bubble").forEach(b => b.remove());
      // Show empty state
      if (emptyState) emptyState.style.display = "";
      // Hide editorial overlay
      const overlay = document.getElementById("tut-editor-overlay");
      if (overlay) overlay.classList.add("hidden");
      // Reset controls
      prevBtn.disabled = true;
      nextBtn.textContent = "Start Tutorial \u2192";
      nextBtn.disabled = false;
      stepLabel.textContent = `Step 0 / ${steps.length}`;
      progressFill.style.width = "0%";
      if (restartBtn) restartBtn.classList.add("hidden");
      // Reset editorial phases
      const edSub = document.getElementById("tut-editor-submitted");
      const edSel = document.getElementById("tut-editor-select-reviewers");
      const edDec = document.getElementById("tut-editor-decision");
      const edDone = document.getElementById("tut-editor-done");
      [edSub, edSel, edDec, edDone].forEach(el => { if (el) el.style.display = "none"; });
      if (edSub) edSub.style.display = "";
      // Reset inventory and paper progress
      ["data", "scripts", "figures", "results", "paper"].forEach(k => {
        const el = document.getElementById("tut-inv-" + k);
        if (el) el.textContent = "0";
      });
      const paperFill = document.getElementById("tut-paper-fill");
      const paperLabel = document.getElementById("tut-paper-label");
      if (paperFill) paperFill.style.width = "0%";
      if (paperLabel) paperLabel.textContent = "0%";
      const secIds = ["tut-sec-abstract", "tut-sec-intro", "tut-sec-methods", "tut-sec-results", "tut-sec-discussion"];
      secIds.forEach(id => { const el = document.getElementById(id); if (el) { el.textContent = "--"; el.classList.remove("has-content"); } });
      // Reset character card states
      const piCard = document.getElementById("tut-pi-card");
      const traineeCard = document.getElementById("tut-trainee-card");
      if (piCard) piCard.className = "lf-char-card";
      if (traineeCard) traineeCard.className = "lf-char-card";
      // Reset iteration display
      const iterEl = document.getElementById("tut-iteration");
      if (iterEl) iterEl.textContent = "1";
      const roleEl = document.getElementById("tut-role");
      if (roleEl) roleEl.textContent = DOMAINS[activeDomain].senior.short;
    }

    // Expose resetTutorial for the domain switcher
    window._resetTutorial = resetTutorial;

    nextBtn.addEventListener("click", () => {
      if (currentStep === -1) {
        nextBtn.textContent = "Next Step \u2192";
        // Stop auto-rotation while tutorial is active
        _tutorialActive = true;
        _stopAutoLoop();
      }
      updateStep(1);
    });
    prevBtn.addEventListener("click", () => updateStep(-1));
    if (restartBtn) restartBtn.addEventListener("click", resetTutorial);

    // Editorial overlay — phase transitions
    const edPanels = {
      submitted: document.getElementById("tut-editor-submitted"),
      select: document.getElementById("tut-editor-select-reviewers"),
      decision: document.getElementById("tut-editor-decision"),
      done: document.getElementById("tut-editor-done"),
    };

    function showEdPhase(phase) {
      Object.values(edPanels).forEach(el => { if (el) el.style.display = "none"; });
      if (edPanels[phase]) edPanels[phase].style.display = "";
    }

    // "Send to Reviewers" → show reviewer selection panel
    const sendBtn = document.getElementById("tut-btn-send-reviewers");
    if (sendBtn) sendBtn.addEventListener("click", () => {
      showEdPhase("select");
      // Render reviewer sprites
      document.querySelectorAll(".tut-reviewer-sprite").forEach(c => {
        const av = c.getAttribute("data-avatar");
        if (av && window.LabSprites) window.LabSprites.renderExpert(c, av);
      });
    });

    // "Back" from reviewer selection → show submission
    const backBtn = document.getElementById("tut-btn-back-submit");
    if (backBtn) backBtn.addEventListener("click", () => showEdPhase("submitted"));

    // "Confirm & Send" → advance to next tutorial step (reviewer)
    const confirmBtn = document.getElementById("tut-btn-confirm-reviewers");
    if (confirmBtn) confirmBtn.addEventListener("click", () => updateStep(1));

    // "Minor Revision" → show done state (this is the last step)
    const minorBtn = document.getElementById("tut-btn-minor");
    if (minorBtn) minorBtn.addEventListener("click", () => {
      showEdPhase("done");
      // Mark tutorial as done in the controls
      nextBtn.textContent = "Done \u2714";
      nextBtn.disabled = true;
      stepLabel.textContent = `Step ${steps.length} / ${steps.length}`;
      progressFill.style.width = "100%";
    });
  }

  // ============================================================
  // MARKETPLACE — GitHub-repo-based character registry
  //
  // Characters link to their source GitHub repos.
  // Users fork the template, create their character,
  // and list their repo. Ranked by GitHub stars.
  // Official characters from albert-ying get a badge.
  // ============================================================
  const CHARACTERS = [
    {
      id: "maria-chen",
      name: "Dr. Maria Chen",
      role: "pi",
      avatar: "neuroscientist",
      title: "Computational Biology PI",
      expertise: "Single-cell genomics, machine learning, and multi-modal data integration",
      goal: "Discover cell-type-specific gene regulatory programs using multi-modal single-cell data",
      skills: ["scanpy", "scvi-tools", "pytorch-lightning", "scientific-writing", "scientific-visualization", "statistical-analysis"],
      personality: [
        "Visionary: identifies novel biological questions from data patterns",
        "Rigorous: demands reproducible computational pipelines with version control",
        "Collaborative: bridges wet lab and dry lab teams effectively"
      ],
      github: "albert-ying/autolab-char-compbio-pi",
      stars: 342,
      official: true
    },
    {
      id: "alex-kumar",
      name: "Alex Kumar",
      role: "trainee",
      avatar: "bioinformatician",
      title: "Bioinformatics Postdoc",
      expertise: "NGS data analysis, pipeline development, and statistical genomics",
      goal: "Build clean, reproducible analysis pipelines and generate publication-quality figures",
      skills: ["scanpy", "pydeseq2", "pysam", "matplotlib", "seaborn", "scikit-learn", "deeptools"],
      personality: [
        "Dedicated: completes tasks thoroughly with comprehensive documentation",
        "Technical: writes self-contained, reproducible code with proper testing",
        "Proactive: identifies additional analyses that strengthen the narrative"
      ],
      github: "albert-ying/autolab-char-bioinfo-trainee",
      stars: 256,
      official: true
    },
    {
      id: "sarah-oconnor",
      name: "Dr. Sarah O'Connor",
      role: "pi",
      avatar: "chemist",
      title: "Medicinal Chemistry PI",
      expertise: "Drug discovery, QSAR modeling, and lead optimization",
      goal: "Identify and optimize novel small molecule inhibitors through computational screening",
      skills: ["rdkit", "datamol", "deepchem", "pytdc", "medchem", "pubchem-database"],
      personality: [
        "Strategic: prioritizes compounds with drug-like properties early",
        "Data-driven: demands SAR analysis before advancing any lead series",
        "Publication-savvy: structures work for high-impact medicinal chemistry journals"
      ],
      github: "albert-ying/autolab-char-medchem-pi",
      stars: 189,
      official: false
    },
    {
      id: "james-park",
      name: "James Park",
      role: "trainee",
      avatar: "ml_engineer",
      title: "ML Research Engineer",
      expertise: "Deep learning, model training, and inference optimization",
      goal: "Implement and benchmark state-of-the-art models with clean, efficient code",
      skills: ["pytorch-lightning", "transformers", "accelerate", "weights-and-biases", "vllm", "flash-attention"],
      personality: [
        "Efficient: writes highly optimized code with proper GPU utilization",
        "Systematic: benchmarks every change with rigorous ablation studies",
        "Clear communicator: documents architecture decisions and trade-offs"
      ],
      github: "albert-ying/autolab-char-ml-engineer",
      stars: 198,
      official: false
    },
    {
      id: "elena-vasquez",
      name: "Dr. Elena Vasquez",
      role: "collaborator",
      avatar: "epidemiologist",
      title: "Clinical Epidemiologist",
      expertise: "Clinical trial design, survival analysis, and real-world evidence",
      goal: "Ensure robust clinical study designs and proper statistical interpretation",
      skills: ["scikit-survival", "statistical-analysis", "statsmodels", "clinical-reports", "clinicaltrials-database"],
      personality: [
        "Methodical: insists on pre-registered analysis plans",
        "Critical: identifies confounders and biases in study designs",
        "Translational: connects statistical findings to clinical implications"
      ],
      github: "albert-ying/autolab-char-clinical-epi",
      stars: 145,
      official: false
    },
    {
      id: "wei-zhang",
      name: "Dr. Wei Zhang",
      role: "collaborator",
      avatar: "statistician",
      title: "Biostatistician",
      expertise: "Bayesian modeling, causal inference, and high-dimensional statistics",
      goal: "Provide rigorous statistical frameworks and validate analytical approaches",
      skills: ["pymc", "statistical-analysis", "statsmodels", "scikit-learn", "shap", "scientific-visualization"],
      personality: [
        "Precise: never allows hand-waving about statistical assumptions",
        "Educational: explains complex methods in accessible terms",
        "Conservative: prefers well-validated methods over trendy approaches"
      ],
      github: "albert-ying/autolab-char-biostatistician",
      stars: 267,
      official: true
    },
    {
      id: "cleon-garzol",
      name: "Cleon Garzol",
      role: "pi",
      avatar: "engineer",
      title: "Founder & Interdisciplinary Research Director",
      expertise: "Biomedical engineering, systems biology, longevity science, geroscience, startup strategy, translational research, signal detection in noisy datasets",
      goal: "Drive high-impact research at the intersection of engineering and biology, with a focus on longevity interventions. Form rapid hypotheses from limited data, stress-test them rigorously, and pivot when evidence demands it.",
      skills: ["scientific-writing", "hypothesis-generation", "scientific-brainstorming", "literature-review", "statistical-analysis", "exploratory-data-analysis", "peer-review", "scientific-visualization"],
      personality: [
        "Decisive: forms strong initial opinions from sparse data — treats early signals as actionable hypotheses",
        "Intellectually Honest: genuinely willing to abandon a position when confronted with compelling counter-evidence",
        "Signal Hunter: exceptional at extracting meaningful patterns from noisy, ambiguous, or contradictory information",
        "Builder's Mindset: thinks in terms of systems, leverage, and scalability — shaped by founding companies",
        "Cross-Pollinator: draws unexpected connections between engineering principles, biological mechanisms, and business strategy",
        "Impatient with Theater: low tolerance for academic posturing or unfalsifiable claims",
        "Contrarian Instinct: gravitates toward overlooked or unpopular hypotheses when consensus feels under-examined"
      ],
      github: "mayi12345/autolab-char-cleon-garzol",
      stars: 0,
      official: false
    },
    {
      id: "michael-florea",
      name: "Michael Florea",
      role: "pi",
      avatar: "geneticist",
      title: "CEO & Scientist — Gene Delivery & Animal Automation",
      expertise: "Whole-body gene delivery systems, animal automation for drug discovery, systems biology, synthetic biology, AI safety through biological approaches, multilingual (Estonian, English, German, Chinese, Japanese)",
      goal: "Eliminate the animal experimentation bottleneck in drug discovery through automation, while pursuing biological strategies that could outpace AI. Build world-class teams of high-agency individuals who need zero hand-holding.",
      skills: ["scientific-writing", "hypothesis-generation", "scientific-brainstorming", "literature-review", "statistical-analysis", "exploratory-data-analysis", "peer-review", "scientific-visualization"],
      personality: [
        "A-Player Filter: only works with high-agency people who take ownership — zero patience for passengers",
        "Founder's Urgency: treats every project like a startup with runway burning — moves fast, ships fast",
        "Deep Scientist: PhD-level rigor in biology — won't cut corners on experimental design",
        "AI Safety Realist: deeply concerned about AI alignment — pursues biological computation to outcompete silicon AI",
        "Systems Thinker: sees gene delivery, animal models, and drug pipelines as interconnected systems",
        "Caring but Direct: genuinely warm and loyal — delivers hard truths without sugarcoating",
        "Polyglot Perspective: draws on Estonian, German, Chinese, and Japanese cultural frameworks for unique angles",
        "Relentless Optimizer: gym rat discipline applied to everything — consistent effort compounds"
      ],
      github: "mayi12345/autolab-char-michael-florea",
      stars: 0,
      official: false
    },
  ];

  function initMarketplace() {
    const grid = document.getElementById("mp-grid");
    const filters = document.querySelectorAll(".mp-filter");
    const searchInput = document.getElementById("mp-search-input");
    const sortSelect = document.getElementById("mp-sort");

    function renderCards(chars) {
      grid.innerHTML = "";
      if (chars.length === 0) {
        grid.innerHTML = '<div style="grid-column:1/-1;text-align:center;padding:40px;color:var(--text-muted);">No characters match your search.</div>';
        return;
      }
      chars.forEach(ch => {
        const card = document.createElement("div");
        card.className = "mp-card";
        card.dataset.role = ch.role;
        const officialBadge = ch.official ? '<span class="mp-official-badge" title="Official character by albert-ying">&#x2705; Official</span>' : '';
        card.innerHTML = `
          <div class="mp-card-top">
            <div class="mp-card-avatar"><canvas width="40" height="64" data-expert="${ch.avatar}"></canvas></div>
            <div class="mp-card-info">
              <div class="mp-card-name-row">
                <span class="mp-card-name">${ch.name}</span>
                ${officialBadge}
              </div>
              <div class="mp-card-title">${ch.title}</div>
              <span class="mp-card-role role-${ch.role}">${ch.role.toUpperCase()}</span>
            </div>
          </div>
          <div class="mp-card-skills">
            ${ch.skills.slice(0, 5).map(s => `<span class="mp-card-skill">${s}</span>`).join("")}
            ${ch.skills.length > 5 ? `<span class="mp-card-skill">+${ch.skills.length - 5}</span>` : ""}
          </div>
          <div class="mp-card-footer">
            <a href="https://github.com/${ch.github}" class="mp-card-repo" target="_blank" rel="noopener" title="${ch.github}">
              <svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
            </a>
            <span class="mp-card-stars" title="${ch.stars} GitHub stars">&#x2B50; ${ch.stars}</span>
            <button class="mp-card-download" data-id="${ch.id}">&#x1F4E5; Download</button>
          </div>
        `;

        // Render sprite
        const canvas = card.querySelector("canvas");
        window.LabSprites.renderExpert(canvas, ch.avatar);

        // Click for detail
        card.addEventListener("click", (e) => {
          if (e.target.closest(".mp-card-download") || e.target.closest(".mp-card-repo")) return;
          showCharacterDetail(ch);
        });

        // Download button
        card.querySelector(".mp-card-download").addEventListener("click", (e) => {
          e.stopPropagation();
          downloadCharacterYAML(ch);
        });

        grid.appendChild(card);
      });
    }

    function filterSortAndRender() {
      const activeFilter = document.querySelector(".mp-filter.active").dataset.filter;
      const query = searchInput.value.toLowerCase().trim();
      const sortBy = sortSelect ? sortSelect.value : "stars";

      let filtered = CHARACTERS;
      if (activeFilter !== "all") {
        filtered = filtered.filter(c => c.role === activeFilter);
      }
      if (query) {
        filtered = filtered.filter(c =>
          c.name.toLowerCase().includes(query) ||
          c.title.toLowerCase().includes(query) ||
          c.expertise.toLowerCase().includes(query) ||
          c.skills.some(s => s.toLowerCase().includes(query))
        );
      }

      // Sort
      filtered = [...filtered];
      if (sortBy === "stars") {
        filtered.sort((a, b) => b.stars - a.stars);
      } else if (sortBy === "name") {
        filtered.sort((a, b) => a.name.localeCompare(b.name));
      }

      renderCards(filtered);
    }

    filters.forEach(f => f.addEventListener("click", () => {
      filters.forEach(ff => ff.classList.remove("active"));
      f.classList.add("active");
      filterSortAndRender();
    }));

    searchInput.addEventListener("input", filterSortAndRender);
    if (sortSelect) sortSelect.addEventListener("change", filterSortAndRender);

    // Initial render
    filterSortAndRender();
  }

  // ============================================================
  // CHARACTER DETAIL MODAL
  // ============================================================
  function showCharacterDetail(ch) {
    const modal = document.getElementById("char-detail-modal");
    const body = document.getElementById("detail-body");
    const title = document.getElementById("detail-title");
    const icon = document.getElementById("detail-icon");

    title.textContent = ch.name;
    icon.innerHTML = ch.role === "pi" ? "&#x1F52C;" : ch.role === "trainee" ? "&#x1F9EA;" : "&#x1F465;";

    const skillCategories = {
      analysis: ["scanpy", "scvi-tools", "pydeseq2", "pysam", "deeptools", "rdkit", "datamol", "deepchem", "medchem"],
      ml: ["pytorch-lightning", "transformers", "accelerate", "scikit-learn", "flash-attention", "vllm", "shap"],
      writing: ["scientific-writing", "clinical-reports"],
      viz: ["scientific-visualization", "matplotlib", "seaborn", "plotly"],
      stats: ["statistical-analysis", "statsmodels", "pymc", "scikit-survival"],
      infra: ["weights-and-biases", "pytdc", "pubchem-database", "clinicaltrials-database"],
    };

    function skillClass(skill) {
      for (const [cat, skills] of Object.entries(skillCategories)) {
        if (skills.includes(skill)) return "skill-" + cat;
      }
      return "skill-infra";
    }

    const officialHTML = ch.official ? '<span class="mp-official-badge" style="font-size:10px;">&#x2705; Official</span>' : '';

    body.innerHTML = `
      <div class="detail-grid">
        <div class="detail-avatar">
          <canvas width="40" height="64" id="detail-sprite"></canvas>
          <span class="mp-card-role role-${ch.role}" style="font-size:10px;">${ch.role.toUpperCase()}</span>
          ${officialHTML}
          <div class="detail-stars">&#x2B50; ${ch.stars} stars</div>
        </div>
        <div class="detail-content">
          <h3>&#x1F3AF; EXPERTISE</h3>
          <p>${ch.expertise}</p>

          <h3>&#x1F4CB; GOAL</h3>
          <p>${ch.goal}</p>

          <h3>&#x1F9EC; SKILLS (Cursor Skills)</h3>
          <p class="detail-skill-note">Each skill maps to a <code>SKILL.md</code> file that the AI agent reads when acting as this character.</p>
          <div class="char-explainer-skills" style="justify-content:flex-start;gap:5px;">
            ${ch.skills.map(s => `<span class="skill-tag ${skillClass(s)}">${s}</span>`).join("")}
          </div>

          <h3>&#x1F3AD; PERSONALITY</h3>
          <ul>${ch.personality.map(p => `<li>${p}</li>`).join("")}</ul>

          <h3>&#x1F517; SOURCE REPO</h3>
          <p><a href="https://github.com/${ch.github}" target="_blank" rel="noopener">github.com/${ch.github}</a></p>

          <div class="detail-actions">
            <button class="btn-primary" onclick="window._downloadYAML('${ch.id}')">&#x1F4E5; Download YAML</button>
            <button class="btn-secondary" onclick="window._copyYAML('${ch.id}')">&#x1F4CB; Copy YAML</button>
            <a href="https://github.com/${ch.github}" target="_blank" rel="noopener" class="btn-github">
              <svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
              View Repo
            </a>
          </div>
        </div>
      </div>
    `;

    // Render sprite
    const canvas = document.getElementById("detail-sprite");
    canvas.style.imageRendering = "pixelated";
    canvas.style.width = "80px";
    canvas.style.height = "128px";
    window.LabSprites.renderExpert(canvas, ch.avatar);

    modal.classList.remove("hidden");

    document.getElementById("detail-modal-close").onclick = () => modal.classList.add("hidden");
    modal.addEventListener("click", (e) => {
      if (e.target === modal) modal.classList.add("hidden");
    });
  }

  // ============================================================
  // YAML GENERATION & DOWNLOAD
  // ============================================================
  function characterToYAML(ch) {
    let yaml = `# ${ch.name || "Custom Character"} — ${ch.title || "Role"}\n`;
    yaml += `# Role: ${ch.role || "pi"}\n`;
    if (ch.github) yaml += `# Source: https://github.com/${ch.github}\n`;
    yaml += `\n`;
    yaml += `title: ${ch.title}\n`;
    yaml += `expertise: ${ch.expertise}\n`;
    yaml += `goal: ${ch.goal}\n`;
    yaml += `skills:\n`;
    (ch.skills || []).forEach(s => { yaml += `  - ${s}\n`; });
    yaml += `personality:\n`;
    (ch.personality || []).forEach(p => { yaml += `  - "${p}"\n`; });
    return yaml;
  }

  function downloadCharacterYAML(ch) {
    const yaml = characterToYAML(ch);
    const blob = new Blob([yaml], { type: "text/yaml" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = (ch.role || "character") + ".yaml";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  // Global helpers for detail modal buttons
  window._downloadYAML = (id) => {
    const ch = CHARACTERS.find(c => c.id === id);
    if (ch) downloadCharacterYAML(ch);
  };
  window._copyYAML = (id) => {
    const ch = CHARACTERS.find(c => c.id === id);
    if (ch) {
      navigator.clipboard.writeText(characterToYAML(ch)).then(() => {
        alert("YAML copied to clipboard!");
      });
    }
  };

  // ============================================================
  // LINK REPO MODAL (replaces old upload modal)
  // ============================================================
  function initLinkRepoModal() {
    const modal = document.getElementById("link-repo-modal");
    if (!modal) return;
    const openBtn = document.getElementById("btn-list-char");
    const closeBtn = document.getElementById("link-repo-modal-close");

    openBtn.addEventListener("click", () => modal.classList.remove("hidden"));
    closeBtn.addEventListener("click", () => modal.classList.add("hidden"));
    modal.addEventListener("click", (e) => {
      if (e.target === modal) modal.classList.add("hidden");
    });
  }

  function escapeHTML(str) {
    return str.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

})();
