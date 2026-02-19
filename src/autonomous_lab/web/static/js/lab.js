/* ==============================================================
   AUTONOMOUS LAB — Game UI Frontend
   Pixel art characters, live polling, inventory, conversation,
   file preview modal, character thought bubbles, markdown,
   expert consultants, biomedical toolkit integration
   ============================================================== */

(function () {
  "use strict";

  // ============================================================
  // PIXEL ART SPRITE DATA (10w x 16h)
  // Color indices: 0=transparent, 1=skin, 2=hair, 3=eyes,
  // 4=coat/clothing, 5=shirt, 6=pants, 7=shoes, 8=mouth, 9=accessory
  // ============================================================

  const PI_PALETTE = {
    0:null, 1:"#f4c7a3", 2:"#8a8a8a", 3:"#333",
    4:"#f0f0f0", 5:"#4a6b9c", 6:"#3a3a5a", 7:"#6b5040",
    8:"#d46b6b", 9:"#5588bb",
  };
  const TRAINEE_PALETTE = {
    0:null, 1:"#f4c7a3", 2:"#7a5c3a", 3:"#333",
    4:"#4a78bc", 5:"#4a78bc", 6:"#4a5a8a", 7:"#e0d8c8",
    8:"#d46b6b", 9:"#4a78bc",
  };

  const FRAME_GLASSES = [
    [0,0,0,2,2,2,2,0,0,0],[0,0,2,2,2,2,2,2,0,0],
    [0,0,2,1,1,1,1,2,0,0],[0,0,1,9,1,1,9,1,0,0],
    [0,0,1,1,1,1,1,1,0,0],[0,0,0,1,8,8,1,0,0,0],
    [0,0,0,0,1,1,0,0,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,4,4,4,5,5,4,4,4,0],[0,4,4,4,5,5,4,4,4,0],
    [0,0,4,4,4,4,4,4,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,0,0,6,6,6,6,0,0,0],[0,0,0,6,0,0,6,0,0,0],
    [0,0,0,6,0,0,6,0,0,0],[0,0,7,7,0,0,7,7,0,0],
  ];
  const FRAME_NO_GLASSES = [
    [0,0,0,2,2,2,2,0,0,0],[0,0,2,2,2,2,2,2,0,0],
    [0,0,2,1,1,1,1,2,0,0],[0,0,1,3,1,1,3,1,0,0],
    [0,0,1,1,1,1,1,1,0,0],[0,0,0,1,8,8,1,0,0,0],
    [0,0,0,0,1,1,0,0,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,4,4,4,4,4,4,4,4,0],[0,4,4,4,4,4,4,4,4,0],
    [0,0,4,4,4,4,4,4,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,0,0,6,6,6,6,0,0,0],[0,0,0,6,0,0,6,0,0,0],
    [0,0,0,6,0,0,6,0,0,0],[0,0,7,7,0,0,7,7,0,0],
  ];

  const PI_FRAME_1 = FRAME_GLASSES;
  const PI_FRAME_2 = FRAME_GLASSES.map((r,y)=>y===3?[0,0,1,1,1,1,1,1,0,0]:r);
  const TRAINEE_FRAME_1 = FRAME_NO_GLASSES;
  const TRAINEE_FRAME_2 = FRAME_NO_GLASSES.map((r,y)=>y===3?[0,0,1,1,1,1,1,1,0,0]:r);

  // ============================================================
  // EXPERT CHARACTER PALETTES — 25+ interdisciplinary roles
  // Each: label, palette {0-9}, glasses bool
  // ============================================================
  const EXPERT_DEFS = {
    // --- Review & Meta ---
    reviewer:          {label:"Critical Reviewer",   glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#5a2020",3:"#333",4:"#2a2a3a",5:"#e0e0e0",6:"#2a2a3a",7:"#3a2a20",8:"#c06060",9:"#881818"}},
    bioethicist:       {label:"Bioethicist",          glasses:true,  palette:{0:null,1:"#d4a070",2:"#4a3020",3:"#333",4:"#404060",5:"#8080a0",6:"#3a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#606080"}},
    science_writer:    {label:"Science Writer",       glasses:false, palette:{0:null,1:"#f4c7a3",2:"#c06030",3:"#333",4:"#3a5a6a",5:"#5a8a9a",6:"#4a5a6a",7:"#6a5a4a",8:"#d46b6b",9:null}},
    grant_reviewer:    {label:"Grant Reviewer",       glasses:true,  palette:{0:null,1:"#e8b888",2:"#6a6a7a",3:"#333",4:"#3a3a4a",5:"#b0b0c0",6:"#3a3a5a",7:"#5a4a3a",8:"#c06060",9:"#5a5a7a"}},
    // --- Life Sciences ---
    immunologist:      {label:"Immunologist",         glasses:false, palette:{0:null,1:"#d4a878",2:"#1a1a30",3:"#333",4:"#2a8a4a",5:"#50b060",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    oncologist:        {label:"Oncologist",           glasses:false, palette:{0:null,1:"#f4c7a3",2:"#3a2a1a",3:"#333",4:"#7a2a4a",5:"#b04070",6:"#4a3a5a",7:"#6a5040",8:"#d46b6b",9:null}},
    neuroscientist:    {label:"Neuroscientist",       glasses:true,  palette:{0:null,1:"#e8b888",2:"#2a2040",3:"#333",4:"#5050a0",5:"#7070c0",6:"#3a3a5a",7:"#e0d8c8",8:"#d46b6b",9:"#4040a0"}},
    geneticist:        {label:"Geneticist",           glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#4a3020",3:"#333",4:"#2a6a8a",5:"#4a90b0",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:"#2a6080"}},
    cell_biologist:    {label:"Cell Biologist",       glasses:false, palette:{0:null,1:"#d4a070",2:"#6a4a2a",3:"#333",4:"#3a7a5a",5:"#5aa07a",6:"#3a5a4a",7:"#6a5a4a",8:"#d46b6b",9:null}},
    microbiologist:    {label:"Microbiologist",       glasses:false, palette:{0:null,1:"#f4c7a3",2:"#8a5a2a",3:"#333",4:"#6a8a3a",5:"#90b050",6:"#4a5a3a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    pathologist:       {label:"Pathologist",          glasses:true,  palette:{0:null,1:"#e8b888",2:"#3a2020",3:"#333",4:"#f0f0f0",5:"#6a3050",6:"#3a3a5a",7:"#6b5040",8:"#c06060",9:"#5a3040"}},
    pharmacologist:    {label:"Pharmacologist",       glasses:false, palette:{0:null,1:"#f4c7a3",2:"#5a4030",3:"#333",4:"#5a3a7a",5:"#8060a0",6:"#4a3a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    structural_bio:    {label:"Structural Biologist", glasses:true,  palette:{0:null,1:"#d4a878",2:"#2a2a2a",3:"#333",4:"#4a6a8a",5:"#6a90b0",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:"#3a6080"}},
    systems_biologist: {label:"Systems Biologist",    glasses:false, palette:{0:null,1:"#f4c7a3",2:"#7a4a20",3:"#333",4:"#3a5a3a",5:"#5a8a5a",6:"#4a5a4a",7:"#6a5040",8:"#d46b6b",9:null}},
    epidemiologist:    {label:"Epidemiologist",       glasses:true,  palette:{0:null,1:"#e8b888",2:"#5a3a2a",3:"#333",4:"#4a5a7a",5:"#6a80a0",6:"#3a4a6a",7:"#5a4a3a",8:"#d46b6b",9:"#3a5080"}},
    // --- Computational ---
    statistician:      {label:"Statistician",         glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#c8a848",3:"#333",4:"#8a7050",5:"#e0e0e0",6:"#4a4a5a",7:"#6a5040",8:"#d46b6b",9:"#886830"}},
    bioinformatician:  {label:"Bioinformatician",     glasses:false, palette:{0:null,1:"#e8b888",2:"#2a2a40",3:"#333",4:"#5060a0",5:"#6070b0",6:"#3a3a5a",7:"#e0d8c8",8:"#d46b6b",9:"#70a0c0"}},
    data_scientist:    {label:"Data Scientist",       glasses:false, palette:{0:null,1:"#f4c7a3",2:"#3a3a50",3:"#333",4:"#3a6a6a",5:"#5a9a9a",6:"#3a5a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    ml_engineer:       {label:"ML Engineer",          glasses:false, palette:{0:null,1:"#d4a070",2:"#2a2a30",3:"#333",4:"#5a5a6a",5:"#8a8a9a",6:"#3a3a4a",7:"#e0d8c8",8:"#d46b6b",9:null}},
    comp_biologist:    {label:"Comp. Biologist",      glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#6a4a2a",3:"#333",4:"#4a4a7a",5:"#6a6aa0",6:"#3a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#4a4a80"}},
    // --- Clinical ---
    clinician:         {label:"Clinician",            glasses:false, palette:{0:null,1:"#d4a070",2:"#3a2020",3:"#333",4:"#f0f0f0",5:"#4090b0",6:"#3a3a5a",7:"#6b5040",8:"#d46b6b",9:null}},
    radiologist:       {label:"Radiologist",          glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#4a4a4a",3:"#333",4:"#2a3a4a",5:"#4a6a8a",6:"#2a3a4a",7:"#5a4a3a",8:"#d46b6b",9:"#3a5a7a"}},
    surgeon:           {label:"Surgeon",              glasses:false, palette:{0:null,1:"#e8b888",2:"#2a2a3a",3:"#333",4:"#2a8a7a",5:"#50b0a0",6:"#2a5a5a",7:"#f0f0f0",8:"#d46b6b",9:null}},
    // --- Interdisciplinary ---
    chemist:           {label:"Chemist",              glasses:true,  palette:{0:null,1:"#f4c7a3",2:"#9a6a2a",3:"#333",4:"#f0f0f0",5:"#c0a020",6:"#4a4a3a",7:"#5a4a3a",8:"#d46b6b",9:"#8a7020"}},
    physicist:         {label:"Physicist",            glasses:true,  palette:{0:null,1:"#e8b888",2:"#5a5a6a",3:"#333",4:"#2a2a3a",5:"#4a4a6a",6:"#2a2a3a",7:"#5a4a3a",8:"#d46b6b",9:"#3a3a5a"}},
    engineer:          {label:"Biomedical Engineer",  glasses:false, palette:{0:null,1:"#f4c7a3",2:"#5a3a20",3:"#333",4:"#8a6a20",5:"#c0a040",6:"#5a5a4a",7:"#6a5040",8:"#d46b6b",9:null}},
    psychologist:      {label:"Psychologist",         glasses:true,  palette:{0:null,1:"#d4a878",2:"#8a4a3a",3:"#333",4:"#6a5a8a",5:"#8a7ab0",6:"#4a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#5a4a80"}},
    ecologist:         {label:"Ecologist",            glasses:false, palette:{0:null,1:"#e8b888",2:"#7a6a3a",3:"#333",4:"#4a7a3a",5:"#6aa04a",6:"#5a6a4a",7:"#7a6a4a",8:"#d46b6b",9:null}},
    // --- Generic fallback ---
    generic:           {label:"Consultant",           glasses:false, palette:{0:null,1:"#f4c7a3",2:"#8a6a40",3:"#333",4:"#707070",5:"#909090",6:"#4a4a4a",7:"#5a4a3a",8:"#d46b6b",9:null}},
  };

  // ============================================================
  // SPRITE RENDERER
  // ============================================================
  function renderSprite(canvas, spriteData, palette) {
    const ctx = canvas.getContext("2d");
    const cols = spriteData[0].length;
    const rows = spriteData.length;
    const px = Math.floor(canvas.width / cols);
    const py = Math.floor(canvas.height / rows);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (let y = 0; y < rows; y++) {
      for (let x = 0; x < cols; x++) {
        const ci = spriteData[y][x];
        const color = palette[ci];
        if (!color) continue;
        ctx.fillStyle = color;
        ctx.fillRect(x * px, y * py, px, py);
      }
    }
  }

  function renderExpertSprite(canvas, expertType) {
    const def = EXPERT_DEFS[expertType] || EXPERT_DEFS.generic;
    const frame = def.glasses ? FRAME_GLASSES : FRAME_NO_GLASSES;
    renderSprite(canvas, frame, def.palette);
  }

  // ============================================================
  // STATE
  // ============================================================
  let currentState = null;
  let currentTurns = [];
  let lastTurnCount = 0;
  let cachedFiles = {};
  let labStartTime = null; // ISO timestamp from state.created_at

  // Domain-configurable labels (overridden from API state.domain_config)
  let domainLabels = {
    senior_label: "Professor",
    senior_short: "PI",
    junior_label: "Trainee",
    junior_short: "Trainee",
    overseer_label: "Editor",
    reviewer_label: "Reviewer",
    artifact: "Paper",
    meeting_log_name: "Meeting Log",
    consultant_label: "Consultant",
  };

  /** Apply domain labels from API response to the UI */
  function applyDomainLabels(dcfg) {
    if (!dcfg) return;
    Object.assign(domainLabels, dcfg);

    // Character card names (dynamic roster handles its own labels via renderTeamRoster)

    // Meeting log panel header
    const convHeader = document.querySelector("#conversation-panel .panel-header h2");
    if (convHeader) convHeader.textContent = domainLabels.meeting_log_name;

    // Paper panel header
    const paperHeader = document.querySelector(".paper-header h2");
    if (paperHeader) paperHeader.textContent = domainLabels.artifact;

    // Feedback dropdown
    const piOpt = document.querySelector('#feedback-target option[value="pi"]');
    if (piOpt) piOpt.textContent = domainLabels.senior_label;
    const trainOpt = document.querySelector('#feedback-target option[value="trainee"]');
    if (trainOpt) trainOpt.textContent = domainLabels.junior_label;

    // Editorial office text
    const editorTitle = document.getElementById("editor-title");
    if (editorTitle) editorTitle.textContent = domainLabels.overseer_label + "'s Desk";
    const editorialBtn = document.getElementById("btn-open-editorial");
    if (editorialBtn) editorialBtn.innerHTML = "&#x1F4EC; " + domainLabels.overseer_label + "'s Office";

    // Empty editorial state hint
    const emptyHint = document.querySelector("#editor-phase-empty .empty-hint");
    if (emptyHint) {
      emptyHint.textContent =
        `The ${domainLabels.overseer_label.toLowerCase()}'s office will open when the ` +
        `${domainLabels.senior_label} submits a ${domainLabels.artifact.toLowerCase()} for review.`;
    }

    // Experts/consultants section divider
    const expertsDivider = document.querySelector("#experts-section .experts-divider span");
    if (expertsDivider) expertsDivider.textContent = domainLabels.consultant_label + "s";
  }

  // Project-scoped URL builder (reads ?project= from location)
  const _urlProject = new URLSearchParams(window.location.search).get("project") || "";
  function apiUrl(path, extraParams) {
    const u = new URL(path, window.location.origin);
    if (_urlProject) u.searchParams.set("project", _urlProject);
    if (extraParams) {
      for (const [k, v] of Object.entries(extraParams)) u.searchParams.set(k, v);
    }
    return u.toString();
  }

  // ============================================================
  // SPRITE ANIMATION LOOP
  // Sprites are now created dynamically by renderTeamRoster().
  // The animation loop looks up canvases by ID each tick so it
  // works even after the DOM is rebuilt by the roster renderer.
  // ============================================================
  function animateSprites() {
    setInterval(() => {
      const piC = document.getElementById("pi-sprite");
      if (piC) {
        renderSprite(piC, PI_FRAME_2, PI_PALETTE);
        setTimeout(() => {
          const piC2 = document.getElementById("pi-sprite");
          if (piC2) renderSprite(piC2, PI_FRAME_1, PI_PALETTE);
        }, 200);
      }
    }, 3000);
    setInterval(() => {
      const trC = document.getElementById("trainee-sprite");
      if (trC) {
        renderSprite(trC, TRAINEE_FRAME_2, TRAINEE_PALETTE);
        setTimeout(() => {
          const trC2 = document.getElementById("trainee-sprite");
          if (trC2) renderSprite(trC2, TRAINEE_FRAME_1, TRAINEE_PALETTE);
        }, 200);
      }
    }, 4200);

    // Initial render (canvases may not exist yet; renderTeamRoster handles initial draw)
    const piC = document.getElementById("pi-sprite");
    if (piC) renderSprite(piC, PI_FRAME_1, PI_PALETTE);
    const trC = document.getElementById("trainee-sprite");
    if (trC) renderSprite(trC, TRAINEE_FRAME_1, TRAINEE_PALETTE);
  }

  // ============================================================
  // POLLING
  // ============================================================
  async function fetchState() {
    try {
      const resp = await fetch(apiUrl("/api/autolab/state"));
      const data = await resp.json();
      if (data.active) {
        currentState = data;
        if (data.files) cachedFiles = data.files;
        // Apply domain labels from config (idempotent)
        if (data.domain_config) applyDomainLabels(data.domain_config);
        // Set lab start time from state (only once)
        if (!labStartTime && data.created_at) {
          labStartTime = new Date(data.created_at).getTime();
          updateLabTimer();
        }
        updateStatusBar(data);
        renderTeamRoster(data);
        updateCharacterStatus(data);
        updateInventory(data);
        updatePaperProgress(data);
        updateProgressBar(data);
        updateExperts(data);
        updateEditorDesk(data);
        renderRecruitingDashboard(data);
        renderSkillsTab(data);
      }
    } catch (e) { console.warn("State fetch error:", e); }
  }

  async function fetchMeetingLog() {
    try {
      const resp = await fetch(apiUrl("/api/autolab/meeting-log"));
      const data = await resp.json();
      if (data.turns) {
        currentTurns = data.turns;
        if (data.turns.length !== lastTurnCount) {
          renderConversation(data.turns);
          lastTurnCount = data.turns.length;
        }
      }
    } catch (e) { console.warn("Meeting log fetch error:", e); }
  }

  // ============================================================
  // LAB TIMER
  // ============================================================
  function updateLabTimer() {
    const el = document.getElementById("lab-timer");
    if (!labStartTime || !el) return;
    const elapsed = Math.max(0, Math.floor((Date.now() - labStartTime) / 1000));
    const h = Math.floor(elapsed / 3600);
    const m = Math.floor((elapsed % 3600) / 60);
    const s = elapsed % 60;
    const pad = (n) => String(n).padStart(2, "0");
    el.textContent = h > 0
      ? `\u23F1 ${h}:${pad(m)}:${pad(s)}`
      : `\u23F1 ${pad(m)}:${pad(s)}`;
  }

  function startPolling() {
    fetchState();
    fetchMeetingLog();
    setInterval(fetchState, 2500);
    setInterval(fetchMeetingLog, 3000);
    setInterval(updateLabTimer, 1000);
  }

  // ============================================================
  // UI UPDATERS
  // ============================================================
  function updateStatusBar(state) {
    document.getElementById("iteration-num").textContent = state.iteration;
    // Use domain labels for the next-role badge
    const roleLabel = state.next_role === "pi"
      ? domainLabels.senior_short
      : (state.next_role === "trainee" ? domainLabels.junior_short || domainLabels.junior_label : state.next_role);
    document.getElementById("next-role").textContent = roleLabel.toUpperCase();
    document.getElementById("status-text").textContent =
      state.status === "active" ? "Running" :
      state.status === "ready_for_review" ? "Review" : state.status;
    const dot = document.getElementById("status-dot");
    dot.className = "status-dot";
    if (state.status === "active") dot.classList.add("working");

    // Project name in tab title and header
    if (state.project_name) {
      document.title = state.project_name + " — Autonomous Lab";
      const h = document.getElementById("project-name-header");
      if (h) h.textContent = state.project_name.toUpperCase();
    }
  }

  // Palette variants for dynamic trainees (cycling through distinct colors)
  const TRAINEE_PALETTE_VARIANTS = [
    // Original trainee blue
    {0:null,1:"#f4c7a3",2:"#7a5c3a",3:"#333",4:"#4a78bc",5:"#4a78bc",6:"#4a5a8a",7:"#e0d8c8",8:"#d46b6b",9:"#4a78bc"},
    // Teal
    {0:null,1:"#f4c7a3",2:"#3a5a5a",3:"#333",4:"#2a8a7a",5:"#50b0a0",6:"#2a5a5a",7:"#e0d8c8",8:"#d46b6b",9:"#2a8a7a"},
    // Purple
    {0:null,1:"#e8b888",2:"#4a2a5a",3:"#333",4:"#6a4a8a",5:"#8a6ab0",6:"#4a3a6a",7:"#e0d8c8",8:"#d46b6b",9:"#6a4a8a"},
    // Orange
    {0:null,1:"#f4c7a3",2:"#6a3a1a",3:"#333",4:"#aa6a2a",5:"#cc8a3a",6:"#6a4a2a",7:"#e0d8c8",8:"#d46b6b",9:"#aa6a2a"},
  ];

  let _lastTraineeKey = "";  // fingerprint to avoid re-rendering unchanged trainee cards
  let _dynamicSpriteIntervals = [];  // track intervals for cleanup

  function updateCharacterStatus(state) {
    // Character status is now primarily handled by renderTeamRoster().
    // This function retains multi-agent (orchestration=multi) trainee support
    // via the #trainees-container for backward compatibility.
    const traineesContainer = document.getElementById("trainees-container");

    // ── Multi-agent mode with trainees list ──
    const isMulti = state.orchestration === "multi" && state.trainees && state.trainees.length > 0;

    if (isMulti && traineesContainer) {
      traineesContainer.style.display = "block";

      // Build a key from trainee names + statuses to detect changes
      const traineeKey = state.trainees.map(t => `${t.name}:${t.status||"pending"}`).join("|");
      if (traineeKey === _lastTraineeKey) return;
      _lastTraineeKey = traineeKey;

      // Clear old sprite intervals
      _dynamicSpriteIntervals.forEach(id => clearInterval(id));
      _dynamicSpriteIntervals = [];

      // Rebuild trainee cards
      traineesContainer.innerHTML = "";

      // Divider header
      const divider = document.createElement("div");
      divider.className = "trainees-divider";
      divider.innerHTML = `<span>${domainLabels.junior_label}s (${state.trainees.length})</span>`;
      traineesContainer.appendChild(divider);

      state.trainees.forEach((t, idx) => {
        const card = document.createElement("div");
        const st = t.status || "pending";
        card.className = `trainee-dynamic-card ${st}`;
        card.dataset.name = t.name;

        // Canvas for sprite
        const canvas = document.createElement("canvas");
        canvas.width = 40;
        canvas.height = 64;
        canvas.className = "trainee-dynamic-sprite";
        canvas.style.imageRendering = "pixelated";
        const palette = TRAINEE_PALETTE_VARIANTS[idx % TRAINEE_PALETTE_VARIANTS.length];
        const frame = FRAME_NO_GLASSES;
        renderSprite(canvas, frame, palette);

        // Animate blink with staggered timing
        const blinkInterval = setInterval(() => {
          const blinkFrame = frame.map((r,y) => y === 3 ? [0,0,1,1,1,1,1,1,0,0] : r);
          renderSprite(canvas, blinkFrame, palette);
          setTimeout(() => renderSprite(canvas, frame, palette), 200);
        }, 3500 + idx * 700);
        _dynamicSpriteIntervals.push(blinkInterval);

        // Info section
        const info = document.createElement("div");
        info.className = "trainee-dynamic-info";

        const nameEl = document.createElement("div");
        nameEl.className = "trainee-dynamic-name";
        nameEl.textContent = t.name;
        info.appendChild(nameEl);

        if (t.focus) {
          const focusEl = document.createElement("div");
          focusEl.className = "trainee-focus";
          focusEl.textContent = t.focus;
          focusEl.title = t.focus;
          info.appendChild(focusEl);
        }

        // Status indicator
        const statusEl = document.createElement("div");
        statusEl.className = "trainee-dynamic-status";
        const statusMap = {
          pending: { cls: "idle", text: "Pending" },
          working: { cls: "working", text: "Working..." },
          done: { cls: "done", text: "Done" },
          failed: { cls: "failed", text: "Failed" },
        };
        const sm = statusMap[st] || statusMap.pending;
        statusEl.innerHTML = `<span class="status-indicator ${sm.cls}"></span><span class="status-label">${sm.text}</span>`;
        info.appendChild(statusEl);

        card.appendChild(canvas);
        card.appendChild(info);

        // Click to show thought bubble
        card.addEventListener("click", () => {
          showFixedBubble(card, t.focus || `${t.name} is ${sm.text.toLowerCase()}...`);
        });

        traineesContainer.appendChild(card);
      });

    } else if (traineesContainer) {
      traineesContainer.style.display = "none";
      _lastTraineeKey = "";
    }
  }

  // ============================================================
  // DYNAMIC TEAM ROSTER — renders PI + recruited characters
  // ============================================================
  let _lastRosterKey = "";  // fingerprint to avoid unnecessary re-renders

  function renderTeamRoster(state) {
    const roster = document.getElementById("team-roster");
    if (!roster) return;

    const recruitment = state.recruitment || { characters: [] };
    const characters = recruitment.characters || [];
    const nextRole = state.next_role || "";

    // Build a fingerprint to avoid re-rendering unchanged roster
    const rosterKey = JSON.stringify({
      nr: nextRole,
      chars: characters.map(c => `${c.slug}:${c.ready}:${JSON.stringify(c.skills||{})}`),
      dl: domainLabels.senior_short,
    });
    if (rosterKey === _lastRosterKey) return;
    _lastRosterKey = rosterKey;

    // Always show PI card
    let html = `
        <div class="char-card ${nextRole === 'pi' ? 'active active-pi' : ''}" data-role="pi">
            <div class="char-header">
                <span class="char-role-icon">&#x1F52C;</span>
                <span class="char-name">${escapeHtmlSafe(domainLabels.senior_label || 'PI')}</span>
            </div>
            <div class="char-sprite-wrap">
                <canvas id="pi-sprite" width="40" height="64"></canvas>
            </div>
            <div class="char-status">
                <span class="status-indicator ${nextRole === 'pi' ? 'working' : 'idle'}"></span>
                <span class="status-label">${nextRole === 'pi' ? 'Thinking...' : 'Idle'}</span>
            </div>
        </div>
    `;

    // Render each recruited character
    for (const char of characters) {
      const isActive = nextRole === char.slug;
      const skills = char.skills || {};
      const badgesHtml = Object.entries(skills).map(([name, info]) => {
        const status = (typeof info === "object" ? info.status : info) || "unknown";
        const icons = {certified: "\u2713", testing: "\u27F3", learning: "\u25CC", failed: "\u2717", queued: "\u00B7"};
        return `<span class="skill-badge ${escapeHtmlSafe(status)}"><span class="skill-icon">${icons[status] || '?'}</span>${escapeHtmlSafe(name)}</span>`;
      }).join("");

      html += `
          <div class="char-card ${isActive ? 'active active-trainee' : ''}" data-role="${escapeHtmlSafe(char.slug || 'trainee')}">
              <div class="char-header">
                  <span class="char-role-icon">&#x1F9EA;</span>
                  <span class="char-name">${escapeHtmlSafe(char.name || char.slug || 'Trainee')}</span>
              </div>
              <div class="char-sprite-wrap">
                  <canvas class="trainee-sprite" data-avatar="${escapeHtmlSafe(char.avatar || 'generic')}" width="40" height="64"></canvas>
              </div>
              <div class="char-status">
                  <span class="status-indicator ${isActive ? 'working' : 'idle'}"></span>
                  <span class="status-label">${isActive ? 'Working...' : (char.ready ? 'Ready' : 'Training')}</span>
              </div>
              ${badgesHtml ? '<div class="skill-badges">' + badgesHtml + '</div>' : ''}
          </div>
      `;
    }

    // If no characters recruited yet, show default trainee card
    if (characters.length === 0) {
      html += `
          <div class="char-card ${nextRole === 'trainee' ? 'active active-trainee' : ''}" data-role="trainee">
              <div class="char-header">
                  <span class="char-role-icon">&#x1F9EA;</span>
                  <span class="char-name">${escapeHtmlSafe(domainLabels.junior_label || 'Trainee')}</span>
              </div>
              <div class="char-sprite-wrap">
                  <canvas id="trainee-sprite" width="40" height="64"></canvas>
              </div>
              <div class="char-status">
                  <span class="status-indicator ${nextRole === 'trainee' ? 'working' : 'idle'}"></span>
                  <span class="status-label">${nextRole === 'trainee' ? 'Working...' : 'Idle'}</span>
              </div>
          </div>
      `;
    }

    roster.innerHTML = html;

    // Render sprites on the newly created canvases
    const piCanvas = document.getElementById("pi-sprite");
    if (piCanvas) renderSprite(piCanvas, PI_FRAME_1, PI_PALETTE);

    const defaultTrainee = document.getElementById("trainee-sprite");
    if (defaultTrainee) renderSprite(defaultTrainee, TRAINEE_FRAME_1, TRAINEE_PALETTE);

    // Render expert-type sprites for recruited characters
    roster.querySelectorAll("canvas.trainee-sprite").forEach((canvas) => {
      const avatarType = canvas.dataset.avatar || "generic";
      const def = EXPERT_DEFS[avatarType];
      if (def) {
        const frame = def.glasses ? FRAME_GLASSES : FRAME_NO_GLASSES;
        renderSprite(canvas, frame, def.palette);
      } else {
        renderSprite(canvas, TRAINEE_FRAME_1, TRAINEE_PALETTE);
      }
    });

    // Attach click handlers for thought bubbles on all cards
    roster.querySelectorAll(".char-card").forEach((card) => {
      card.addEventListener("click", () => {
        const role = card.dataset.role || "trainee";
        showThought(role);
      });
    });
  }

  // ============================================================
  // RECRUITING DASHBOARD — shown in center column during recruiting phase
  // ============================================================
  function renderRecruitingDashboard(state) {
    const panel = document.getElementById("conversation-content");
    if (!panel) return;
    if (state.phase !== "recruiting") return;

    const recruitment = state.recruitment || {};
    const characters = recruitment.characters || [];
    if (characters.length === 0) return;

    let html = '<div class="recruiting-dashboard">';
    html += '<div class="recruiting-title">&#x1F4CB; RECRUITING PHASE</div>';
    html += '<p style="color:var(--text-dim);margin-bottom:16px;">The ' + escapeHtmlSafe(domainLabels.senior_label) + ' is assembling a research team...</p>';

    for (const char of characters) {
      html += '<div class="recruiting-card">';
      html += '<div class="recruiting-card-header">';
      html += '<span class="recruiting-card-name">' + escapeHtmlSafe(char.name || char.slug || "Character") + '</span>';
      html += '<span class="recruiting-card-source">' + escapeHtmlSafe(char.source || "fresh") + '</span>';
      html += '</div>';

      for (const [skillName, info] of Object.entries(char.skills || {})) {
        const status = (typeof info === "object" ? info.status : info) || "queued";
        const widths = {certified: "100%", testing: "70%", learning: "40%", queued: "10%", failed: "100%"};
        const colors = {certified: "green", failed: "red", testing: "gold", learning: "gold", queued: "gold"};
        html += '<div class="skill-progress-row">';
        html += '<span class="skill-progress-name">' + escapeHtmlSafe(skillName) + '</span>';
        html += '<div class="skill-progress-bar"><div class="skill-progress-fill ' + escapeHtmlSafe(status) + '" style="width:' + (widths[status] || "10%") + '"></div></div>';
        html += '<span class="skill-progress-status" style="color:var(--' + (colors[status] || "gold") + ')">' + escapeHtmlSafe(status) + '</span>';
        html += '</div>';
      }

      html += '</div>';
    }

    html += '</div>';
    panel.innerHTML = html;
  }

  // ============================================================
  // PROGRESS BAR — PI-determined
  // ============================================================
  function updateProgressBar(state) {
    const pct = Math.max(0, Math.min(100, state.progress || 0));
    const fill = document.getElementById("paper-bar-fill");
    const label = document.getElementById("paper-bar-label");
    if (fill) fill.style.width = pct + "%";
    if (label) label.textContent = pct + "%";
  }

  function updatePaperProgress(state) {
    const progress = state.paper_progress || {};
    const sections = ["abstract", "introduction", "methods", "results", "discussion"];
    sections.forEach((sec) => {
      const info = progress[sec] || { exists: false, words: 0 };
      const el = document.getElementById(`sec-${sec}`);
      if (el) {
        if (info.exists && info.words > 0) {
          el.textContent = `${info.words}w`;
          el.className = "sec-words has-content";
          el.title = `Click to view ${sec}.tex`;
          el.onclick = () => openFilePreview(`paper/sections/${sec}.tex`);
        } else {
          el.textContent = "--";
          el.className = "sec-words";
          el.style.cursor = "default";
          el.onclick = null;
        }
      }
    });
  }

  // ============================================================
  // EXPERTS
  // ============================================================
  function updateExperts(state) {
    const experts = state.experts || [];
    const section = document.getElementById("experts-section");
    const list = document.getElementById("experts-list");
    if (!section || !list) return;

    if (experts.length === 0) {
      section.style.display = "none";
      return;
    }
    section.style.display = "block";

    if (list.dataset.count === String(experts.length)) return;
    list.dataset.count = String(experts.length);
    list.innerHTML = "";

    experts.forEach((exp) => {
      const card = document.createElement("div");
      card.className = "expert-card";

      const canvas = document.createElement("canvas");
      canvas.width = 40;
      canvas.height = 64;
      canvas.style.imageRendering = "pixelated";
      renderExpertSprite(canvas, exp.avatar || "generic");

      const info = document.createElement("div");
      info.className = "expert-info";
      info.innerHTML = `<div class="expert-name">${escapeHtmlSafe(exp.name || "Expert")}</div><div class="expert-role">${escapeHtmlSafe(exp.role || "")}</div>`;

      card.appendChild(canvas);
      card.appendChild(info);

      card.addEventListener("click", () => showFixedBubble(card, exp.thought || `${exp.role || "Consulting"}...`));

      list.appendChild(card);
    });
  }

  // ============================================================
  // INVENTORY — Expandable folders with file lists
  // ============================================================
  const FILE_ICONS = {
    ".py":"\uD83D\uDC0D",".tex":"\uD83D\uDCC4",".bib":"\uD83D\uDCDA",
    ".csv":"\uD83D\uDCCA",".tsv":"\uD83D\uDCCA",".json":"\u2699\uFE0F",
    ".png":"\uD83D\uDDBC\uFE0F",".jpg":"\uD83D\uDDBC\uFE0F",".jpeg":"\uD83D\uDDBC\uFE0F",
    ".pdf":"\uD83D\uDCC4",".svg":"\uD83C\uDFA8",".md":"\uD83D\uDCDD",
    ".txt":"\uD83D\uDCDD",".sh":"\uD83D\uDCBB",".log":"\uD83D\uDCDC",
    ".parquet":"\uD83D\uDCCA",".h5ad":"\uD83E\uDDEC",".h5":"\uD83E\uDDEC",
  };

  function getFileIcon(filename) {
    const ext = "." + filename.split(".").pop().toLowerCase();
    return FILE_ICONS[ext] || "\uD83D\uDCC1";
  }

  function isImageFile(filename) {
    const ext = filename.split(".").pop().toLowerCase();
    return ["png","jpg","jpeg","gif","webp","svg","bmp","tiff","tif"].includes(ext);
  }

  function updateInventory(state) {
    const counts = state.file_counts || {};
    const files = state.files || {};
    const folders = ["data","scripts","figures","results","paper"];

    folders.forEach((f) => {
      const countEl = document.querySelector(`#inv-${f} .inv-count`);
      if (countEl) countEl.textContent = counts[f] || 0;

      let listEl = document.getElementById(`inv-${f}-files`);
      if (!listEl) {
        const section = document.getElementById(`inv-${f}`);
        if (section) {
          listEl = document.createElement("div");
          listEl.className = "inv-file-list";
          listEl.id = `inv-${f}-files`;
          section.after(listEl);
        }
      }

      if (listEl) {
        const folderFiles = files[f] || [];
        if (listEl.dataset.count !== String(folderFiles.length)) {
          listEl.dataset.count = String(folderFiles.length);
          listEl.innerHTML = "";
          folderFiles.forEach((filePath) => {
            const name = filePath.split("/").pop();
            const item = document.createElement("div");
            item.className = "inv-file";
            item.innerHTML = `<span class="file-icon">${getFileIcon(name)}</span><span class="file-name" title="${escapeHtmlSafe(filePath)}">${escapeHtmlSafe(name)}</span>`;
            item.addEventListener("click", (e) => {
              e.stopPropagation();
              openFilePreview(filePath);
            });
            listEl.appendChild(item);
          });
        }
      }
    });
  }

  function setupInventory() {
    ["data","scripts","figures","results","paper"].forEach((f) => {
      const el = document.getElementById(`inv-${f}`);
      if (!el) return;
      el.classList.add("inv-folder");
      const row = el.querySelector(".inv-row");
      if (row) {
        const arrow = document.createElement("span");
        arrow.className = "folder-arrow";
        arrow.textContent = "\u25B6";
        row.insertBefore(arrow, row.firstChild);
      }
      el.addEventListener("click", () => {
        el.classList.toggle("expanded");
        const listEl = document.getElementById(`inv-${f}-files`);
        if (listEl) listEl.style.display = el.classList.contains("expanded") ? "block" : "none";
      });
    });
  }

  // ============================================================
  // FILE PREVIEW MODAL
  // ============================================================
  function openFilePreview(filePath) {
    const overlay = document.getElementById("file-modal-overlay");
    const title = document.getElementById("modal-title");
    const body = document.getElementById("modal-body");

    const fileName = filePath.split("/").pop();
    title.textContent = fileName;
    body.innerHTML = '<div class="loading-text">Loading...</div>';
    overlay.classList.remove("hidden");

    if (isImageFile(fileName)) {
      const img = document.createElement("img");
      img.className = "preview-image";
      img.src = apiUrl("/api/autolab/file", { path: filePath, t: Date.now() });
      img.alt = fileName;
      img.onload = () => { body.innerHTML = ""; body.appendChild(img); };
      img.onerror = () => {
        body.innerHTML = '<div class="loading-text">Could not load image.<br>The file may not exist yet.</div>';
      };
    } else {
      fetch(apiUrl("/api/autolab/file", { path: filePath }))
        .then((r) => {
          if (!r.ok) throw new Error(`HTTP ${r.status}`);
          return r.json();
        })
        .then((data) => {
          if (data.error) {
            body.innerHTML = `<div class="loading-text">${escapeHtmlSafe(data.error)}</div>`;
            return;
          }
          const lang = data.language || "text";
          body.innerHTML = `
            <span class="preview-lang">${lang.toUpperCase()}</span>
            <pre class="preview-code">${escapeHtmlSafe(data.content || "")}</pre>`;
        })
        .catch((err) => {
          body.innerHTML = `<div class="loading-text">Failed to load file: ${escapeHtmlSafe(err.message)}</div>`;
        });
    }
  }

  function closeModal() {
    document.getElementById("file-modal-overlay").classList.add("hidden");
  }

  function setupModal() {
    document.getElementById("modal-close").addEventListener("click", closeModal);
    document.getElementById("file-modal-overlay").addEventListener("click", (e) => {
      if (e.target === e.currentTarget) closeModal();
    });
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") closeModal();
    });
  }

  // ============================================================
  // THOUGHT BUBBLES — position:fixed to avoid clipping
  // ============================================================
  const IDLE_THOUGHTS = {
    pi: [
      "Need higher impact...",
      "Is this novel enough?",
      "Check figure quality.",
      "What's the story here?",
      "Hmm, interesting data...",
      "Nature-worthy figures?",
      "Invite a reviewer?",
      "What about controls?",
      "Need a statistician...",
      "Check biomedical toolkit?",
    ],
    trainee: [
      "Running the pipeline...",
      "Let me check PubMed...",
      "Debugging this script...",
      "Almost done with plots!",
      "Reading the paper...",
      "Updating the methods...",
      "Checking the stats...",
      "Hmm, p < 0.05?",
      "Normalizing the data...",
      "Trying a new model...",
    ],
  };

  let thoughtTimers = { pi: null, trainee: null };

  function getLastThought(role) {
    for (let i = currentTurns.length - 1; i >= 0; i--) {
      if (currentTurns[i].role === role && currentTurns[i].summary) {
        const first = currentTurns[i].summary.split(/[.!?\n]/)[0].trim();
        return first.length > 50 ? first.slice(0, 47) + "..." : first;
      }
    }
    const pool = IDLE_THOUGHTS[role];
    return pool[Math.floor(Math.random() * pool.length)];
  }

  /** Show a fixed-position pixel bubble near a DOM element */
  function showFixedBubble(anchorEl, text) {
    // Remove any existing fixed bubble
    const old = document.getElementById("fixed-thought-bubble");
    if (old) old.remove();

    const rect = anchorEl.getBoundingClientRect();
    const bubble = document.createElement("div");
    bubble.id = "fixed-thought-bubble";
    bubble.className = "thought-bubble visible";
    bubble.textContent = text;
    bubble.style.position = "fixed";
    bubble.style.left = (rect.right + 10) + "px";
    bubble.style.top = (rect.top + rect.height / 2) + "px";
    bubble.style.transform = "translateY(-50%)";
    bubble.style.zIndex = "100";
    document.body.appendChild(bubble);

    setTimeout(() => { if (bubble.parentNode) bubble.remove(); }, 3500);
  }

  function showThought(role) {
    // Dynamic roster: find card by data-role attribute
    const card = document.querySelector(`.char-card[data-role="${role}"]`);
    if (!card) return;
    const text = getLastThought(role);
    showFixedBubble(card, text);
  }

  function setupThoughtBubbles() {
    // Click handlers are attached by renderTeamRoster() since cards are dynamic.
    // Schedule random thought bubbles.
    function scheduleRandom() {
      const delay = 12000 + Math.random() * 13000;
      setTimeout(() => {
        showThought(Math.random() < 0.5 ? "pi" : "trainee");
        scheduleRandom();
      }, delay);
    }
    setTimeout(scheduleRandom, 8000);
  }

  // ============================================================
  // MARKDOWN RENDERER — lightweight
  // ============================================================
  function renderMarkdown(text) {
    let html = text
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");

    // Headings
    html = html.replace(/^### (.+)$/gm, '<div class="md-h3">$1</div>');
    html = html.replace(/^## (.+)$/gm, '<div class="md-h2">$1</div>');

    // Bold then italic
    html = html.replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>");
    html = html.replace(/\*(.+?)\*/g, "<em>$1</em>");

    // Lists
    html = html.replace(/^- (.+)$/gm, '<div class="md-li">\u2022 $1</div>');
    html = html.replace(/^(\d+)\. (.+)$/gm, '<div class="md-li">$1. $2</div>');

    // Paragraph breaks & line breaks
    html = html.replace(/\n\n/g, '<div class="md-break"></div>');
    html = html.replace(/\n/g, "<br>");

    return html;
  }

  // ============================================================
  // CONVERSATION RENDERER
  // ============================================================
  function renderConversation(turns) {
    const container = document.getElementById("conversation-content");
    if (!turns.length) {
      container.innerHTML =
        `<div class="empty-state"><p>Waiting for the first turn...</p><p class="hint">The ${domainLabels.senior_label} will speak first.</p></div>`;
      return;
    }

    container.innerHTML = "";
    turns.forEach((turn, idx) => {
      const bubble = document.createElement("div");
      bubble.className = `turn-bubble ${turn.role}`;

      const roleIcon = turn.role === "pi" ? "\uD83D\uDD2C" : "\uD83E\uDDEA";
      // Multi-trainee turns have summaries with " | " separators
      const isMultiTrainee = turn.role === "trainee" && turn.summary && turn.summary.includes(" | ");
      const roleName = turn.role === "pi"
        ? domainLabels.senior_label
        : (isMultiTrainee ? domainLabels.junior_label + "s (team)" : domainLabels.junior_label);
      const summaryText = turn.summary || "(no summary)";
      const hasDetails = turn.content && turn.content.length > 0;
      const detailId = `detail-${idx}`;

      bubble.innerHTML = `
        <div class="bubble-header">
          <span class="bubble-role ${turn.role}"><span>${roleIcon}</span> ${roleName}</span>
          <span class="bubble-meta">Iter ${turn.iteration} &middot; ${turn.date || ""}</span>
        </div>
        <div class="bubble-summary">${renderMarkdown(summaryText)}</div>
        ${hasDetails
          ? `<div class="bubble-details collapsed" id="${detailId}">${renderMarkdown(turn.content)}</div>
             <button class="toggle-details" data-target="${detailId}">Show details</button>`
          : ""}`;
      container.appendChild(bubble);
    });

    container.querySelectorAll(".toggle-details").forEach((btn) => {
      btn.addEventListener("click", () => {
        const target = document.getElementById(btn.dataset.target);
        if (target) {
          const collapsed = target.classList.toggle("collapsed");
          btn.textContent = collapsed ? "Show details" : "Hide details";
        }
      });
    });

    const scroll = document.getElementById("conversation-scroll");
    if (lastTurnCount > 0) {
      scroll.scrollTop = scroll.scrollHeight;
    } else {
      scroll.scrollTop = 0;
    }
  }

  // ============================================================
  // RIGHT COLUMN TABS
  // ============================================================
  function setupRightTabs() {
    const tabs = document.querySelectorAll('.right-tab');
    tabs.forEach(tab => {
      tab.addEventListener('click', () => {
        tabs.forEach(t => t.classList.remove('active'));
        document.querySelectorAll('.right-tab-content').forEach(c => c.classList.remove('active'));
        tab.classList.add('active');
        const targetId = 'tab-' + tab.dataset.tab;
        const target = document.getElementById(targetId);
        if (target) target.classList.add('active');
      });
    });
  }

  function renderSkillsTab(state) {
    const container = document.getElementById('skills-content');
    if (!container) return;

    const chars = ((state.recruitment || {}).characters) || [];
    if (!chars.length) {
      container.innerHTML = '<p style="font-size:11px;color:var(--text-dim);">No characters recruited yet.</p>';
      return;
    }

    let html = '<div style="font-family:\'Press Start 2P\',monospace;font-size:8px;color:var(--gold);margin-bottom:8px;">TEAM SKILLS</div>';

    for (const char of chars) {
      for (const [name, info] of Object.entries(char.skills || {})) {
        const status = info.status || 'unknown';
        const colorMap = {certified: 'var(--green)', testing: 'var(--gold)', learning: 'var(--text-dim)', failed: 'var(--red)', queued: 'var(--text-muted)'};
        const color = colorMap[status] || 'var(--text-dim)';
        html += '<div style="background:var(--bg-panel-light, var(--bg-panel));border:1px solid var(--border);border-radius:4px;padding:8px;margin-bottom:6px;">';
        html += '<div style="display:flex;justify-content:space-between;align-items:center;">';
        html += '<span style="font-size:11px;color:var(--text);font-weight:600;">' + escapeHtmlSafe(name) + '</span>';
        html += '<span style="font-family:\'Press Start 2P\',monospace;font-size:7px;color:' + color + ';">' + escapeHtmlSafe(status) + '</span>';
        html += '</div>';
        html += '<div style="font-size:10px;color:var(--text-dim);margin-top:4px;">Owner: ' + escapeHtmlSafe(char.slug || char.name || '?') + '</div>';
        if (info.source) html += '<div style="font-size:10px;color:var(--text-dim);">Source: ' + escapeHtmlSafe(info.source) + '</div>';
        if (info.tests_passed !== undefined) html += '<div style="font-size:10px;color:var(--text-dim);">Tests: ' + info.tests_passed + '/' + (info.tests_total || '?') + ' pass</div>';
        html += '</div>';
      }
    }

    let certified = 0, testing = 0, queued = 0;
    chars.forEach(function(c) {
      Object.values(c.skills || {}).forEach(function(i) {
        if (i.status === 'certified') certified++;
        else if (i.status === 'testing' || i.status === 'learning') testing++;
        else queued++;
      });
    });
    html += '<div style="font-size:10px;color:var(--text-dim);margin-top:8px;padding-top:8px;border-top:1px solid var(--border);">';
    html += certified + ' certified &nbsp; ' + testing + ' in-progress &nbsp; ' + queued + ' queued</div>';

    container.innerHTML = html;
  }

  // ============================================================
  // HELPERS
  // ============================================================
  function escapeHtmlSafe(text) {
    const div = document.createElement("div");
    div.textContent = text;
    return div.innerHTML;
  }

  // ============================================================
  // EDITOR'S DESK — Journal-style editorial workflow
  // ============================================================

  // Available reviewer pool (subset of EXPERT_DEFS useful for peer review)
  const REVIEWER_POOL = [
    {name:"Dr. Alvarez",  role:"Immunologist",         avatar:"immunologist"},
    {name:"Dr. Yamamoto", role:"Statistician",          avatar:"statistician"},
    {name:"Dr. Okonkwo",  role:"Oncologist",            avatar:"oncologist"},
    {name:"Dr. Chen",     role:"Bioinformatician",      avatar:"bioinformatician"},
    {name:"Dr. Berg",     role:"Neuroscientist",        avatar:"neuroscientist"},
    {name:"Dr. Patel",    role:"Geneticist",            avatar:"geneticist"},
    {name:"Dr. Kowalski", role:"Cell Biologist",        avatar:"cell_biologist"},
    {name:"Dr. Singh",    role:"Pharmacologist",        avatar:"pharmacologist"},
    {name:"Dr. Müller",   role:"Structural Biologist",  avatar:"structural_bio"},
    {name:"Dr. Kim",      role:"Epidemiologist",        avatar:"epidemiologist"},
    {name:"Dr. Santos",   role:"Data Scientist",        avatar:"data_scientist"},
    {name:"Dr. Zhao",     role:"Systems Biologist",     avatar:"systems_biologist"},
    {name:"Dr. Petrov",   role:"Pathologist",           avatar:"pathologist"},
    {name:"Dr. Harper",   role:"Clinician",             avatar:"clinician"},
    {name:"Dr. Rossi",    role:"Microbiologist",        avatar:"microbiologist"},
    {name:"Dr. Tanaka",   role:"Comp. Biologist",       avatar:"comp_biologist"},
    {name:"Dr. Johansson",role:"ML Engineer",           avatar:"ml_engineer"},
    {name:"Dr. Dubois",   role:"Chemist",               avatar:"chemist"},
    {name:"Dr. Fischer",  role:"Bioethicist",           avatar:"bioethicist"},
    {name:"Dr. Lee",      role:"Radiologist",           avatar:"radiologist"},
    {name:"Dr. Thompson", role:"Biomedical Engineer",   avatar:"engineer"},
    {name:"Dr. Vasquez",  role:"Psychologist",          avatar:"psychologist"},
    {name:"Dr. Nguyen",   role:"Critical Reviewer",     avatar:"reviewer"},
    {name:"Dr. Brown",    role:"Science Writer",        avatar:"science_writer"},
  ];

  let selectedReviewers = [];
  let editorPhase = "none";

  let _selectingReviewers = false; // flag to prevent poll reset during reviewer selection
  let _manualOpen = false;          // flag: user manually opened the overlay
  let _editorialBusy = false;       // flag: an editorial API call is in flight
  let _userDismissedPhase = null;   // phase the user dismissed — don't re-show until phase changes
  let _lastEditorialHash = "";      // fingerprint of last-rendered editorial state

  function openEditorDesk() {
    const overlay = document.getElementById("editor-desk-overlay");
    if (!overlay) return;
    _manualOpen = true;
    _userDismissedPhase = null; // clear dismiss memory when manually opening
    overlay.classList.remove("hidden");
    // Populate with latest editorial state
    const editorial = (currentState && currentState.editorial) || { phase: "none" };
    const roundEl = document.getElementById("editor-round");
    if (roundEl) roundEl.textContent = editorial.round ? `Round ${editorial.round}` : "";
    document.querySelectorAll(".editor-phase").forEach(el => el.style.display = "none");

    const phase = editorial.phase || "none";
    if (phase === "submitted") {
      showPhaseSubmitted(editorial);
    } else if (phase === "reviewers_invited" || phase === "under_review") {
      showPhaseReviewing(editorial);
    } else if (phase === "reviews_complete") {
      showPhaseDecision(editorial);
    } else if (phase === "decision_made") {
      showPhaseDone(editorial);
    } else {
      // Nothing submitted yet — show empty state
      document.getElementById("editor-phase-empty").style.display = "block";
    }
  }

  function closeEditorDesk() {
    const overlay = document.getElementById("editor-desk-overlay");
    if (overlay) overlay.classList.add("hidden");
    _manualOpen = false;
    _selectingReviewers = false;
    // Remember which phase we dismissed so the poll doesn't re-open it
    _userDismissedPhase = editorPhase || null;
  }

  function updateEditorDesk(state) {
    const editorial = state.editorial || { phase: "none" };
    const overlay = document.getElementById("editor-desk-overlay");
    if (!overlay) return;

    const phase = editorial.phase || "none";
    editorPhase = phase;

    // If an editorial API call is in flight, skip ALL poll updates to prevent UI reset
    if (_editorialBusy) return;

    // If user is actively selecting reviewers, don't let the poll override the invite view
    if (_selectingReviewers && phase === "submitted") return;

    // If user manually opened the overlay and phase is "none", keep it open
    if (_manualOpen && phase === "none") return;

    // ── Auto-show logic ──
    // Only auto-pop the overlay when the EDITOR needs to act:
    //   "submitted"        → editor decides: desk reject or invite reviewers
    //   "reviews_complete"  → editor makes final decision (accept/minor/major/reject)
    //
    // Do NOT auto-pop during:
    //   "reviewers_invited" / "under_review" → AI is playing reviewers, no editor action
    //   "decision_made" → decision was already handled, PI revision loop is starting
    //   "none" → no editorial workflow
    //
    // The user can ALWAYS manually view any phase via the "Editorial Office" button.
    const actionPhases = ["submitted", "reviews_complete"];

    if (actionPhases.includes(phase)) {
      // Editor needs to act — auto-show (unless user dismissed this exact phase)
      if (_userDismissedPhase === phase) return;
      if (_userDismissedPhase && _userDismissedPhase !== phase) _userDismissedPhase = null;
      overlay.classList.remove("hidden");
      _manualOpen = false;
    } else if (_manualOpen) {
      // User manually opened — keep open, update content below
    } else {
      // No action needed and not manually opened — hide
      overlay.classList.add("hidden");
      _selectingReviewers = false;
      return;
    }

    // Update round
    const roundEl = document.getElementById("editor-round");
    if (roundEl) roundEl.textContent = `Round ${editorial.round || 1}`;

    // Update timeout countdown badge
    const badge = document.getElementById("editor-timeout-badge");
    if (badge) {
      const timeoutMin = state.editor_timeout_minutes || 0;
      const waitingSince = editorial.waiting_since || "";
      const needsAction = ["submitted", "reviews_complete"].includes(phase);

      if (timeoutMin > 0 && waitingSince && needsAction) {
        const started = new Date(waitingSince);
        const now = new Date();
        const elapsedMs = now - started;
        const timeoutMs = timeoutMin * 60 * 1000;
        const remainMs = Math.max(0, timeoutMs - elapsedMs);
        const remainMin = Math.floor(remainMs / 60000);
        const remainSec = Math.floor((remainMs % 60000) / 1000);

        if (remainMs <= 0) {
          badge.textContent = `AI ${domainLabels.overseer_label.toUpperCase()} ACTIVE`;
          badge.classList.add("urgent");
        } else {
          const padSec = String(remainSec).padStart(2, "0");
          badge.textContent = `AI in ${remainMin}:${padSec}`;
          badge.classList.toggle("urgent", remainMin < 5);
        }
      } else if (timeoutMin === 0 && needsAction) {
        badge.textContent = "NO TIMEOUT";
        badge.classList.remove("urgent");
      } else {
        badge.textContent = "";
        badge.classList.remove("urgent");
      }
    }

    // ── Fingerprint the editorial state to avoid redundant re-renders ──
    // Re-rendering destroys DOM elements and resets scroll positions.
    const editorialHash = JSON.stringify({
      p: phase,
      r: editorial.round,
      rv: Object.keys(editorial.reviews || {}),
      d: editorial.decision,
    });

    if (editorialHash === _lastEditorialHash) {
      // Nothing changed — skip the expensive DOM rebuild to preserve scroll
      return;
    }
    _lastEditorialHash = editorialHash;

    // Hide all phases
    document.querySelectorAll(".editor-phase").forEach(el => el.style.display = "none");

    if (phase === "submitted") {
      showPhaseSubmitted(editorial);
    } else if (phase === "reviewers_invited" || phase === "under_review") {
      _selectingReviewers = false;
      showPhaseReviewing(editorial);
    } else if (phase === "reviews_complete") {
      _selectingReviewers = false;
      showPhaseDecision(editorial);
    } else if (phase === "decision_made") {
      _selectingReviewers = false;
      showPhaseDone(editorial);
    } else if (phase === "none") {
      // Show empty state for manual open
      const emptyEl = document.getElementById("editor-phase-empty");
      if (emptyEl) emptyEl.style.display = "block";
    }
  }

  function showPhaseSubmitted(editorial) {
    const panel = document.getElementById("editor-phase-submitted");
    panel.style.display = "block";
    document.getElementById("editor-cover-letter").innerHTML =
      renderMarkdown(editorial.cover_letter || "(No cover letter provided)");
    // Show manuscript info from paper progress
    const info = document.getElementById("editor-manuscript-info");
    if (currentState && currentState.paper_progress) {
      const pp = currentState.paper_progress;
      const lines = Object.entries(pp).map(([sec, d]) =>
        `<div><strong>${sec}</strong>: ${d.words > 0 ? d.words + " words" : "not written"}</div>`
      );
      info.innerHTML = lines.join("");
    } else {
      info.textContent = "Manuscript sections will be shown here.";
    }
  }

  function showPhaseReviewing(editorial) {
    const panel = document.getElementById("editor-phase-reviewing");
    panel.style.display = "block";
    const list = document.getElementById("review-progress-list");
    list.innerHTML = "";
    (editorial.reviewers || []).forEach(rev => {
      const done = editorial.reviews && editorial.reviews[rev.id];
      const isNext = currentState && currentState.next_role === rev.id;
      const item = document.createElement("div");
      item.className = "review-progress-item";
      const canvas = document.createElement("canvas");
      canvas.width = 24; canvas.height = 38;
      canvas.style.imageRendering = "pixelated";
      renderExpertSprite(canvas, rev.avatar || "generic");
      item.appendChild(canvas);
      const info = document.createElement("div");
      info.innerHTML = `<div style="font-family:'Press Start 2P';font-size:7px;color:var(--gold);">${escapeHtmlSafe(rev.name)}</div><div style="font-size:9px;color:var(--text-dim);">${escapeHtmlSafe(rev.role)}</div>`;
      item.appendChild(info);
      const status = document.createElement("span");
      status.className = "rpi-status " + (done ? "rpi-done" : isNext ? "rpi-active" : "rpi-pending");
      status.textContent = done ? "DONE" : isNext ? "REVIEWING..." : "PENDING";
      item.appendChild(status);
      list.appendChild(item);
    });
  }

  function showPhaseDecision(editorial) {
    const panel = document.getElementById("editor-phase-decision");
    panel.style.display = "block";
    const reports = document.getElementById("review-reports");
    reports.innerHTML = "";
    (editorial.reviewers || []).forEach(rev => {
      const review = (editorial.reviews || {})[rev.id];
      if (!review) return;
      const card = document.createElement("div");
      card.className = "review-report";
      const hdr = document.createElement("div");
      hdr.className = "review-report-header";
      const canvas = document.createElement("canvas");
      canvas.width = 24; canvas.height = 38;
      canvas.style.imageRendering = "pixelated";
      renderExpertSprite(canvas, rev.avatar || "generic");
      hdr.appendChild(canvas);
      const nameEl = document.createElement("span");
      nameEl.className = "rev-header-name";
      nameEl.textContent = `${rev.name} (${rev.role})`;
      hdr.appendChild(nameEl);
      const recBadge = document.createElement("span");
      const rec = (review.recommendation || "").replace(/ /g, "_");
      recBadge.className = `rev-header-rec rec-${rec}`;
      recBadge.textContent = (review.recommendation || "N/A").toUpperCase().replace(/_/g, " ");
      hdr.appendChild(recBadge);
      card.appendChild(hdr);
      const body = document.createElement("div");
      body.className = "review-report-body";
      body.innerHTML = renderMarkdown(review.report || "(No report)");
      card.appendChild(body);
      reports.appendChild(card);
    });
  }

  function showPhaseDone(editorial) {
    const panel = document.getElementById("editor-phase-done");
    panel.style.display = "block";
    const result = document.getElementById("decision-result");
    const dec = (editorial.decision || "").toUpperCase().replace(/_/g, " ");
    result.innerHTML = `Decision: ${dec}<br><br>${escapeHtmlSafe(editorial.decision_feedback || "")}<br><br>Returning to ${domainLabels.senior_label}...`;
    // Auto-hide after 5 seconds only if not user-opened
    if (!_manualOpen) {
      setTimeout(() => {
        document.getElementById("editor-desk-overlay").classList.add("hidden");
      }, 5000);
    }
  }

  function setupEditorDesk() {
    // Open editorial office button (top bar)
    document.getElementById("btn-open-editorial").addEventListener("click", openEditorDesk);

    // Close button (X in overlay header)
    document.getElementById("btn-close-editorial").addEventListener("click", closeEditorDesk);

    // Also close on clicking the backdrop (outside the desk)
    document.getElementById("editor-desk-overlay").addEventListener("click", (e) => {
      if (e.target.id === "editor-desk-overlay") closeEditorDesk();
    });

    // Helper: editorial API call with lock, error handling, and project_dir
    async function editorialFetch(url, body) {
      _editorialBusy = true;
      // Attach project_dir from current state for multi-instance safety
      if (currentState && currentState.project_dir) {
        body.project_dir = currentState.project_dir;
      }
      try {
        const resp = await fetch(url, {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify(body),
        });
        if (!resp.ok) {
          const err = await resp.json().catch(() => ({}));
          const msg = err.error || `Server error (${resp.status})`;
          console.error("Editorial API error:", msg);
          showEditorToast("Error: " + msg);
          return false;
        }
        return true;
      } catch (e) {
        console.error("Editorial fetch failed:", e);
        showEditorToast("Network error - is the server running?");
        return false;
      } finally {
        // Brief delay so the next poll picks up the new state
        await new Promise(r => setTimeout(r, 600));
        _editorialBusy = false;
        fetchState();
      }
    }

    // Desk reject — show inline feedback panel
    document.getElementById("btn-desk-reject").addEventListener("click", () => {
      // Hide the action buttons, show the feedback panel
      const actionsDiv = document.getElementById("btn-desk-reject").closest(".editor-actions");
      if (actionsDiv) actionsDiv.style.display = "none";
      const panel = document.getElementById("desk-reject-panel");
      panel.style.display = "block";
      const textarea = document.getElementById("desk-reject-feedback");
      textarea.value = "";
      textarea.focus();
    });

    // Cancel desk reject — go back to action buttons
    document.getElementById("btn-cancel-reject").addEventListener("click", () => {
      document.getElementById("desk-reject-panel").style.display = "none";
      // Re-show the action buttons
      const actionsDiv = document.querySelector("#editor-phase-submitted .editor-actions");
      if (actionsDiv) actionsDiv.style.display = "";
    });

    // Confirm desk reject — send feedback
    document.getElementById("btn-confirm-reject").addEventListener("click", async () => {
      const feedback = document.getElementById("desk-reject-feedback").value.trim();
      if (!feedback) {
        showEditorToast("Please provide feedback to the authors before rejecting.");
        return;
      }
      await editorialFetch("/api/autolab/editorial/desk-reject", { feedback });
    });

    // Invite reviewers
    document.getElementById("btn-invite-reviewers").addEventListener("click", () => {
      _selectingReviewers = true;
      document.getElementById("editor-phase-submitted").style.display = "none";
      showReviewerSelection();
    });

    // Back button
    document.getElementById("btn-back-to-submission").addEventListener("click", () => {
      _selectingReviewers = false;
      document.getElementById("editor-phase-invite").style.display = "none";
      document.getElementById("editor-phase-submitted").style.display = "block";
    });

    // Confirm reviewers
    document.getElementById("btn-confirm-reviewers").addEventListener("click", async () => {
      if (selectedReviewers.length === 0) return;
      _selectingReviewers = false;
      await editorialFetch("/api/autolab/editorial/invite-reviewers", {
        reviewers: selectedReviewers,
      });
      selectedReviewers = [];
    });

    // Decision buttons
    document.querySelectorAll(".decision-buttons .editor-btn").forEach(btn => {
      btn.addEventListener("click", async () => {
        const decision = btn.dataset.decision;
        if (!decision) return;
        const feedback = document.getElementById("editor-feedback").value.trim();
        await editorialFetch("/api/autolab/editorial/decision", {
          decision, feedback,
        });
      });
    });
  }

  function showReviewerSelection() {
    const panel = document.getElementById("editor-phase-invite");
    panel.style.display = "block";
    selectedReviewers = [];
    updateSelectedDisplay();

    const pool = document.getElementById("reviewer-pool");
    pool.innerHTML = "";
    REVIEWER_POOL.forEach(rev => {
      const opt = document.createElement("div");
      opt.className = "reviewer-option";
      const canvas = document.createElement("canvas");
      canvas.width = 24; canvas.height = 38;
      canvas.style.imageRendering = "pixelated";
      renderExpertSprite(canvas, rev.avatar);
      opt.appendChild(canvas);
      const info = document.createElement("div");
      info.className = "rev-info";
      info.innerHTML = `<div class="rev-name">${escapeHtmlSafe(rev.name)}</div><div class="rev-role">${escapeHtmlSafe(rev.role)}</div>`;
      opt.appendChild(info);
      opt.addEventListener("click", () => {
        if (opt.classList.contains("selected")) {
          opt.classList.remove("selected");
          selectedReviewers = selectedReviewers.filter(r => r.name !== rev.name);
        } else if (selectedReviewers.length < 5) {
          opt.classList.add("selected");
          selectedReviewers.push({...rev});
        }
        updateSelectedDisplay();
      });
      pool.appendChild(opt);
    });
  }

  function updateSelectedDisplay() {
    document.getElementById("selected-count").textContent = selectedReviewers.length;
    const confirmBtn = document.getElementById("btn-confirm-reviewers");
    confirmBtn.disabled = selectedReviewers.length === 0;
  }

  // ============================================================
  // FEEDBACK SUBMISSION
  // ============================================================
  function setupFeedback() {
    const input = document.getElementById("feedback-input");
    const send = document.getElementById("feedback-send");
    const status = document.getElementById("feedback-status");

    async function submitFeedback() {
      const text = input.value.trim();
      if (!text) return;
      const target = document.getElementById("feedback-target").value;
      try {
        const resp = await fetch(apiUrl("/api/autolab/feedback"), {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({ text, target }),
        });
        const data = await resp.json();
        if (data.status === "ok") {
          input.value = "";
          // Show queue count so user knows messages accumulate
          const qc = currentState && currentState.feedback_queue_count
            ? currentState.feedback_queue_count + 1
            : 1;
          showStatus(`Queued (${qc})`, status);
        } else {
          showStatus("Empty", status);
        }
      } catch (e) {
        showStatus("Queued ✓", status);  // Optimistic — don't show "Error" to user
        console.error("Feedback send issue (non-fatal):", e);
      }
    }

    send.addEventListener("click", submitFeedback);
    input.addEventListener("keydown", (e) => {
      if (e.key === "Enter") submitFeedback();
    });
  }

  function showStatus(msg, el) {
    el.textContent = msg;
    el.classList.add("show");
    setTimeout(() => el.classList.remove("show"), 2000);
  }

  // Toast notification inside editor desk overlay
  function showEditorToast(msg) {
    let toast = document.getElementById("editor-toast");
    if (!toast) {
      toast = document.createElement("div");
      toast.id = "editor-toast";
      toast.style.cssText = "position:fixed;top:20px;left:50%;transform:translateX(-50%);z-index:99999;background:#cc3333;color:#fff;font-family:'Press Start 2P',monospace;font-size:8px;padding:10px 20px;border-radius:4px;box-shadow:0 2px 10px rgba(0,0,0,0.5);transition:opacity 0.3s;pointer-events:none;";
      document.body.appendChild(toast);
    }
    toast.textContent = msg;
    toast.style.opacity = "1";
    setTimeout(() => { toast.style.opacity = "0"; }, 4000);
  }

  // ============================================================
  // INIT
  // ============================================================
  document.addEventListener("DOMContentLoaded", () => {
    animateSprites();
    startPolling();
    setupFeedback();
    setupInventory();
    setupModal();
    setupThoughtBubbles();
    setupRightTabs();
    setupEditorDesk();
  });
})();
