/* ==============================================================
   AUTONOMOUS LAB — Pixel Art Sprite System
   Ported from the monitoring UI (lab.js)
   10w x 16h pixel characters with indexed palettes
   ============================================================== */

(function () {
  "use strict";

  // ============================================================
  // PALETTES
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

  // ============================================================
  // SPRITE FRAMES (10x16)
  // ============================================================
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
  // Long-hair variants — hair flows down sides of the face
  const FRAME_LONG_HAIR_GLASSES = [
    [0,0,0,2,2,2,2,0,0,0],[0,0,2,2,2,2,2,2,0,0],
    [0,2,2,1,1,1,1,2,2,0],[0,2,1,9,1,1,9,1,2,0],
    [0,2,1,1,1,1,1,1,2,0],[0,0,0,1,8,8,1,0,0,0],
    [0,0,0,0,1,1,0,0,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,4,4,4,5,5,4,4,4,0],[0,4,4,4,5,5,4,4,4,0],
    [0,0,4,4,4,4,4,4,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,0,0,6,6,6,6,0,0,0],[0,0,0,6,0,0,6,0,0,0],
    [0,0,0,6,0,0,6,0,0,0],[0,0,7,7,0,0,7,7,0,0],
  ];
  const FRAME_LONG_HAIR = [
    [0,0,0,2,2,2,2,0,0,0],[0,0,2,2,2,2,2,2,0,0],
    [0,2,2,1,1,1,1,2,2,0],[0,2,1,3,1,1,3,1,2,0],
    [0,2,1,1,1,1,1,1,2,0],[0,0,0,1,8,8,1,0,0,0],
    [0,0,0,0,1,1,0,0,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,4,4,4,4,4,4,4,4,0],[0,4,4,4,4,4,4,4,4,0],
    [0,0,4,4,4,4,4,4,0,0],[0,0,4,4,4,4,4,4,0,0],
    [0,0,0,6,6,6,6,0,0,0],[0,0,0,6,0,0,6,0,0,0],
    [0,0,0,6,0,0,6,0,0,0],[0,0,7,7,0,0,7,7,0,0],
  ];

  // ============================================================
  // EXPERT DEFINITIONS — 29 interdisciplinary roles
  // ============================================================
  const EXPERT_DEFS = {
    reviewer:          {label:"Critical Reviewer",   glasses:true,  longHair:false, palette:{0:null,1:"#f4c7a3",2:"#5a2020",3:"#333",4:"#2a2a3a",5:"#e0e0e0",6:"#2a2a3a",7:"#3a2a20",8:"#c06060",9:"#881818"}},
    bioethicist:       {label:"Bioethicist",          glasses:true,  longHair:true,  palette:{0:null,1:"#d4a070",2:"#4a3020",3:"#333",4:"#404060",5:"#8080a0",6:"#3a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#606080"}},
    science_writer:    {label:"Science Writer",       glasses:false, longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#c06030",3:"#333",4:"#3a5a6a",5:"#5a8a9a",6:"#4a5a6a",7:"#6a5a4a",8:"#d46b6b",9:null}},
    grant_reviewer:    {label:"Grant Reviewer",       glasses:true,  longHair:false, palette:{0:null,1:"#e8b888",2:"#6a6a7a",3:"#333",4:"#3a3a4a",5:"#b0b0c0",6:"#3a3a5a",7:"#5a4a3a",8:"#c06060",9:"#5a5a7a"}},
    immunologist:      {label:"Immunologist",         glasses:false, longHair:true,  palette:{0:null,1:"#d4a878",2:"#1a1a30",3:"#333",4:"#2a8a4a",5:"#50b060",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    oncologist:        {label:"Oncologist",           glasses:false, longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#3a2a1a",3:"#333",4:"#7a2a4a",5:"#b04070",6:"#4a3a5a",7:"#6a5040",8:"#d46b6b",9:null}},
    neuroscientist:    {label:"Neuroscientist",       glasses:true,  longHair:true,  palette:{0:null,1:"#e8b888",2:"#2a2040",3:"#333",4:"#5050a0",5:"#7070c0",6:"#3a3a5a",7:"#e0d8c8",8:"#d46b6b",9:"#4040a0"}},
    geneticist:        {label:"Geneticist",           glasses:true,  longHair:false, palette:{0:null,1:"#f4c7a3",2:"#4a3020",3:"#333",4:"#2a6a8a",5:"#4a90b0",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:"#2a6080"}},
    cell_biologist:    {label:"Cell Biologist",       glasses:false, longHair:true,  palette:{0:null,1:"#d4a070",2:"#6a4a2a",3:"#333",4:"#3a7a5a",5:"#5aa07a",6:"#3a5a4a",7:"#6a5a4a",8:"#d46b6b",9:null}},
    microbiologist:    {label:"Microbiologist",       glasses:false, longHair:false, palette:{0:null,1:"#f4c7a3",2:"#8a5a2a",3:"#333",4:"#6a8a3a",5:"#90b050",6:"#4a5a3a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    pathologist:       {label:"Pathologist",          glasses:true,  longHair:true,  palette:{0:null,1:"#e8b888",2:"#3a2020",3:"#333",4:"#f0f0f0",5:"#6a3050",6:"#3a3a5a",7:"#6b5040",8:"#c06060",9:"#5a3040"}},
    pharmacologist:    {label:"Pharmacologist",       glasses:false, longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#5a4030",3:"#333",4:"#5a3a7a",5:"#8060a0",6:"#4a3a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    structural_bio:    {label:"Structural Biologist", glasses:true,  longHair:false, palette:{0:null,1:"#d4a878",2:"#2a2a2a",3:"#333",4:"#4a6a8a",5:"#6a90b0",6:"#3a4a5a",7:"#5a4a3a",8:"#d46b6b",9:"#3a6080"}},
    systems_biologist: {label:"Systems Biologist",    glasses:false, longHair:false, palette:{0:null,1:"#f4c7a3",2:"#7a4a20",3:"#333",4:"#3a5a3a",5:"#5a8a5a",6:"#4a5a4a",7:"#6a5040",8:"#d46b6b",9:null}},
    epidemiologist:    {label:"Epidemiologist",       glasses:true,  longHair:true,  palette:{0:null,1:"#e8b888",2:"#5a3a2a",3:"#333",4:"#4a5a7a",5:"#6a80a0",6:"#3a4a6a",7:"#5a4a3a",8:"#d46b6b",9:"#3a5080"}},
    statistician:      {label:"Statistician",         glasses:true,  longHair:false, palette:{0:null,1:"#f4c7a3",2:"#c8a848",3:"#333",4:"#8a7050",5:"#e0e0e0",6:"#4a4a5a",7:"#6a5040",8:"#d46b6b",9:"#886830"}},
    bioinformatician:  {label:"Bioinformatician",     glasses:false, longHair:false, palette:{0:null,1:"#e8b888",2:"#2a2a40",3:"#333",4:"#5060a0",5:"#6070b0",6:"#3a3a5a",7:"#e0d8c8",8:"#d46b6b",9:"#70a0c0"}},
    data_scientist:    {label:"Data Scientist",       glasses:false, longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#3a3a50",3:"#333",4:"#3a6a6a",5:"#5a9a9a",6:"#3a5a5a",7:"#5a4a3a",8:"#d46b6b",9:null}},
    ml_engineer:       {label:"ML Engineer",          glasses:false, longHair:false, palette:{0:null,1:"#d4a070",2:"#2a2a30",3:"#333",4:"#5a5a6a",5:"#8a8a9a",6:"#3a3a4a",7:"#e0d8c8",8:"#d46b6b",9:null}},
    comp_biologist:    {label:"Comp. Biologist",      glasses:true,  longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#6a4a2a",3:"#333",4:"#4a4a7a",5:"#6a6aa0",6:"#3a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#4a4a80"}},
    clinician:         {label:"Clinician",            glasses:false, longHair:true,  palette:{0:null,1:"#d4a070",2:"#3a2020",3:"#333",4:"#f0f0f0",5:"#4090b0",6:"#3a3a5a",7:"#6b5040",8:"#d46b6b",9:null}},
    radiologist:       {label:"Radiologist",          glasses:true,  longHair:false, palette:{0:null,1:"#f4c7a3",2:"#4a4a4a",3:"#333",4:"#2a3a4a",5:"#4a6a8a",6:"#2a3a4a",7:"#5a4a3a",8:"#d46b6b",9:"#3a5a7a"}},
    surgeon:           {label:"Surgeon",              glasses:false, longHair:false, palette:{0:null,1:"#e8b888",2:"#2a2a3a",3:"#333",4:"#2a8a7a",5:"#50b0a0",6:"#2a5a5a",7:"#f0f0f0",8:"#d46b6b",9:null}},
    chemist:           {label:"Chemist",              glasses:true,  longHair:true,  palette:{0:null,1:"#f4c7a3",2:"#9a6a2a",3:"#333",4:"#f0f0f0",5:"#c0a020",6:"#4a4a3a",7:"#5a4a3a",8:"#d46b6b",9:"#8a7020"}},
    physicist:         {label:"Physicist",            glasses:true,  longHair:false, palette:{0:null,1:"#e8b888",2:"#5a5a6a",3:"#333",4:"#2a2a3a",5:"#4a4a6a",6:"#2a2a3a",7:"#5a4a3a",8:"#d46b6b",9:"#3a3a5a"}},
    engineer:          {label:"Biomedical Engineer",  glasses:false, longHair:false, palette:{0:null,1:"#f4c7a3",2:"#5a3a20",3:"#333",4:"#8a6a20",5:"#c0a040",6:"#5a5a4a",7:"#6a5040",8:"#d46b6b",9:null}},
    psychologist:      {label:"Psychologist",         glasses:true,  longHair:true,  palette:{0:null,1:"#d4a878",2:"#8a4a3a",3:"#333",4:"#6a5a8a",5:"#8a7ab0",6:"#4a3a5a",7:"#5a4a3a",8:"#d46b6b",9:"#5a4a80"}},
    ecologist:         {label:"Ecologist",            glasses:false, longHair:true,  palette:{0:null,1:"#e8b888",2:"#7a6a3a",3:"#333",4:"#4a7a3a",5:"#6aa04a",6:"#5a6a4a",7:"#7a6a4a",8:"#d46b6b",9:null}},
    generic:           {label:"Consultant",           glasses:false, longHair:false, palette:{0:null,1:"#f4c7a3",2:"#8a6a40",3:"#333",4:"#707070",5:"#909090",6:"#4a4a4a",7:"#5a4a3a",8:"#d46b6b",9:null}},
  };

  // ============================================================
  // RENDERER
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

  function renderPI(canvas) {
    renderSprite(canvas, FRAME_GLASSES, PI_PALETTE);
  }

  function renderTrainee(canvas) {
    renderSprite(canvas, FRAME_NO_GLASSES, TRAINEE_PALETTE);
  }

  function renderExpert(canvas, type) {
    const def = EXPERT_DEFS[type] || EXPERT_DEFS.generic;
    const frame = def.longHair
      ? (def.glasses ? FRAME_LONG_HAIR_GLASSES : FRAME_LONG_HAIR)
      : (def.glasses ? FRAME_GLASSES : FRAME_NO_GLASSES);
    renderSprite(canvas, frame, def.palette);
  }

  // Export globally
  window.LabSprites = {
    renderPI,
    renderTrainee,
    renderExpert,
    EXPERT_DEFS,
    PI_PALETTE,
    TRAINEE_PALETTE,
  };
})();
