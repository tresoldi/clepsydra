/*
 * clepsydra browser interface.
 *
 * Two engines' worth of work with two different shapes:
 *   Pair — a fraction of a second, so it runs on the main thread and updates
 *          as the visitor changes a control.
 *   Tree — tens of seconds of MCMC, so it runs in a worker with a progress bar
 *          and a Cancel that terminates it.
 *
 * Everything the page shows comes from src/gd_analysis.c through
 * web/clepsydra_wasm.c — the same code the CLI runs.
 */
"use strict";

const $ = (id) => document.getElementById(id);

/* ---- state ----------------------------------------------------------- */

let engine = null;             /* main-thread module, for the pairwise model */
let worker = null;             /* the chronology sampler */
let treeRunning = false;
let lastTreeRequest = null;    /* re-sent, varied, for the sensitivity re-run */
let lastPairSvg = null;        /* kept so the figures can be downloaded */
let lastTreeSvg = null;
/* The baseline chronology, kept while the widened re-run happens, so the two
 * can be compared node by node. */
let sensitivityBase = null;
let sensitivityPending = false;
const CAL_WIDEN = 2.0;
const datasetCache = new Map();  /* id -> {tree, forms, calibrations, languages} */

/* Data the visitor supplied. It never leaves the tab: there is no server to
 * send it to, which is the point of compiling the engine to WebAssembly. */
const CUSTOM = "__custom";
let customPair = null;         /* {forms, languages} */
let customTree = null;         /* {tree, forms, calibrations, calibrated, languages} */

/* ---- wasm helpers ---------------------------------------------------- */

/*
 * Arguments go on the heap rather than through ccall's "string" type, which
 * copies onto the 256 KB wasm stack; a real forms file overflows it and fails
 * as an opaque out-of-bounds trap. Same rule as web/worker.js.
 */
function callJson(mod, name, args) {
  const ptrs = args.map((s) => mod.stringToNewUTF8(s === null ? "" : s));
  try {
    const out = mod.ccall(name, "number", ptrs.map(() => "number"), ptrs);
    const json = mod.UTF8ToString(out);
    mod.ccall("clepsydra_free", null, ["number"], [out]);
    return JSON.parse(json);
  } finally {
    for (const p of ptrs) mod._free(p);
  }
}

/*
 * The engine's own defaults, read once at load through clepsydra_defaults_json.
 * Nothing on this page restates a default that lives in C: the numbers used to
 * appear here, in index.html's value= attributes, in the CLI and in the
 * library, and they had already drifted -- chains was 4 in the library and 2
 * here. `defaults.tree` and `defaults.pair` are keyed by the same option names
 * the engine parses back.
 */
let defaults = { tree: {}, pair: {} };

function applyDefaultsToControls() {
  const set = (id, v) => { const el = $(id); if (el && v !== undefined) el.value = String(v); };
  /* Only the controls nothing else owns. tree-warmup, tree-samples and
   * tree-chains are set per dataset by onTreeDataset, because each dataset
   * brings the chain length it actually needs -- a single global default
   * would fail the page's own convergence check on the hard ones. Writing
   * the library default into them here would either be overwritten a moment
   * later or, worse, win. */
  set("tree-seed", defaults.tree.seed);
  set("pair-tlo", defaults.pair.t_lo);
  set("pair-thi", defaults.pair.t_hi);
}

/* ---- small helpers --------------------------------------------------- */

function esc(s) {
  return String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

/* Years BP with thousands separators; dates are read, not computed, from here. */
function bp(v) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  return Math.round(v).toLocaleString("en-US");
}

function showError(el, what, detail) {
  el.innerHTML = '<div class="error"><span class="what">' + esc(what) + "</span>" +
    (detail ? '<div class="detail">' + esc(detail) + "</div>" : "") + "</div>";
}

function clearError(el) { el.innerHTML = ""; }

function datasetById(id) {
  if (id === CUSTOM) {
    return {
      id: CUSTOM,
      name: "Your own data",
      note: "Your file, read in this tab and sent nowhere.",
      source: "supplied by you",
      rate_logmean: Number($("tree-logmean").value) || -8.8,
      rate_logsd: Number($("tree-logsd").value) || 0.8,
      warmup: 4000, samples: 4000, chains: 2,
    };
  }
  return DATASETS.find((d) => d.id === id);
}

/* Adds (or refreshes) the "your own data" entry and selects it. */
function selectCustom(sel) {
  if (![...sel.options].some((o) => o.value === CUSTOM)) {
    const opt = document.createElement("option");
    opt.value = CUSTOM;
    opt.textContent = "Your own data";
    sel.appendChild(opt);
  }
  sel.value = CUSTOM;
}

function readFile(input) {
  return new Promise((resolve, reject) => {
    const file = input.files && input.files[0];
    if (!file) { resolve(null); return; }
    const reader = new FileReader();
    reader.onload = () => resolve(String(reader.result));
    reader.onerror = () => reject(new Error("could not read " + file.name));
    reader.readAsText(file);
  });
}

/* A forms table with no usable language column is the common way a wrong file
 * gets picked; saying so beats an engine error about a missing column. */
function checkForms(text, statusEl) {
  const langs = languagesOf(text);
  if (langs.length < 2) {
    statusEl.textContent = "no language_id column, or fewer than two languages";
    statusEl.className = "file-status err";
    return null;
  }
  statusEl.textContent = langs.length + " languages read";
  statusEl.className = "file-status";
  return langs;
}

function optionLines(obj) {
  return Object.entries(obj)
    .filter(([, v]) => v !== null && v !== undefined && v !== "")
    .map(([k, v]) => k + "=" + v)
    .join("\n");
}

/* ---- datasets -------------------------------------------------------- */

/* Fetched on demand rather than inlined: the four bundled datasets are ~670 KB
 * together, and a visitor who only opens the Pair tab should not pay for the
 * ones they never select. */
async function loadDataset(id) {
  if (id === CUSTOM) return customTree || customPair;
  if (datasetCache.has(id)) return datasetCache.get(id);
  const ds = datasetById(id);
  const [tree, forms, calibrations] = await Promise.all([
    fetch(ds.dir + "/tree.nwk").then((r) => r.text()),
    fetch(ds.dir + "/forms.tsv").then((r) => r.text()),
    fetch(ds.dir + "/calibrations.tsv").then((r) => r.text()),
  ]);
  const entry = {
    tree, forms, calibrations,
    languages: languagesOf(forms),
    calibrated: calibratedTargets(calibrations),
  };
  datasetCache.set(id, entry);
  return entry;
}

/* The language column of a forms table, in first-seen order. Reading it here
 * keeps the pickers honest: they offer exactly what the data contains. */
function languagesOf(formsTsv) {
  const lines = formsTsv.split("\n");
  if (lines.length < 2) return [];
  const header = lines[0].split("\t");
  const col = header.indexOf("language_id");
  if (col < 0) return [];
  const seen = new Set();
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue;
    const f = lines[i].split("\t");
    if (f.length > col && f[col]) seen.add(f[col]);
  }
  return Array.from(seen).sort();
}

/* Which nodes carry a calibration. The result payload says what each age is,
 * not where it came from, so the figure reads that from the calibration file —
 * otherwise an age the model was handed would be drawn as one it worked out. */
function calibratedTargets(calTsv) {
  const out = new Set();
  const lines = calTsv.split("\n");
  if (lines.length < 2) return out;
  const col = lines[0].split("\t").indexOf("target_id");
  if (col < 0) return out;
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue;
    const f = lines[i].split("\t");
    if (f.length > col && f[col]) out.add(f[col]);
  }
  return out;
}

function fillSelect(sel, values, keep) {
  sel.innerHTML = "";
  for (const v of values) {
    const opt = document.createElement("option");
    opt.value = typeof v === "string" ? v : v.id;
    opt.textContent = typeof v === "string" ? v : v.name;
    sel.appendChild(opt);
  }
  if (keep && values.some((v) => (typeof v === "string" ? v : v.id) === keep)) {
    sel.value = keep;
  }
}

function noteFor(ds, withTiming) {
  let html = esc(ds.note);
  if (withTiming && ds.expect_s) {
    html += " <b>Roughly " + (ds.expect_s >= 60
      ? Math.round(ds.expect_s / 60) + " minute" + (ds.expect_s >= 90 ? "s" : "")
      : ds.expect_s + " seconds") + "</b> at these settings.";
  }
  return html + ' <span class="source">' + esc(ds.source) + "</span>";
}

/* ---- pair tab -------------------------------------------------------- */

async function onPairDataset() {
  const ds = datasetById($("pair-dataset").value);
  $("pair-note").innerHTML = noteFor(ds);
  const data = await loadDataset(ds.id);
  const prevA = $("pair-a").value, prevB = $("pair-b").value;
  fillSelect($("pair-a"), data.languages, prevA);
  fillSelect($("pair-b"), data.languages, prevB);
  if ($("pair-a").value === $("pair-b").value && data.languages.length > 1) {
    $("pair-b").value = data.languages[1];
  }
  runPair();
}

async function runPair() {
  const errEl = $("pair-error");
  const ds = datasetById($("pair-dataset").value);
  const data = await loadDataset(ds.id);
  const a = $("pair-a").value, b = $("pair-b").value;
  if (!a || !b) return;
  if (a === b) {
    $("pair-results").classList.remove("active");
    showError(errEl, "pick two different languages",
              "A language has diverged from itself zero years ago; there is nothing to estimate.");
    return;
  }

  const opts = optionLines({
    t_lo: Number($("pair-tlo").value) || defaults.pair.t_lo,
    t_hi: Number($("pair-thi").value) || 8000,
  });

  const r = callJson(engine, "clepsydra_date_pair_json",
                     [data.forms, a, b, "", opts]);
  if (r.ok !== true) {
    $("pair-results").classList.remove("active");
    showError(errEl, r.error || "the pair could not be dated", r.detail);
    return;
  }
  clearError(errEl);
  renderPair(r, ds);
}

function renderPair(r, ds) {
  $("pair-median").textContent = bp(r.median);
  $("pair-interval").textContent =
    "90% interval " + bp(r.et90_lo) + "–" + bp(r.et90_hi) + " BP";

  const retention = (100 * r.retention).toFixed(1);
  $("pair-summary").innerHTML = [
    ["shared concepts", r.n_shared],
    ["cognate pairs", r.n_cognate],
    ["retention", retention + "%"],
    ["mean", bp(r.mean) + " BP"],
    ["sd", bp(r.sd)],
    ["model", r.model],
  ].map(([k, v]) => "<span>" + esc(k) + " <b>" + esc(v) + "</b></span>").join("");

  /* Prior dominance is the number that says whether the data spoke at all. It
   * belongs beside the date, not buried in a diagnostics panel. */
  const dom = r.prior_dominance;
  let caveat = "";
  if (dom >= 0.5) {
    caveat = "<b>The prior is doing most of the work here</b> (prior dominance " +
      dom.toFixed(2) + "): this pair carries little information about its own " +
      "date, and the interval above is close to what you would get before " +
      "looking at the data.";
  } else {
    caveat = "Prior dominance " + dom.toFixed(2) + " — the data is moving the " +
      "estimate away from the prior.";
  }
  if (!r.absolute_anchored) {
    caveat += " No tip age was supplied, so the scale comes from the rate " +
      "prior rather than from an anchored calendar date.";
  }
  caveat += " Data: " + esc(ds.source) + ".";
  $("pair-caveat").innerHTML = caveat;

  lastPairSvg = drawPairDensity($("pair-figure"), r);
  $("pair-results").classList.add("active");
}

/* ---- your own data --------------------------------------------------- */

async function onPairFile() {
  const status = $("pair-file-status");
  try {
    const text = await readFile($("pair-file"));
    if (text === null) return;
    const langs = checkForms(text, status);
    if (langs === null) return;
    customPair = { forms: text, languages: langs, calibrated: new Set() };
    selectCustom($("pair-dataset"));
    await onPairDataset();
  } catch (e) {
    status.textContent = e.message;
    status.className = "file-status err";
  }
}

async function onTreeFiles() {
  const status = $("tree-file-status");
  try {
    const [tree, forms, cals] = await Promise.all([
      readFile($("tree-file-nwk")),
      readFile($("tree-file-forms")),
      readFile($("tree-file-cal")),
    ]);
    /* All three or none: a chronology needs a tree, data and an anchor in
     * time, and guessing a missing one would be inventing an assumption. */
    if (!tree || !forms || !cals) {
      status.textContent = "waiting for all three files";
      status.className = "file-status err";
      return;
    }
    const langs = checkForms(forms, status);
    if (langs === null) return;
    customTree = {
      tree, forms, calibrations: cals,
      languages: langs,
      calibrated: calibratedTargets(cals),
    };
    selectCustom($("tree-dataset"));
    await onTreeDataset();
    status.textContent = langs.length + " languages read; press Run";
    status.className = "file-status";
  } catch (e) {
    status.textContent = e.message;
    status.className = "file-status err";
  }
}

/* ---- decay tab ------------------------------------------------------- */

/* The default model first, then the alternatives. `cls` picks the stroke the
 * stylesheet gives each curve. */
const DECAY_LAWS = [
  { id: "", label: "P1 (constant rate, default)", cls: "law-p1" },
  { id: "swadesh", label: "Swadesh (constant rate)", cls: "law-swadesh" },
  { id: "starostin", label: "Starostin (accelerating)", cls: "law-starostin" },
  { id: "vasilyev", label: "Vasilyev (stable core)", cls: "law-vasilyev" },
  { id: "two-component", label: "Two-component", cls: "law-two-component" },
];

async function onDecayDataset() {
  const ds = datasetById($("decay-dataset").value);
  $("decay-note").innerHTML = noteFor(ds);
  const data = await loadDataset(ds.id);
  const prevA = $("decay-a").value, prevB = $("decay-b").value;
  fillSelect($("decay-a"), data.languages, prevA);
  fillSelect($("decay-b"), data.languages, prevB);
  if ($("decay-a").value === $("decay-b").value && data.languages.length > 1) {
    $("decay-b").value = data.languages[1];
  }
  runDecay();
}

async function runDecay() {
  const errEl = $("decay-error");
  const ds = datasetById($("decay-dataset").value);
  const data = await loadDataset(ds.id);
  const a = $("decay-a").value, b = $("decay-b").value;
  if (!a || !b) return;
  if (a === b) {
    $("decay-results").classList.remove("active");
    showError(errEl, "pick two different languages",
              "There is nothing to date between a language and itself.");
    return;
  }

  /* Every law over the same range, so the curves are comparable rather than
   * each drawn on its own convenient axis. */
  const opts = optionLines({ t_lo: defaults.pair.t_lo, t_hi: defaults.pair.t_hi });
  const series = [];
  for (const law of DECAY_LAWS) {
    const r = callJson(engine, "clepsydra_date_pair_json",
                       [data.forms, a, b, law.id, opts]);
    if (r.ok !== true) {
      $("decay-results").classList.remove("active");
      showError(errEl, r.error || "the comparison failed", r.detail);
      return;
    }
    series.push({ law, result: r, cls: law.cls });
  }
  clearError(errEl);
  renderDecay(series, ds);
}

function renderDecay(series, ds) {
  drawDensityOverlay($("decay-figure"), series);

  $("decay-legend").innerHTML = series.map((s) =>
    '<span><i class="legend-line ' + s.cls + '"></i> ' + esc(s.law.label) +
    "</span>").join("");

  let html = "<thead><tr><th>law</th><th>median BP</th><th>90% interval</th>" +
             "<th>width</th></tr></thead><tbody>";
  for (const s of series) {
    html += "<tr><td>" + esc(s.law.label) + "</td>" +
      "<td>" + bp(s.result.median) + "</td>" +
      "<td>" + bp(s.result.et90_lo) + "–" + bp(s.result.et90_hi) + "</td>" +
      "<td>" + bp(s.result.et90_hi - s.result.et90_lo) + "</td></tr>";
  }
  html += "</tbody>";
  $("decay-table").innerHTML = html;

  /* The spread between laws is the honest headline: it is the part of the
   * uncertainty that no single model's interval reports. */
  const medians = series.map((s) => s.result.median);
  const spread = Math.max(...medians) - Math.min(...medians);
  $("decay-caveat").innerHTML =
    "The five laws put this divergence between <b>" + bp(Math.min(...medians)) +
    "</b> and <b>" + bp(Math.max(...medians)) + "</b> years BP — a spread of " +
    bp(spread) + " years that no single model's interval reports, because each " +
    "interval is computed as if its own law were true. Picking the law that " +
    "gives the date you expected would be the least defensible way to resolve " +
    "it; comparing held-out predictions is the defensible one. Every law was " +
    "given the same 100–8,000 BP range to search; the figure is trimmed to " +
    "where the curves actually are. Data: " + esc(ds.source) + ".";

  $("decay-results").classList.add("active");
}

/* ---- validation tab -------------------------------------------------- */

/* Parses the committed held-out-calibration report. It is the output of
 * tools/validation/calibration_benchmark.py, shipped as-is rather than
 * retyped, so the page cannot drift from the result it claims. */
function parseLoco(tsv) {
  const lines = tsv.trim().split("\n");
  const header = lines[0].split("\t");
  return lines.slice(1).map((line) => {
    const f = line.split("\t");
    const row = {};
    header.forEach((h, i) => {
      const v = f[i];
      row[h] = /^[-0-9.]+$/.test(v) ? parseFloat(v) : v;
    });
    return row;
  });
}

let validationLoaded = false;

async function loadValidation() {
  if (validationLoaded) return;
  validationLoaded = true;

  const rows = parseLoco(await fetch("data/loco.tsv").then((r) => r.text()));
  drawLoco($("loco-figure"), rows);

  const byClock = {};
  for (const r of rows) {
    const c = (byClock[r.clock] = byClock[r.clock] ||
               { n: 0, hit: 0, err: 0, is: 0, isPrior: 0 });
    c.n++;
    c.hit += r.covered ? 1 : 0;
    c.err += r.abs_err;
    c.is += r.interval_score;
    c.isPrior += r.prior_interval_score;
  }
  $("loco-summary").innerHTML = Object.entries(byClock).map(([clock, c]) =>
    "<span>" + esc(clock) + " clock <b>" + c.hit + "/" + c.n + " covered</b>" +
    " (mean error " + bp(c.err / c.n) + " yr; interval score " +
    bp(c.is / c.n) + " against " + bp(c.isPrior / c.n) +
    " for the prior alone)</span>").join("");

  $("loco-caveat").innerHTML =
    "The relaxed clock covers the withheld date far more often than the strict " +
    "one, with smaller errors. The strict clock is overconfident: narrow " +
    "intervals that miss. The <b>interval score</b> (width plus 20× any miss " +
    "distance, lower is better) adds a harder test: does the posterior beat " +
    "the same run's interval with the cognates switched off? Overall, neither " +
    "clock does. The relaxed clock beats it on Sino-Tibetan and ties on " +
    "Indo-European, so the other dated tips already say most of what the " +
    "cognates can add. This is a small-N check — nine held-out points — not a " +
    "coverage guarantee. Romance-Latin is missed under both clocks for a " +
    "structural reason: freeing the family's only ancient tip leaves no " +
    "tip-age signal anywhere in the tree.";

  let html = "<thead><tr><th>clock</th><th>family</th><th>held-out tip</th>" +
             "<th>attested BP</th><th>inferred</th><th>90% interval</th>" +
             "<th>covered</th></tr></thead><tbody>";
  for (const r of rows) {
    html += "<tr><td>" + esc(r.clock) + "</td><td>" + esc(r.family) + "</td>" +
      "<td>" + esc(r.target) + "</td><td>" + bp(r.true_bp) + "</td>" +
      "<td>" + bp(r.median_bp) + "</td>" +
      "<td>" + bp(r.hpd90_lo) + "–" + bp(r.hpd90_hi) + "</td>" +
      '<td class="' + (r.covered ? "ok" : "flag") + '">' +
      (r.covered ? "yes" : "no") + "</td></tr>";
  }
  html += "</tbody>";
  $("loco-table").innerHTML = html;

  /* The simulation evidence, stated as claim + what it cannot show. Numbers
   * come from the project's own test suite, named so they can be re-run. */
  const claims = [
    ["T1 likelihood", "validity",
     "Agrees with an independent brute-force enumeration on small trees",
     "Says nothing about whether the model fits real languages"],
    ["T1 sampler", "validity",
     "Simulation-based calibration: posterior ranks are uniform",
     "Only shows the sampler recovers what the model itself generates"],
    ["Browser vs CLI", "validity",
     "Identical seeded results, checked on the Romance chronology",
     "Agreement between two builds, not evidence either is right"],
    ["Well-specified fit", "adequacy",
     "Posterior-predictive p ≈ 0.56 — the data look like the model's own",
     "A single dataset under a model chosen to be correct"],
    ["Rate-heterogeneous data", "adequacy",
     "Posterior-predictive p ≈ 0 — the misfit is detected, not hidden",
     "The canonical validity-without-adequacy failure"],
    ["Held-out real dates", "adequacy",
     "Relaxed clock 8/9 covered vs strict 3/9 across three families",
     "Nine points; attestation midpoints carry their own windows"],
  ];
  let vhtml = "<thead><tr><th>check</th><th>kind</th><th>what it shows</th>" +
              "<th>what it cannot show</th></tr></thead><tbody>";
  for (const [what, kind, shows, cannot] of claims) {
    vhtml += "<tr><td>" + esc(what) + "</td>" +
      '<td class="' + (kind === "validity" ? "ok" : "") + '">' + esc(kind) + "</td>" +
      "<td>" + esc(shows) + "</td><td>" + esc(cannot) + "</td></tr>";
  }
  vhtml += "</tbody>";
  $("validation-table").innerHTML = vhtml;
}

/* ---- guide ----------------------------------------------------------- */

function renderGuide() {
  $("guide-body").innerHTML = GUIDE_SECTIONS.map(
    (s) => "<h3>" + esc(s.title) + "</h3>" + s.html).join("");
}

function toggleGuide(open) {
  $("guide").classList.toggle("open", open);
  $("guide").setAttribute("aria-hidden", open ? "false" : "true");
}

/* ---- tree tab -------------------------------------------------------- */

function startWorker() {
  worker = new Worker("worker.js");
  worker.onmessage = onWorkerMessage;
  worker.onerror = (e) => {
    treeDone();
    showError($("tree-error"), "the engine could not start",
              e && e.message ? e.message : String(e));
  };
}

function onWorkerMessage(event) {
  const msg = event.data;
  if (msg.type === "ready") return;
  if (msg.type === "progress") {
    const chains = Number($("tree-chains").value) || 1;
    const done = (msg.chain * msg.total + msg.iter) / (chains * msg.total);
    $("tree-progress-fill").style.width = (100 * done).toFixed(1) + "%";
    $("tree-progress-label").textContent =
      "chain " + (msg.chain + 1) + " of " + chains + " — " +
      msg.iter.toLocaleString("en-US") + " of " +
      msg.total.toLocaleString("en-US") + " iterations";
    return;
  }
  if (msg.type === "result") {
    let r;
    try {
      r = JSON.parse(msg.json);
    } catch (e) {
      treeDone();
      showError($("tree-error"), "the result could not be read", String(e));
      return;
    }
    if (r.ok !== true) {
      treeDone();
      showError($("tree-error"), r.error || "the run failed", r.detail);
      return;
    }
    clearError($("tree-error"));

    if (msg.tag === "widened") {
      treeDone();
      renderSensitivity(sensitivityBase, r);
      sensitivityBase = null;
      return;
    }

    renderTree(r, msg.elapsedMs);
    if (sensitivityPending) {
      /* The same run again with every soft calibration window loosened: the
       * difference between the two is the part of the chronology the
       * calibrations are holding up. */
      sensitivityPending = false;
      sensitivityBase = r;
      $("tree-progress-label").textContent =
        "re-running with calibrations loosened ×" + CAL_WIDEN + "…";
      worker.postMessage(Object.assign({}, lastTreeRequest, {
        tag: "widened",
        opts: lastTreeRequest.opts + "\ncal_widen=" + CAL_WIDEN,
      }));
      return;
    }
    treeDone();
    return;
  }
  if (msg.type === "fatal") {
    treeDone();
    showError($("tree-error"), "the run failed", msg.message);
  }
}

function treeStarted() {
  treeRunning = true;
  $("tree-run").hidden = true;
  $("tree-cancel").hidden = false;
  $("tree-progress").classList.add("active");
  $("tree-progress-fill").style.width = "0%";
  $("tree-progress-label").textContent = "starting…";
  $("tree-results").classList.remove("active");
  clearError($("tree-error"));
}

function treeDone() {
  treeRunning = false;
  $("tree-run").hidden = false;
  $("tree-cancel").hidden = true;
  $("tree-progress").classList.remove("active");
}

/* Each dataset brings the chain length it actually needs (see datasets.js):
 * one global default would either fail its own convergence check on the hard
 * datasets or make the easy ones needlessly slow. */
async function onTreeDataset() {
  const ds = datasetById($("tree-dataset").value);
  if (ds.id === CUSTOM) {
    /* No measured setting exists for data nobody has run before, so the
     * controls are left where they are and the convergence check — not a
     * promise from this page — decides whether the run was long enough. */
    $("tree-note").innerHTML = noteFor(ds);
    return;
  }
  $("tree-note").innerHTML = noteFor(ds, true);
  $("tree-warmup").value = ds.warmup;
  $("tree-samples").value = ds.samples;
  $("tree-chains").value = ds.chains;
}

async function runTree() {
  if (treeRunning) return;
  const ds = datasetById($("tree-dataset").value);
  treeStarted();
  const data = await loadDataset(ds.id);

  if (worker === null) startWorker();
  const opts = optionLines({
    warmup: Math.max(1, Number($("tree-warmup").value) || defaults.tree.n_warmup),
    samples: Math.max(1, Number($("tree-samples").value) || defaults.tree.n_samples),
    chains: Math.max(1, Number($("tree-chains").value) || defaults.tree.n_chains),
    seed: Number($("tree-seed").value) || 0,
    rate_logmean: ds.rate_logmean,
    rate_logsd: ds.rate_logsd,
    relaxed: $("tree-relaxed").checked ? 1 : 0,
    /* Report about 200 times per chain: often enough that the bar moves,
     * rarely enough that postMessage is not the bottleneck. */
    /* Enough thinned draws for a readable trace, few enough that the payload
     * stays small. */
    trace_draws: 120,
    progress_every: Math.max(1, Math.floor(
      ((Number($("tree-warmup").value) || defaults.tree.n_warmup) +
       (Number($("tree-samples").value) || defaults.tree.n_samples)) / 200)),
  });

  sensitivityPending = $("tree-sensitivity").checked;
  sensitivityBase = null;
  $("tree-sens-panel").hidden = true;
  lastTreeRequest = {
    type: "tree",
    newick: data.tree,
    forms: data.forms,
    calibrations: data.calibrations,
    opts,
  };
  worker.postMessage(lastTreeRequest);
}

/* Baseline against loosened calibrations, node by node. */
function renderSensitivity(base, widened) {
  if (!base || !widened) return;
  const byNode = new Map();
  for (const n of widened.nodes) byNode.set(n.node, n);

  let worst = 0, worstLabel = "";
  let html = "<thead><tr><th>node</th><th>median BP</th>" +
             "<th>with windows ×" + CAL_WIDEN + "</th><th>shift</th>" +
             "<th>interval width</th><th>widened width</th></tr></thead><tbody>";
  for (const b of base.nodes) {
    const w = byNode.get(b.node);
    if (!w) continue;
    const shift = Math.abs(w.median - b.median);
    if (shift > worst) { worst = shift; worstLabel = b.label || "(node)"; }
    html += "<tr><td>" + esc(b.label || "(unlabelled)") + "</td>" +
      "<td>" + bp(b.median) + "</td><td>" + bp(w.median) + "</td>" +
      '<td class="' + (shift > 200 ? "flag" : "ok") + '">' + bp(shift) + "</td>" +
      "<td>" + bp(b.hpd_hi - b.hpd_lo) + "</td>" +
      "<td>" + bp(w.hpd_hi - w.hpd_lo) + "</td></tr>";
  }
  html += "</tbody>";
  $("tree-sens-table").innerHTML = html;
  $("tree-sens-hint").innerHTML =
    "Every soft calibration window was widened by a factor of " + CAL_WIDEN +
    " and the chronology re-estimated. The largest node-age shift is <b>" +
    bp(worst) + " years</b>" + (worstLabel ? " (" + esc(worstLabel) + ")" : "") +
    ". Fixed tips and hard bounds have no window to loosen and are unchanged.";
  $("tree-sens-panel").hidden = false;
}

/* Cancelling means terminating the worker: it is inside one synchronous call
 * into WebAssembly and cannot hear a message until that call returns. A fresh
 * worker starts immediately, so the next run does not wait for this one. */
function cancelTree() {
  if (worker !== null) worker.terminate();
  worker = null;
  treeDone();
  showError($("tree-error"), "stopped",
            "The run was cancelled, so there is no posterior to show — a " +
            "half-finished chain is not an answer.");
  startWorker();
}

function renderTree(r, elapsedMs) {
  $("tree-summary").innerHTML = [
    ["languages", r.n_languages],
    ["nodes", r.n_nodes],
    ["concepts", r.n_concepts],
    ["chains", r.n_chains + " × (" + r.n_warmup + " + " + r.n_samples + ")"],
    ["clock", r.relaxed ? "relaxed" : "strict"],
    ["took", (elapsedMs / 1000).toFixed(1) + " s"],
  ].map(([k, v]) => "<span>" + esc(k) + " <b>" + esc(v) + "</b></span>").join("");

  /* Node ages, oldest first: the reading order of a chronology. */
  const nodes = r.nodes.slice().sort((a, b2) => b2.median - a.median);
  let html = "<thead><tr><th>node</th><th>median BP</th><th>90% interval</th>" +
             "<th>mean</th><th>sd</th><th>R̂</th><th>ESS</th></tr></thead><tbody>";
  for (const n of nodes) {
    html += "<tr><td>" + esc(n.label || "(unlabelled)") + "</td>" +
      "<td>" + bp(n.median) + "</td>" +
      "<td>" + bp(n.hpd_lo) + "–" + bp(n.hpd_hi) + "</td>" +
      "<td>" + bp(n.mean) + "</td>" +
      "<td>" + bp(n.sd) + "</td>" +
      '<td class="' + (n.rhat > 1.01 ? "flag" : "ok") + '">' + n.rhat.toFixed(3) + "</td>" +
      '<td class="' + (n.ess < 200 ? "flag" : "ok") + '">' + Math.round(n.ess) + "</td></tr>";
  }
  html += "</tbody>";
  $("tree-table").innerHTML = html;

  /* Convergence, stated in words before it is stated in numbers. */
  const worstRhat = Math.max(...r.params.map((p) => p.rhat || 0));
  const worstEss = Math.min(...r.params.map((p) => (p.ess === null ? Infinity : p.ess)));
  const converged = worstRhat <= 1.01 && worstEss >= 200;
  $("tree-diag-hint").innerHTML = converged
    ? "Every parameter has R̂ ≤ 1.01 and ESS ≥ 200: the chains agree and the " +
      "sample is large enough to summarise."
    : '<b style="color:var(--warn)">This run has not converged.</b> ' +
      "Worst R̂ " + worstRhat.toFixed(3) + ", smallest ESS " +
      Math.round(worstEss) + ". Raise warm-up and samples, or add chains — " +
      "the ages above should not be quoted as they stand.";

  let dhtml = "<thead><tr><th>parameter</th><th>type</th><th>mean</th><th>sd</th>" +
              "<th>R̂</th><th>ESS</th><th>acceptance</th></tr></thead><tbody>";
  for (const p of r.params) {
    const name = p.type === "age" && p.label ? "age:" + p.label : p.type;
    dhtml += "<tr><td>" + esc(name) + "</td><td>" + esc(p.type) + "</td>" +
      "<td>" + fmt(p.mean) + "</td><td>" + fmt(p.sd) + "</td>" +
      '<td class="' + (p.rhat > 1.01 ? "flag" : "ok") + '">' + fmt(p.rhat, 3) + "</td>" +
      '<td class="' + (p.ess < 200 ? "flag" : "ok") + '">' + fmt(p.ess, 0) + "</td>" +
      "<td>" + fmt(p.acceptance, 3) + "</td></tr>";
  }
  dhtml += "</tbody>";
  $("tree-diagnostics").innerHTML = dhtml;

  const ds = datasetById($("tree-dataset").value);
  const data = datasetCache.get(ds.id);
  lastTreeSvg = drawDatedTree($("tree-figure"), r,
                              data ? data.calibrated : new Set());
  drawDiagnostics($("tree-diag-figure"), r);

  let caveat = "These ages depend on the rate prior (log-mean " +
    ds.rate_logmean + ", log-sd " + ds.rate_logsd + ") and on the calibrations " +
    "supplied with the dataset; they are not independent measurements of time. ";
  if (r.branch_rate_confounded) {
    caveat += "<b>Branch rates and node ages are confounded in this run</b> " +
      "(most negative correlation " + fmt(r.min_branch_duration_corr, 2) +
      "): a faster branch and an older node explain the data equally well, so " +
      "the ages are less identified than their intervals suggest. ";
  }
  if (r.n_skipped_rows > 0) {
    caveat += r.n_skipped_rows + " form rows named a language that is not a " +
      "tip of this tree and were skipped. ";
  }
  caveat += "Data: " + esc(ds.source) + ".";
  $("tree-caveat").innerHTML = caveat;

  $("tree-results").classList.add("active");
}

function fmt(v, places) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  if (places === undefined) {
    return Math.abs(v) >= 1000 || Math.abs(v) < 0.01
      ? v.toExponential(3) : v.toFixed(4);
  }
  return v.toFixed(places);
}

/* ---- init ------------------------------------------------------------ */

(async function init() {
  fillSelect($("pair-dataset"), DATASETS);
  fillSelect($("tree-dataset"), DATASETS);
  fillSelect($("decay-dataset"), DATASETS);
  renderGuide();

  for (const btn of document.querySelectorAll(".tab")) {
    btn.addEventListener("click", () => {
      document.querySelector(".tab.on").classList.remove("on");
      btn.classList.add("on");
      document.querySelector(".tab-content:not([hidden])").hidden = true;
      $("tab-" + btn.dataset.tab).hidden = false;
      /* The validation evidence is static and a little heavy; it is fetched
       * the first time someone looks at it, not on page load. */
      if (btn.dataset.tab === "validation") loadValidation();
    });
  }

  try {
    engine = await createClepsydra();
  } catch (e) {
    showError($("pair-error"), "the engine could not load",
              e && e.message ? e.message : String(e));
    return;
  }

  try {
    const d = callJson(engine, "clepsydra_defaults_json", []);
    if (d && d.ok) { defaults = { tree: d.tree || {}, pair: d.pair || {} }; }
  } catch (e) {
    /* Leave the markup's values in place; they are a reasonable last resort. */
  }
  applyDefaultsToControls();

  const version = engine.ccall("clepsydra_version", "string", [], []);
  $("version").textContent = "v" + version;
  fetch("BUILD_INFO", { cache: "no-cache" })
    .then((r) => (r.ok ? r.text() : ""))
    .then((text) => {
      const m = text.match(/^git_commit\s+(\S+)/m);
      if (m) $("version").textContent += " (" + m[1] + ")";
    })
    .catch(() => {});

  $("pair-dataset").addEventListener("change", onPairDataset);
  $("pair-a").addEventListener("change", runPair);
  $("pair-b").addEventListener("change", runPair);
  $("pair-tlo").addEventListener("change", runPair);
  $("pair-thi").addEventListener("change", runPair);

  $("tree-dataset").addEventListener("change", onTreeDataset);
  $("tree-run").addEventListener("click", runTree);
  $("tree-cancel").addEventListener("click", cancelTree);
  $("decay-dataset").addEventListener("change", onDecayDataset);
  $("decay-a").addEventListener("change", runDecay);
  $("decay-b").addEventListener("change", runDecay);
  $("pair-file").addEventListener("change", onPairFile);
  for (const id of ["tree-file-nwk", "tree-file-forms", "tree-file-cal"]) {
    $(id).addEventListener("change", onTreeFiles);
  }
  $("guide-open").addEventListener("click", () => toggleGuide(true));
  $("guide-close").addEventListener("click", () => toggleGuide(false));
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") toggleGuide(false);
  });
  $("pair-download").addEventListener("click",
    () => downloadSvg(lastPairSvg, "clepsydra-divergence-posterior.svg"));
  $("tree-download").addEventListener("click",
    () => downloadSvg(lastTreeSvg, "clepsydra-dated-tree.svg"));

  await onPairDataset();
  await onTreeDataset();
  await onDecayDataset();
  startWorker();
})();
