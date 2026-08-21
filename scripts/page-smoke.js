/*
 * page-smoke.js — drives the published page in a real browser and checks that
 * it does what it claims.
 *
 *   node scripts/page-smoke.js
 *
 * The page is the one part of this project that no C test can reach: a working
 * engine plus a broken worker, a mis-wired control or a caveat that never
 * renders all look identical from the command line. This starts a static server
 * and headless Chrome, then asserts the things that would embarrass us if they
 * broke — the pairwise date matches what the CLI computes, a chronology runs to
 * completion off the main thread, an under-converged run says so, and Cancel
 * really stops it.
 *
 * Needs: node (>= 18, for global fetch and WebSocket), python3, and one of
 * google-chrome / chromium. Skips with exit code 0 if the browser is absent, so
 * it can run anywhere without becoming a false alarm.
 */
"use strict";

const { spawn, spawnSync } = require("child_process");
const path = require("path");
const fs = require("fs");
const os = require("os");

const repo = path.resolve(__dirname, "..");
const PORT = 8731 + (process.pid % 200);
const DEBUG_PORT = 9222 + (process.pid % 200);

function findBrowser() {
  for (const b of ["google-chrome", "chromium", "chromium-browser"]) {
    const r = spawnSync("command", ["-v", b], { shell: true, encoding: "utf8" });
    if (r.status === 0 && r.stdout.trim()) return r.stdout.trim();
  }
  return null;
}

/* The expected pairwise date comes from the CLI, not from a number pasted in
 * here: the point of the check is that the two agree, so hard-coding one would
 * test nothing. */
function cliPairMedian(a, b) {
  const out = fs.mkdtempSync(path.join(os.tmpdir(), "clepsydra-smoke-"));
  const r = spawnSync(path.join(repo, "build", "clepsydra"),
    ["date-pair", "--forms", path.join(repo, "web/data/romance/forms.tsv"),
     "--lang-a", a, "--lang-b", b, "--out", out, "--t-lo", "100", "--t-hi", "8000"],
    { encoding: "utf8" });
  if (r.status !== 0) throw new Error("the CLI failed: " + (r.stderr || ""));
  const tsv = fs.readFileSync(path.join(out, "posterior_pair.tsv"), "utf8");
  fs.rmSync(out, { recursive: true, force: true });
  const m = tsv.match(/^divergence_median_bp\t(\S+)/m);
  return Math.round(parseFloat(m[1]));
}

async function connect(browserPath) {
  const chrome = spawn(browserPath, [
    "--headless", "--disable-gpu", "--no-sandbox",
    "--remote-debugging-port=" + DEBUG_PORT,
    "--user-data-dir=" + fs.mkdtempSync(path.join(os.tmpdir(), "clepsydra-chrome-")),
    "about:blank",
  ], { stdio: "ignore" });

  for (let i = 0; i < 60; i++) {
    await new Promise((r) => setTimeout(r, 250));
    try {
      const res = await fetch(`http://127.0.0.1:${DEBUG_PORT}/json/version`);
      if (res.ok) return chrome;
    } catch { /* not up yet */ }
  }
  throw new Error("the browser never opened its debugging port");
}

async function openPage(url) {
  const target = await (await fetch(
    `http://127.0.0.1:${DEBUG_PORT}/json/new?${encodeURIComponent(url)}`,
    { method: "PUT" })).json();
  const ws = new WebSocket(target.webSocketDebuggerUrl);
  const waiting = new Map();
  let id = 0;
  ws.onmessage = (ev) => {
    const m = JSON.parse(ev.data);
    if (m.id && waiting.has(m.id)) { waiting.get(m.id)(m); waiting.delete(m.id); }
  };
  await new Promise((res) => { ws.onopen = res; });

  const send = (method, params) => new Promise((res) => {
    const myId = ++id;
    waiting.set(myId, res);
    ws.send(JSON.stringify({ id: myId, method, params }));
  });

  return {
    async eval(expression) {
      const r = await send("Runtime.evaluate",
                           { expression, returnByValue: true, awaitPromise: true });
      if (r.result && r.result.exceptionDetails) {
        throw new Error(r.result.exceptionDetails.text);
      }
      return r.result && r.result.result ? r.result.result.value : undefined;
    },
    /* Polls until the expression stops returning `pending`. Everything here is
     * real work in wasm, so virtual time cannot help; this waits for it. */
    async until(expression, pending, timeoutMs) {
      const started = Date.now();
      let value;
      while (Date.now() - started < timeoutMs) {
        await new Promise((r) => setTimeout(r, 400));
        try { value = await this.eval(expression); } catch { /* still loading */ }
        if (value !== undefined && value !== null && String(value) !== pending) return value;
      }
      throw new Error("timed out waiting for: " + expression);
    },
    close() { ws.close(); },
  };
}

const checks = [];
function check(name, ok, detail) {
  checks.push({ name, ok, detail });
  console.log((ok ? "  ok   " : "  FAIL ") + name + (detail ? "  — " + detail : ""));
}

async function main() {
  const browser = findBrowser();
  if (!browser) {
    console.log("page-smoke: no chrome/chromium found; skipping");
    return 0;
  }
  if (!fs.existsSync(path.join(repo, "web", "clepsydra.wasm"))) {
    console.log("page-smoke: web/clepsydra.wasm is missing; run web/build-wasm.sh");
    return 2;
  }
  if (!fs.existsSync(path.join(repo, "build", "clepsydra"))) {
    console.log("page-smoke: build/clepsydra is missing; run make native");
    return 2;
  }

  const server = spawn("python3", ["-m", "http.server", String(PORT)],
                       { cwd: path.join(repo, "web"), stdio: "ignore" });
  let chrome = null, page = null;
  try {
    await new Promise((r) => setTimeout(r, 700));
    chrome = await connect(browser);
    page = await openPage(`http://localhost:${PORT}/index.html`);

    /* ---- the page loads and the engine starts ---- */
    const version = await page.until(
      "document.getElementById('version').textContent", "", 60000);
    check("the engine loads and reports its version", /^v\d/.test(version), version);

    /* ---- the pairwise date matches the CLI ---- */
    const shown = await page.until(
      "document.getElementById('pair-median').textContent", "—", 60000);
    const expected = cliPairMedian("Catalan", "French");
    check("the pairwise date matches the CLI",
          parseInt(shown.replace(/,/g, ""), 10) === expected,
          "page " + shown + " vs CLI " + expected);

    /* ---- the posterior is drawn, not just summarised ---- */
    const density = await page.eval(`JSON.stringify({
      curve: document.querySelectorAll('#pair-figure .curve').length,
      band: document.querySelectorAll('#pair-figure .band').length,
      median: document.querySelectorAll('#pair-figure .median').length,
    })`);
    const dfig = JSON.parse(density);
    check("the pairwise posterior is drawn with its interval",
          dfig.curve === 1 && dfig.band === 1 && dfig.median === 1,
          JSON.stringify(dfig));

    /* ---- a chronology runs in the worker, and an unconverged run says so ---- */
    const treeExpr = `(() => {
      if (!window.__smokeRun) {
        document.querySelector('[data-tab=tree]').click();
        document.getElementById('tree-warmup').value = 300;
        document.getElementById('tree-samples').value = 300;
        document.getElementById('tree-run').click();
        window.__smokeRun = 1;
      }
      if (!document.getElementById('tree-results').classList.contains('active')) return 'running';
      return JSON.stringify({
        rows: document.querySelectorAll('#tree-table tbody tr').length,
        diag: document.getElementById('tree-diag-hint').textContent,
        caveat: document.getElementById('tree-caveat').textContent,
      });
    })()`;
    const treeOut = JSON.parse(await page.until(treeExpr, "running", 180000));
    check("a chronology runs off the main thread", treeOut.rows > 0,
          treeOut.rows + " node ages");
    check("a deliberately short run is reported as unconverged",
          /has not converged/.test(treeOut.diag),
          treeOut.diag.slice(0, 60));
    check("the result carries its rate-prior caveat",
          /rate prior/.test(treeOut.caveat));

    /* ---- the chronology is drawn, on a real time axis ----
     * Latin is calibrated to 2050 BP and the living lects to 0, so on a
     * calendar-time axis Latin's tip must sit measurably left of theirs. A
     * cladogram that ignored the ages would line them all up, and every other
     * check here would still pass. */
    const treeFig = JSON.parse(await page.eval(`(() => {
      const f = document.getElementById('tree-figure');
      const tips = [...f.querySelectorAll('.tip')].map(t => ({
        label: t.textContent, x: parseFloat(t.getAttribute('x')) }));
      const latin = tips.find(t => t.label === 'Latin');
      const spanish = tips.find(t => t.label === 'Spanish');
      return JSON.stringify({
        branches: f.querySelectorAll('.branch').length,
        hpd: f.querySelectorAll('.hpd').length,
        calibrated: f.querySelectorAll('.node-cal').length,
        estimated: f.querySelectorAll('.node-est').length,
        latinLeftOfSpanish: latin && spanish ? latin.x < spanish.x - 20 : null,
      });
    })()`));
    check("the dated tree is drawn with interval bars",
          treeFig.branches > 0 && treeFig.hpd > 0, JSON.stringify(treeFig));
    check("ancient tips sit at their own ages on the time axis",
          treeFig.latinLeftOfSpanish === true);
    check("calibrated and estimated ages are drawn differently",
          treeFig.calibrated > 0 && treeFig.estimated > 0,
          treeFig.calibrated + " given, " + treeFig.estimated + " inferred");

    /* ---- the diagnostics figure, including a per-chain trace ---- */
    const diag = JSON.parse(await page.eval(`JSON.stringify({
      bars: document.querySelectorAll('#tree-diag-figure .track').length,
      thresholds: document.querySelectorAll('#tree-diag-figure .threshold').length,
      traceA: document.querySelectorAll('#tree-diag-figure .trace-a').length,
      traceB: document.querySelectorAll('#tree-diag-figure .trace-b').length,
      /* Every drawn trace names the parameter it belongs to, and that name
       * must be one of the rows. The trace used to be a parallel array paired
       * with params by index -- a reordering would have mislabelled all of
       * them silently, and nothing here could have seen it. */
      tracesNamed: (() => {
        const labels = new Set(
          Array.from(document.querySelectorAll('#tree-diag-figure .row-label'))
               .map((n) => n.textContent));
        const paths = Array.from(
          document.querySelectorAll('#tree-diag-figure .trace-a, #tree-diag-figure .trace-b'));
        if (!paths.length) return 'no traces';
        const bad = paths.filter((p) => !labels.has(p.getAttribute('data-param')));
        return bad.length === 0 ? 'all named' : (bad.length + ' unmatched');
      })(),
    })`));
    check("convergence is drawn against its thresholds, with a trace per chain",
          diag.bars > 0 && diag.thresholds > 0 && diag.traceA > 0 && diag.traceB > 0,
          JSON.stringify(diag));
    check("every trace names the parameter it belongs to",
          diag.tracesNamed === "all named", diag.tracesNamed);

    /* ---- the shipped defaults converge ----
     * datasets.js claims each dataset's default chain length reaches R̂ <= 1.01
     * and ESS >= 200. That claim is exactly the kind that rots silently when a
     * prior or a move changes, so it is checked rather than asserted. */
    const defaultExpr = `(() => {
      if (!window.__smokeDefault) {
        document.getElementById('tree-dataset').dispatchEvent(new Event('change'));
        window.__smokeDefault = 1;
        setTimeout(() => document.getElementById('tree-run').click(), 200);
        return 'running';
      }
      if (!document.getElementById('tree-results').classList.contains('active')) return 'running';
      return document.getElementById('tree-diag-hint').textContent;
    })()`;
    await page.eval("document.getElementById('tree-results').classList.remove('active')");
    const defaultDiag = await page.until(defaultExpr, "running", 240000);
    check("the dataset's default settings converge",
          !/has not converged/.test(defaultDiag), defaultDiag.slice(0, 70));

    /* ---- the decay tab compares laws rather than picking one ---- */
    const decay = JSON.parse(await page.until(`(() => {
      document.querySelector('[data-tab=decay]').click();
      if (!document.getElementById('decay-results').classList.contains('active')) return 'running';
      const medians = [...document.querySelectorAll('#decay-table tbody tr')]
        .map(tr => parseInt(tr.children[1].textContent.replace(/,/g, ''), 10));
      return JSON.stringify({
        curves: document.querySelectorAll('#decay-figure .curve').length,
        rows: medians.length,
        distinct: new Set(medians).size,
        caveat: document.getElementById('decay-caveat').textContent,
      });
    })()`, "running", 60000));
    check("all five decay laws are drawn together",
          decay.curves === 5 && decay.rows === 5,
          decay.curves + " curves, " + decay.rows + " rows");
    /* If the laws all agreed there would be nothing to compare, and a silent
     * fallback to one law would look exactly like agreement. */
    check("the laws give visibly different dates", decay.distinct >= 4,
          decay.distinct + " distinct medians of " + decay.rows);
    check("the spread between laws is stated as uncertainty",
          /spread of/.test(decay.caveat));

    /* ---- the validation tab shows the held-out result, both clocks ---- */
    const valid = JSON.parse(await page.until(`(() => {
      document.querySelector('[data-tab=validation]').click();
      if (document.querySelectorAll('#loco-table tbody tr').length === 0) return 'running';
      return JSON.stringify({
        rows: document.querySelectorAll('#loco-table tbody tr').length,
        strict: document.querySelectorAll('#loco-figure .point.strict').length,
        relaxed: document.querySelectorAll('#loco-figure .point.relaxed').length,
        summary: document.getElementById('loco-summary').textContent,
        claims: document.querySelectorAll('#validation-table tbody tr').length,
      });
    })()`, "running", 60000));
    check("the held-out benchmark is drawn for both clocks",
          valid.strict === 9 && valid.relaxed === 9,
          valid.strict + " strict, " + valid.relaxed + " relaxed");
    check("its coverage is reported per clock",
          /3\/9 covered/.test(valid.summary) && /8\/9 covered/.test(valid.summary),
          valid.summary.trim().slice(0, 70));
    check("validity and adequacy claims are listed apart",
          valid.claims >= 6, valid.claims + " claims");

    /* ---- calibration sensitivity runs a second, varied chronology ----
     * The claim is that loosening the calibration windows moves the dates; if
     * the second run were silently identical to the first, the panel would look
     * fine and mean nothing. */
    const sens = JSON.parse(await page.until(`(() => {
      if (!window.__smokeSens) {
        window.__smokeSens = 1;
        document.querySelector('[data-tab=tree]').click();
        document.getElementById('tree-warmup').value = 400;
        document.getElementById('tree-samples').value = 400;
        document.getElementById('tree-sensitivity').checked = true;
        document.getElementById('tree-run').click();
        return 'running';
      }
      if (document.getElementById('tree-sens-panel').hidden) return 'running';
      const rows = [...document.querySelectorAll('#tree-sens-table tbody tr')];
      return JSON.stringify({
        rows: rows.length,
        hint: document.getElementById('tree-sens-hint').textContent,
        anyShift: rows.some(tr =>
          parseInt(tr.children[3].textContent.replace(/,/g, ''), 10) > 0),
      });
    })()`, "running", 240000));
    check("loosening the calibrations re-dates the tree",
          sens.rows > 0 && sens.anyShift === true,
          sens.rows + " nodes compared");
    check("the largest calibration-driven shift is reported",
          /largest node-age shift/.test(sens.hint));
    await page.eval("document.getElementById('tree-sensitivity').checked = false");

    /* ---- a visitor's own data is accepted and used ----
     * Driven through the DataTransfer API so it exercises the real file input
     * and the real change handler, not a shortcut into the state. */
    const custom = JSON.parse(await page.until(`(async () => {
      if (!window.__smokeFile) {
        window.__smokeFile = 1;
        document.querySelector('[data-tab=pair]').click();
        const tsv = [
          'language_id\\tconcept_id\\tcognate_id',
          'Alpha\\tc1\\tk1', 'Beta\\tc1\\tk1',
          'Alpha\\tc2\\tk2', 'Beta\\tc2\\tk3',
          'Alpha\\tc3\\tk4', 'Beta\\tc3\\tk4',
        ].join('\\n');
        const dt = new DataTransfer();
        dt.items.add(new File([tsv], 'mine.tsv', { type: 'text/tab-separated-values' }));
        const input = document.getElementById('pair-file');
        input.files = dt.files;
        input.dispatchEvent(new Event('change'));
        return 'running';
      }
      const sel = document.getElementById('pair-dataset');
      if (sel.value !== '__custom') return 'running';
      const langs = [...document.getElementById('pair-a').options].map(o => o.value);
      return JSON.stringify({
        status: document.getElementById('pair-file-status').textContent,
        languages: langs,
        median: document.getElementById('pair-median').textContent,
      });
    })()`, "running", 60000));
    check("a supplied forms file is read and dated",
          custom.languages.join(",") === "Alpha,Beta" && custom.median !== "—",
          custom.status + " → " + custom.median + " BP");

    /* ---- the guide opens ---- */
    const guide = JSON.parse(await page.eval(`(() => {
      document.getElementById('guide-open').click();
      return JSON.stringify({
        open: document.getElementById('guide').classList.contains('open'),
        sections: document.querySelectorAll('#guide-body h3').length,
      });
    })()`));
    check("the guide opens with its sections", guide.open && guide.sections >= 8,
          guide.sections + " sections");
    await page.eval("document.getElementById('guide-close').click()");

    /* ---- cancel really stops it ---- */
    const cancelExpr = `(() => {
      if (!window.__smokeCancel) {
        document.getElementById('tree-warmup').value = 20000;
        document.getElementById('tree-samples').value = 20000;
        document.getElementById('tree-run').click();
        window.__smokeCancel = 1;
        setTimeout(() => document.getElementById('tree-cancel').click(), 2000);
        return 'running';
      }
      const e = document.getElementById('tree-error').textContent;
      if (!e) return 'running';
      return JSON.stringify({ text: e, idle: !document.getElementById('tree-run').hidden });
    })()`;
    const cancelOut = JSON.parse(await page.until(cancelExpr, "running", 90000));
    check("cancel stops the run and offers no partial answer",
          /stopped/.test(cancelOut.text) && cancelOut.idle,
          cancelOut.text.slice(0, 50));
  } finally {
    if (page) page.close();
    if (chrome) chrome.kill();
    server.kill();
  }

  const failed = checks.filter((c) => !c.ok).length;
  console.log(`${checks.length - failed}/${checks.length} page checks passed`);
  return failed === 0 ? 0 : 1;
}

main().then((code) => process.exit(code)).catch((e) => {
  console.error("page-smoke: " + (e && e.stack ? e.stack : e));
  process.exit(1);
});
