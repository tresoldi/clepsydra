/*
 * viz.js — the three figures the page draws.
 *
 * Hand-rolled SVG, no libraries: the page has no runtime dependencies, works
 * offline, and cannot be broken by a CDN. Colours come from the stylesheet's
 * custom properties through CSS classes, so every figure follows the visitor's
 * light or dark theme without being redrawn.
 *
 *   drawDatedTree   — a chronology on a calendar-time axis, with an interval
 *                     bar on every estimated node. The centrepiece.
 *   drawPairDensity — the whole posterior over a divergence time.
 *   drawDiagnostics — R-hat and ESS against their thresholds, plus a trace.
 *
 * Every figure draws uncertainty as a first-class object rather than an
 * annotation: the width of a bar is the result, not decoration on it.
 */
"use strict";

const SVG_NS = "http://www.w3.org/2000/svg";

function el(name, attrs, parent) {
  const node = document.createElementNS(SVG_NS, name);
  for (const k in attrs) {
    if (attrs[k] !== null && attrs[k] !== undefined) {
      node.setAttribute(k, String(attrs[k]));
    }
  }
  if (parent) parent.appendChild(node);
  return node;
}

function svgRoot(container, width, height, cls) {
  container.innerHTML = "";
  const svg = el("svg", {
    class: cls,
    width: width,
    height: height,
    viewBox: "0 0 " + width + " " + height,
    xmlns: SVG_NS,
  }, container);
  return svg;
}

/* Round tick steps: 1, 2, 5 x 10^k. A time axis labelled 1,732-year steps
 * would be arithmetically fine and unreadable. */
function niceStep(span, targetTicks) {
  const raw = span / Math.max(1, targetTicks);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) {
    if (m * mag >= raw) return m * mag;
  }
  return 10 * mag;
}

function ticksTo(maxValue, targetTicks) {
  const step = niceStep(maxValue, targetTicks);
  const out = [];
  for (let v = 0; v <= maxValue + 1e-9; v += step) out.push(v);
  return out;
}

function thousands(v) {
  return Math.round(v).toLocaleString("en-US");
}

/* ------------------------------------------------------------ dated tree ---- */

/*
 * A time-scaled tree: the x axis is calendar time in years before present, the
 * present at the right, so a node's horizontal position IS its age. Tips sit at
 * their own ages — an ancient language like Latin stands where its attestation
 * puts it, not level with the living languages — which is what makes this a
 * chronology rather than a cladogram.
 *
 *   result      the date-tree payload from the engine
 *   calibrated  Set of node labels that carry a calibration, so an age the
 *               model was told apart from one it inferred
 */
function drawDatedTree(container, result, calibrated) {
  const byId = new Map();
  for (const n of result.tree.nodes) byId.set(n.id, n);

  /* Posterior summaries for the sampled nodes; fixed nodes keep the age the
   * calibration pinned them to. */
  const post = new Map();
  for (const n of result.nodes) post.set(n.node, n);

  const ageOf = (id) => (post.has(id) ? post.get(id).median : byId.get(id).age);

  /* Tip order from a depth-first walk, so sister groups stay adjacent. */
  const children = new Map();
  for (const n of result.tree.nodes) children.set(n.id, []);
  for (const n of result.tree.nodes) {
    if (n.parent !== null && n.parent !== undefined) children.get(n.parent).push(n.id);
  }
  const tips = [];
  (function walk(id) {
    const kids = children.get(id);
    if (kids.length === 0) { tips.push(id); return; }
    for (const k of kids) walk(k);
  })(result.tree.root);

  const y = new Map();
  tips.forEach((id, i) => y.set(id, i));
  (function assign(id) {
    const kids = children.get(id);
    if (kids.length === 0) return y.get(id);
    const ys = kids.map(assign);
    const mid = (Math.min(...ys) + Math.max(...ys)) / 2;
    y.set(id, mid);
    return mid;
  })(result.tree.root);

  /* The axis has to hold the oldest interval bound, not just the oldest
   * estimate: an interval clipped by the frame would understate the doubt. */
  let maxAge = 0;
  for (const n of result.tree.nodes) maxAge = Math.max(maxAge, ageOf(n.id));
  for (const n of result.nodes) maxAge = Math.max(maxAge, n.hpd_hi);
  if (!(maxAge > 0)) maxAge = 1;

  const labelWidth = Math.min(190,
    14 + 7.2 * Math.max(...tips.map((id) => (byId.get(id).label || "").length)));
  const margin = { top: 16, right: labelWidth + 12, bottom: 44, left: 14 };
  const rowHeight = 26;
  const plotWidth = Math.max(280,
    Math.min(760, (container.clientWidth || 900) - margin.left - margin.right));
  const plotHeight = tips.length * rowHeight;
  const width = margin.left + plotWidth + margin.right;
  const height = margin.top + plotHeight + margin.bottom;

  /* Present at the right: reading right-to-left is reading back in time. */
  const X = (age) => margin.left + (maxAge - age) / maxAge * plotWidth;
  const Y = (row) => margin.top + row * rowHeight + rowHeight / 2;

  const svg = svgRoot(container, width, height, "fig fig-tree");

  /* Axis first, so the tree is drawn over it. */
  const axisY = margin.top + plotHeight + 10;
  el("line", { class: "axis", x1: margin.left, x2: margin.left + plotWidth,
               y1: axisY, y2: axisY }, svg);
  for (const t of ticksTo(maxAge, 6)) {
    const x = X(t);
    if (x < margin.left - 1 || x > margin.left + plotWidth + 1) continue;
    el("line", { class: "grid", x1: x, x2: x, y1: margin.top, y2: axisY }, svg);
    el("line", { class: "axis", x1: x, x2: x, y1: axisY, y2: axisY + 4 }, svg);
    el("text", { class: "tick", x: x, y: axisY + 16, "text-anchor": "middle" }, svg)
      .textContent = thousands(t);
  }
  el("text", { class: "axis-label", x: margin.left + plotWidth / 2,
               y: axisY + 33, "text-anchor": "middle" }, svg)
    .textContent = "years before present";

  /* Interval bars under the branches: the estimate's width is the point. */
  for (const n of result.nodes) {
    const row = y.get(n.node);
    if (row === undefined) continue;
    const x1 = X(n.hpd_hi), x2 = X(n.hpd_lo);
    el("rect", { class: "hpd", x: Math.min(x1, x2), y: Y(row) - 6,
                 width: Math.max(2, Math.abs(x2 - x1)), height: 12,
                 rx: 2 }, svg);
  }

  /* Branches, drawn as a rectangular cladogram in time. */
  (function branches(id) {
    const kids = children.get(id);
    if (kids.length === 0) return;
    const px = X(ageOf(id));
    const kidYs = kids.map((k) => Y(y.get(k)));
    el("line", { class: "branch", x1: px, x2: px,
                 y1: Math.min(...kidYs), y2: Math.max(...kidYs) }, svg);
    for (const k of kids) {
      el("line", { class: "branch", x1: px, x2: X(ageOf(k)),
                   y1: Y(y.get(k)), y2: Y(y.get(k)) }, svg);
      branches(k);
    }
  })(result.tree.root);

  /* Nodes. A calibrated age is drawn as a square and an inferred one as a
   * circle, so what the model was told never looks like what it worked out. */
  for (const n of result.tree.nodes) {
    const row = y.get(n.id);
    const cx = X(ageOf(n.id)), cy = Y(row);
    const isCal = calibrated && calibrated.has(n.label);
    if (isCal) {
      el("rect", { class: "node-cal", x: cx - 3.5, y: cy - 3.5, width: 7, height: 7 }, svg);
    } else if (post.has(n.id)) {
      el("circle", { class: "node-est", cx: cx, cy: cy, r: 3 }, svg);
    }
    if (n.is_leaf) {
      /* Beside its own tip, not at the right margin: an ancient language sits
       * where its age puts it, and a label parked at the present would undo
       * exactly the thing this axis exists to show. */
      el("text", { class: "tip", x: cx + 7, y: cy + 4 }, svg)
        .textContent = n.label || "";
    } else if (n.label && post.has(n.id)) {
      /* The age above the node, on a panel of background so it stays legible
       * where it crosses a branch or an interval bar. */
      const text = n.label + " " + thousands(post.get(n.id).median);
      const w = text.length * 5.3 + 6;
      el("rect", { class: "label-bg", x: cx - w / 2, y: cy - 18,
                   width: w, height: 12, rx: 2 }, svg);
      el("text", { class: "node-label", x: cx, y: cy - 9,
                   "text-anchor": "middle" }, svg).textContent = text;
    }
  }

  return svg;
}

/* --------------------------------------------------------- pair density ---- */

/*
 * The posterior over a divergence time, drawn whole. The page reports a median
 * because people ask for one, but the median is a summary of this — so this is
 * what gets the space.
 */
function drawPairDensity(container, result) {
  const t = result.grid.t, d = result.grid.density;
  if (!t || t.length < 2) { container.innerHTML = ""; return null; }

  const margin = { top: 14, right: 16, bottom: 44, left: 16 };
  const width = Math.max(320, Math.min(760, container.clientWidth || 720));
  const height = 210;
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  const tMax = t[t.length - 1], tMin = t[0];
  const dMax = Math.max(...d);
  const X = (v) => margin.left + (v - tMin) / (tMax - tMin) * plotWidth;
  const Y = (v) => margin.top + plotHeight - (v / dMax) * plotHeight;

  const svg = svgRoot(container, width, height, "fig fig-density");

  const axisY = margin.top + plotHeight;
  for (const tick of ticksTo(tMax, 6)) {
    if (tick < tMin) continue;
    const x = X(tick);
    el("line", { class: "grid", x1: x, x2: x, y1: margin.top, y2: axisY }, svg);
    el("text", { class: "tick", x: x, y: axisY + 16, "text-anchor": "middle" }, svg)
      .textContent = thousands(tick);
  }

  /* The central 90% as a filled band, the rest of the curve outlined: the
   * interval is part of the shape, not a caption under it. */
  const inBand = (v) => v >= result.et90_lo && v <= result.et90_hi;
  let band = "";
  for (let i = 0; i < t.length; i++) {
    if (!inBand(t[i])) continue;
    band += (band === "" ? "M" + X(t[i]) + "," + axisY + " L" : " L") +
            X(t[i]) + "," + Y(d[i]);
  }
  if (band !== "") {
    const last = t.filter(inBand).slice(-1)[0];
    band += " L" + X(last) + "," + axisY + " Z";
    el("path", { class: "band", d: band }, svg);
  }

  let line = "";
  for (let i = 0; i < t.length; i++) {
    line += (i === 0 ? "M" : " L") + X(t[i]) + "," + Y(d[i]);
  }
  el("path", { class: "curve", d: line }, svg);
  el("line", { class: "axis", x1: margin.left, x2: margin.left + plotWidth,
               y1: axisY, y2: axisY }, svg);

  const mx = X(result.median);
  el("line", { class: "median", x1: mx, x2: mx, y1: margin.top, y2: axisY }, svg);
  el("text", { class: "node-label", x: mx, y: margin.top - 2,
               "text-anchor": "middle" }, svg)
    .textContent = thousands(result.median);

  el("text", { class: "axis-label", x: margin.left + plotWidth / 2,
               y: axisY + 33, "text-anchor": "middle" }, svg)
    .textContent = "years before present";

  return svg;
}

/* ---------------------------------------------------------- diagnostics ---- */

/*
 * Convergence at a glance: R-hat and ESS as bars against the thresholds the
 * page judges by, and the thinned trace beside them. A stuck chain, or two
 * chains that never met, is obvious in the trace and merely a number in a
 * table — which is why the trace is here.
 */
function drawDiagnostics(container, result) {
  const params = result.params;
  if (!params || params.length === 0) { container.innerHTML = ""; return null; }

  const rowHeight = 22;
  const nameWidth = 170, barWidth = 110, valueWidth = 46, traceWidth = 230;
  const gap = 12;
  const width = nameWidth + (barWidth + valueWidth) * 2 + traceWidth + gap * 3 + 16;
  const height = 26 + params.length * rowHeight + 10;
  const svg = svgRoot(container, width, height, "fig fig-diag");

  const xRhat = nameWidth + gap;
  const xEss = xRhat + barWidth + valueWidth + gap;
  const xTrace = xEss + barWidth + valueWidth + gap;

  el("text", { class: "col-head", x: xRhat, y: 12 }, svg).textContent = "R̂ (want ≤ 1.01)";
  el("text", { class: "col-head", x: xEss, y: 12 }, svg).textContent = "ESS (want ≥ 200)";
  el("text", { class: "col-head", x: xTrace, y: 12 }, svg).textContent = "trace, by chain";

  /* R-hat is read against 1: a bar from 1.00 to 1.10 spends its whole length on
   * the part of the scale that matters. */
  const rhatSpan = 0.10;
  const essCap = 1000;

  params.forEach((p, i) => {
    const top = 26 + i * rowHeight;
    const mid = top + rowHeight / 2;
    const name = p.type === "age" && p.label ? "age:" + p.label : p.type;
    el("text", { class: "row-label", x: 0, y: mid + 4 }, svg).textContent = name;

    const rhatBad = p.rhat > 1.01;
    const rhatLen = Math.max(1, Math.min(1, (p.rhat - 1) / rhatSpan) * barWidth);
    el("rect", { class: "track", x: xRhat, y: mid - 5, width: barWidth, height: 10 }, svg);
    el("rect", { class: rhatBad ? "bar-bad" : "bar-ok", x: xRhat, y: mid - 5,
                 width: rhatLen, height: 10 }, svg);
    /* Where the threshold falls, so the bar is read against it, not guessed. */
    el("line", { class: "threshold", x1: xRhat + barWidth * (0.01 / rhatSpan),
                 x2: xRhat + barWidth * (0.01 / rhatSpan),
                 y1: mid - 8, y2: mid + 8 }, svg);
    el("text", { class: "row-value", x: xRhat + barWidth + 4, y: mid + 4 }, svg)
      .textContent = p.rhat === null ? "—" : p.rhat.toFixed(3);

    const essBad = p.ess < 200;
    const essLen = Math.max(1, Math.min(1, (p.ess || 0) / essCap) * barWidth);
    el("rect", { class: "track", x: xEss, y: mid - 5, width: barWidth, height: 10 }, svg);
    el("rect", { class: essBad ? "bar-bad" : "bar-ok", x: xEss, y: mid - 5,
                 width: essLen, height: 10 }, svg);
    el("line", { class: "threshold", x1: xEss + barWidth * (200 / essCap),
                 x2: xEss + barWidth * (200 / essCap), y1: mid - 8, y2: mid + 8 }, svg);
    el("text", { class: "row-value", x: xEss + barWidth + 4, y: mid + 4 }, svg)
      .textContent = p.ess === null ? "—" : Math.round(p.ess);

    /* The trace comes with the parameter, so it cannot be paired with the
     * wrong one. It used to be result.trace[i], matched by position against
     * params[i]. */
    if (p.trace && p.trace.length) {
      const chains = p.trace;
      let lo = Infinity, hi = -Infinity;
      for (const ch of chains) for (const v of ch) {
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
      const span = hi - lo || 1;
      chains.forEach((ch, ci) => {
        let dpath = "";
        ch.forEach((v, k) => {
          const x = xTrace + (k / Math.max(1, ch.length - 1)) * traceWidth;
          const yy = top + rowHeight - 3 - ((v - lo) / span) * (rowHeight - 6);
          dpath += (k === 0 ? "M" : " L") + x.toFixed(1) + "," + yy.toFixed(1);
        });
        /* Tagged with the parameter it belongs to, so the pairing is
         * observable from outside rather than trusted. */
        el("path", { class: ci % 2 === 0 ? "trace-a" : "trace-b",
                     "data-param": name, d: dpath }, svg);
      });
    }
  });

  return svg;
}

/* ------------------------------------------------------------- download ---- */

/*
 * Saves a figure as a standalone SVG. The page's colours live in CSS custom
 * properties, which mean nothing outside the page, so they are resolved and
 * inlined here — otherwise the downloaded file renders as invisible black on
 * black for half of its viewers.
 */
function downloadSvg(svg, filename) {
  if (!svg) return;
  const clone = svg.cloneNode(true);
  const source = svg.ownerDocument.defaultView;
  const walk = (from, to) => {
    const style = source.getComputedStyle(from);
    const wanted = ["fill", "stroke", "stroke-width", "stroke-dasharray",
                    "font-family", "font-size", "opacity"];
    let inline = "";
    for (const prop of wanted) {
      const v = style.getPropertyValue(prop);
      if (v) inline += prop + ":" + v + ";";
    }
    to.setAttribute("style", inline);
    for (let i = 0; i < from.children.length; i++) {
      walk(from.children[i], to.children[i]);
    }
  };
  walk(svg, clone);
  clone.setAttribute("xmlns", SVG_NS);

  const blob = new Blob(
    ['<?xml version="1.0" encoding="UTF-8"?>\n', clone.outerHTML],
    { type: "image/svg+xml" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}

/* -------------------------------------------------- competing decay laws ---- */

/*
 * Several posteriors over the same divergence time, overlaid.
 *
 * The point of drawing them together is that they disagree: alternative
 * retention laws fit the same cognate counts and put the date in different
 * centuries. Showing one curve at a time would let a visitor read whichever law
 * they opened as the answer; showing them together makes the choice visible and
 * unresolved, which is what it is until held-out prediction resolves it.
 */
function drawDensityOverlay(container, series) {
  const usable = series.filter((s) => s.result && s.result.grid &&
                                      s.result.grid.t.length > 1);
  if (usable.length === 0) { container.innerHTML = ""; return null; }

  const margin = { top: 16, right: 16, bottom: 44, left: 16 };
  const width = Math.max(320, Math.min(760, container.clientWidth || 720));
  const height = 240;
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  /* The axis stops where the last curve does, not where the search range does:
   * five posteriors squeezed into the left fifth of a mostly empty frame are
   * hard to compare, and comparing them is the whole point of this figure. The
   * range actually searched is stated in the caption instead. */
  let tMin = Infinity, tMax = -Infinity, dMax = 0, gridMax = 0;
  for (const s of usable) {
    const t = s.result.grid.t, d = s.result.grid.density;
    tMin = Math.min(tMin, t[0]);
    gridMax = Math.max(gridMax, t[t.length - 1]);
    tMax = Math.max(tMax, s.result.et90_hi * 1.35);
    dMax = Math.max(dMax, Math.max(...d));
  }
  tMax = Math.min(tMax, gridMax);
  const X = (v) => margin.left + (v - tMin) / (tMax - tMin) * plotWidth;
  const Y = (v) => margin.top + plotHeight - (v / dMax) * plotHeight;

  const svg = svgRoot(container, width, height, "fig fig-density fig-overlay");
  const axisY = margin.top + plotHeight;
  for (const tick of ticksTo(tMax, 6)) {
    if (tick < tMin) continue;
    const x = X(tick);
    el("line", { class: "grid", x1: x, x2: x, y1: margin.top, y2: axisY }, svg);
    el("text", { class: "tick", x: x, y: axisY + 16, "text-anchor": "middle" }, svg)
      .textContent = thousands(tick);
  }

  usable.forEach((s) => {
    const t = s.result.grid.t, d = s.result.grid.density;
    let line = "";
    for (let i = 0; i < t.length; i++) {
      if (t[i] > tMax) break;
      line += (i === 0 ? "M" : " L") + X(t[i]) + "," + Y(d[i]);
    }
    el("path", { class: "curve " + s.cls, d: line }, svg);
    const mx = X(s.result.median);
    if (s.result.median <= tMax) {
      el("line", { class: "median-tick " + s.cls, x1: mx, x2: mx,
                   y1: axisY - 6, y2: axisY + 4 }, svg);
    }
  });

  el("line", { class: "axis", x1: margin.left, x2: margin.left + plotWidth,
               y1: axisY, y2: axisY }, svg);
  el("text", { class: "axis-label", x: margin.left + plotWidth / 2,
               y: axisY + 33, "text-anchor": "middle" }, svg)
    .textContent = "years before present";
  return svg;
}

/* ------------------------------------------- held-out calibration (LOCO) ---- */

/*
 * Attested age against inferred age, with the credible interval as a whisker
 * and the diagonal as perfect agreement. A point on the diagonal is a date the
 * model recovered without being shown it; a whisker that fails to cross the
 * diagonal is a miss the model was confident about, which is the failure worth
 * seeing. Two clocks are drawn together because the comparison between them is
 * the result.
 */
function drawLoco(container, rows) {
  if (!rows || rows.length === 0) { container.innerHTML = ""; return null; }

  const margin = { top: 16, right: 18, bottom: 46, left: 54 };
  const width = Math.max(340, Math.min(680, container.clientWidth || 640));
  const height = 400;
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  let hi = 0;
  for (const r of rows) hi = Math.max(hi, r.true_bp, r.hpd90_hi);
  hi = Math.ceil(hi / 1000) * 1000;
  const X = (v) => margin.left + (v / hi) * plotWidth;
  const Y = (v) => margin.top + plotHeight - (v / hi) * plotHeight;

  const svg = svgRoot(container, width, height, "fig fig-loco");

  for (const tick of ticksTo(hi, 5)) {
    el("line", { class: "grid", x1: X(tick), x2: X(tick),
                 y1: margin.top, y2: margin.top + plotHeight }, svg);
    el("line", { class: "grid", x1: margin.left, x2: margin.left + plotWidth,
                 y1: Y(tick), y2: Y(tick) }, svg);
    el("text", { class: "tick", x: X(tick), y: margin.top + plotHeight + 16,
                 "text-anchor": "middle" }, svg).textContent = thousands(tick);
    el("text", { class: "tick", x: margin.left - 6, y: Y(tick) + 4,
                 "text-anchor": "end" }, svg).textContent = thousands(tick);
  }

  /* Perfect agreement. Everything else on this figure is read against it. */
  el("line", { class: "diagonal", x1: X(0), y1: Y(0), x2: X(hi), y2: Y(hi) }, svg);

  /* Both clocks were given the same held-out target, so they share an x. Nudged
   * apart by a few pixels, otherwise whichever is drawn second hides the other
   * entirely — and the comparison between them is the whole result. */
  const DODGE = 4;
  for (const r of rows) {
    const relaxed = r.clock === "relaxed";
    const cls = relaxed ? "relaxed" : "strict";
    const x = X(r.true_bp) + (relaxed ? DODGE : -DODGE);
    el("line", { class: "whisker " + cls, x1: x, x2: x,
                 y1: Y(r.hpd90_lo), y2: Y(r.hpd90_hi) }, svg);
    if (relaxed) {
      el("circle", { class: "point relaxed", cx: x, cy: Y(r.median_bp), r: 4 }, svg);
    } else {
      el("rect", { class: "point strict", x: x - 3.5, y: Y(r.median_bp) - 3.5,
                   width: 7, height: 7 }, svg);
    }
  }

  el("text", { class: "axis-label", x: margin.left + plotWidth / 2,
               y: margin.top + plotHeight + 34, "text-anchor": "middle" }, svg)
    .textContent = "attested age (years BP)";
  el("text", { class: "axis-label", x: 12, y: margin.top + plotHeight / 2,
               "text-anchor": "middle",
               transform: "rotate(-90 12 " + (margin.top + plotHeight / 2) + ")" }, svg)
    .textContent = "inferred age (years BP)";
  return svg;
}
