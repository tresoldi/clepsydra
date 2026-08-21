/*
 * The guide panel's text.
 *
 * Written for someone who knows historical linguistics but not this tool, and
 * who is entitled to know what the numbers do and do not mean before quoting
 * one. Kept separate from app.js so the prose can be edited without touching
 * the machinery.
 */
"use strict";

const GUIDE_SECTIONS = [
  {
    title: "What this page is",
    html: `
      <p>Clepsydra estimates <em>when</em> languages diverged from
      basic-vocabulary cognate data. Give it which languages share a cognate
      for each concept, and — for a whole family — a tree and some age
      calibrations, and it returns a posterior distribution over each
      divergence date.</p>
      <p>The whole engine is compiled to WebAssembly and runs in this tab.
      Nothing is uploaded. It is the same C code the command-line tool runs, so
      a date here is the date the CLI computes, not an approximation of it.</p>`,
  },
  {
    title: "Reading a result",
    html: `
      <p>Every estimate is an interval with a point inside it. The point is a
      median because people ask for one; the interval is the actual answer. A
      wide interval means the data does not pin the date down — that is
      information, not a malfunction, and narrowing it by choosing a tighter
      prior would be inventing certainty.</p>
      <p>Dates are in <strong>years BP</strong> — years before the present.
      2,000 BP is roughly the turn of the first millennium.</p>`,
  },
  {
    title: "The pairwise model",
    html: `
      <p>Two languages, no tree. Given how many of the shared concepts are
      cognate, and a prior on how fast words are replaced, it integrates over
      the rate to get a posterior over the divergence time. This is the
      hierarchical constant-rate model, <code>P1</code>.</p>
      <p>Watch <strong>prior dominance</strong>. It reports how much of the
      answer came from the prior rather than the data. Above about 0.5 the pair
      is telling you very little about its own date, and the interval you see
      is close to what you would have got without looking at the data at all.</p>`,
  },
  {
    title: "The chronology model",
    html: `
      <p>A whole family at once. Each concept evolves down the tree: along a
      branch a language either keeps the class it inherited or replaces it with
      one that has never occurred anywhere. This is the strict novel-state
      model, <code>T1</code>, and it is fitted by MCMC over every node age at
      once, subject to the calibrations.</p>
      <p>The tree topology is an <em>input</em>. Clepsydra does not infer it and
      does not check it; if the tree is wrong, the dates are answers to the
      wrong question.</p>`,
  },
  {
    title: "Calibrations, and where absolute time comes from",
    html: `
      <p>Cognate data alone gives relative time — how much change separates two
      points, in units of rate × time. Turning that into calendar years takes
      either a rate prior (an assumption about how fast vocabulary changes) or a
      calibration: a node or tip whose age is known independently, such as an
      attested ancient language.</p>
      <p>On the dated tree, calibrated ages are drawn as open squares and
      inferred ages as filled circles. The distinction matters: a square is
      something the model was told, and it will happily reproduce it.</p>`,
  },
  {
    title: "Convergence, and when not to believe a run",
    html: `
      <p>MCMC gives a correct answer only if it has run long enough. Two
      numbers say whether it has. <strong>R̂</strong> compares chains that
      started in different places; at 1.00 they agree, and above about 1.01 they
      have not met yet. <strong>ESS</strong> is how many effectively independent
      draws the chain produced; below about 200 the summaries are too noisy to
      quote.</p>
      <p>The convergence panel shows both against those thresholds, plus the
      trace of each chain. If the chains wander apart, or one gets stuck, it is
      visible there long before it is obvious in a table. When the page says a
      run has not converged, the dates in it should not be quoted — raise the
      warm-up and sample counts, or add chains.</p>`,
  },
  {
    title: "Strict and relaxed clocks",
    html: `
      <p>The strict clock assumes every branch changes at the same rate. The
      relaxed clock gives each branch its own multiplier, shrunk toward one.
      Real lineages are rate-heterogeneous, so the relaxed clock is often the
      better description — but it is also more flexible, and more flexible is
      not the same as more correct.</p>
      <p>Which is better is an empirical question, answered by held-out
      prediction rather than preference. On the Validation tab you can see it
      answered both ways: on simulated strict-clock data the relaxed clock does
      not help, and on real data across three families it clearly does.</p>
      <p>Watch for the confounding warning. A faster branch and an older node
      explain the same data equally well, so when branch rates are free the ages
      can be less identified than their intervals suggest.</p>`,
  },
  {
    title: "Decay laws",
    html: `
      <p>The default model assumes vocabulary is replaced at a constant rate, so
      the share of retained cognates decays exponentially. That is an
      assumption, and there are others: Starostin's accelerating law, Vasilyev's
      two-component law with a stable core, and a configurable
      two-component form.</p>
      <p>The Decay tab runs all of them on one pair so you can see how far apart
      they land. They are offered as competing models to compare — not as
      settings to pick the one you like. Choosing between them honestly means
      comparing held-out predictions in a shared observable space.</p>`,
  },
  {
    title: "What this cannot do",
    html: `
      <p>It does not infer trees, judge cognate coding, detect loanwords for
      you, or reconstruct proto-forms. It cannot tell you that a date is right —
      only what the model and data imply, given assumptions you can inspect.</p>
      <p>Deep-time estimates are the weakest case. The further back you go, the
      more of the answer comes from the rate prior rather than from the data;
      the bundled deep Indo-European example is included precisely because it
      shows that honestly, with an interval wide enough to span competing
      scholarly hypotheses.</p>
      <p>The defaults here are also shorter than a publication run would be, to
      keep the page responsive. Real analyses belong on the command-line tool,
      where the run can be as long as it needs to be and the whole thing is
      reproducible from a manifest.</p>`,
  },
];
