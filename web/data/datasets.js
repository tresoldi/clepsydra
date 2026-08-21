/*
 * The datasets the page ships with, and the settings each example uses.
 *
 * These are the same files and the same rate priors as the committed examples
 * under examples/<name>/run.sh, so a run here reproduces a run there. Every
 * dataset is real cognate data under CC-BY-4.0; the citation travels with it
 * onto the page rather than living only in a repository file.
 */
"use strict";

const DATASETS = [
  {
    id: "romance",
    name: "Romance (IE-CoR)",
    dir: "data/romance",
    note: "Eight Romance lects with Latin as an ancient tip. Converges cleanly; " +
          "the friendliest place to start.",
    source: "Heggarty et al. 2023, IE-CoR (lexibank/iecor v1.0), CC-BY-4.0",
    rate_logmean: -8.80,
    rate_logsd: 0.80,
    warmup: 4000, samples: 4000, chains: 2, expect_s: 25,
  },
  {
    id: "polynesian",
    name: "Polynesian (Walworth 2018)",
    dir: "data/polynesian",
    note: "A shallow, cleanly dated family with a node-calibrated " +
          "Proto-Polynesian around 2900 BP.",
    source: "Walworth 2018 (lexibank/walworthpolynesian), CC-BY-4.0",
    rate_logmean: -8.80,
    rate_logsd: 1.20,
    warmup: 4000, samples: 4000, chains: 2, expect_s: 28,
  },
  {
    id: "sinotibetan",
    name: "Sino-Tibetan (Sagart 2019)",
    dir: "data/sinotibetan",
    note: "A second family, with three ancient tips. Subgroup ages are well " +
          "estimated; the deep root is honestly weakly identified.",
    source: "Sagart et al. 2019 (lexibank/sagartst), CC-BY-4.0",
    rate_logmean: -8.80,
    rate_logsd: 1.00,
    warmup: 8000, samples: 8000, chains: 2, expect_s: 100,
  },
  {
    id: "five_ie",
    name: "Five ancient Indo-European (IE-CoR)",
    dir: "data/five_ie",
    note: "Deep-time tip dating, and deliberately the hard case: the root age " +
          "is rate-prior-sensitive, the interval is wide on purpose, and it " +
          "needs longer chains than the others to converge.",
    source: "Heggarty et al. 2023, IE-CoR (lexibank/iecor v1.0), CC-BY-4.0",
    rate_logmean: -8.80,
    rate_logsd: 0.80,
    warmup: 12000, samples: 12000, chains: 4, expect_s: 95,
  },
];

/*
 * The per-dataset chain lengths above were measured, not guessed: each is a
 * setting that reaches R̂ <= 1.01 and ESS >= 200 on that dataset — the same
 * thresholds the page judges a run by. They differ per dataset because one
 * global default would either fail its own convergence check on Sino-Tibetan
 * and deep Indo-European, or make Romance needlessly slow; and a default that
 * quietly fails its own check would teach visitors to ignore the warning.
 *
 * `expect_s` is roughly how long that run takes in a browser (about twice the
 * native time), so nobody starts a two-minute run expecting five seconds.
 */
