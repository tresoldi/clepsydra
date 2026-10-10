# Changelog

Notable changes per release. Dates are release dates; the reasoning behind each
structural change is in the commit that made it, and in `docs/deepening-plan.md`
and `docs/c99-alignment-plan.md`.

## Unreleased

Work on `docs/roadmap.md`. Breaking changes to defaults, the CLI, output
columns and example numbers are expected throughout 0.x and are listed here.

### Fixed

- **Missing cognate ids were counted as matches** (roadmap E0.1). A `.`, `?`
  or empty `cognate_id`, or a row with `status` `missing`, `excluded` or
  `absent`, was interned as an ordinary cognate class. Two languages that both
  lacked a concept therefore "shared a cognate" in the pairwise models, and
  formed a spurious class in T1. Such rows are now skipped and counted:
  `missing_cognate_rows` in `run_manifest.tsv`, `n_missing_cognates` in both
  result structs.
- **Unknown or disjoint languages no longer "date"** (E0.8). `date-pair` with
  a language the forms file never mentions, or with two languages sharing no
  observed concept, is an error instead of returning the prior. `date-tree`
  fails when a tree tip has no observed form, naming the tip.
- **Language count and tip resolution** (E0.10). `n_languages` counts tips
  instead of assuming a binary tree (`(nnodes + 1) / 2`). A forms language
  that matches an internal node label is an error.

## 0.2.0 — 2026-08-21

First release intended for use outside the repository. The engine's numbers are
unchanged: every one of the 22 worked-example outputs is byte-identical to
0.1.0-m1 apart from the version line, and the four alternative decay laws whose
moments moved are noted below.

### Added

- **Installation.** `make install` places the static and shared libraries, the
  public headers under `include/clepsydra/`, the `clepsydra` command and a
  pkg-config file. `DESTDIR` and `PREFIX` are honoured; `make uninstall`
  removes exactly what was placed. The soname is `libclepsydra.so.0.2` —
  major *and* minor, because a 0.x minor release may still break ABI.
- **`gd_options.h`** — every run option described once: name, CLI flag, type,
  location and bounds. The CLI, the WebAssembly shim and the browser page read
  the same table, so all three now accept the same options with the same
  validation.
- **`gd_report.h`** — the reported tables and the renderers over them. Output
  goes to a caller-supplied sink, so the library still does no I/O.
- **`gd_analysis_sensitivity_score`** and `gd_analysis_divergence_ceiling` —
  the calibration-sensitivity verdict, previously unreachable inside the CLI.
- **`gd_analysis_trace_draw`** and result accessors, so callers no longer
  reproduce the sampler's index arithmetic.
- **`make benchmark`**, with `docs/benchmark-baseline.md` stating how small an
  effect it can resolve (about 3%).
- T2, the binary cognate-set likelihood, can now be attached to the sampler
  through `gd_mcmc_t2_create`. It cannot yet be *sampled*: that needs a second
  rate parameter, which is a modelling decision.

### Changed

- **`gd_mcmc_state` carries a rate vector** rather than one scalar log-rate.
- **The parsers are buffer-only.** `gd_tsv_parse` and `gd_newick_parse` no
  longer take a path; nothing in `src/` opens a file.
- **One reading of the forms table.** The "highest `primary_weight` wins" rule
  had three implementations that disagreed; an unreadable weight is now taken
  as 1.0 everywhere, and counted in `n_unreadable_weights`. The pairwise
  pipeline used to refuse such a file.
- **Alternative decay laws report slightly different moments.** Sharing one
  posterior-grid summariser moved `mean` by at most 1.0e-2 years and `sd` by
  2.2e-2, on values of 3000–7000 years; medians and every equal-tailed interval
  are unchanged. P1, the default, did not move at all.
- `run_manifest.tsv` reports `unreadable_primary_weights` when there were any.
- The browser payload carries each parameter's trace inside the parameter.

### Fixed

- A heap-buffer-overflow reading a TSV field containing an embedded NUL,
  reachable from a pasted file in the browser. Found by libFuzzer.
- `run_manifest.tsv` could record a commit the binary was not built from.
- The Makefile tracked no header dependencies, so editing a header produced a
  mixed build rather than a rebuilt one.
- The WebAssembly shim truncated error messages into a 192-byte buffer.
- `gd.h` did not include the whole public API it claims to.

## 0.1.0-m1

Milestones 0–12: the pairwise and tree likelihoods, the fixed-tree sampler with
relaxed clocks and robust observations, calibrations, convergence diagnostics,
the validation suite and model cards, the CLI, the reproducible worked examples
on IE-CoR data, and the WebAssembly browser page.
