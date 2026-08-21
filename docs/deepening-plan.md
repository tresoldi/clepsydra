# Deepening Plan (D1–D6)

Follow-on to `docs/c99-alignment-plan.md`. That plan settled *where the API is*
and *what the conventions are*; this one settles *how much each module hides*.

The assessment behind it is an architecture review of `master` at `3b43160`,
weighted toward the hot spots of the last 60 commits — the `gd_analysis`
runner, the CLI, and the M10 browser layer. Six candidates were rated Strong.
Two further candidates (draw addressing, the likelihood seam) were rated Worth
exploring and are recorded under **Deferred work** below.

---

## Decisions settled before the work starts

| | |
|---|---|
| **Behaviour** | Every step lands **behaviour-preserving**, verified by the worktree/diff procedure. A moved number during a deepening step is a bug, not a result. Each divergence the work surfaces is then fixed in its **own** commit, with the output diff as the evidence. |
| **Public surface** | Two new public headers — `include/gd_options.h` and `include/gd_report.h` — added to the `gd.h` umbrella. `gd_forms` and `gd_grid` stay internal (`src/`), because only the library needs them. |
| **P0** | `gd_pairwise.c` stays test-only. D2 extracts the shared summariser out of it; whether P0 is wired into `gd_analysis_pair_run` or retired is a scientific call, made separately and at no cost. |
| **Report scope** | The five tabular outputs only. The shim's topology and trace emitters stay bespoke — they have no TSV counterpart. |
| **Option set** | One registry row set, reachable from **both** surfaces, validated by **both**. The CLI gains the six options only the shim reaches today; the browser gains the four bounds only the CLI checks today. |
| **Sensitivity** | Library scoring plus CLI. Page UI is a feature, not a refactoring, and is left to M10 proper. |
| **Cadence** | One branch, one commit per step, in the order D1…D6. Full check set on every commit. Stop for review after **D4**, the largest and most user-visible. |

---

## The check set, run on every commit

```sh
make && make test && make analyze && make release && make shim-check
source ~/emsdk/emsdk_env.sh && web/build-wasm.sh   # rewrites web/BUILD_INFO
scripts/wasm-provenance.sh --check
scripts/wasm-verify.sh
node scripts/page-smoke.js                         # google-chrome is present, so it will not self-skip
```

Behaviour preservation, per `STYLE.md` and the end of `c99-alignment-plan.md`
— note the second glob, which is what picks up `calib_sensitivity/run.sh`:

```sh
git worktree add -q --detach /tmp/base <base-sha>
(cd /tmp/base && make native && for e in examples/*/run.sh examples/*/*/run.sh; do "$e"; done)
make native && for e in examples/*/run.sh examples/*/*/run.sh; do "$e"; done
for f in /tmp/base/examples/*/output/*; do diff "$f" "${f#/tmp/base/}"; done
git status --porcelain examples/   # loco.tsv, widen.tsv and friends are TRACKED — must stay clean
```

`make fuzz-run` additionally on D1 (it changes both fuzzed entry points) and on
D3 (it changes what they parse into).

**Every new check gets broken on purpose once and watched go red** before the
commit lands. Four times during C1–C8 a green check turned out not to be
checking; the cheap version of the habit is to ask what the check would print
if the thing it tests were broken.

---

## D1 — Remove the parsers' path branch

*Review candidate 8. Already on `c99-alignment-plan.md`'s "what to do next".*

`gd_tsv_parse` and `gd_newick_parse` take a path first and a buffer second. All
44 call sites in `src/`, `tests/` and `fuzz/` pass `NULL`; the real file reading
lives in `cli/commands.c:189-201`. The branch is unreachable, unfuzzed, and the
less careful of the two — `gd_newick.c:34-36` calls `fseek` twice unchecked and
never bounds the size.

**Do:** drop the parameter from both signatures, delete `read_file_to_buffer`
(`gd_tsv.c:18-63`) and `read_file` (`gd_newick.c:28-51`) and the first branch of
each three-way `if`. `owns_buf` becomes a constant.

**Then:** `grep -rn 'fopen(' src/` finds nothing, which is what `STYLE.md`
already claims. Strike the `gd_tsv`/`gd_newick` item from
`c99-alignment-plan.md` §Deferred and the "one exception" sentence from
`STYLE.md` §"The library does no I/O".

The plan's stated worry — "touching it means rewriting every test's fixture
loading" — no longer applies: every fixture is already a C string literal and
`tests/fixtures/` holds only a `.gitkeep`.

**Risk:** none numerical. Mechanical across 46 files.

---

## D2 — One posterior grid

*Review candidate 4.*

"Grid → normalize → CDF → moments, quantiles, HPD" exists three times:
`gd_pairwise.c:77-341`, `gd_pairwise_p1.c:267-377` and `gd_decay.c:97-189`.
`inv_cdf_p1` is line-for-line `inv_cdf`. The one complete copy — the only one
that fills `gd_pair_result.hpd_*` (`gd_pairwise.h:32-34`) — is in the module no
pipeline calls, so those fields are populated in tests and zero in production.

**New module** `src/gd_grid.{c,h}`:

```c
typedef struct gd_grid_summary {
    double mean, sd, median, mode;
    double et_lo_90, et_hi_90;
    double hpd_lo_80, hpd_hi_80, hpd_lo_90, hpd_hi_90, hpd_lo_95, hpd_hi_95;
} gd_grid_summary;

/* t and logdens are n points of an unnormalised log-density over t. */
int gd_grid_summarize(const double *t, const double *logdens, size_t n,
                      gd_grid_summary *out, gd_error *err);
```

**Order matters, because this is the step that can move numbers.** P1 is the
shipped default, so seed the module with P1's algorithm verbatim and switch in
this order, diffing the examples after each: **P1 first** (must be
bit-identical), **then the four decay laws**, **then P0**. If decay's
`quantile_from_cdf` turns out to differ numerically from P1's, that is a finding
to report — not something to absorb quietly.

**Tests:** new `tests/unit/test_grid.c` — normalisation, CDF monotonicity, HPD
narrower than the equal-tailed interval on a skewed density, analytic cases
where the answer is known, and the degenerate single-point grid.

**Note:** the public `gd_analysis_pair_result` exposes `et_lo_90`/`et_hi_90`,
not HPD, so no shipped output changes. Filling the internal HPD fields is safe.

---

## D3 — One reading of the forms table

*Review candidate 1.*

The rule "highest `primary_weight` wins, first-seen breaks ties" is written
three times — `gd_analysis.c:525-567` (tree), `gd_analysis.c:946-990` (pair),
`tools/analysis/concept_stability.py:63-83` — and the two C copies **already
disagree**: `:544` discards a malformed weight, `:958-963` fails the run on it.
Calibration rows are likewise decoded twice (`:274-318` and `:568-628`), which
`include/gd_analysis.h:60-64` explicitly warns against.

**New module** `src/gd_forms.{c,h}` returning one language-indexed matrix:

```c
typedef struct gd_forms {
    size_t   nconcept, nlang;
    char   **concepts;       /* first-seen order */
    char   **languages;      /* first-seen order */
    size_t  *rows_per_lang;  /* so a caller can count rows it discards */
    int     *labels;         /* nconcept * nlang; interned class, -1 = missing */
} gd_forms;
```

Two invariants are load-bearing for byte-identical example outputs and must be
verified against the current code before the switch, not assumed:

1. **Concepts are interned over every row**, including rows whose language is
   not in the tree — that is what pass one at `gd_analysis.c:494` does, and
   `gd_analysis.c:371` calls first-seen concept order load-bearing for the
   committed outputs.
2. **`n_skipped_rows` counts rows, not languages** (`gd_analysis.c:529`). With a
   tree-agnostic intake it is computed at the mapping step as the sum of
   `rows_per_lang` over languages with no matching tree label.

**Carrying the divergence honestly.** `gd_forms_read` takes a policy for the
malformed-weight case, so in D3 the tree path keeps ignoring and the pair path
keeps failing — behaviour preserved, and the disagreement becomes an explicit,
named, tested parameter instead of an accident. The follow-up fix commit
collapses the policy to one value and shows the diff.

**Also:** `tree_prepare` starts calling `gd_analysis_calibrations_read` instead
of decoding the calibration table a second time.

**Tests:** new `tests/unit/test_forms.c` — weight selection and tie-breaking,
first-seen ordering, synonymy, missing cells, both weight policies, and a
language present in the forms file but not the tree.

---

## D4 — One registry for run options — *stop for review here*

*Review candidate 2. The top recommendation.*

Filling `gd_analysis_tree_config` is spelled four times with different names,
different coverage and different validation. `chains` is already **4** in
`gd_analysis.c:349`, **"4"** in `commands.c:291` and **2** in `app.js:667` and
`index.html:139`. The shim's 18 slots are type-erased through `void *` and once
wrote four bytes into a `bool`; no compiler saw it, and `make shim-check` is
`-fsyntax-only`, so it never could.

**New public header** `include/gd_options.h`:

```c
typedef enum gd_opt_type { GD_OPT_SIZE, GD_OPT_DOUBLE, GD_OPT_U64, GD_OPT_BOOL } gd_opt_type;

typedef struct gd_option {
    const char *name;      /* "rate_logmean" — canonical key, what the shim takes */
    const char *flag;      /* "--rate-logmean" — CLI spelling */
    gd_opt_type type;
    size_t      offset;    /* offsetof into the config struct */
    double      min, max;  /* bounds; min_open distinguishes > from >= */
    bool        min_open;
    const char *help;
} gd_option;

size_t           gd_options_tree_count(void);
const gd_option *gd_options_tree(size_t i);
size_t           gd_options_pair_count(void);
const gd_option *gd_options_pair(size_t i);

int gd_option_set(const gd_option *opts, size_t n, void *config,
                  const char *name, size_t name_len,
                  const char *value, size_t value_len, gd_error *err);
int gd_option_format(const gd_option *o, const void *config,
                     char *buf, size_t bufsz, gd_error *err);
```

**Defaults do not go in the table.** They stay in the existing
`gd_analysis_tree_config_default` / `_pair_config_default`, which are already
the single correct source. Every surface *reads* them from a filled default
struct through `gd_option_format` rather than retyping them — which is what
retires the string literals in `commands.c:289-296`, the usage prose at `:333`
and `:501`, the JS fallbacks at `app.js:84-86,665-679`, and the `value=`
attributes in `index.html`. Keeping defaults out of the table also avoids
narrowing a `uint64_t` seed through a `double`.

**The `void *` hazard goes away without C11.** Rows are emitted by a macro that
pairs `offsetof` with the C99 negative-array-size assert on
`sizeof(((T *)0)->field)`. That is exactly the check that would have caught the
`OPT_INT`-into-`bool` incident, and it fires at compile time in the ordinary
`make`, not only under a browser smoke test.

**Surfaces after this step:** the CLI gains `--hpd-level`,
`--block-age-scaling`, `--trace-draws`, `--cal-drop`, `--cal-widen` and
`--progress-every`; the shim gains the bounds, so a page can no longer send
`rate_logsd=0` or `epsilon=2` into the sampler. A new
`clepsydra_defaults_json()` entry point lets the page populate its controls
instead of hardcoding them.

**Three things the plan did not anticipate, settled during the work:**

- `progress_every` could not take `--progress`. That flag is a 0|1 toggle for
  the stderr reporter, not an interval, so reusing it would have made
  `--progress 1` mean "report every iteration". The interval is
  `--progress-every`; `--progress` stays what it was.
- The canonical option name is stated per row, not derived from the struct
  field. The shim's key set (`warmup`, not `n_warmup`) is a published
  interface that `app.js`, `worker.js` and `wasm_verify.js` all spell.
- Integer options are parsed as integers. Routing them through a `double`, as
  the old shim did, cannot represent `GD_ANALYSIS_KEEP_ALL` or a seed past
  2^53 — the round-trip test caught it immediately.

**Tests:** new `tests/unit/test_options.c` — every row round-trips through
`set` then `format`; unknown names and out-of-bounds values are refused with the
right `GD_E_*`; a coverage test asserts each registry row's `offset` lands
inside its config struct; and one test asserts the tree registry and
`gd_analysis_tree_config` have not drifted apart in count.

**Docs:** `docs/cli-manual.md` for the new flags, `docs/web.md` for the defaults
endpoint, and the `opt_slot` entry in `STYLE.md` §"Things to leave alone" is
deleted, because the thing it warns about no longer exists.

---

## D5 — A result reporting module, two adapters

*Review candidate 3.*

Six independent serializers over the same results: five `out_path → fopen →
fprintf ×n → fclose` blocks in `cli/commands.c` (`:406, :419, :446, :626, :807`,
47 `fprintf` calls between them) and the shim's `jbuf` emitters. The output
schema lives inside command bodies and is re-parsed by name in
`tools/validation/calibration_benchmark.py:121-143`. None of it is testable:
`cli/commands.c` has no C test, and CI only checks the file has ≥1 row.

**New public header** `include/gd_report.h`, built on a table value rather than
a byte sink — so tests assert on cells, not on parsed strings:

```c
typedef struct gd_table gd_table;   /* columns carry name, type and precision */

int  gd_report_nodes      (const gd_analysis_tree_result *r, gd_table **out, gd_error *err);
int  gd_report_diagnostics(const gd_analysis_tree_result *r, gd_table **out, gd_error *err);
int  gd_report_manifest   (const gd_analysis_manifest    *m, gd_table **out, gd_error *err);
int  gd_report_pair       (const gd_analysis_pair_result *r, gd_table **out, gd_error *err);
int  gd_report_sensitivity(const gd_sensitivity_row *rows, size_t n, gd_table **out, gd_error *err);
void gd_table_destroy(gd_table *t);

/* Renderers. The library still does no I/O: the sink takes bytes. */
typedef struct gd_sink { void *ctx; int (*write)(void *, const char *, size_t); } gd_sink;
int gd_table_write_tsv (const gd_table *t, const gd_sink *sink, gd_error *err);
int gd_table_write_json(const gd_table *t, const gd_sink *sink, gd_error *err);
```

Column names, order and per-column precision (`%.1f`, `%.3f`, `%.4f`) come
straight from the current `fprintf` calls, so the TSV is byte-identical. The
CLI's sink wraps `FILE *`, the shim's wraps `jbuf`, and tests use a growable
buffer — three adapters, and the seam is real.

This also retires the `"T1_novel_state"` literal being written independently at
`commands.c:454` and `clepsydra_wasm.c:362` with nothing comparing them.

**Tests:** new `tests/unit/test_report.c` — a golden table per output, plus a
byte-compare of the TSV renderer against a fixture lifted from a current
example output. `cmd_date_pair`'s four hand-rolled
`result_destroy(); return 1;` tails collapse into the one `goto done` shape
`STYLE.md` mandates.

---

## D6 — Move the sensitivity verdict behind the seam

*Review candidate 5.*

Leaving `calib-sensitivity` out of the public header was right — it composes
from `drop_index` and `widen`. But the **scoring** went with it: `node_median`
and `max_other_shift` (`commands.c:525-547`), the 100× divergence ceiling
(`:606-618`) and the covered / unbounded / `.` decision (`:651-670`) are
scientific policy in `static` CLI functions, testable only by shelling out.

**Add to `include/gd_analysis.h`:**

```c
typedef struct gd_sensitivity_row {
    double baseline_median, variant_median, shift, max_other_shift;
    bool   diverged, covered, unbounded;
} gd_sensitivity_row;

int gd_analysis_sensitivity_score(const gd_analysis_tree_result *baseline,
                                  const gd_analysis_tree_result *variant,
                                  const gd_analysis_calibration *withheld,
                                  gd_sensitivity_row *out, gd_error *err);
```

The CLI keeps the loop and the output; the row goes through `gd_report`
from D5.

**This step has an unusually sharp check:** `examples/romance_tree/calib_sensitivity/loco.tsv`
and `widen.tsv` are **tracked**, and the CI example loop regenerates them. If
the scoring moved, `git status --porcelain examples/` goes dirty on its own.

**Tests:** new `tests/unit/test_sensitivity.c` — the ceiling fires and does not
fire on either side of 100×, each of the three row outcomes, and
`max_other_shift` over a hand-built pair of results.

---

**D1-D6 are complete.** What actually landed, and where it differed from the
plan above, is in the commit messages; the short version is that every step
was verified against all 22 example outputs, and three of them moved numbers
only where a commit says so.

---

## Closing out

- **`CONTEXT.md`** — create it. This plan names four concepts the repository has
  never written down as terms: *forms table*, *posterior grid*, *run option*,
  *sensitivity verdict*. Four new modules is the moment a glossary starts
  paying for itself.
- **`STYLE.md`** — update the dependency graph with `gd_forms`, `gd_grid`,
  `gd_options`, `gd_report`; delete the `opt_slot` and `gd_tsv`-path entries
  that D4 and D1 make untrue.
- **`CONTRIBUTING.md`** — same graph, plus a D-series note alongside the M- and
  C-series ones.
- **`docs/c99-alignment-plan.md`** §Deferred — strike the `gd_tsv`/`gd_newick`
  item, which D1 completes.

---

## Part two — the deferred candidates (E1–E4)

The two candidates the review rated *Worth exploring*, planned after D1–D6
landed. Both were deferred for stated reasons; those reasons are the first
thing this part deals with.

### Decisions settled before the work starts

| | |
|---|---|
| **The benchmark gate** | `c99-alignment-plan.md` §Deferred refuses a vtable in the sampler's inner loop "with no benchmark behind it", and `make benchmark` is still `@echo`. E1 makes it real, and E3 is gated on it. This is also that plan's own first recommendation, so it is owed regardless. |
| **Seam shape** | `gd_mcmc_state` carries a small rate vector, not one scalar and not an opaque handle. The inner loop still dereferences a plain array. |
| **T2 scope** | An adapter plus SBC, no user-facing surface. That makes the seam real and proves T2 samples correctly without committing to how binary cognate-set data would be presented — which is a data-format decision, not a refactoring one. |
| **Addressing** | Internal accessors, plus one public accessor *alongside* the existing fields. The fields stay readable, so this does not reopen §Deferred's "value structs stay non-opaque". |
| **Cadence** | One commit per step, same check set as D1–D6, same behaviour rule: a moved number during a refactor step is a bug. |

---

### E1 — Make `make benchmark` real

*Prerequisite for E3. Independently owed since C8.*

**Do:** a reproducible workload with pinned seeds — the Romance and five-IE
examples at their committed chain lengths — reporting wall-clock per phase and
likelihood calls per second, one machine-readable line per measurement so two
runs can be diffed. Record a baseline in `docs/` next to the plan.

**Not** part of `make test`: it takes minutes, and a slow suite is a suite
people stop running.

**Risk:** a benchmark that is not reproducible is worse than none, because it
invites conclusions. Report the spread across at least three runs, and state
the machine. If the spread is wider than the effect E3 is looking for, say so
rather than quoting a number.

**Done.** `docs/benchmark-baseline.md`. Two deviations from the sketch above,
both recorded there: CPU seconds rather than wall clock, and iterations per
second rather than likelihood calls — counting calls needs a wrapper around
the model callback, which would insert exactly the indirection E3 exists to
measure. The instrument resolves about 3%; the first two runs are discarded
because without that the reported spread was 47%, and it was the processor's
frequency ramp rather than the code.

---

### E2 — Address a draw, don't compute its offset

*Review candidate 6.*

`gd_mcmc_result.samples[(chain * n_draws + draw) * n_param + p]` is re-derived
at seven sites: `gd_analysis.c:691,708` and `gd_mcmc.c:725,930,933,936,1032`.
Two of those are the sampler decoding its own output in `compute_confounding`,
which re-derives a node→column map `gd_mcmc_run` already had in hand. The
public `gd_analysis_tree_result.trace` formula is documented at
`include/gd_analysis.h:200-205` and re-derived once more, in
`web/clepsydra_wasm.c:433`.

**Do:**

```c
/* internal, gd_mcmc.h */
size_t gd_mcmc_result_find(const gd_mcmc_result *r, gd_mcmc_ptype t, size_t node);
double gd_mcmc_result_draw(const gd_mcmc_result *r, size_t p, size_t chain, size_t k);

/* public, gd_analysis.h — beside the fields, not instead of them */
double gd_analysis_trace_draw(const gd_analysis_tree_result *r,
                              size_t param, size_t chain, size_t k);
```

**Risk: low, and checked rather than assumed.** `gd_mcmc_chain_sample`
(`gd_mcmc.c:722`) is already a bounds-checked public accessor, not an inner-loop
call, and the other sites are post-processing. E1's baseline confirms it.

**Not in scope:** keying the JSON trace by parameter. `viz.js:339` pairs
`result.trace[i]` with `params[i]` positionally with no key, so a reordering
would mislabel every trace silently — real, but it is a page change with its
own page-smoke checks, and it is recorded below rather than smuggled in here.

---

### E3 — Widen the likelihood seam to a rate vector

*Review candidate 7, first half. Gated on E1.*

`gd_mcmc_model` is a function-pointer interface with exactly one
implementation, `t1_model_loglik`, defined **inside the sampler** at
`gd_mcmc.c:1055`. `src/gd_mcmc.h:34` includes `gd_model_novel_state.h` for no
other reason — that include is the compile-time dependency on T1 that
§Deferred names.

**Do:** replace `gd_mcmc_state.logr` with a rate vector, and move
`gd_mcmc_t1_create` / `_destroy` / `_model` and `t1_model_loglik` into
`gd_model_novel_state.{c,h}`, so a likelihood owns its own adapter and the
sampler stops implementing its own interface. `src/gd_mcmc.h:34` goes with
them; `gd_analysis.c` and four test files include the model header instead.

**Risk: this is the step that can move every tree number.** The likelihood
feeds MCMC accept/reject, so a single changed bit diverges a chain. Point
`rates[0]` at what `logr` was and leave the arithmetic untouched — the
committed outputs are the check, and they must be byte-identical. Then run E1
and report the inner-loop cost against its baseline, whichever way it comes
out.

---

### E4 — The second adapter

*Review candidate 7, second half.*

`gd_model_binary.c` is 426 lines of brute-force-validated T2 likelihood that
production reaches from nowhere: `grep gd_t2_` finds only the module and its
own test. With E3's rate vector it fits.

**Do:** `gd_mcmc_t2_create` beside T1's adapter, mapping the vector to
`lambda_gain`/`lambda_loss`, plus `tests/simulation/test_mcmc_t2.c` running
SBC the way `test_mcmc_tree.c` does for T1 — using `gd_t2_simulate_class`,
which already exists.

**Why this is the step that pays:** one adapter is a hypothetical seam. Two
make it real, and only then has the indirection E3 introduced been earned.

**Risk:** SBC is a statistical check and can fail for tuning reasons rather
than correctness ones. Budget for that, and do not weaken the test to make it
pass — a rank histogram that needs a wider tolerance to look uniform is
reporting something.

---

### Recorded, not planned

- **Keying the JSON trace by parameter.** `viz.js:339`'s positional pairing is
  a real latent mislabel. It wants a payload change, a `viz.js`/`app.js` edit
  and a page-smoke check that a trace belongs to the parameter it is drawn
  under — a page task, not an engine one.
- **A user-facing T2 surface.** T2 wants differently-shaped input than the
  forms table. Deciding how that is presented is a data-format question worth
  its own discussion, and E4 deliberately stops before it.
- **A second sampled rate.** Until `gd_mcmc` grows a rate-pair parameter block,
  T2's adapter can be called but not sampled. That is the next step for anyone
  who wants binary-cognate chronologies, and it is a modelling decision (what
  prior does a loss rate get?) as much as a coding one.

---

## Deferred work

Attractive, and deliberately **not** in this plan:

- **Draw addressing** (review candidate 6). The sampler's output layout is
  public and its index formula is re-derived at five call sites, including
  `gd_mcmc.c:904` decoding its own output. Real, but adjacent to a decision
  `c99-alignment-plan.md` §Deferred already took — value structs stay
  non-opaque — so it wants its own argument rather than a ride on this plan.
- **The likelihood seam** (review candidate 7). `gd_mcmc_model` has exactly one
  adapter, defined inside the sampler; T2 is 426 brute-force-validated lines
  that cannot reach it because `gd_mcmc_state` carries a single scalar `logr`.
  The existing deferral names its own trigger — "revisit if a second likelihood
  needs the sampler" — and that condition is now met. But the same entry warns
  that a vtable in the inner loop is a performance decision with no benchmark
  behind it, and `make benchmark` is still a stub. **Benchmarks first.**
- **Splitting `gd_mcmc.c` or `gd_analysis.c`.** Unchanged from
  `c99-alignment-plan.md`. D3 and D5 shrink `gd_analysis.c` as a side effect,
  which is the right way for that to happen.
