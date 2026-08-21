# The browser page (M10)

`web/` is a static page that runs the whole inference engine in the visitor's
browser via WebAssembly. It is deployed from the repository as-is: the engine is
a committed artifact, so the deploy compiles nothing.

The page is at `clepsydra.tresoldi.org` once the repository is public; until
then, serve it locally:

```
cd web && python3 -m http.server 8080
```

## Why it exists, and what it is not

A dating tool is easy to misread as an oracle. The page exists to let someone
try the models on real data without installing anything, and — just as much — to
put the models' limits in front of them while they do: the interval rather than
the point, the convergence check, the prior-dominance number, the held-out
evidence, and a plain statement of what the tool cannot do.

It is **not** a replacement for the CLI. Its default chain lengths are chosen so
a run finishes while someone waits, and a result meant for publication belongs
in a `clepsydra date-tree` run with a manifest beside it.

## Architecture

```
index.html + app.js + viz.js
        │                        pairwise: straight to wasm on the main thread
        ├──> clepsydra.js/.wasm  (0.3 s — updates as controls change)
        │
        └──> worker.js ────────> clepsydra.js/.wasm
                                 chronology: MCMC off the main thread,
                                 with progress and cancel
                                        │
                                 web/clepsydra_wasm.c   (JSON in/out)
                                        │
                                 src/gd_analysis.c      ← the CLI calls this too
```

The last line is the important one. `src/gd_analysis.{c,h}` is a shared runner:
buffers in, result structs out, no filesystem, no `stdio`. `cli/commands.c`
wraps it with file reading and TSV formatting; `web/clepsydra_wasm.c` wraps it
with JSON. Neither contains any inference logic, so the page and the command
line cannot compute different numbers — a property that is checked, not assumed
(see *Verification*).

## Building the engine

```
source ~/emsdk/emsdk_env.sh
./web/build-wasm.sh
```

`emcc` compiles `src/*.c` plus the shim directly — there is no CMake step, since
this is a plain Makefile project. The build uses the same warning set as the
native build plus `-Werror`: the browser build is not held to looser standards.
It links `-sFILESYSTEM=0`; the engine reads only the buffers it is handed.

`web/clepsydra.js` and `web/clepsydra.wasm` (~135 KB) are **committed**. That
keeps the deploy trivial but introduces one risk — shipping an engine older than
the sources it claims to be — which `scripts/wasm-provenance.sh` closes:

```
./scripts/wasm-provenance.sh --write    # done automatically by build-wasm.sh
./scripts/wasm-provenance.sh --check    # gate for CI and the deploy
```

It hashes every source that goes into the artifact (`src/`, `include/`,
`web/clepsydra_wasm.c`) into `web/BUILD_INFO`. Editing a source without
rebuilding is then a failing check rather than a silently stale page. `cli/` is
deliberately **not** hashed: it is not part of the browser build, and including
it would fail the check for changes that cannot reach the page.

Sources are hashed rather than the artifact because `emcc` is not reproducible
byte-for-byte across toolchain versions, and because the check must run
anywhere `bash` and `git` are — no emsdk, no compiler.

## The JavaScript ↔ C surface

Five entry points, all returning caller-owned strings that JavaScript hands back
to `clepsydra_free`:

| function | returns |
|---|---|
| `clepsydra_date_tree_json(newick, forms, calibrations, opts)` | node ages, parameters, topology, optional trace |
| — with `cal_drop=` / `cal_widen=` | the same, under a stated variation of the calibrations |
| `clepsydra_date_pair_json(forms, langA, langB, decay, opts)` | posterior summary **and the whole grid** |
| `clepsydra_defaults_json()` | every option's default, as the engine holds it |
| `clepsydra_version()`, `clepsydra_git_commit()` | build identity |
| `clepsydra_free(ptr)` | releases a returned string |

**Options are `key=value` lines, not JSON.** Every option the engine takes is a
number, so a key/value list needs no quoting, no escaping and no parser worth
testing; the arguments that could need escaping (language names, the decay law)
are passed separately. An unknown key is an error rather than a silent default —
a misspelt `sample=4000` quietly running 2000 would be a bad way to learn about
typos.

**The keys come from the option registry** (`include/gd_options.h`), which the
CLI reads too, so the page and the command line accept the same set with the
same bounds. Before D4 they did not: the page could send a `rate_logsd` of 0 or
an `epsilon` of 1 that the CLI would have refused, and six options were
reachable only from here. The registry also replaced a hand-maintained table of
`void *` slots whose declared types nothing checked — see `CONTRIBUTING.md` for the
one time that went wrong.

**The trace travels inside its parameter.** `params[i].trace` is
`[[chain0 draws], [chain1 draws], …]`, and `trace_draws` says how many. It was
a top-level `trace` array indexed by parameter until the parallel arrays were
retired: `viz.js` paired `result.trace[i]` with `params[i]` by position, so any
reordering of `params` would have mislabelled every trace with nothing able to
notice. Each drawn path now carries `data-param`, and `page-smoke.js` checks
that name against the row labels.

**Defaults are read, not restated.** `clepsydra_defaults_json()` returns
`{ok, tree:{…}, pair:{…}}` keyed by the same option names, formatted by the
same code that parses them back. `app.js` calls it once at load and fills its
controls from the result, so a default lives in exactly one place — C — instead
of being retyped in `app.js`, `index.html` and the CLI, where `chains` had
already drifted to 2 against the library's 4. The per-dataset chain lengths in
`data/datasets.js` are deliberately *not* defaults: each dataset brings the
length it needs to pass the page's own convergence check.

**Pass datasets as heap pointers, never as `ccall("string")`.** This is the one
non-obvious rule in the whole page. Emscripten's `"string"` argument type copies
the value onto the *wasm stack*; a 100 KB forms file overflows it, and the
failure surfaces as `RuntimeError: memory access out of bounds` from a wasm
frame with nothing to connect it to the cause. Use `stringToNewUTF8`, pass
`"number"`, and `_free` afterwards — see the `call` helper in `web/app.js`,
`web/worker.js` and `scripts/wasm_verify.js`. `STACK_SIZE` is raised to 256 KB
for headroom only; it is not the fix.

**Progress and cancellation.** `gd_mcmc_config` carries an optional
`on_progress` callback (default `NULL`, so the draws are bit-identical without
it — there is a test for exactly that). The shim bridges it to
`Module.onProgress` via `EM_JS`; a truthy return cancels and the run reports
`GD_E_CANCELLED` rather than a partial chain. In the page, cancelling actually
*terminates the worker*: it sits inside one synchronous wasm call and cannot
receive a message until that call returns. A fresh worker starts immediately.

## Verification

Three layers, none of which is "it looked fine":

1. **`./scripts/wasm-verify.sh`** — runs the Romance chronology and the
   Spanish–Portuguese pair through both engines with the same seed and compares
   every reported quantity. The tolerance is relative and documented; the
   observed difference is **zero**. This is M10's central acceptance criterion.
2. **`node scripts/page-smoke.js`** — drives headless Chrome over the DevTools
   Protocol and checks 20 properties of the live page: the pairwise date matches
   a freshly-run CLI (not a pasted constant), a chronology completes in the
   worker, a deliberately short run is *reported as unconverged*, the shipped
   defaults *do* converge, ancient tips sit left of modern ones on the time
   axis, five decay laws give distinct dates, both clocks render in the held-out
   plot, a supplied file is read and dated, and Cancel leaves no partial answer.
3. **`make test`** — includes `shim-check`, which compiles the adapter's
   non-Emscripten branch with the native toolchain so the browser surface cannot
   rot on a machine without emsdk.

## Cost

Measured under node; browsers land in the same range. Wasm is roughly 2× native.

| dataset | run | time | memory |
|---|---|---|---|
| Romance (8 tips, 170 concepts) | 2 × (4000+4000) | ~25 s | ~19 MB |
| Sino-Tibetan (250 concepts) | 2 × (8000+8000) | ~100 s | ~17 MB |
| Deep Indo-European | 4 × (12000+12000) | ~95 s | ~14 MB |
| any pairwise date | — | 0.3 s | negligible |

## Data

Four datasets ship under `web/data/`, copied from the committed examples so a
run on the page reproduces a run in `examples/`. They are fetched on demand
rather than inlined: together they are ~670 KB, and a visitor who only opens the
Pair tab should not pay for the rest. Every one is real published cognate data
under CC-BY-4.0, and the citation is shown on the page beside the selector.

Each dataset carries its own chain length in `web/data/datasets.js`. These were
**measured**, not guessed: each is a setting that reaches R̂ ≤ 1.01 and ESS ≥ 200
on that data, which are the same thresholds the page judges a run by. A global
default would either fail its own convergence check on the harder datasets or
make Romance needlessly slow — and a default that quietly fails its own check
would teach visitors to ignore the warning.

Visitors can also supply their own files ("Use your own data"), which are read
in the tab and go nowhere: there is no server to send them to.

## The figures

`web/viz.js` draws everything as hand-rolled SVG. **No runtime dependencies** —
no CDN, nothing to break when a third party changes; the page works offline and
from any static host. Colours come from the stylesheet's custom properties via
CSS classes, so figures follow the visitor's light/dark theme; `downloadSvg()`
resolves them to literal values on the way out, or a downloaded file would
render black-on-black for half its viewers.

- **Dated tree** — the centrepiece, and the reason the other sibling projects'
  tree renderers were not reused: this axis is *calendar time*, and the tree is
  heterochronous. The present is at the right, so a node's horizontal position
  is its age, and an ancient tip stands where its attestation puts it rather
  than level with the living languages. Every estimated node carries a 90%
  interval bar; calibrated ages are open squares, inferred ages filled circles,
  so what the model was told never looks like what it worked out.
- **Pairwise posterior** — the whole distribution, with the central 90% filled
  and the median dashed. The page reports a median because people ask for one;
  the median is a summary of this, so this gets the space.
- **Convergence** — R̂ and ESS as bars against their thresholds (R̂ scaled
  1.00–1.10, so the whole bar is the part that matters), plus the thinned trace
  of each chain. Chains are drawn separately on purpose: two chains that never
  met is exactly what R̂ measures, and a merged trace would hide it.
- **Held-out calibration** — attested against inferred age with interval
  whiskers and a diagonal for perfect agreement. Both clocks are dodged a few
  pixels apart, because they share an x and whichever drew second would
  otherwise erase the other — destroying the one comparison the figure exists
  to make.

## Deployment

`.github/workflows/pages.yml` publishes `web/` on a push to `main`. Nothing is
compiled: the artifact is committed, so the job's work is proving the artifact
deserves to be published. Four gates:

1. **`wasm-provenance.sh --check`** — the engine was built from the sources now
   in the tree.
2. **Every asset the page loads exists and is non-empty** — the list is
   extracted from `index.html`'s `src`/`href` attributes and from
   `new Worker(...)` in `app.js`, rather than duplicated in the workflow, so
   adding or renaming a file cannot bypass the check.
3. **Every bundled dataset is complete** — read from `data/datasets.js`. These
   are fetched at runtime, so a missing one is invisible until a visitor picks
   that dataset.
4. **`web/data/loco.tsv` is byte-identical to
   `examples/loco_benchmark/loco_multi_report.tsv`** — the benchmark shown on
   the Validation tab must be the one the repository actually produced.

The workflow triggers on `src/**` and `include/**` as well as `web/**`. Those
are not published, but a change there without a rebuild leaves a stale engine on
the page — and if the workflow watched only `web/`, that push would simply not
run the provenance check and the old engine would keep serving. Triggering on
the sources turns that case into a loud failure that says "rebuild the engine".

`web/CNAME` sets the custom domain to `clepsydra.tresoldi.org`; the deploy
checks it is present, since losing it silently reverts the site to the
`github.io` address. Two things live outside this repository and have to match
it: the repository must be **public** (GitHub Pages does not publish from a
private repository on a free plan), and DNS for `clepsydra.tresoldi.org` must
point at `tresoldi.github.io` (a `CNAME` record; an apex domain would need A
records instead). Until both are in place the deploy will run its checks and
then fail at the publish step.

The page always runs locally with `cd web && python3 -m http.server 8080`.

## Not on the page

- **No CLDF or JSON parsing.** The C engine reads TSV and Newick; conversion
  from CLDF is Python-side (`tools/datasets/`), as on the CLI.
- **No tree inference.** The topology is an input here as everywhere else.
- **No full leave-one-calibration-out.** The Tree tab can loosen every
  calibration window and re-date (two runs, and it says so before starting).
  Withholding each calibration in turn needs N+1 chronologies, which is a
  command-line job: `clepsydra calib-sensitivity --mode loco`. The engine
  supports both — `gd_analysis_tree_config.calibrations` — so the difference is
  a judgement about what is reasonable to ask of a browser tab, not a
  capability gap.
- **No T2 (binary cognate-set) tab.** T2 exists in the engine and is tested, but
  comparing it with T1 requires a shared observable space, and offering both
  side by side would invite exactly the raw likelihood comparison that
  `docs/model-cards.md` rules out.
- **No multi-threading.** Chains run sequentially; `SharedArrayBuffer` would
  need COOP/COEP headers that GitHub Pages does not set.
