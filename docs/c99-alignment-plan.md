# C99 Alignment Plan

Assessment of `clepsydra` against *Modern C99 Refactoring Guide for Coding
Agents*, and a milestone plan to close the gaps. Structure follows the guide's
§38.

Milestones here are numbered **C1–C8** so they cannot be confused with the
project's own scientific milestones M0–M12 in `PLAN.md`. They are engineering
work, not science, and none of them is meant to change a single posterior.

Baseline measured on `1c604a0` (2026-08-20). The assessment below describes
that baseline; where a milestone has since landed, its section says so.

**Status: C1–C8 complete.** The plan is finished; see the closing section.

---

## Repository assessment

### Size and shape

| Layer | Files | Lines |
|---|---|---|
| Library core (`src/`) | 17 modules (`.c` + `.h`) | 8,933 |
| Public headers (`include/`) | 3 | 280 |
| CLI (`cli/`) | 2 + header | 867 |
| WebAssembly shim (`web/clepsydra_wasm.c`) | 1 | 501 |
| Tests (`tests/`) | 16 C files | 6,179 |
| Python toolchain (`tools/`) | 21 modules | 3,053 |

About 10,600 lines of first-party C, plus 6,200 of C tests.

### Build system

GNU Make, single `Makefile`, no generated sources, no vendored code, no
external dependencies beyond libm. Retain it — the guide (§2.1) says prefer the
existing build system when it works, and this one does. `CC` and `SANFLAGS` are
overridable, which is how CI drives its matrix.

Compile flags are already **exactly** the guide's §16 baseline:

```
-std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
-Wstrict-prototypes -Wmissing-prototypes
```

Verified: a clean `make native` on gcc 15 emits **zero warnings**. CI enforces
`-Werror` on gcc and clang separately, which is the guide's recommended posture
(strict in CI, not imposed on downstream builds).

Gap: only one configuration exists (`-O2 -g`). There is no `-O0 -g` debug build
and no `-O2 -DNDEBUG` release build (§17). Low-cost to add, low urgency, since
`NDEBUG` currently changes nothing — see "Errors" below.

### C dialect actually used

Strict C99, C standard library and libm only. No compiler extensions, no VLAs
on user-controlled sizes, no `__attribute__`, no anonymous unions. The claim in
`CONTRIBUTING.md` holds up under inspection. `-Wpedantic` passing on both gcc and
clang is real evidence, not an aspiration.

### Public API — the central finding

**The repository has two APIs, and the declared one is not the real one.**

`include/gd.h` declares 35 functions and presents itself as the public surface
("frozen API surface"). Of those, **12 return `GD_E_NOT_IMPLEMENTED`** —
`gd_load_forms_tsv`, `gd_load_priors_tsv`, `gd_load_calibrations_tsv`,
`gd_load_analysis_tsv`, `gd_pairwise_integrate`, and the entire `gd_chain_*`
family. `src/gd_context.c` (560 lines) exists almost entirely to back this
surface, and **no shipped consumer calls any of it**: outside its own
`tests/unit/test_context.c`, `gd_context` appears nowhere in `cli/`, `web/`, or
the rest of `src/`.

The API that both real consumers use is `src/gd_analysis.h` — an *internal*
header. `cli/commands.c` and `web/clepsydra_wasm.c` are both compiled with
`-Isrc` and include it directly.

This is not an accident of sloppiness; it is history. `gd.h` is the M0/M1
skeleton, and `gd_analysis` is what M10 built when the browser needed a
buffer-only pipeline. But the result is that the guide's §2.3 question — "which
symbols are intentionally public?" — has no answer a reader can trust:

- 162 symbols are externally linked from `libclepsydra.a`;
- 38 are declared in `include/`;
- **124 are internal-header symbols with external linkage.**

The 124 are not accidental leakage in the guide's sense — every one is declared
in some `src/gd_*.h`, and the test suite deliberately links against them
(verified: no symbol is exported without a declaration anywhere). But nothing
distinguishes "internal, exposed for tests" from "internal, exposed by
oversight", and no linker script or visibility attribute marks the boundary.

Ironically, `src/gd_analysis.h` is the *best-documented header in the
repository* — it states ownership, nullability, defaults, buffer-length
conventions and lifetime rules for nearly every field, and would satisfy the
guide's §30 on its own. It is a public header living in a private directory.

Dead public surface: `gd_model_vtable` is declared in `include/gd_types.h` and
referenced nowhere in the entire tree.

### Major modules and dependencies

17 modules in `src/`, each a `.c` + `.h` pair with a narrow responsibility. No
`utils.c`, `helpers.c`, `misc.c` or `common.c` — the guide's §24 anti-pattern is
absent. The one deviation is `gd_pairwise_p1.c`, a second `.c` behind
`gd_pairwise.h`, which is a deliberate split of a large module rather than a
dumping ground.

The dependency graph was extracted and is **acyclic**, pointing inward exactly
as §23 prescribes:

```
                    gd_analysis
                         │
        ┌────────────────┼────────────────┐
        │                │                │
     gd_mcmc      gd_model_novel_state  gd_pairwise(_p1)
        │           gd_model_binary      gd_decay
        │                │                │
   gd_calibration   gd_newick        gd_diagnostics
        │                │
     gd_tree ─────── gd_tsv ─────── gd_rng ─────── gd_math
        └────────────────┴────────────────┘
                         │
                    gd_memory
                         │
                     gd_error
                         │
                     gd_types
```

One coupling worth naming: `gd_mcmc.h` includes `gd_model_novel_state.h` and
declares `gd_mcmc_t1_create`, so the sampler has a compile-time dependency on
one specific likelihood. The `gd_model_vtable` in `gd_types.h` was evidently
meant to break this and never was.

`gd_context` sits outside this graph entirely, depending on
`gd_memory`/`gd_tsv`/`gd_newick`/`gd_tree` and depended upon by nothing.

### Test coverage

**176 assertions across 16 binaries, all passing, all AddressSanitizer-clean**
(measured, not quoted; 168 across 15 after C1 retired the skeleton's own test). Structured by kind, which is what §27 asks for:

- `tests/unit/` — module-level (math, rng, tsv, newick, calibration, diagnostics, context)
- `tests/likelihood/` — T1 and T2 validated against independent brute-force enumeration
- `tests/simulation/` — SBC rank uniformity, recovery, R̂, posterior-predictive adequacy
- `tests/regression/` — pairwise, decay laws, and the analysis pipeline

CI adds clang/UBSan (with `integer` checks, unsigned wrap deliberately
excluded), a `-Werror` build on both compilers, every worked example from a
clean checkout, 46 Python tests, WebAssembly-vs-CLI numerical equivalence, and a
headless-Chrome page smoke test.

This is a stronger baseline than the guide's Milestones 0, 1, 2, 10 and 18
assume. **The characterization work the guide front-loads is already done.**

Note: `README.md` claimed "173/173"; the measured figure was 176, and is 168
after C1. Corrected in the README.

### Memory-management patterns

- **Zero mutable global state.** Verified across all of `src/`: the only
  file-scope datum is `static const char GD_VERSION_STRING[]`. The guide's §10
  and Milestone 8 are already satisfied; §10's own caveat ("do not eliminate
  harmless immutable constants") covers the one survivor.
- Ownership is generally explicit, with matched `_create`/`_destroy` and
  `_result_destroy` pairs and a documented handoff rule in `gd_analysis.h` ("on
  a successful run the tree moves into the result").
- **The arena is dead.** `gd_arena_create` was called at exactly one site
  (`gd_context.c:154`) and `gd_arena_alloc` from nowhere: a working, documented,
  block-bump allocator that nothing allocates from. Since C1 retired
  `gd_context.c`, it has **no callers at all**, which makes C3's decision about
  it easier rather than harder.
- **Allocation discipline is split** (settled by C3). `gd_malloc_checked` /
  `gd_calloc_checked` exist and do checked multiplication, but are used in only
  4 of 17 modules (23 call sites). Elsewhere there are ~128 raw `malloc`/`calloc` calls, of
  which **~25 are `malloc(n * sizeof(T))` with unchecked multiplication** on an
  `n` derived from parsed input (concept counts, node counts, draw counts) —
  `gd_diagnostics.c`, `gd_mcmc.c`, `gd_pairwise_p1.c`, `gd_model_binary.c`,
  `gd_calibration.c`. The `calloc` sites are fine; the C library checks those.

### Error-handling patterns

`gd_error.h` defines a clean 18-value `GD_E_*` code space and `gd.h` documents
the contract as "returns 0 on success, nonzero `gd_error` code on failure".
**Four different `int` conventions actually coexist:**

1. **`GD_E_*` code returned directly** — the documented contract, used in
   `gd_tsv`, `gd_newick`, `gd_calibration`.
2. **Bare `0`/`1` success flag** — `gd_analysis_*` returns literal `1` on
   failure and puts the real code in `err`. Since `GD_E_NOMEM == 1`, a caller
   following the documented contract would read every analysis failure as an
   allocation failure.
3. **`0`/`-1`** — `gd_rng_serialize_impl`, `gd_rng_restore_impl`,
   `rank_normalize`, `build_grid`, `cell_double`.
4. **Predicate `int`** — `gd_error_failed`, `gd_is_finite`,
   `gd_cal_check_order`; plus three-valued `gd_sign`.

**This is a contract inconsistency, not a live defect.** Both real callers
(`cli/commands.c`, `web/clepsydra_wasm.c`) were checked: they treat the return
as a boolean and read `err.code` for the real classification, which is correct
under convention 2. Nothing is currently mis-reporting. But the header says one
thing and the code does another, and the next caller has a one-in-four chance
of guessing wrong.

Other observations:

- `assert()` is used **nowhere** — the library validates inputs and returns
  errors, which is the right half of §12's guidance, but it also means there is
  no expression of programmer invariants, and `NDEBUG` is inert.
- `exit()`/`abort()` appear nowhere in `src/`. Verified.
- `gd_error.message` is a **borrowed `const char *` with no owned storage**.
  Runtime-composed messages need somewhere to live, and two different
  workarounds exist: the caller-supplied `gd_analysis_diag.msgbuf` scratch, and
  a `char msg[192]` buried inside `gd_model_novel_state.c`. With `msgbuf ==
  NULL` a caller silently gets a less specific message.

### Buffer and string representations

Good: the parse boundary is `(const char *buf, size_t len)` throughout, no
`strlen()` on caller data, no `strcpy`/`strcat`/`sprintf`/`gets`/`strtok`
anywhere in the tree, and `snprintf` used consistently. One `atof` survives at
`gd_analysis.c:876` (a form weight), where the codebase elsewhere uses the
error-returning `gd_tsv_parse_double`.

Gaps:

- **`len == 0` means "NUL-terminated"** in `gd_tsv_parse`, `gd_newick_parse` and
  `gd_analysis_*_input`. Documented and convenient, but it makes a genuinely
  empty buffer inexpressible and is precisely the sentinel-dependence §8 warns
  about.
- **Fixed-size buffers with unchecked truncation.** The CLI builds output paths
  into `char path[1200]` with 9 `snprintf` calls, **none of which check the
  return value**. There is no overflow — `snprintf` is bounded — but a long
  `--out` directory truncates silently and the run writes to the wrong file.
  Also `char model[64]` in `gd_analysis_pair_result` and `char msgbuf[256]` in
  the CLI.

### Type semantics

- **`stdbool.h` is included nowhere in the entire repository.** About 30 struct
  fields carry boolean semantics as `int`: `is_hard`, `relaxed`, `is_leaf`,
  `feasible`, `excluded`, `resolved`, `lca_built`, `compatible`,
  `absolute_anchored`, `branch_rate_confounded`, `block_age_scaling`,
  `store_pointwise`, `sample_rate`, `sample_branch`, and so on.
- **Magic-value `int`s where an enum is meant.** `gd_options.branch_lengths`
  (`/* 0=ignore,1=start,2=duration */`) and
  `gd_analysis_calibration.type` (holds a `gd_cal_type`, typed as `int`).
- **Stringly-typed configuration.** `gd_options` has 8 `const char *` fields
  whose valid values are enumerated in comments
  (`"strict"|"family"|"relaxed"`, `"none"|"global"|"concept"`, …), and
  `gd_analysis_pair_config.decay` likewise. The compiler checks none of it.
- `size_t` is used correctly for sizes and counts throughout; fixed-width types
  appear where width is semantic (`uint64_t` seeds and PCG32 state). The guide's
  warning against mechanical `int` conversion is well respected.
- `const` is used generously — 159 const-qualified pointer parameters and
  fields across the headers.

### Separation of core from I/O

Largely correct and deliberately so. `gd_analysis` is buffer-only by design
("this module does the science"), which is exactly §22, and it is *why* the same
engine runs in the CLI and in the browser and agrees to the last digit.

One exception: **`gd_tsv.c` and `gd_newick.c` both `fopen`/`fread`/`fclose`**.
Their signatures are `(path, buf, len)` with path-or-buffer semantics, so the
two lowest-level parsers carry filesystem I/O into the library core. It works
because the WebAssembly build simply passes `path == NULL`, but it means the
core is not I/O-free, and a fuzz harness or a binding pulls stdio along with it.

### Tooling status

| Tool | Status |
|---|---|
| gcc build | ✅ warning-free, `-Werror` in CI |
| clang build | ✅ warning-free, `-Werror` in CI |
| AddressSanitizer | ✅ default for `make test`, CI-enforced |
| UndefinedBehaviorSanitizer | ✅ CI (clang, `+integer`) |
| LeakSanitizer | ✅ implied by ASan default |
| Static analysis | ❌ **none** |
| Fuzzing | ❌ **none** |
| Benchmarks | ❌ `make benchmark` is `@echo "Benchmark target (M2+)"` |
| CI | ✅ 4 jobs, 6 configurations |

No `LICENSE`, no `STYLE.md`, no `CONTRIBUTING.md`. The project notes' conventions
section is doing the job of a style guide, informally.

### Architectural risks

1. **The two-API split** is the one that compounds. Every future contributor
   must independently rediscover that `include/` is decorative and `src/` is
   the API.
2. **Unchecked allocation-size arithmetic** on input-derived counts is the only
   finding with a plausible path to memory corruption, and it is the one the
   existing sanitizer coverage is least likely to surface (it needs a hostile or
   absurdly large dataset, not a normal one).
3. **The parsers are unfuzzed.** `gd_tsv_parse` and `gd_newick_parse` consume
   untrusted text, already take `(buf, len)`, and are the exact shape §21 calls
   a strong fuzzing candidate. The browser page makes this concrete: a visitor
   pastes arbitrary Newick into a page that runs C in their own tab.
4. **Dead scaffolding** (arena, vtable, `gd_context` loaders) reads as
   infrastructure a maintainer must respect rather than debris they may delete.

---

## Proposed target architecture

The end state is not far from where the repository already stands. Four
structural changes:

**One API, honestly located.** `libclepsydra`'s public surface becomes what its
consumers actually use. `gd_analysis.h` moves to `include/clepsydra/` (or stays
put and is declared internal, with a thin public façade over it — see C1 for
the decision). `include/gd.h` stops advertising 12 functions that return
`GD_E_NOT_IMPLEMENTED`. `cli/` and `web/` stop needing `-Isrc`.

**One error convention.** Every `int`-returning function in the tree means the
same thing by its return value, and the header says which thing. Predicates
return `bool` and are visibly not status functions.

**One allocation convention.** Every allocation whose size is a product goes
through a checked helper. Either the arena earns its place or it goes.

**Analysable and fuzzable core.** Static analysis runs in CI with a documented
baseline; the two parsers have fuzz harnesses calling production code; the
library core does no filesystem I/O.

Everything else the guide asks for — module boundaries, acyclic dependencies,
no global state, sanitizers, dual compilers, `const` discipline, pointer+length
buffers, `goto cleanup`, minimal preprocessor use — this repository already
does, and the plan should not disturb it.

---

## Milestone plan

### C1 — Resolve the two-API split — **DONE**

Resolved in favour of the second reading: `gd_analysis` is where the library
actually landed, so it was promoted and the skeleton retired. What changed:

- `src/gd_analysis.h` → **`include/gd_analysis.h`**. It no longer includes the
  internal `gd_mcmc.h`; the one thing it needed from it, `gd_mcmc_ptype`, moved
  to `gd_types.h`, where it belongs — a result labels its parameters with it.
- `src/gd_tree.h` split. The read surface (the structs, `GD_NODE_NONE`,
  `gd_tree_find_label`, `gd_tree_mrca`, `gd_tree_depth`, `gd_tree_is_leaf`) is
  now **`include/gd_tree.h`**; construction (`_create`, `_destroy`,
  `_add_child`, `_build_traversals`, `_build_lca`) moved to the internal
  **`src/gd_tree_build.h`**. A caller receives a built tree and destroys it
  through the result that owns it, never directly.
- **`include/gd.h` is now an umbrella** over the four public headers, plus
  `gd_version` / `gd_git_commit`. Every function it declares is implemented.
  It documents the ownership vocabulary and states that `src/gd_*.h` are not
  part of the API.
- `include/gd_types.h` **trimmed from 116 lines to the vocabulary that is
  genuinely shared**: `gd_error`, `gd_mcmc_ptype`, the `gd_rng` forward
  declaration. Deleted with the skeleton: `gd_context`, `gd_chain`, `gd_sample`,
  `gd_progress`, `gd_validation_report`, `gd_model_id`, `gd_options`, and the
  never-referenced `gd_model_vtable`. `gd_pair_result` was not public at all —
  no consumer touches it — so it moved to `src/gd_pairwise.h` beside the engines
  that fill it.
- **`src/gd_context.c` (560 lines) deleted**, along with
  `tests/unit/test_context.c`. `gd_version` / `gd_git_commit` were the only live
  things in it and are now `src/gd_version.c`. The version string is unchanged
  on purpose: every committed example manifest carries it.
- **The `gd_rng_*_impl` duplication is gone.** Nine thunks existed only to give
  the old `gd.h` a set of public names distinct from the internal ones; with
  that header retired they were dead weight, so the two name sets collapsed into
  one. `gd_rng_next` had zero callers and went with them.
- **`-Isrc` dropped from the CLI and shim builds.** Both now compile against
  `include/` alone, enforced by the build rules. Tests keep it, deliberately and
  with a comment saying why. Test sources also stopped reaching into `src/` via
  `../../src/x.h` relative paths, which is what had left the unit-test rule
  inconsistent with the other three.

Measured outcome:

| | before | after |
|---|---|---|
| Exported symbols | 162 | **132** |
| Declared in `include/` | 38 | **25** |
| Returning `GD_E_NOT_IMPLEMENTED` | 12 | **0** |
| Consumers needing `-Isrc` | 2 | **0** |

Verified: gcc and clang both build warning-free under `-Werror`; all 18 headers
still compile standalone; **168/168 tests pass** under ASan (176 before, minus
the retired skeleton's own 8); the rebuilt WebAssembly engine agrees with the
CLI with **largest relative difference 0**; **22/22 page checks** pass; and
`examples/pair_simple` and `examples/romance_tree` reproduce their committed
artifacts **byte for byte**.

The original milestone specification follows.

**Objective.** Make `libclepsydra`'s public surface be the one its consumers use.

**Current evidence.** 12 of 35 `gd.h` functions return `GD_E_NOT_IMPLEMENTED`;
`gd_context` has no consumer outside its own test; both real consumers include
`src/gd_analysis.h` via `-Isrc`; 124 of 162 exported symbols are
internal-header symbols; `gd_model_vtable` is referenced nowhere.

**Concrete changes.**
1. Decide the public surface. Recommended: promote `gd_analysis.h` to
   `include/` as the public API (it is already documented to that standard) and
   demote `gd.h` to the frozen-but-unimplemented status it actually has —
   either deleting the stub declarations or moving them to a clearly labelled
   `include/gd_planned.h`. This is a decision for the maintainer, not a default.
2. Delete `gd_model_vtable` from `include/gd_types.h`, or implement it and use
   it to break `gd_mcmc.h`'s dependency on `gd_model_novel_state.h`. Not both.
3. Retire `src/gd_context.c`'s 12 stub bodies. Either implement the loaders on
   top of `gd_analysis`, or remove the declarations so the library does not ship
   functions that cannot work.
4. Drop `-Isrc` from the `cli/` and `web/` build rules once the promoted header
   is reachable from `include/`.
5. Document the tests' privileged access to internal headers explicitly — they
   legitimately link against internals, and that should be a stated policy, not
   an inference.

**Affected files.** `include/gd.h`, `include/gd_types.h`, `src/gd_context.c`,
`src/gd_analysis.h`, `tests/unit/test_context.c`, `Makefile`,
`web/build-wasm.sh`, `CONTRIBUTING.md`, `README.md`.

**Tests.** No behavioural change intended. `tests/unit/test_context.c` shrinks
or is deleted with `gd_context`. Existing suite must stay green apart
from the skeleton's own tests; the WebAssembly-vs-CLI equivalence check is the
real guard here, since it exercises both consumers against the same engine.

**Tooling.** `make test`, `-Werror` build on gcc and clang, `./web/build-wasm.sh`
then `./scripts/wasm-verify.sh` and `node scripts/page-smoke.js` (the shim's
include path changes, so the engine must be rebuilt and re-verified — this is
already the documented rule in `CONTRIBUTING.md`).

**Risks.** Touches the WebAssembly build. `scripts/wasm-provenance.sh` hashes
`src include web/clepsydra_wasm.c`, so moving a header between those directories
invalidates provenance and forces a rebuild — expected, but it means this
milestone cannot be done on a machine without emsdk.

**Compatibility.** Breaking, by intent, and the right moment for it: the
repository was private until the M10 page shipped, has two in-tree consumers and
no known external users. Per §29 this is a deliberate, isolated, documented
break — not an accident during an internal refactor.

**Acceptance criteria.** No shipped function returns `GD_E_NOT_IMPLEMENTED`;
public headers compile standalone (all 18 already do — keep it that way);
`cli/` and `web/` build without `-Isrc`; a reader can name the public API from
`include/` alone.

**Depends on.** Nothing.

---

### C2 — Unify the error model — **DONE**

What changed:

- **The convention is written down** at the top of `include/gd_error.h` and is
  now true: an `int`-returning function returns a `GD_E_*` code, and that code
  equals `err->code`. `gd_error_set` and the new `gd_error_setf` **return the
  code**, so a failure branch reads
  `return gd_error_set(err, GD_E_PARSE, "...")` and cannot drift from what it
  recorded.
- **`gd_error` owns its message.** `const char *message` became
  `char message[GD_ERROR_MSG_MAX]`, and `gd_error_setf` formats a message that
  names a runtime value straight into it. This is what retired **both**
  workarounds: `gd_analysis_diag` (with its `msgbuf`/`msgbuf_size` and the
  degraded wording when a caller passed none) is **deleted from the public
  API**, and so is the `char msg[192]` that `gd_t1_ws` carried. Two callers and
  a test no longer pass scratch buffers. Side effect: `gd_model_novel_state.c`
  no longer needs `<stdio.h>`.
- **The bare `1` is gone.** Twelve `return 1` sites in `gd_analysis.c` — where
  `1` silently read as `GD_E_NOMEM` — now return the real code. The five `0`/`-1`
  sites (`gd_rng_serialize`/`_restore`, `rank_normalize`, `build_grid`,
  `cell_double`) return codes too; `gd_rng`'s magic `16` became
  `GD_RNG_SERIAL_SIZE`.
- **The 44 `goto done` paths** keep their shape. A `fail_code()` helper resolves
  the status from the error at the label, so no failure branch has to restate a
  code it just set, and none of them can disagree with `err`.
- **Predicates return `bool`**: `gd_error_failed`, `gd_is_finite`,
  `gd_tsv_is_missing`, `gd_tsv_field_eq`, `gd_cal_check_order`,
  `gd_t1_is_compatible`. Two internal helpers turned up returning **`1` for
  success** — `t1_extend` and `branch_logs`, the exact inverted-polarity trap
  the convention exists to prevent — and are now `bool` as well. `gd_sign`
  stays `int` (three-valued) and the qsort comparator stays `int` (that is its
  contract).
- **`assert()` entered the codebase**, four of them, for programmer invariants
  only: null-handle preconditions on internal statics in `gd_mcmc.c` and
  `gd_error.c`, each with a comment saying why it is an invariant and not
  input validation. `gd_diagnostics.h` stayed dependency-free — only its `.c`
  gained `gd_error.h`, so a caller of the diagnostics still needs nothing but
  `<stddef.h>`.
- **142 positional `gd_error err = {0, NULL, 0, NULL}` initialisers** became
  `GD_ERROR_INIT`. The two dead fields nothing ever read or wrote — `line` and
  `file` — went with them.
- A `GD_PRINTF_LIKE` shim gives gcc and clang format-string checking on
  `gd_error_setf`, isolated in one place so no call site sees an extension.

New test: **`tests/unit/test_error.c`, 11 cases** covering the contract itself —
that the returned code equals `err.code` across seven codes, that a NULL error
still returns the code, that a message survives the buffer it was built from
(the regression the old design could not pass), that an overlong message
truncates rather than overflows, and that a real parser agrees with its own
return value.

Verified: gcc and clang `-Werror` clean; **179/179 tests pass** under ASan
(168 before, plus the 11 new); Wasm engine vs CLI **difference 0**; 22/22 page
checks; `pair_simple` and `romance_tree` still byte-identical.

Also fixed in passing: a pre-existing `-Wmaybe-uninitialized` warning in
`tests/likelihood/test_t2_likelihood.c` (present at `1c604a0`, unrelated to
this work). The brute-force oracle allocated its state array with an unchecked
`malloc`; it is now a checked `calloc`.

The original milestone specification follows.

**Objective.** One meaning per `int` return, stated in the header.

**Current evidence.** Four conventions coexist (`GD_E_*` code, bare `0`/`1`,
`0`/`-1`, predicate). `gd_analysis_*` returns literal `1` while `gd.h` documents
"nonzero `gd_error` code" and `GD_E_NOMEM == 1`. Both current callers happen to
be correct; the contract is not.

**Concrete changes.**
1. Choose one status convention and write it into `gd_error.h`. Recommended:
   `int` returns are always a `GD_E_*` code, `0 == GD_E_OK`.
2. Convert `gd_analysis_*`'s bare `return 1` sites to return `err.code`.
3. Convert the `0`/`-1` sites (`gd_rng_serialize_impl`, `gd_rng_restore_impl`,
   `rank_normalize`, `build_grid`, `cell_double`) to the same convention.
4. Give predicates `bool` returns so they are visibly not status functions
   (`gd_error_failed`, `gd_is_finite`, `gd_tsv_is_missing`,
   `gd_cal_check_order`). `gd_sign` is three-valued and stays `int`.
5. Give `gd_error` owned message storage — a fixed inline `char msg[N]` in the
   struct, or a documented owned-string field — and retire both workarounds
   (`gd_analysis_diag.msgbuf` and the `char msg[192]` in
   `gd_model_novel_state.c`). This is the change that makes the whole error
   model honest, and it is the largest single edit in this milestone.
6. Introduce `assert()` for genuine programmer invariants (null-handle
   preconditions on internal helpers), giving `NDEBUG` a meaning. Never for
   input validation.

**Affected files.** `include/gd_error.h`, `src/gd_error.c`, `src/gd_analysis.{c,h}`,
`src/gd_rng.{c,h}`, `src/gd_diagnostics.c`, `src/gd_pairwise.c`,
`src/gd_model_novel_state.c`, `cli/commands.c`, `web/clepsydra_wasm.c`.

**Tests.** Add a unit test per error class asserting the returned code equals
`err.code` — the invariant the current split violates. Add a test that a
runtime-composed message survives past the call that set it (the `msgbuf`
problem, tested directly for the first time).

**Tooling.** `make test`; UBSan job; rebuild and re-verify the WebAssembly
engine, since `clepsydra_wasm.c` reads `err.code` into its JSON.

**Risks.** The `gd_error` message-storage change touches every module. Doing it
as a mechanical, behaviour-preserving pass, separately from any other change, is
what keeps it reviewable.

**Compatibility.** Source-breaking for anyone reading `gd_error.message`
directly. In-tree that is only the CLI and the shim.

**Acceptance criteria.** One documented convention; no bare `1`/`-1` status
returns; predicates return `bool`; error messages needing runtime data have
owned storage and no caller-supplied scratch; `assert()` used only for
invariants.

**Depends on.** C1 (the public surface should be settled before its error
contract is rewritten).

---

### C3 — Allocation discipline and checked size arithmetic — **DONE**

What changed:

- **`grep -n 'malloc(' src/*.c` now finds nothing outside `gd_memory.c`.** That
  is the rule, and it is greppable on purpose: every array allocation in the
  library goes through a checked helper that takes count and element size
  separately, so the multiplication happens where it can be checked.
- Two helpers were missing and are now there: **`gd_malloc_array_checked`**
  (the uninitialised counterpart to `gd_calloc_checked`, so a caller does not
  have to pay for zeroing to get a checked size) and
  **`gd_realloc_array_checked`** (which preserves realloc's contract: on
  failure the original block is still valid and still the caller's).
- **`gd_size_mul`** handles the counts that are themselves products. Eight
  sites needed it — `chains x draws x parameters` in the sampler, `concepts x
  nodes` in the analysis loader, `draws x parameters` in a chain's draw
  storage, `m x n` in the diagnostics — where checking only the final multiply
  would have left the first one free to wrap.
- Converted: **67** `malloc(n * sizeof(T))` sites, **12** reallocs (7 of which
  multiplied), and **20** `calloc(n, size)` calls. The callocs were already
  overflow-safe through the C library; routing them through the helper is what
  makes the rule absolute rather than something a reader has to reason about
  case by case.
- **Three sites escaped the first mechanical pass** because their count was a
  ternary — `malloc((n ? n : 1) * sizeof(T))` — and carried exactly the defect
  the milestone targets. They are why the audit is a grep for the *allocator*
  rather than for the pattern.
- Two local string-duplication helpers in `gd_analysis.c` now call
  `gd_strdup_n`, which already did the `len + 1` check they were open-coding
  without one.
- **The arena is gone** — 118 lines from `gd_memory.c` and 12 from its header.
  C1 left it with zero callers; the guide is explicit that an allocator
  abstraction has to be earned by a real lifetime problem, and this one never
  was. Every module owns its allocations through the object that outlives them,
  which is the simpler model. `gd_memory.c` went from 190 lines to 108.
- **The `atof` is gone.** A malformed `primary_weight` used to read as `0.0`,
  which silently dropped a form out of the pairwise comparison and moved the
  date with nothing said. It is now `GD_E_BAD_NUMERIC`, and the message names
  the language and concept, because a bad cell in a large forms file is
  otherwise very hard to find. A `.` in the same column still means "absent".

New tests:

- **`tests/unit/test_memory.c`, 10 cases** — that `gd_size_mul` computes and
  rejects at the boundary and leaves its output alone when it fails; that each
  checked allocator refuses a wrapping size with `GD_E_NOMEM`; that realloc
  keeps the old block on failure (getting that wrong leaks on every growth
  failure); that a zero count is not an error, since an empty dataset is a
  legitimate input. Every case is rejected by arithmetic before any allocator
  is asked for anything, so the suite stays fast.
- A regression in `test_analysis.c` for the weight: malformed is refused and
  named, `.` still means absent.

Verified: gcc and clang `-Werror` clean; **189/189 tests pass** under ASan;
Wasm engine vs CLI **difference 0**; 22/22 page checks; examples still
byte-identical.

The original milestone specification follows.

**Objective.** No allocation size derived from input can overflow silently.

**Current evidence.** ~25 `malloc(n * sizeof(T))` sites with unchecked
multiplication where `n` comes from parsed data; `gd_malloc_checked` /
`gd_calloc_checked` exist but are used in only 4 of 17 modules; the arena is
allocated from nowhere.

**Concrete changes.**
1. Route every product-sized allocation through `gd_calloc_checked` /
   `gd_malloc_checked`. Priority order by input-derived `n`:
   `gd_diagnostics.c` (9), `gd_mcmc.c` (28 allocation sites), `gd_pairwise_p1.c`
   (8), `gd_model_binary.c`, `gd_model_novel_state.c`, `gd_calibration.c`.
2. Decide the arena's fate. Either use it for the analysis-lifetime allocations
   it was built for, or delete `gd_memory`'s arena half. §25 is explicit that an
   allocator abstraction must be justified by a real lifetime model; an
   allocator with no callers is not.
3. Replace the surviving `atof` at `gd_analysis.c:876` with
   `gd_tsv_parse_double`, so a malformed form weight is an error rather than
   `0.0`.

**Affected files.** `src/gd_memory.{c,h}`, `src/gd_diagnostics.c`,
`src/gd_mcmc.c`, `src/gd_pairwise_p1.c`, `src/gd_model_binary.c`,
`src/gd_model_novel_state.c`, `src/gd_calibration.c`, `src/gd_analysis.c`.

**Tests.** A test that a deliberately huge declared count fails cleanly with
`GD_E_NOMEM` rather than allocating a wrapped-around buffer. A regression test
that a malformed weight field is rejected.

**Tooling.** `make test` under ASan; the UBSan job (its `integer` checks are
what would catch the overflow if a test reaches one).

**Risks.** Low. Mechanical and locally verifiable. The one real risk is
behavioural drift in `gd_analysis.c:876`: a dataset with a malformed weight
currently runs and silently gets `0.0`, and would start failing. That is a
correctness improvement, but it is a behaviour change and belongs in the commit
message.

**Compatibility.** Internal only, except the `atof` fix.

**Acceptance criteria.** No unchecked multiplication in an allocation size;
arena either used or gone; every numeric field parsed through the error-returning
parser.

**Depends on.** C2 (uses the unified status convention).

---

### C4 — Boolean and enumerated semantics — **DONE**

What changed:

- **32 boolean fields across seven headers became `bool`**, along with the
  locals and parameters that carry them, and their `= 1` / `= 0` assignments
  became `true` / `false` (68 lines). The three `int`s left that could be
  mistaken for booleans are deliberate and now visibly so: `gd_error.code` is a
  status, and the two T1 conflict class ids carry a `-1` sentinel.
- **`gd_cal_type` is public** (`gd_types.h`), so the enum a caller reads back
  from `gd_analysis_calibrations_read` is the one the engine wrote.
  `gd_analysis_calibration.type` is typed as it rather than as `int`.
- **`gd_analysis_pair_config.decay` (a `const char *`) became
  `gd_analysis_pair_config.model` (a `gd_pair_model`).** `GD_PAIR_P1` is
  deliberately 0, so a zeroed config asks for the model the project actually
  recommends, and the four decay laws are named members rather than strings
  the compiler cannot check.
- **`gd_pair_model_from_string` is now the single place a model name is
  spelled**, with `gd_pair_model_name` for the manifest. Both the CLI and the
  browser parse through it, so they cannot drift into accepting different
  spellings, and an unknown name is `GD_E_INVALID_ARG` rather than a silent
  fall back to P1.
- **The second parser is gone.** `gd_decay_model_from_string` had no
  production caller left once the boundary parsed to the public enum — two
  parsers with different accepted sets is exactly the drift this milestone
  exists to prevent, so it and its test were retired and the header says why
  there is no replacement.
- **`gd_cal_type_from_string` was reshaped** from `int` with a `-1` sentinel to
  a `bool` predicate with an out-parameter, and its nine-way if-chain became a
  table. clang found this one: `-Wsign-conversion` rejected assigning the
  `int` result to the now-enum field, and a cast would have papered over the
  sentinel rather than removing it. It now leaves `*out` untouched when it
  refuses, so a caller that ignores the result cannot end up holding a
  plausible-looking type.

Note on what C4 did *not* have to do: the plan called for giving
`gd_options.branch_lengths` a real enum and replacing eight stringly-typed
`gd_options` fields. C1 deleted `gd_options` entirely with the rest of the
skeleton, so that work evaporated — the only stringly-typed configuration left
in the shipped API was the one decay field.

**A bug this milestone introduced, and what caught it.** The WebAssembly
shim parses its options through a table of
`{ name, kind, void *target }` slots. Two of its entries said `OPT_INT` and
pointed at `gd_analysis_tree_config.relaxed` and `.block_age_scaling` — which
C4 had just turned into `bool`. `*(int *)slot->target = v` then wrote four
bytes into a one-byte field on every browser run, corrupting whatever the
compiler had laid out after it. The draws came out different.

Nothing in the C build could see it: the target is a `void *`, so the type is
erased and `-Werror` on both compilers had nothing to complain about. The C
test suite could not see it either — it does not go through the shim. What
found it was **`scripts/page-smoke.js`**, whose "the dataset's default settings
converge" check went from passing to reporting ESS 181. Diagnosis was by
elimination: the CLI's output was byte-identical to `HEAD`'s at the page's own
settings, so the engine was fine and the shim was not.

The fix is an `OPT_BOOL` kind, plus a comment on the slot table saying in
as many words that the compiler cannot check it and it must be updated in the
same commit as any field whose type changes. Worth recording as the shape of
the risk rather than the incident: a hand-maintained `void *` dispatch table is
the one place in this codebase where a type change is invisible to the
toolchain.

Verified: gcc and clang `-Werror` clean; **188/188 tests pass** under ASan (189
before, minus the retired duplicate-parser test); Wasm engine vs CLI
**difference 0**; **22/22 page checks**, with the browser's diagnostics
matching `HEAD`'s to the digit again; examples still byte-identical, and the
manifest still records `P1_hierarchical_constant_rate`.

The original milestone specification follows.

**Objective.** Let the compiler check what comments currently assert.

**Current evidence.** `stdbool.h` included nowhere; ~30 boolean-semantic `int`
fields; `branch_lengths` and `gd_analysis_calibration.type` are magic-value
`int`s; 9 stringly-typed configuration fields with valid values listed only in
comments.

**Concrete changes.**
1. Include `<stdbool.h>` and convert the ~30 genuinely boolean fields and
   parameters to `bool`. Only the genuinely boolean ones — §34 explicitly warns
   against mechanical conversion, and `int code`, `int type`,
   `int time_origin_year` and `gd_sign`'s result must not be touched.
2. Give `gd_options.branch_lengths` a real enum
   (`GD_BRANCH_IGNORE`/`_START`/`_DURATION`). Type
   `gd_analysis_calibration.type` as `gd_cal_type`.
3. Replace the stringly-typed `gd_options` fields and
   `gd_analysis_pair_config.decay` with enums, keeping string→enum parsing at
   the CLI and shim boundary where strings actually arrive. The parse functions
   already exist in spirit (`gd_cal_type_from_string`); this generalises them.

**Affected files.** `include/gd_types.h`, `src/gd_analysis.{c,h}`,
`src/gd_calibration.{c,h}`, `src/gd_mcmc.{c,h}`, `src/gd_tree.h`,
`src/gd_decay.{c,h}`, `cli/commands.c`, `web/clepsydra_wasm.c`.

**Tests.** Existing suite is the guard — this is a pure typing change. Add a
test that an unrecognised `--decay` or clock string is rejected with
`GD_E_INVALID_ARG` at the boundary rather than silently falling through to a
default.

**Tooling.** `-Werror` on both compilers matters most here: `-Wconversion` and
`-Wswitch` do the real work of proving the conversion is complete.

**Risks.** `bool` in a public struct changes its layout and size. Since C1 has
already established that this API breaks deliberately, sequence C4 after it.
The WebAssembly ABI is affected — the shim passes structs by pointer, so a
rebuild and re-verify is mandatory, not optional.

**Compatibility.** ABI-breaking. Bundle with C1's break if both land in one
release.

**Acceptance criteria.** `bool` used for boolean semantics; no magic-value
`int` where an enum exists; string parsing confined to the CLI and shim;
`-Wswitch` clean.

**Depends on.** C1.

---

### C5 — Buffer and truncation audit — **DONE**

What changed:

- **Every output path in the CLI goes through one checked helper.** Five
  unchecked `snprintf` calls became `out_path()`, which refuses rather than
  truncates. This was **not** cosmetic: at `1c604a0`, `date-pair --out <a
  1191-character directory>` — a legitimate path, well inside `PATH_MAX` —
  **reported success and wrote a file named `posteri`** instead of
  `posterior_pair.tsv`. Verified by running the old binary from a scratch
  worktree; the same invocation now exits 1 and writes nothing. Refusing is
  the only honest answer, since nothing in the CLI can shorten a path the user
  chose.
- The parameter-name buffer in `diagnostics.tsv` is checked too, and grew from
  128 to 512. A truncated name there is not cosmetic either: two long node
  labels can truncate to the same string, and the table would show one
  parameter measured twice rather than two parameters.
- The checksum buffer's size stopped being the literal `24` repeated at three
  sites.
- **`gd_analysis_pair_result.model` stopped being `char[64]`** and became a
  borrowed `const char *` into the static name table C4 introduced, alongside a
  `gd_pair_model model_id` for a caller that wants to branch rather than
  print. No public struct carries a fixed-size string buffer any more.
- **`GD_LEN_NUL_TERMINATED` replaced the `len == 0` overload.** A length of 0
  now means an empty buffer, which is a real input; the sentinel asks the
  library to measure. The old spelling made an empty buffer inexpressible and
  let a zeroed input struct quietly mean "call `strlen` on whatever this
  pointer reaches" — the shim relied on exactly that, and now says so. Callers
  updated: the shim (three buffers from `stringToNewUTF8`), the tests, and
  `gd_analysis_calibrations_read`. The CLI never relied on it: it reads real
  lengths from files.

New test: `test_analysis_empty_buffer_is_empty_not_measured` — that a 0 length
is honoured as empty, that the same buffer measured is a working dataset (which
is what makes the first assertion about the length rather than about the data),
that a length shorter than the buffer is not rounded up to the terminator, and
that an empty calibration table reads as zero calibrations.

Verified: gcc and clang `-Werror` clean; **189/189 tests pass** under ASan;
Wasm engine vs CLI **difference 0**; 22/22 page checks; examples still
byte-identical.

The original milestone specification follows.

**Objective.** No silent truncation; no inexpressible input.

**Current evidence.** 9 unchecked `snprintf` calls building `char path[1200]`
in the CLI; `char model[64]` in a result struct; `len == 0` overloaded to mean
"NUL-terminated" across the parse boundary.

**Concrete changes.**
1. Check every `snprintf` return in `cli/commands.c`; a truncated output path
   must fail loudly, not write to the wrong file.
2. Replace `gd_analysis_pair_result.model[64]` with a borrowed `const char *`
   pointing at a static model-name table — the values are a closed set, so a
   copy buys nothing.
3. Resolve the `len == 0` overload. Either document it as a hard convention in
   one place and keep it, or add an explicit
   `GD_LEN_NUL_TERMINATED` sentinel so that a genuinely empty buffer is
   expressible and the intent is visible at the call site. Recommended: the
   sentinel.

**Affected files.** `cli/commands.c`, `src/gd_analysis.{c,h}`, `src/gd_tsv.{c,h}`,
`src/gd_newick.{c,h}`, `web/clepsydra_wasm.c`.

**Tests.** A CLI test with an over-long `--out` path asserting a clean failure.
A parser test passing a genuinely empty buffer.

**Tooling.** `make test`; the examples job (it exercises real CLI paths).

**Risks.** Low, but the `len == 0` change touches every parse call site.

**Compatibility.** Source-compatible if the sentinel is `0`-valued by default.

**Acceptance criteria.** No unchecked `snprintf`; no fixed-size string buffer
in a public result struct; buffer-length semantics stated once and expressible
for the empty case.

**Depends on.** C1, C2.

---

### C6 — Static analysis — **DONE**

`scan-build` is not installed on this machine, but it is only a wrapper:
`clang --analyze` is the same analyzer. Together with `gcc -fanalyzer` that
gives **two** analyzers, both already part of the toolchains CI pins, so the
zero-dependency posture is intact and there is nothing to install.

- **`make analyze`**, plus `analyze-gcc` and `analyze-clang` separately. They
  find different things, and every finding below came from exactly one of them.
  Not part of `make test`: analysis is slower than compilation.
- **Coverage includes `web/clepsydra_wasm.c`**, not just `src/` and `cli/`.
  That was a deliberate addition — the shim is the one file no C test executes,
  and it is what runs in a visitor's browser. It is also where two of the six
  findings were.
- **A CI job** gates on it, separately from the `-Werror` build because it
  answers a different question: warnings are about a line, the analyzers are
  about a path.
- **`docs/static-analysis.md`** records the baseline, each finding, and the
  suppression policy.

**The gate is zero findings with no suppressions anywhere in the tree.** That
is only a defensible gate because the baseline was small enough to clear
completely — six findings, four fixed on their merits and two the same analyzer
limitation, resolved by writing the code more plainly rather than by silencing
the tool.

**Two real defects, both on paths no test reaches** — which is the whole
argument for running these at all, given a 189-assertion sanitizer-clean suite:

1. **A null dereference on an out-of-memory cleanup path** in
   `gd_analysis.c`. `class_names` and `class_count` are allocated and checked
   together; if `class_count` was the one that failed, the cleanup ran with
   `class_names` live and `class_count` NULL and then indexed it.
2. **A leak on a partial input load** in `cli/commands.c`. `load_tree_inputs`
   read three files and returned failure without freeing the ones that had
   loaded, while its callers return straight out without calling
   `free_tree_inputs`. Its contract is now absolute: on failure it owns nothing.

The other two:

3. **A dead store** in the checkpoint reader — a true positive, but deleting
   the increment would have been the wrong fix. It is now
   `assert((size_t)(q - buf) == body)`: the invariant the symmetric walk was
   implicitly asserting, and one that would catch a serializer/deserializer
   layout disagreement instead of restoring a chain of shifted values.
4. **Two uninitialised-read reports in the JSON writer**, one root cause, both
   analyzer limitations rather than defects — and both fixed rather than
   suppressed, because in each case the code the analyzer could not follow did
   not need to be written that way. `jb_printf`'s `char tmp[512]` is now
   initialised at its declaration (one store, against a suppression that would
   also have hidden a genuine uninitialised read later), and while there,
   truncation became a failure rather than something to emit — a clipped number
   is valid JSON carrying the wrong value. `jb_str`'s two-element string for
   every ordinary character became a `jb_putc`.

A note on the tooling itself: **neither analyzer exits non-zero on a finding**,
so the targets grep their output. The first version of `analyze-gcc` did not,
and cheerfully reported "clean" while gcc was printing a warning — recorded in
`docs/static-analysis.md` because it is the failure mode any future analyzer
target will have too.

Verified: gcc and clang `-Werror` clean; `make analyze` clean; **189/189 tests
pass** under ASan; Wasm engine vs CLI **difference 0**; 22/22 page checks;
examples still byte-identical.

The original milestone specification follows.

**Objective.** Catch path-sensitive defects the compiler and sanitizers miss.

**Current evidence.** No static analyser is configured or referenced anywhere in
the repository.

**Concrete changes.**
1. Add `scan-build` (Clang Static Analyzer) and gcc's `-fanalyzer` as separate
   `make` targets. Both are already-installed toolchain features, which fits the
   project's zero-dependency posture.
2. Run both, record the full baseline in `docs/`, and triage per §20 — do not
   gate CI on day one.
3. Fix high-confidence findings. Expect overlap with C3's allocation work,
   which is a reason to sequence C6 after it.
4. Add a CI job once the baseline is clean, with any suppression narrow and
   commented at the site.

**Affected files.** `Makefile`, `.github/workflows/ci.yml`, a new
`docs/static-analysis.md`, plus whatever the findings touch.

**Tests.** None new; this is a tooling milestone.

**Risks.** False positives on the numerical code, particularly around
log-space guards where the analyser cannot know a value is finite. Triage
honestly rather than suppressing broadly.

**Compatibility.** None.

**Acceptance criteria.** Analyser invocation is a documented `make` target;
baseline is recorded; high-confidence findings fixed; suppressions narrow and
documented.

**Depends on.** C3.

---

### C7 — Fuzz the parsers — **DONE**

Three libFuzzer harnesses in `fuzz/`, a `make fuzz` / `make fuzz-run` pair, a
CI job, committed corpora (~630 files, 91 KB), and `docs/fuzzing.md`. All
three call the production parsers directly. They were twenty lines each rather
than a rewrite, because the parsers already took `(buf, len)` with
`path == NULL` — **the separation the WebAssembly build needed is the same one
fuzzing needs**, which is the guide's §22 paying off years after the fact.

`fuzz_analysis` is the one that mattered. The individual parsers can each be
correct while an index derived from one and used against another is not: the
forms table is interned against tip labels, and calibrations are resolved onto
node indices. It splits one byte string into three inputs with a length prefix,
so the fuzzer mutates the boundaries as well as the contents.

**Two flags turned out to be load-bearing, and one of them nearly slipped
through.** `-fsanitize=...,undefined` alone lets UBSan *print* a runtime error
while the fuzzer carries on and exits 0 — a finding would scroll past a green
run. I caught it by building a harness with a deliberate signed overflow and
watching it "pass"; the fix is `-fno-sanitize-recover`, which the CI UBSan job
already set. This is the third time in this plan that a green check turned out
not to be checking (see also C6's `analyze` target and C4's shim). `make
fuzz-run` was then verified to fail on a real finding — a harness with a
deliberate heap overflow — rather than assumed to.

**The finding: `a,(b);` — a forest read as a tree.**

Not a memory-safety slip in the parser but a **missing structural invariant**.
`a,(b);` has two top-level nodes and no enclosing parentheses; the parser
accepted it and returned a structure whose second node had
`parent == GD_NODE_NONE`. Downstream, **five** places read a node's parent
whenever the node is not the root — the calibration propagator, the calibration
sampler, both likelihoods, and the confounding diagnostic — and each indexed an
array at `SIZE_MAX`. UBSan caught the pointer arithmetic in one of them; the
other four were the same bug waiting for a different input.

Fixed at the boundary *and* at the sites, deliberately:
`gd_tree_build_traversals` now rejects a disconnected tree (every consumer
calls it after construction, whether parsed or hand-built), and all five
indexing sites test the parent rather than testing "is this the root", because
`gd_tree` is public and a caller can hand one over.

Per the plan, the reproducer became **an ordinary test** —
`test_newick_rejects_a_forest` — rather than only a corpus entry, and it also
asserts the rooted forms still parse, so the check rejects forests rather than
commas.

Roughly **16 million executions** in total; the ~11 million after the fix found
nothing further. The thin crash surface is not a surprise: C3 audited the
allocations and C5 the truncations, and the parsers never used `strlen` on
caller data. The value here is mostly ahead — a gate on future parser changes.

Verified: gcc and clang `-Werror` clean; `make analyze` clean; **190/190 tests
pass** under ASan; Wasm engine vs CLI **difference 0**; 22/22 page checks;
examples still byte-identical.

The original milestone specification follows.

**Objective.** Exercise the untrusted-input boundary automatically.

**Current evidence.** `gd_tsv_parse` and `gd_newick_parse` consume arbitrary
text and already take `(buf, len)` — the §21 shape. The browser page lets a
visitor paste arbitrary Newick into a page running this C in their own tab. No
fuzzing exists.

**Concrete changes.**
1. Add `fuzz/fuzz_tsv.c` and `fuzz/fuzz_newick.c` as libFuzzer harnesses
   calling production parsers directly — no duplicated parsing logic, per §12's
   acceptance criteria.
2. Add a third harness over `gd_analysis_tree_prepare`, the whole
   dataset-ingestion path, once C1 has settled where it lives.
3. Seed corpora from `tests/fixtures/` and the committed `examples/` datasets.
4. Add a `make fuzz` target building with `-fsanitize=fuzzer,address,undefined`.
5. Run a short bounded fuzz in CI (seconds, not hours) so regressions surface;
   longer runs stay manual.
6. Any crash found becomes a committed regression test under
   `tests/regression/`.

**Affected files.** New `fuzz/`, `Makefile`, `.github/workflows/ci.yml`,
`tests/fixtures/`.

**Tests.** The harnesses themselves, plus whatever regression cases fall out.

**Tooling.** clang with libFuzzer; ASan and UBSan together.

**Risks.** The parsers will very likely produce findings — that is the point,
but it means C7's scope is genuinely open-ended. Timebox the triage and let
unfixed findings be recorded rather than blocking.

**Compatibility.** None.

**Acceptance criteria.** Both parsers have harnesses over production code;
seed corpora committed; a bounded CI run exists; discovered crashes are
regression tests.

**Depends on.** C3, C5 (fuzzing before the allocation and buffer work is done
mostly rediscovers what C3 and C5 already fix).

---

### C8 — Contributor contract — **DONE**

- **`LICENSE`** — MIT, `Copyright (c) 2026 Tiago Tresoldi`. Not a free choice
  on my part: the sibling projects this one is explicitly modelled on
  (`merkmal`, `regulae`) carry the identical text under the same author, so the
  house default was already established and the repository was simply missing
  it. It ships CC-BY-4.0 data whose terms are stated in
  `examples/DATA_SOURCE.md`; the code's own were not stated anywhere.
- **`STYLE.md`** — the conventions C2–C7 settled, in one place: the error
  model, the ownership vocabulary, the allocation rule, the buffer and length
  rules, the type rules, the `goto done` shape, and the library's no-I/O
  property. It says what the rule is; this plan says why.
- **`make debug` and `make release`** (the guide's §17 gap). Each builds into
  its own object directory, so switching cannot mix objects built under
  different assertion policies. Verified: all three configurations produce
  identical numbers, and `NDEBUG` genuinely removes the assertions —
  `__assert_fail` is present in the default binary and absent from the release
  one. That check is what gives C2's `assert()` real meaning: the release
  configuration is where an assertion with a side effect would go missing, and
  a CI job now builds it on both compilers.
- **One `goto` label.** There were three — `done` (76), `cleanup` (10),
  `out_fail` (10). The library was already uniform on `done` and only
  `cli/commands.c` mixed the other two, in separate functions, so the
  difference had never meant anything. All `done` now.
- **The reading order** in `CONTRIBUTING.md` starts at `STYLE.md`, and `README.md`
  lists the full check set.

A note for whoever reads `STYLE.md` first: its "verify it actually fails"
section is the one rule this plan produced by being wrong three times.
`make analyze` reported clean while gcc printed warnings; the fuzz targets let
UBSan report and still exited 0; and the shim's option table silently corrupted
memory because no compiler could see through its `void *`. All three were
green. The habit that catches them is breaking the check on purpose once and
watching it go red.

Verified: gcc and clang `-Werror` clean; `make release` builds on both;
`make analyze` clean; `make fuzz-run` clean; **190/190 tests pass** under ASan;
Wasm engine vs CLI **difference 0**; 22/22 page checks; examples still
byte-identical.

The original milestone specification follows.

**Objective.** Put the conventions in writing where the guide expects them.

**Current evidence.** No `LICENSE`, no `STYLE.md`, no `CONTRIBUTING.md`. The
project notes' conventions section is doing this job informally, and it is good
— but it is interleaved with milestone history, and the data licences in
`examples/DATA_SOURCE.md` (CC-BY-4.0, attribution required) have no
corresponding statement for the code itself.

**Concrete changes.**
1. Add `LICENSE`. The repository is public and ships CC-BY-4.0 data; the code's
   own terms are currently unstated.
2. Add `STYLE.md`: the C99 house style, the error convention settled in C2, the
   ownership vocabulary (`_create`/`_destroy`/`_result_destroy`), the allocation
   rule from C3, the buffer-length convention from C5, the `goto` label
   vocabulary (currently three names — `done` ×74, `cleanup` ×9, `out_fail` ×6 —
   pick one).
3. Document the build configurations: debug, release, ASan, UBSan, analyser,
   fuzz. Add the missing `-O0 -g` debug and `-O2 -DNDEBUG` release targets to
   the `Makefile` (§17).
5. State the tests' privileged access to internal headers as policy.

**Affected files.** New `LICENSE`, new `STYLE.md`, `Makefile`, `README.md`,
`CONTRIBUTING.md`.

**Tests.** None.

**Tooling.** None.

**Risks.** None.

**Compatibility.** None.

**Acceptance criteria.** A new contributor can determine from the repository
alone how to build, test, sanitize, analyse and fuzz; where the public API
lives; how ownership and errors work; and under what licence.

**Depends on.** C2, C3, C5 (it documents their outcomes).

---

## Deferred work

Attractive, and deliberately **not** in this plan:

- **Replacing Make with CMake or Meson.** §2.1 says prefer the existing system
  when it works. This one is 120 lines, has no generated sources, drives a
  6-configuration CI matrix, and cross-compiles to WebAssembly through a
  separate script. Changing it buys nothing this plan needs.
- **Splitting `gd_analysis.c` (999 lines) or `gd_mcmc.c` (1,069 lines).** Both
  are large but internally sectioned and single-purpose. `gd_analysis.c`'s
  private TSV-table layer is a plausible extraction, but §24 is explicit: do not
  promote a helper to a shared module until a second module needs it. Nothing
  else does.
- **Making value structs opaque.** `gd_analysis_tree_result`,
  `gd_pair_result` and the summary structs are value-like by design and both
  consumers read their fields directly. §34 warns specifically against forcing
  opacity here.
- **ThreadSanitizer / MemorySanitizer.** The library is single-threaded and has
  no uninitialised-read history under ASan and UBSan. §19 warns against enabling
  sanitizer combinations blindly.

---

## First recommended milestone

**C1 — Resolve the two-API split.** *(Done — see the C1 section for what landed.)*

It is the smallest change with the largest effect on the question the guide
cares most about: *which functions form the public API?* Right now the
repository answers that question wrongly, in writing, in a file called
`gd.h`. Every subsequent milestone — the error contract, the type changes, the
fuzz targets, the style guide — has to name a public surface, and each one
would otherwise re-litigate this.

It is also the milestone with the least behavioural risk: no numerical code
changes, and the WebAssembly-vs-CLI equivalence check plus the passing suite
make a regression loud.

One decision is genuinely the maintainer's and should be settled before the work
starts: **is `gd.h`'s frozen surface a commitment the project still intends to
honour, or is it superseded scaffolding?** If the `gd_context`/`gd_chain` API is
still the intended shape for v1 — with `gd_analysis` a temporary M10 expedient —
then C1 means implementing those 12 stubs on top of `gd_analysis`, and the
milestone is much larger. If `gd_analysis` is where the library actually landed,
C1 means promoting it and retiring the skeleton, which is a week's careful work
at most.

The assessment above supports the second reading — `gd_analysis.h` is better
documented, buffer-only, browser-proven and consumed by everything that ships —
but the intent behind `PLAN.md`'s frozen API is not something the code can
settle.

**Resolved 2026-08-20: the second reading.** `gd_analysis.h` was promoted and
the skeleton retired, and C2 followed it: the error model is unified and
`gd_error` now owns its message, which removed the last workaround the old
design forced on callers. C3 followed: every array allocation in the library
now goes through a checked helper, and the arena and the `atof` are gone. C4
made booleans `bool` and closed sets enums; C5 audited the buffers and found a
CLI that reported success while writing to a truncated filename. C6 put two
analyzers behind a zero-findings CI gate and cleared the six findings they
raised. C7 fuzzed the three ingestion entry points and found a forest being
read as a tree. C8 wrote the contract down. **The plan is complete.**


---

## Outcome

All eight milestones are done. Measured against the guide's §37 checklist,
every box is now true except the ones this plan deliberately deferred:
benchmarks (§37 asks for them before optimisation; there is no optimisation
work, and `make benchmark` is still a stub), and removing the last filesystem
I/O from `gd_tsv` and `gd_newick`.

### What changed, in numbers

| | before | after |
|---|---|---|
| Exported symbols | 162 | 132 |
| Declared in `include/` | 38 | 25 |
| Functions returning `GD_E_NOT_IMPLEMENTED` | 12 | 0 |
| Consumers needing `-Isrc` | 2 | 0 |
| `int` return conventions in use | 4 | 1 |
| Raw allocator calls in `src/` outside `gd_memory.c` | 128 | 0 |
| Unchecked `malloc(n * sizeof(T))` | ~66 | 0 |
| `stdbool.h` users | 0 | 12 headers |
| Static analyzers | 0 | 2, gating CI at zero |
| Fuzz harnesses | 0 | 3, gating CI |
| Tests | 176 | 190 |

Net about −400 lines, with three new documents, a licence and a style guide.

### The eight defects this found

The refactor was supposed to be behaviour-preserving, and every example output
is still byte-for-byte identical. It nonetheless turned up eight real defects,
none of which the 176-assertion sanitizer-clean suite reached:

1. **`atof` reading a malformed weight as `0.0`** (C3) — silently dropping a
   form from a pairwise comparison and moving the date.
2. **~66 unchecked allocation-size multiplications** (C3) on counts taken from
   parsed input.
3. **A type-confused option table in the shim** (C4) — writing four bytes into
   a one-byte field on every browser run, changing the draws.
4. **A CLI reporting success while writing to a truncated filename** (C5) —
   `posteri` instead of `posterior_pair.tsv`, exit code 0.
5. **A null dereference on an out-of-memory cleanup path** (C6).
6. **A leak on a partial input load** (C6).
7. **A dead store hiding an unstated checkpoint-layout invariant** (C6).
8. **A forest parsed as a tree** (C7) — five places then indexing an array at
   `SIZE_MAX`.

The pattern is worth naming: **not one of them was in the numerical code**.
The likelihoods, the sampler and the diagnostics — the parts that are hard, and
that the test suite was built around — came through untouched. Every defect was
at a boundary: reading input, building a path, passing an option, cleaning up
after a failure. That is where the tests were not looking, and it is what the
guide's emphasis on explicit contracts at boundaries is for.

### What to do next

The **Deferred work** section above still stands, less one item. Of what remains,
the one worth doing first:

- **Benchmarks**, because `make benchmark` is a stub and MCMC runtime is real.
  Not to optimise — to have a number before anyone tries.

Two more of the items above have since been done, in `docs/deepening-plan.md`:
**benchmarks** in E1 — `make benchmark` is real and
`docs/benchmark-baseline.md` states how small an effect it can resolve — and
**breaking `gd_mcmc`'s compile-time dependency on T1** in E3, whose stated
trigger ("revisit if a second likelihood needs the sampler") had been met since
M7 by a fully validated T2 that could not reach it. The performance objection
recorded against it was answered rather than waived: E1 came first, and E3
reports its measurement.

**Removing the path branch from `gd_tsv` and `gd_newick`** — listed here as the
other first candidate — was done in **D1** of `docs/deepening-plan.md`. The
worry recorded above, that it meant rewriting every test's fixture loading,
turned out not to apply: every fixture was already a C string literal, all 42
call sites already passed `NULL`, and the change was mechanical.

**A habit worth keeping.** Four times in this plan a green check turned out not
to be checking:

- `make analyze` reported clean while gcc printed warnings (its exit status was
  never going to be non-zero);
- the fuzz targets let UBSan report a finding and still exited 0;
- the shim's `void *` option table corrupted memory invisibly to every
  compiler;
- and **the behaviour-preservation check used during C1–C8 was vacuous.**
  Example outputs were verified with `git status --porcelain examples/`, but
  `examples/*/output/` is in `.gitignore` and no example output is tracked, so
  that command could never have reported a difference. It was re-done properly
  before the commit — building `HEAD` in a scratch worktree, running every
  example under both binaries, and diffing: **13 output files across all seven
  examples, byte-for-byte identical.** The conclusion held; the evidence for it
  had not.

The last one is the instructive one, because the check was mine and I ran it
eight times. `STYLE.md` says to break a check on purpose and watch it go red.
The cheaper version, which would have caught this one: ask what the check would
print if the thing it tests were broken.

**Verifying behaviour preservation, properly.** For a change that should not
move the numbers:

```sh
git worktree add -q --detach /tmp/base <base-commit>
(cd /tmp/base && make native && for e in $(find examples -name run.sh); do sh "$e"; done)
make native && for e in $(find examples -name run.sh); do sh "$e"; done

# Every output file, however deeply nested, and say how many.
find examples -path '*/output/*' -type f | sort > /tmp/outputs.list
n=0; bad=0
while read -r f; do
    n=$((n + 1))
    cmp -s "/tmp/base/$f" "$f" || { bad=$((bad + 1)); echo "DIFFERS: $f"; }
done < /tmp/outputs.list
echo "compared $n files, $bad differ"

# The generated files under examples/ that ARE tracked must stay clean too.
git status --porcelain examples/
```

Why it is spelled this way, from what the D-series ran into:

- **`find`, not `examples/*/run.sh`.** That glob misses
  `examples/romance_tree/calib_sensitivity/run.sh`, which is why CI spells two
  globs. `find` needs no second guess.
- **`find` and `-type f`, not `examples/*/output/*`.** The old glob does work:
  nine of the 22 outputs are one level deeper, under
  `families_pairwise/<family>/output/`, and the glob yields those *directories*,
  which `diff` then compares entry by entry. That is a load-bearing accident.
  A reader who tightens the loop to regular files — which looks like a strict
  improvement — silently drops those nine, and that is exactly what happened
  during D1: the run compared 13 files, reported no differences, and read like
  a pass. Enumerating files explicitly removes the trap.
- **Print the count.** A comparison that matched nothing at all also prints no
  differences. `compared 22 files, 0 differ` is a result; silence is not, and
  the count is what makes the previous trap visible instead of invisible.

And run it against the commit you are comparing to, not `HEAD` — `HEAD` is the
commit you are about to change.
