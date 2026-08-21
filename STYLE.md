# Style

The conventions this codebase actually follows. Most were settled during the
C99 alignment work; `docs/c99-alignment-plan.md` records why each one is what
it is, and is worth reading before proposing a change to any of them.

The short version: **boring, explicit C99**. A maintainer should be able to look
at a function and see what goes in, what comes out, who owns each object, what
may fail, and how failure is reported — without reading anything else.

---

## Language and toolchain

Strict **C99**. C standard library and libm only — no third-party dependency in
`src/`, and none in the Python toolchain either. No compiler extensions except
`GD_PRINTF_LIKE` in `gd_error.h`, which is guarded on `__GNUC__` and isolated so
no call site sees it.

Everything compiles under this warning set, on gcc **and** clang:

```
-std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
-Wstrict-prototypes -Wmissing-prototypes
```

Not `-Werror` by default — a newer compiler should not break a user's build. CI
turns warnings into errors instead, where the compiler version is pinned.

### Configurations

| | Command | Optimisation | Assertions |
|---|---|---|---|
| default | `make` | `-O2 -g` | on |
| debug | `make debug` | `-O0 -g` | on |
| release | `make release` | `-O2 -DNDEBUG` | **off** |
| sanitized | `make test` | `-O2 -g` + ASan | on |

The default build keeps assertions: an assertion that never runs is not a
check. `make release` exists so the one configuration that removes them gets
built and exercised occasionally — it is where an assertion with side effects
would go missing. There are none, and that is how it stays true.

### The full check set

```
make            # build
make test       # 251 assertions, ASan-clean
make analyze    # gcc -fanalyzer + clang --analyze; gate is zero findings
make fuzz-run   # libFuzzer over the parsers, ASan + UBSan
make shim-check # the WebAssembly adapter compiles natively
make release    # the NDEBUG configuration still builds
make benchmark  # iterations/sec on two examples; NOT part of make test
```

See `docs/static-analysis.md`, `docs/fuzzing.md` and
`docs/benchmark-baseline.md` — the last of which states how small an effect
the benchmark can resolve (about 3%), because a performance number without
that is not evidence.

---

## Modules

A `.c` and its header form a module with one responsibility. There is no
`utils.c`, and there should not be: if a helper is used by one module, keep it
`static` in that module.

**`include/` is the public API and `src/` is not.** `gd.h` is an umbrella over
`gd_analysis.h` (the pipelines), `gd_tree.h` (the topology a result is reported
against, read-only), `gd_types.h` and `gd_error.h`. The CLI and the WebAssembly
shim build against `include/` alone — no `-Isrc` — and the build rules enforce
it. If either ever needs a header from `src/`, the fix is to make that header
public, not to add the flag back.

**Tests are the one consumer allowed into `src/`.** That is deliberate: a unit
test exists to exercise a module directly, not only through the pipelines that
happen to call it.

Dependencies point inward and the graph is acyclic:

```
gd_types → gd_error → gd_memory → {gd_tree, gd_tsv, gd_math, gd_rng}
        → {gd_newick, gd_calibration, gd_forms, gd_grid}
        → gd_mcmc_model → {models, gd_mcmc} → gd_analysis
        → {gd_sensitivity, gd_report}
gd_options sits beside gd_analysis: it describes that module's configs.
```

`gd_mcmc_model.h` is the seam between the sampler and a likelihood, and it is
its own header for a reason: **both sides include it and neither includes the
other.** The sampler must not know which likelihoods exist, and a likelihood
must not depend on the sampler. Until E3 those types lived in `gd_mcmc.h`,
which then had to include `gd_model_novel_state.h` to declare the one adapter
— and that is most of why the interface had exactly one implementation.

The D-series added five: `gd_forms` (one reading of the forms table),
`gd_grid` (one posterior-grid summariser), `gd_options` (one option registry),
`gd_report` (the reported tables) and `gd_sensitivity` (the calibration
verdict). `docs/deepening-plan.md` records why each exists.

---

## Errors

**One convention, no exceptions.** The full statement is at the top of
`include/gd_error.h`; the rules:

- A function returning `int` returns a **`GD_E_*` code**, and that code equals
  `err->code`. Not "nonzero means it went wrong" — the value *is* the
  classification.
- Never a bare `1` or `-1`. `1` is `GD_E_NOMEM`, and `-1` is not a code at all.
- A function that **answers a question** returns `bool`, so it is visibly not a
  status. `gd_error_failed`, `gd_tsv_is_missing`, `gd_cal_check_order` are
  predicates; `gd_tsv_parse` is not. Watch the polarity: a predicate's `true`
  is usually the opposite of a status's zero.
- A function returning a **pointer** fails with `NULL` and puts the
  classification in `err`.
- `gd_error_set` and `gd_error_setf` **return the code**, so a failure branch
  reads `return gd_error_set(err, GD_E_PARSE, "...")` and cannot drift from
  what it recorded.
- Passing `NULL` for the `gd_error *` is always allowed. The code is still
  returned; only the detail is dropped.

**The error owns its message.** `gd_error_setf` formats a message naming a
runtime value straight into the error, so nothing needs a caller-supplied
buffer. `err.message` is an array: never `NULL`, `""` when unset — test it with
`err.message[0]`, not against `NULL`.

Initialise with `gd_error err = GD_ERROR_INIT;`, never a positional brace list.

**`assert()` is for programmer invariants only** — conditions a caller cannot
produce through the API. Bad data from a file, a user or a browser is always a
returned code. The library never calls `exit()` or `abort()`.

---

## Memory

**No raw `malloc`/`calloc`/`realloc` in `src/` outside `gd_memory.c`.**
`grep -n 'malloc(' src/*.c` should find nothing else.

Every array allocation goes through `gd_malloc_array_checked`,
`gd_calloc_checked` or `gd_realloc_array_checked`, which take **count and
element size separately** so the multiplication is checked. The counts come
from parsed input — concepts, nodes, draws — and `malloc(n * sizeof(T))` is the
shape that wraps silently. When the count is itself a product (chains × draws ×
parameters), multiply it out with `gd_size_mul` first.

The strategy is plain heap ownership: each allocation is owned by the object
that outlives it and released in that object's `_destroy`. **No arena, no
allocator callbacks** — there was an arena, and it was removed because nothing
allocated from it. An allocator abstraction has to be earned by a real lifetime
problem.

### Ownership vocabulary

| | |
|---|---|
| `_create` / `_destroy` | an owned handle; destroy exactly once |
| `_prepare` / `_run` | a job is prepared, then run once |
| `_result_destroy` | frees a result and everything it owns |
| `_read` / `_free` | an enumerated array and its matching release |
| `const T *` | borrowed; valid only as long as its owner is |

Anything returned as `const char *` inside a result — a node label, a model
name — is borrowed from that result and dies with it. Say so in the header.

Failure paths free what they allocated: **on failure a function owns nothing.**

---

## Buffers and lengths

A `*_len` is a byte count. `GD_LEN_NUL_TERMINATED` asks the library to measure
it. **`0` means an empty buffer**, which is a real input.

Fields from the TSV parser are **not** NUL-terminated in place — use the
length. No `strlen` on caller data.

**Reported tables go through `gd_report`** (`include/gd_report.h`), not through
`fprintf` in a command body. A report builds a `gd_table`; a renderer turns it
into TSV or JSON and hands the bytes to a caller-supplied sink — a `FILE *` on
the command line, a JSON buffer in the shim, memory in the tests. That is what
makes the output schema testable at all: before D5 the column names lived
inside `cli/commands.c`, which has no C test, and CI could only check that a
file had at least one row.

A column carries **two** names, its TSV spelling and its JSON one, because the
two schemas differ and always did. Do not unify them to tidy the struct: the
JSON is a published payload the page reads by key.

**Check every `snprintf` that builds a path or a name.** It cannot overflow,
but a truncated path names a *different file*: `date-pair` once reported
success while writing `posteri` instead of `posterior_pair.tsv`. Output paths
go through `out_path()` in `cli/commands.c`. No fixed-size string buffer
belongs in a public struct.

---

## Types

`bool` for boolean semantics; `size_t` for sizes, counts and indices;
fixed-width integers where the width is semantic (seeds, PCG32 state, file
formats). Do **not** mechanically convert every `int` — the ones left that look
boolean are deliberate: `gd_error.code` is a status, and the T1 conflict class
ids carry a `-1` sentinel.

A user-facing name from a closed set is parsed to an **enum at the boundary**.
`gd_pair_model_from_string` is the single place a pairwise-model name is
spelled, so the CLI and the browser cannot drift into accepting different ones,
and an unknown name is `GD_E_INVALID_ARG` rather than a silent default.

`const` wherever mutation is not intended. No pointer typedefs.

---

## Control flow

`goto done` for multi-resource cleanup. **One label name — `done`** — for the
whole tree; there were three, and the difference never meant anything. The
shape:

```c
int rc = GD_E_INTERNAL;       /* resolved from err at `done` */
...
if (thing == NULL) { rc = gd_error_set(err, GD_E_NOMEM, "..."); goto done; }
...
rc = GD_E_OK;
done:
    if (rc != GD_E_OK) rc = fail_code(err);
    /* release in reverse order of acquisition */
    return rc;
```

`goto` for anything other than cleanup: no.

**Iterative tree traversal only.** No deep recursion — the browser stack is
256 KB, and a deep tree from a user is not a hypothetical.

Minimal preprocessor. Prefer a `static inline` function to a function-like
macro.

---

## The library does no I/O

`gd_analysis` works from buffers, never paths. The caller does its own I/O and
its own formatting; the library does the science. That is why the page and the
CLI cannot compute different numbers — and, as it turned out, why the fuzz
harnesses are twenty lines each rather than a rewrite.

`gd_tsv` and `gd_newick` are buffer-only too, as of D1. There is no exception
left: `grep -rn 'fopen(' src/` finds nothing, and no `FILE *` appears anywhere
in `src/`. Two modules still include `<stdio.h>` — `gd_error.c` for `vsnprintf`
and `gd_analysis.c` for `snprintf` — which format into a caller-owned buffer
and touch no stream.

Nothing in `src/` calls `printf` or `exit`.

---

## Tests

`make test` must stay at 100%. Tests go next to the thing they test: unit tests
per module, likelihood tests against independent brute force, simulation tests
for SBC and recovery, regression tests for pipelines.

**A bug gets a named test.** A fuzz reproducer goes in the corpus *and* becomes
an ordinary test — the corpus records that an input once broke something; only
a named test records what the bug was.

When adding a check that can fail, **verify it actually fails.** Five separate
times now, a green check turned out not to be checking. Four during the
alignment work: a
`make analyze` target that reported clean while gcc printed warnings; fuzz
targets where UBSan reported and the run still exited 0; a shim option table
that silently corrupted memory because no compiler could see through its
`void *`; and a behaviour-preservation check run with
`git status --porcelain examples/` when `examples/*/output/` is gitignored, so
it could never have shown a difference.

And once during the D-series, in a subtler way. The documented loop globs
`examples/*/output/*` and diffs each entry; nine of the 22 outputs sit one
level deeper, so the glob yields their *directories* and `diff` compares those
entry by entry. It works, by accident. Tightening the loop to regular files —
which looks like a strict improvement — drops those nine silently, and that is
what a D1 run did: 13 files compared, no differences, indistinguishable from a
pass. Enumerating with `find` and printing the count removes both the accident
and the trap.

Break it on purpose once, watch it go red, then fix it. The cheaper version:
ask what the check would print if the thing it tests were broken.

**Example outputs are not tracked** — `examples/*/output/` is gitignored. A few
*generated* files under `examples/` are tracked, though, and the example run
regenerates them: `romance_tree/calib_sensitivity/{loco,widen}.tsv` and the
LOCO benchmark reports. Those show drift as a dirty tree, so
`git status --porcelain examples/` is a real check for them and a vacuous one
for everything else.

To show a change did not move the numbers, build the base commit in a
worktree, run the examples under both binaries and compare every output —
see the end of `docs/c99-alignment-plan.md` for the exact commands. Use them
rather than writing the loop from memory: `examples/*/run.sh` really does miss
`calib_sensitivity/run.sh`, and the output loop only reaches the nine deeper
files under `families_pairwise/<family>/output/` because `diff` happens to
descend into the directories the glob yields. Print the file count, too — a
comparison that matched nothing also reports no differences.

---

## Things to leave alone

- **The version string is load-bearing output.** `src/gd_version.c` holds it,
  and every `run_manifest.tsv` records it, so changing it moves example
  outputs. That is not a reason never to change it — 0.1.0-m1 became 0.2.0 for
  the first public release — but it is a reason to change it deliberately, in
  its own commit, and to check that the manifests are the only thing that
  moved.
- The example outputs. They are byte-for-byte reproducible and that is the
  point; if a change moves them, that is the finding.
- `web/clepsydra.{js,wasm}` are committed artifacts. Rebuild with
  `web/build-wasm.sh` and re-verify after touching `src/`, `include/` or the
  shim — `scripts/wasm-provenance.sh --check` is what stops a stale engine
  shipping.
- **The run options live in one table**, `src/gd_options.c`, behind
  `include/gd_options.h`. The CLI, the WebAssembly shim and the page all read
  it, so an option is named, bounded and typed once. Add an option by adding a
  row; do not add a parser. Defaults are deliberately *not* in the table — they
  stay in `gd_analysis_*_config_default`, and a surface that wants to show one
  reads it back with `gd_option_format`.

  Each row pairs `offsetof` with its C type and asserts the two agree, using
  the C99 negative-array-size trick. This replaced a `void *` table the
  compiler could not check, which was wrong once: two fields became `bool`
  while their slots still said `int`, so the parser wrote four bytes into one
  and corrupted what followed. That is now a build error.
