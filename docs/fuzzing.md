# Fuzzing

Three libFuzzer harnesses over the code that consumes untrusted bytes.

```
make fuzz                              # build them
make fuzz-run                          # short bounded run over the corpora (CI)
make fuzz-run FUZZ_SECONDS=300         # longer

build/fuzz_newick fuzz/corpus/newick   # open-ended, locally
```

clang only — libFuzzer is a clang feature. ASan, UBSan and the integer checks
are all on.

## Why these entry points

The browser page lets a visitor paste a Newick tree into a textarea and upload
a forms file, and that text reaches the C parsers unaltered. It is the
shortest path from a stranger's keyboard to this code.

| Harness | Target | What it is for |
|---|---|---|
| `fuzz_newick` | `gd_newick_parse` | The tree a visitor pastes. |
| `fuzz_tsv` | `gd_tsv_parse` + the field helpers | Every dataset the project reads is TSV. |
| `fuzz_analysis` | `gd_analysis_tree_prepare`, `gd_analysis_calibrations_read` | Where the three inputs stop being independent. |

Each harness calls the production parser **directly**. None pre-validates its
input, because a harness that cleans its input first only proves the cleaner
works.

`fuzz_analysis` is the one worth explaining. The individual parsers can each be
correct while an index derived from one and used against another is not: the
forms table is interned against the tree's tip labels, and the calibrations are
resolved onto node indices. That cross-referencing is what it exercises. It
stops before `gd_analysis_tree_job_run` — fuzzing the sampler would measure the
RNG, not the parsing.

It splits one byte string into three inputs with a length prefix: the first
four bytes are two little-endian 16-bit lengths for the tree and the forms
table, and the rest is the calibration table. The fuzzer therefore mutates the
boundaries as well as the contents.

Nothing here touches the filesystem. All three parsers already took
`(buf, len)` with `path == NULL` meaning "use the buffer", which is why the
harnesses are twenty lines each rather than a rewrite — the separation the
WebAssembly build needed is the same one fuzzing needs.

## Two flags that are not optional

**`-fno-sanitize-recover`.** Without it UBSan *prints* a runtime error and the
fuzzer carries on and exits 0, so a finding scrolls past a green run. The first
version of these targets had exactly that hole; it was caught by building a
harness with a deliberate signed overflow and watching it "pass". The CI UBSan
job already set the flag, and the fuzz targets now match it.

**`-fno-sanitize=unsigned-integer-overflow,unsigned-shift-base`**, for the same
reason CI excludes them: the PCG32 stream relies on unsigned wrap deliberately.

The same care applies to `make fuzz-run` itself, which is verified to exit
non-zero on a finding — checked by adding a harness with a deliberate heap
overflow and confirming the target failed and stopped, rather than trusting
that it would.

## Corpora

`fuzz/corpus/<harness>/` is committed. It holds hand-written seeds — real
trees and tables from `web/data/romance/`, plus the shapes a parser has to
handle (a polytomy, a bare `;`, an empty file, a header-only table, comments
and `.` markers) — merged with a coverage-minimised set from an actual fuzzing
session.

About 510 files and 96 KB in total. `du` will claim ten times that: the units
are tens of bytes each and every one occupies a filesystem block.

`make fuzz-run` writes its discoveries to `build/fuzz-scratch/` (gitignored),
not into the committed corpus — libFuzzer writes new units to the *first*
corpus argument, and a run that quietly adds a few hundred files to the
repository is not a check, it is a mess. The committed corpus is passed second,
read-only, so a run still starts from real parser states.

To grow the committed corpus deliberately, fuzz into a scratch directory and
merge back:

```
build/fuzz_newick /tmp/grow fuzz/corpus/newick -max_total_time=600
build/fuzz_newick -merge=1 -max_len=1024 /tmp/min /tmp/grow fuzz/corpus/newick
cp /tmp/min/* fuzz/corpus/newick/
```

**Check the merge did something before deleting anything.** `-merge=1` against
a missing or stale binary is a silent no-op, and the obvious next step —
clearing the old units first — then throws the corpus away.

## When something is found

libFuzzer writes the reproducer to `build/crash-*` or `build/leak-*`.

1. Confirm it reproduces: `build/fuzz_newick build/crash-<hash>`.
2. Fix the defect.
3. **Turn the input into an ordinary test**, next to the code it broke, so it
   is checked by `make test` rather than only by a fuzzer that happens to
   rediscover it. The corpus is not a substitute for a named test that says
   what the bug was.
4. Add the input to the corpus as well.

## Findings

### `a,(b);` — a forest read as a tree (2026-08-21)

Found by `fuzz_analysis` after about seven minutes, and the reason the flags
above matter: it is a UBSan report, so with `-fsanitize-recover` left at its
default it would have printed and the run would have exited 0.

`a,(b);` is not a tree. It has two top-level nodes and no enclosing
parentheses. The Newick parser accepted it and returned a structure with a
second node whose `parent` was `GD_NODE_NONE`.

Downstream, **five** places read a node's parent whenever the node is not the
root — in the calibration propagator, the calibration sampler, both
likelihoods, and the confounding diagnostic. Each of them indexed an array at
`SIZE_MAX`. UBSan caught the pointer arithmetic in
`gd_calibrations_propagate`; the other four were the same bug waiting for a
different input.

Fixed in two places, deliberately:

- **`gd_tree_build_traversals` now rejects a disconnected tree** — exactly one
  node may lack a parent, and it must be the root. That is the boundary: every
  consumer calls it after construction, whether the tree was parsed or
  assembled by hand in a test.
- **All five indexing sites test the parent** rather than testing "is this the
  root". `gd_tree` is public, so a caller can hand one over; these do not
  depend on the check above having run.

Regression test: `test_newick_rejects_a_forest` in `tests/unit/test_newick.c`,
which also asserts the rooted forms are still accepted — so the check rejects
forests rather than commas. The reproducer is in both corpora.

## Status

At the C7 baseline (2026-08-21), roughly **16 million executions** across the
three harnesses under ASan and UBSan produced **one finding**, above. The
~11 million executions after the fix produced nothing further.

Read the rest for what it is: the parsers were already buffer-and-length based,
had no `strlen` on caller data, and had just been through C3's allocation audit
and C5's truncation audit, so the crash-and-leak surface was thin going in. The
one real finding was not a memory-safety slip in the parser but a **missing
structural invariant** — the parser was safe with its bytes and wrong about its
output — and that is the kind of thing only a fuzzer driving the whole
ingestion path was going to reach.
