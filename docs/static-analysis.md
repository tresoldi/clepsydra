# Static analysis

Two analyzers run over the C sources, both already part of the toolchains CI
pins. Neither adds a dependency; there is no separate tool to install.

```
make analyze          # both
make analyze-gcc      # gcc -fanalyzer
make analyze-clang    # clang --analyze (the Clang Static Analyzer)
```

They are separate targets because they find different things: gcc's is stronger
on allocation and lifetime, clang's on path-sensitive null and dead-value
reasoning. Every finding below was found by exactly one of them.

Coverage is `src/*.c`, `cli/*.c` **and `web/clepsydra_wasm.c`**. The shim is
included deliberately: it is the one file a C test never executes, and it is
what runs in a visitor's browser.

Neither tool exits non-zero on a finding, so the targets check the output
rather than the status. A target that reported "clean" while the tool printed a
warning would be worse than no target at all — which is not hypothetical, it is
what the first version of `analyze-gcc` did.

`make analyze` is **not** part of `make test`: analysis is slower than
compilation and belongs in CI and in a deliberate local run.

---

## Baseline

Established 2026-08-21 at the C6 milestone (`docs/c99-alignment-plan.md`).

| Analyzer | Findings at baseline | Now |
|---|---|---|
| `gcc -fanalyzer` | 2 (both in the shim, one root cause) | 0 |
| `clang --analyze` | 4 | 0 |

**The gate is zero, with no suppressions anywhere in the tree.** That is only a
reasonable gate because the baseline was small enough to clear completely — six
findings, four of which were fixed on their merits and two of which were the
same analyzer limitation, resolved by writing the code more plainly rather than
by silencing the tool.

If a future finding is a genuine false positive that cannot be written around,
suppress it **at the site, with a comment saying why**, and record it here.
There is nothing to record yet.

---

## What the baseline found

Two real defects, both on paths no test reaches. This is the argument for
running the analyzers at all: the suite is 189 assertions and sanitizer-clean,
and it exercises neither of these.

### 1. Null dereference on an out-of-memory cleanup path

`src/gd_analysis.c`, reported by `clang --analyze` as
`core.NullDereference`.

`class_names` and `class_count` are allocated together and checked together.
If `class_names` succeeded and `class_count` was the allocation that failed,
the `goto done` cleanup ran with `class_names` live and `class_count` NULL, and
then read `class_count[c]`.

Fixed by treating an absent `class_count` as a count of zero, which is correct:
nothing has been added to the per-class arrays at the point that allocation
could fail.

### 2. Leak on a partial input load

`cli/commands.c`, reported by `clang --analyze` as `unix.Malloc` (twice — one
site, two pointers).

`load_tree_inputs` read three files and returned `-1` if any was missing,
without freeing the ones that had loaded. Callers return straight out on that
failure without calling `free_tree_inputs`, so a tree and a forms file leaked
whenever the calibrations file was unreadable.

Fixed by making the contract absolute: on failure the function owns nothing.

### 3. Dead store in the checkpoint reader

`src/gd_mcmc.c`, reported by `clang --analyze` as `deadcode.DeadStores`.

The last `q += 8` in `gd_mcmc_chain_restore` advanced a cursor nothing read
afterwards. A true positive, but deleting the increment would have been the
wrong fix: the increment is what makes the walk symmetric, and its real purpose
is an invariant nobody had written down.

It is now `assert((size_t)(q - buf) == body)` — the cursor must land exactly
where `ckpt_body_size` predicted. A serializer and a deserializer that
disagreed about the layout would otherwise restore a chain full of shifted
values and carry on sampling from it.

### 4. Uninitialised read in the JSON writer (analyzer limitation)

`web/clepsydra_wasm.c`, reported by `gcc -fanalyzer` as
`-Wanalyzer-use-of-uninitialized-value`, through two call chains with the same
root.

Neither report was a defect, but both were fixed rather than suppressed,
because in both cases the code the analyzer could not follow was also code that
did not need to be written that way:

- `jb_printf` passed a `char tmp[512]` to `jb_puts` after `vsnprintf`. The
  analyzer cannot model `vsnprintf` writing to the buffer. `tmp` is now
  initialised at its declaration — one store on an error-free path, against a
  suppression that would also have hidden a genuine uninitialised read later.
  While there, truncation became a failure rather than something to emit: a
  clipped number is valid JSON carrying the wrong value.
- `jb_str` built a two-element string for every ordinary character just to
  append one byte. That is now `jb_putc`, which is plainer code and leaves
  nothing for the analyzer to misread.

---

## Notes on interpreting findings here

The numerical code draws false positives that need judgement rather than a
reflex fix. The likelihood works in log space with explicit `-HUGE_VAL`
sentinels and guards that an analyzer cannot always prove reachable or
unreachable. Before "fixing" a finding in `gd_model_*.c`, `gd_math.c` or
`gd_pairwise*.c`, check whether the path it describes can actually occur given
the guards upstream — and if it cannot, prefer restructuring so that is
evident over asserting it in a comment.
