# Benchmark baseline

`make benchmark` was `@echo "Benchmark target (M2+)"` from M2 until E1. This
records what it measures, what it measured on the machine below, and — the
part that matters — **how small an effect it can actually resolve.**

It exists because E3 puts an indirection in the sampler's inner loop, and
`docs/c99-alignment-plan.md` §Deferred refuses that change "with no benchmark
behind it". A number without a stated resolution would not have satisfied that
objection; it would only have looked like it did.

---

## What is measured

Iterations per second through `gd_analysis_tree_job_run`, on two committed
example datasets, seeded so repeats differ only in timing.

**CPU seconds, not wall clock.** `clock()` is C99, needs no POSIX, and reports
processor time — which is the question. A busy scheduler makes wall clock
noisy without making the inner loop slower.

**Iterations, not likelihood calls.** An earlier draft counted calls by
wrapping the model callback. It was dropped: wrapping the callback inserts
exactly the kind of indirection E3 exists to measure, so the instrument would
have changed the reading. Calls per iteration are fixed for a given
configuration, so iterations per second moves with them — and it is measured on
the path that ships.

**The first two runs are discarded.** Without that the numbers are the
processor's frequency ramp, not the code:

| discarded runs | romance spread across 3 repeats |
|---|---|
| 0 | **46.7%** (890 → 1239 → 1469 iters/s, monotonically rising) |
| 1 | 14.1% |
| 2 | **1.0%** |

Whichever dataset runs first pays the ramp. With none discarded, the reported
spread was larger than any effect worth looking for.

---

## The baseline

```
AMD Ryzen 7 4700U (8 threads), Linux 7.1.3-100.fc43.x86_64
gcc 15.2.1 20260123 (Red Hat 15.2.1-7), -O2 -g
make benchmark   # 1500 warmup + 1500 samples x 2 chains, seed 42, 3 repeats
```

| dataset | median iters/s | within-run spread |
|---|---|---|
| romance | ~1460 | 0.4–1.1% |
| five_ie | ~2290 | 0.2–0.4% |

`five_ie` is faster per iteration than `romance` because it carries fewer
concepts (863 form rows against 1375), not because it is a smaller tree.

---

## How small an effect this can resolve

Within one invocation the spread is 0.2–1.1%. **Across invocations it is
larger**: three consecutive unmodified runs gave medians of

- romance: 1474.3, 1455.7, 1473.5 → **1.3%** apart
- five_ie: 2291.0, 2240.7, 2288.7 → **2.2%** apart

So the number to trust is the between-invocation figure, not the tidier spread
printed at the end of a run.

**A single before/after comparison resolves an effect of roughly 3% or more.
Below that, repeat the comparison several times or do not claim it.**

This was checked rather than assumed: a deliberate per-call cost added to
`t1_model_loglik` moved romance by −1.4% and five_ie by −2.3% — visible, but
close enough to the between-run drift that one comparison would not have
settled it. That is the resolution, stated honestly.

---

## Comparing a change

```sh
make benchmark 2>&1 | tee /tmp/before.txt     # on the base commit
# ... make the change ...
make benchmark 2>&1 | tee /tmp/after.txt
```

Compare the medians, not the individual repeats. If the difference is under
about 3%, run both again before drawing a conclusion — and if it stays under
that, the honest report is "no effect this benchmark can resolve", not a
number.

To match a specific question, override the workload:

```sh
make benchmark BENCH_WARMUP=4000 BENCH_SAMPLES=4000 BENCH_CHAINS=4
```

`make benchmark` is deliberately not part of `make test`: it takes minutes, and
a slow suite is a suite people stop running.
