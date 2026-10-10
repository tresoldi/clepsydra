"""M1/M2 at pilot scale: the D2 clock pilot redone with the full engine.

The D2 pilot (`tools/clock_pilot`, `docs/research/clock-pilot.md`) estimated
each calibrated node's log clock rate from pairwise single-rate distances.
This package keeps its nodes, datasets, meta-analysis and LOFO, and replaces
the per-node estimator: each node is dated, without any calendar information,
by `clepsydra date-tree` on a small target-like subtree (roadmap D5), and its
lexical age L (expected replacements per concept) gives beta = ln L - ln T.

Pre-registration: `docs/research/clock-engine/plan.md`. Report:
`docs/research/clock-engine.md`.

Usage (deterministic; hours, almost all of it the sequential dating runs)::

    cd tools
    python3 -m clock_engine.run --clepsydra PATH/TO/clepsydra   # everything
    python3 -m clock_engine.run --stage build     # select nodes, write subtrees
    python3 -m clock_engine.run --stage date      # run clepsydra, one at a time
    python3 -m clock_engine.run --stage analyse   # meta-analysis and LOFO

Dating runs whose outputs exist are not repeated (`--force` redoes them).
Every clepsydra run is wrapped in `systemd-run --user --scope` with a 1 GB
memory cap (`--no-scope` disables the wrapper).

Modules: `subtree` (doculect choice, pruned Glottolog topology, input files),
`dating` (clepsydra runs, convergence handling, reading the target node),
`analysis` (meta-analysis, censored fit, LOFO, comparison with the pilot),
`run` (orchestration).
"""
