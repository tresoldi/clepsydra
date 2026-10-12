"""Blind dating benchmark on GLED and ASJP (roadmap: LOCO expansion).

For every calibrated node of `reference/calibrations/curated_nodes.tsv` that
the 40-concept GLED list covers, clepsydra dates the node from GLED cognates
with no calibration at all, and the ASJP method of Holman et al. (2011) dates
the same node from ASJP word lists with the published constants. Both are
scored against the node's independently known age. A second table compares
family roots with published Bayesian estimates.

Modules: `data` (GLED, ASJP, Glottolog, language sampling), `targets` (node
sets and the reference table), `asjp` (LDND dating), `runner` (clepsydra runs
and scoring). See `docs/research/gled-benchmark.md`.
"""
