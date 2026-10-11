# Indo-European tip-dating (reproduction-style, T1)

Dates a fixed 5-language Indo-European topology from **real IE-CoR cognate data**
(170 meanings) using five *attested ancient* languages whose fixed tip ages
anchor the chronology and, jointly with the cognate data, inform the clock rate
(tip-dating). Reproduction-*style*: it illustrates the method on real data, not a
re-run of any specific published analysis.

```
./run.sh          # ~25 s
```

## Inputs

- `data/` — IE-CoR subset: Vedic (early), Ancient Greek, Latin, Old Irish, Old
  Church Slavonic. CC-BY-4.0; see `../DATA_SOURCE.md`.
- `tree.nwk` — an **illustrative** fixed topology (one hypothesis of the
  inter-branch structure). v1 does not infer topology.
- `calibrations.tsv` — fixed tip ages from the languages' approximate attestation
  (Vedic 3200, Greek 2400, Latin 2050, Old Irish 1250, OCS 1050 BP; illustrative),
  and a broad uniform prior on the PIE root.

## Representative output (seed 42, 4 chains × 6000)

| node | median BP | 90% HPD | lexical age (90% HPD) | R̂ | ESS |
|---|---|---|---|---|---|
| proto_indo_european | ~6150 | 4421–8436 | 1.70 (1.30–2.15) | 1.00 | 1432 |
| core_ie | ~6020 | 4332–8339 | 1.66 (1.27–2.10) | 1.00 | 1433 |
| hellenic_italic | ~5830 | 4184–8050 | 1.60 (1.21–2.03) | 1.00 | 1449 |
| italo_celtic | ~5400 | 3785–7468 | 1.48 (1.12–1.88) | 1.00 | 1493 |

The PIE-root posterior is **wide and rate-prior-sensitive**. Its median sits
near the Steppe range (~5500–6000 BP), but its 90% HPD reaches the lower end
of the Anatolian one (~8000–9500 BP), so five ancient tips cannot decide
between them. That is the honest scientific situation for deep-time lexical dating,
and the point of the example: the tool reports the uncertainty and the
diagnostics rather than a false point date. Changing `--rate-logmean/-logsd`
shifts the root age materially, which the run manifest documents.

Since the proper age prior (lexical root plus uniform ages below it) and the
joint moves of `docs/inference.md`, this deep tree mixes well: every node has
ESS above 1,400. Earlier versions reached ~170.

## Prior only against posterior

![Five-IE (ancient tips): 90% HPD of each node under the prior alone and under the posterior](prior_overlay.svg)

The grey bar is the node's age with the cognate likelihood switched off (the same run's `prior_*` columns: tree prior, calibrations and clock prior only). The blue bar is the posterior. `prior_overlap` in `posterior_nodes.tsv` measures how much the two share.

**The root is mostly prior.** Its overlap is 0.72: the posterior, 4,421–8,436 BP,
is the prior's 3,804–9,039 BP trimmed at both ends. The width of the root is
therefore the width of the clock and root priors, not a finding from the
cognates. The data do act on the inner nodes, which they push much older than
the prior's uniform ratios put them (Italo-Celtic 2,835 → 5,397 BP, overlap
0.20): the ancient tips are lexically far apart, so their common ancestors must
be old.

Regenerate after `./run.sh` with
`python3 ../../tools/plotting/overlay.py output/posterior_nodes.tsv --svg prior_overlay.svg --title "Five-IE (ancient tips)"`, which also prints the numbers as a table.
