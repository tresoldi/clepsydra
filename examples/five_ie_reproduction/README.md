# Indo-European tip-dating (reproduction-style, T1)

Dates a fixed 5-language Indo-European topology from **real IE-CoR cognate data**
(170 meanings) using five *attested ancient* languages whose fixed tip ages
anchor the chronology and, jointly with the cognate data, inform the clock rate
(tip-dating). Reproduction-*style*: it illustrates the method on real data, not a
re-run of any specific published analysis.

```
./run.sh          # ~30 s
```

## Inputs

- `data/` — IE-CoR subset: Vedic (early), Ancient Greek, Latin, Old Irish, Old
  Church Slavonic. CC-BY-4.0; see `../DATA_SOURCE.md`.
- `tree.nwk` — an **illustrative** fixed topology (one hypothesis of the
  inter-branch structure). v1 does not infer topology.
- `calibrations.tsv` — fixed tip ages from the languages' approximate attestation
  (Vedic 3200, Greek 2400, Latin 2050, Old Irish 1250, OCS 1050 BP; illustrative),
  and a broad uniform prior on the PIE root.

## Representative output (seed 42, 4 chains × 6000, block age-scaling on)

| node | median BP | 90% HPD | R̂ | ESS |
|---|---|---|---|---|
| proto_indo_european | ~7870 | 5994–9999 | 1.02 | 168 |

The PIE-root posterior is **wide and rate-prior-sensitive** — its 90% HPD spans
the gap between the Steppe (~5500–6000 BP) and Anatolian (~8000–9500 BP)
hypotheses. That is the honest scientific situation for deep-time lexical dating,
and the point of the example: the tool reports the uncertainty and the
diagnostics rather than a false point date. Changing `--rate-logmean/-logsd`
shifts the root age materially, which the run manifest documents.

This deep tree mixes poorly under single-site moves; the CLI enables the block
age-scaling move (`docs/inference.md` §2.2), which lifts ESS from ~13 to ~170.
