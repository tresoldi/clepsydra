# Cross-family pairwise gallery (P1)

Runs `clepsydra date-pair` for a representative language pair in **nine language
families**, on real cognate data from Lexibank's LexiCore collection. Pairwise
(P1) divergence needs no tree and no calibration, so it works for any family —
this is the breadth demonstration.

```
./run.sh          # ~1 s; writes output/<family>/posterior_pair.tsv for each
```

## Families and pairs (`pairs.tsv`)

| family | pair | source |
|---|---|---|
| Dravidian | Betta Kurumba – Kodava | Kolipakam et al. 2018 (DravLEX) |
| Turkic | Nogai – Yakut | Savelyev & List 2020 |
| Uralic | Komi-Permyak – Udmurt | Syrjänen et al. 2018 (UraLex) |
| Semitic | Amharic – Mesqan | Feleke 2021 |
| Polynesian | Hawaiian – Samoan | Walworth 2018 |
| Tupían | Parintintin – Wayampi | Gerardi et al. 2021 (TuLeD) |
| Arawakan | Achagua – Baniwa | Chacon et al. |
| Alor-Pantar (Papuan) | Kaera – Teiwa | Robinson & Holton 2012 |
| Koreanic | Gyeonggi – Middle Korean | Lee 2015 |

Each pair was selected automatically by `tools/datasets/build_gallery.py` to be
**well-attested and moderately diverged** (many shared concepts, retention near
0.55) — except Koreanic, where the informative comparison is a modern variety
against its attested ancestor (Middle Korean, retention ~0.94). All are
CC-BY-4.0; see `../DATA_SOURCE.md` for full citations and regeneration.

## What the numbers mean (and a caveat)

`run.sh` prints shared concepts, retention fraction, and the P1 divergence
median. The divergences are **lexical time scaled by the rate prior** (per
kiloyear); `absolute_anchored = 0` in every `output/<family>/posterior_pair.tsv`,
because a bare pair carries no tip-age anchor. Do not read these as calendar
dates — that is the honest reporting the tool enforces.

Note too that the retention-matched pairs (~0.55) land at similar divergence
medians (~2000 BP) *by construction*: they were chosen for comparable divergence,
and under a fixed rate prior comparable retention implies comparable lexical
time. The point of the gallery is **breadth** — the same engine and data pipeline
run unchanged across nine unrelated families — not a cross-family date comparison
(which would require a shared observable space and per-family calibration).

Swap in any two languages from a family's `forms.tsv` with `--lang-a/--lang-b`,
or point `date-tree` at a curated subset for a full chronology (see the Romance,
5-IE, and Sino-Tibetan tree examples).
