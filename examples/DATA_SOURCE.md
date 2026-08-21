# Example data: sources, licences, and regeneration

The cognate data under `examples/*/data/` are small subsets of open
**Lexibank / LexiCore** datasets — Lexibank's core collection of expert
cognate-coded wordlists — converted to Clepsydra's strict-TSV input schema
(`docs/input-formats.md`). Two families are used:

- **Indo-European** (`examples/pair_simple`, `romance_tree`,
  `five_ie_reproduction`) from **IE-CoR**, via `tools/datasets/iecor.py`;
- **Sino-Tibetan** (`examples/sinotibetan_tree`) from **Sagart et al. 2019**, via
  the family-agnostic `tools/datasets/lexibank.py`.

The generic converter handles any Lexibank CLDF Wordlist with a `cognates.csv`
(IE-CoR, Sagart Sino-Tibetan, ABVD Austronesian, Semitic, Bantu, …), so more
families plug in the same way.

`examples/families_pairwise/` adds a **nine-family pairwise gallery** built by
`tools/datasets/build_gallery.py` from these LexiCore datasets (each CC-BY-4.0,
cited in `examples/families_pairwise/pairs.tsv`):

| family | dataset | citation |
|---|---|---|
| Dravidian | `lexibank/dravlex` | Kolipakam et al. 2018 |
| Turkic | `lexibank/savelyevturkic` | Savelyev & List 2020 |
| Uralic | `lexibank/uralex` | Syrjänen et al. 2018 |
| Semitic | `lexibank/felekesemitic` | Feleke 2021 |
| Polynesian | `lexibank/walworthpolynesian` | Walworth 2018 |
| Tupían | `lexibank/tuled` | Gerardi et al. 2021 |
| Arawakan | `lexibank/chaconarawakan` | Chacon et al. |
| Alor-Pantar | `lexibank/robinsonap` | Robinson & Holton 2012 |
| Koreanic | `lexibank/leekoreanic` | Lee 2015 |

Regenerate the gallery (clone each dataset under one root, then):

```
cd tools
python3 -m datasets.build_gallery --root /path/with/dataset/clones \
    --out ../examples/families_pairwise
```

## Attribution (required by the licence)

> Heggarty, Paul, Cormac Anderson & Matthew Scarborough. *Indo-European Cognate
> Relationships database (IE-CoR)*, version 1.0, 2019.
> https://github.com/lexibank/iecor

Licensed **CC-BY-4.0** (https://creativecommons.org/licenses/by/4.0/). The CLDF
release used here is `lexibank/iecor` (derived dataset). The subsets committed in
this repository are redistributed under the same CC-BY-4.0 licence, with the
attribution above. IE-CoR underpins Heggarty et al. 2023, *Science* 381:eabg0818.

> Sagart, Laurent, Guillaume Jacques, Yunfan Lai, Robin J. Ryder, Valentin
> Thouzeau, Simon J. Greenhill & Johann-Mattis List. *Sino-Tibetan Database of
> Lexical Cognates*, 2019. https://github.com/lexibank/sagartst — the data behind
> Sagart et al. 2019, *PNAS* 116(21):10317–10322. Licensed **CC-BY-4.0**.

Only the converted subsets are committed (a few hundred to ~2000 forms each); the
full CLDF releases are **not** vendored here.

## Regenerating the subsets

```
# 1. Obtain the IE-CoR CLDF release (either clone the repo or fetch the CSVs):
git clone https://github.com/lexibank/iecor /tmp/iecor
#    (or download cldf/languages.csv, parameters.csv, forms.csv, cognates.csv)

# 2. Regenerate each example's data with the converter (from tools/):
cd tools
python3 -m datasets.iecor --cldf /tmp/iecor/cldf \
    --languages Latin Italian Spanish Portuguese French Romanian Catalan SardinianLogudoro \
    --out ../examples/romance_tree/data
python3 -m datasets.iecor --cldf /tmp/iecor/cldf \
    --languages VedicEarly AncientGreek Latin OldIrish OldChurchSlavonic \
    --out ../examples/five_ie_reproduction/data
python3 -m datasets.iecor --cldf /tmp/iecor/cldf \
    --languages Spanish Portuguese \
    --out ../examples/pair_simple/data

# Sino-Tibetan (Sagart et al. 2019), via the generic Lexibank converter:
git clone https://github.com/lexibank/sagartst /tmp/sagartst
python3 -m datasets.lexibank --cldf /tmp/sagartst/cldf \
    --languages Sinitic_Old_Chinese Sinitic_Beijing Tibetan_Old_Tibetan Tibetan_Lhasa \
                Burmish_Old_Burmese Burmish_Rangoon Kiranti_Limbu rGyalrong_Japhug \
    --out ../examples/sinotibetan_tree/data

# Polynesian (Walworth 2018):
git clone https://github.com/lexibank/walworthpolynesian /tmp/walworthpolynesian
python3 -m datasets.lexibank --cldf /tmp/walworthpolynesian/cldf \
    --languages Tongan Niuean Samoan EastFutuna Rapanui Hawaiian Maori Tahitian \
    --out ../examples/polynesian_tree/data
```

The converter is deterministic, so regeneration reproduces the committed files
byte-for-byte for a given IE-CoR release.

## Bringing your own cognate data (non-CLDF)

If your cognate data is not in CLDF — e.g. a hand-built cognate matrix from a
phylogenetics project — use `tools/datasets/matrix.py` instead of the Lexibank
converters. It reads two plain layouts into the same Clepsydra TSV:

```
# wide: one row per language, one column per meaning, cells = cognate-class labels
python3 -m datasets.matrix --format wide --input cognates.csv --out mydata \
    --family hmong_mien [--synonym-sep ';'] [--missing NA X]

# long: a tidy table with language / concept / cognate columns
python3 -m datasets.matrix --format long --input forms.csv --out mydata \
    --lang-col Language --concept-col Concept --cognate-col Cognateset
```

Missing cells (`?`, `-`, `.`, `NA`, empty) become no form; multiple labels in a
wide cell can be split into weighted synonym candidates. The output plugs
straight into `clepsydra date-pair` and `date-tree`. (A Nexus/BEAST *binary*
presence-absence matrix is a different encoding and is not read by this tool.)

## What the Lexibank converters do

- Restricts IE-CoR to the chosen doculects (matched by `ascii_name`).
- Emits `languages.tsv`, `concepts.tsv`, `forms.tsv`; each form carries its
  `cognate_id` (the IE-CoR cognate set, namespaced by concept), a
  `primary_weight` split over synonyms in a cell, and `loan_probability` /
  `cognacy_confidence` derived from IE-CoR's loan and doubt flags.
- Cognate-set identity is used only for equality within a concept, exactly what
  the T1 likelihood needs.

The trees (`tree.nwk`) and calibrations (`calibrations.tsv`) in each example are
**illustrative** fixed inputs authored for the demonstration, not authoritative
phylogenetic or chronological claims — see each example's README and
`docs/scientific-methods.md`.
