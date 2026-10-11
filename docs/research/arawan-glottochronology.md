# Arawan: clepsydra against classical lexicostatistical dating

How old is Proto-Arawan by the classical methods, at their best, and how do
they compare with the Bayesian result (`arawan-result.md`)? The same KAHD
cognate data (9 lists, 253 concepts) and the same tip dates are used
throughout: Arua at 131 BP, Paumari 1833 at 165–170 BP.

    python3 tools/glottochronology/dating.py \
        --forms forms.tsv --tree tree_labeled.nwk \
        --tip-age arawa=131 --tip-age paumari_1833=167.5 \
        --calibrate paumari_1833 paumari 167.5 \
        --upgma paumari,suruwaha,deni,kulina,banawa,jamamadi,jarawara,arawa

## Methods

| method | list | rate | formula for two lists sharing a share C of cognates |
|---|---|---|---|
| Swadesh (1955) | 100 items (79 in KAHD) | r = 0.86 per millennium | t = ln C / (2 ln r) |
| Lees (1953) | 200 items (120 in KAHD) | r = 0.805 per millennium | t = ln C / (2 ln r) |
| Starostin (2000) | 100 items | λ = 0.05 per millennium per lineage | C = exp(−2λ·C·t²), so t = √(ln C / (−2λC)) |
| Starostin, single λ | 100 items | λ = 0.05 | t = √(ln C / (−λC)), the form often quoted |
| self-calibrated | all 253 | fitted on Paumari 1833 vs Paumari | t = ln C / (2 ln r) |
| ASJP (Holman et al. 2011) | 40 items, forms | published | 1,764 BP for "Arauan", ±29% |

Starostin's method has two further steps that cannot be applied here.
Loans are to be removed before counting, and KAHD carries no loan coding.
The 110-item list adds ten items that KAHD does not have. Keeping loans
lowers C and makes Starostin's dates older. Starostin's formula is taken from
Vasilyev's formalisation (Swadesh Centenary Conference, MPI-EVA 2009): per
lineage N(t) = e^(−λ·N·t²), and for two daughters the exponent doubles.

A node is dated from its cross-clade comparisons: matches and compared
concepts are summed over every pair of tips on either side, the classical
practice. A dated tip has been separate for fewer years, which the formulas
correct for (exactly for Swadesh and Lees, through the mean separation for
Starostin). Intervals are 5–95% of a bootstrap over concepts (2,000
replicates). They cover the sampling error in C only, not the error in the
rate constant, which the classical methods treat as known.

## The check that can be made: Paumari 1833

The 1833 Paumari list and modern Paumari are the same language 167 years
apart. They share 45 of the 68 concepts both record (0.66).

| method | Paumari 1833 / Paumari split, BP | true |
|---|---|---|
| Lees | 1,018 (645–1,519) | ≈ 167 |
| Swadesh | 1,380 (801–2,153) | |
| Starostin | 2,488 (1,754–3,417) | |
| clepsydra (`arawan-result.md` §1, n7) | 496 (166–1,393) | |
| clepsydra, date held out (`arawan-result.md` §4) | list placed at 260 (0–751) | 165–170 |

Every classical method puts the split six to fifteen times too early. A
rate fitted to this pair, the one in-family calibration, is r = 0.085 per
millennium; on it the whole family is 220 years old. The 1833 list differs
from modern Paumari mostly through error in the list, not through time.
The classical methods have no term for that and read it as time.
Clepsydra's list-error prior (D1.3, posterior error ≈ 0.19) absorbs most
of it. It still covers the true value in its 90% interval.

## A control with a known answer: Romance

The same code on the Romance example (IE-CoR, 8 lists, 170 concepts, Latin
at 2,050 BP). Latin is close to the ancestor of the modern languages, so the
root (Latin against the rest) should be near 2,050 BP, and Proto-Romance
around 1,500–1,800 BP. The concept IDs are not mapped to the Swadesh lists,
so every method uses all 170 concepts (a list mismatch for the 100-item
constants).

| method | root (≈ 2,050) | Proto-Romance (≈ 1,500–1,800) | Spanish–Portuguese |
|---|---|---|---|
| Lees | 2,151 (1,978–2,363) | 835 (681–993) | 399 |
| Swadesh | 2,645 | 1,201 | 573 |
| Starostin | 3,847 | 2,281 | 1,434 |
| clepsydra, default | 2,580 (2,282–2,822) | 1,662 (1,383–1,959) | 377 |

Lees recovers the root and puts Proto-Romance 700–1,000 years too young.
Starostin is too old throughout, and places Proto-Romance before Latin.
Clepsydra's Proto-Romance is in the expected range; its root is somewhat
old, and Spanish–Portuguese is too young under every method.

## Arawan

### On the MCC tree

Node dates in BP, with bootstrap 5–95% for the classical methods and 90%
HPD for clepsydra.

| node | Lees | Swadesh | Starostin | clepsydra, tip dates | clepsydra, empirical clock |
|---|---|---|---|---|---|
| root (Paumari vs rest) | 1,749 (1,448–2,111) | 2,453 (1,934–3,062) | 3,926 (3,240–4,871) | 3,092 (671–7,755) | 2,235 (563–5,171) |
| n6 (Suruwaha vs rest) | 2,056 | 3,050 | 4,804 | 2,856 (592–7,197) | 2,072 (521–4,799) |
| n5 | 1,051 | 1,490 | 2,653 | 1,676 (324–4,230) | 1,212 (317–2,826) |
| n4 Deni–Kulina | 808 | 982 | 1,996 | 834 (153–2,125) | 598 (144–1,431) |
| n2 Madi | 190 | 156 | 703 | 185 (28–493) | 133 (26–328) |
| n1 Jamamadi–Jarawara | 103 | 90 | 527 | 124 (18–344) | 90 (15–226) |

Every classical method dates n6 older than the root: pooled over pairs,
Suruwaha is lexically farther from the rest than Paumari is. On this tree
that is impossible. The classical dates are not a chronology of the tree.

### On each method's own tree

Built by UPGMA on the method's pairwise dates over the eight lists with
enough data (Arua included, the 1833 list left out):

| method | own-tree root, BP (bootstrap 5–95%) | outgroup |
|---|---|---|
| Lees | 1,962 (1,574–2,428) | Suruwaha |
| Swadesh | 2,940 (2,292–3,887) | Suruwaha |
| Starostin | 4,653 (3,656–6,214) | Suruwaha |
| Starostin, single λ | 6,576 (5,203–8,714) | Suruwaha |
| ASJP (published) | 1,764 (≈ 1,250–2,280) | — |

All four put Suruwaha, not Paumari, as the first split; the BEAST sample
gives the MCC grouping of everything but Paumari 70% support.

## What the comparison shows

1. **The classical methods do not agree with each other.** The same
   cognates give a root between 1,750 and 6,600 BP depending on the list and
   the formula: a factor of nearly four, from choices the data cannot make.
   Their bootstrap intervals are narrow (±20%) because they leave the rate
   out. With the rate's real spread across families (σ ≈ 0.77 on the log
   scale, `clock-pilot.md`) each would be as wide as clepsydra's.
2. **They fail the only test the data offer.** Paumari 1833 comes out
   1,000–2,500 years old instead of 167. Clepsydra, with a list-error term,
   covers the truth.
3. **They do not give a chronology of the tree.** On the tree they date a
   child older than its parent.
4. **Where they are right, clepsydra agrees and carries the uncertainty.**
   The shallow nodes (Madi, Jamamadi–Jarawara, Deni–Kulina) come out at the
   same ages under Lees, Swadesh and clepsydra. At the root, clepsydra's
   tip-dated 90% interval (671–7,755) contains every classical estimate, and
   its empirical-clock median (2,235) sits between Lees and
   Swadesh.

The best classical estimate is Lees on the 200-item list. It is the only
variant that recovers the Romance root, though it puts Proto-Romance far too
young. On Arawan it gives about 1,750–1,960 BP, close to the ASJP 1,764. Clepsydra's result is not a different number so much as an
honest width: the classical point estimates are inside its interval, and
the width is what the two historical lists can actually support.
