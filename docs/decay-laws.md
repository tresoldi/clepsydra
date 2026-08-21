# Alternative pairwise decay laws

> M12 extension. Ideas ported from the `distglot` toolkit (an earlier
> distributional-glottochronology project). Clepsydra's pairwise likelihood
> assumes a **constant** replacement rate; several published models assume it is
> not constant. These are added as *competing hypotheses to compare*, never as
> defaults — the model-cards principle that flexibility is not correctness
> applies here too.

## The laws (`src/gd_decay.{c,h}`)

`s(T)` is the expected fraction of shared cognates between two languages that
diverged `T` years ago (both lineages have evolved for `T`). With `t = T/1000`
in millennia:

| model | s(T) | idea |
|---|---|---|
| `swadesh` | `r^{2t}` | classical **constant** rate (the baseline; = P0) |
| `starostin` | `exp(-2λt²)` | Starostin 2000: retained words are increasingly stable, so replacement **accelerates** with time |
| `vasilyev` | `f + (1−f)·e^{-λt}(1+λt)` | Vasilyev's coherent + independent two-component form, with an ultra-stable **floor** `f` (published constants `f=0.08`, `λ=0.45`) |
| `two-component` | `w_s(a+b·e^{-kt}) + w_r·r^{t}` | a configurable stable-core + regular-vocabulary mixture |

The key qualitative difference is the **tail**: the constant-rate law decays to
zero, while the two-component and Vasilyev laws keep a stable floor — two very
distant languages still share a core. That changes deep dates a lot.

> **Fixed constants, not fit.** The Vasilyev floor `f` and rate `λ` (and the
> two-component `a,b,k,w_s`) are constants *of the law*, taken as published —
> they are **not** estimated from your data. They live in `gd_decay_params`
> (defaults `vasilyev_floor=0.08`, `vasilyev_rate=0.45`) so they are explicit and
> overridable, never magic numbers buried in the code. The consequence is an
> identifiability caveat: an absolute date under these laws rests on those fixed
> constants as much as on the observed similarity. Near the floor the curve is
> nearly flat, so deep dates are **weakly identified** — a small similarity
> change moves the date a lot. This is why the laws are competing hypotheses to
> compare by held-out prediction (`docs/model-cards.md`), not defaults.

## How they are used

Each law is wrapped in a proper Bayesian posterior over the divergence time
(`gd_decay_posterior`, uniform prior on `T`, grid integration), so an
alternative law yields a full posterior with intervals — not a point estimate,
unlike the original inversion in distglot. From the CLI:

```
clepsydra date-pair --forms forms.tsv --lang-a A --lang-b B --out DIR \
    --decay starostin        # or vasilyev, two-component, swadesh
```

Omitting `--decay` keeps the default hierarchical constant-rate model (P1). The
chosen model is recorded in `posterior_pair.tsv` (`model` field).

## Why compare rather than pick

On real data the laws disagree by design. For example, Spanish–Portuguese from
IE-CoR (`examples/pair_simple`):

| model | divergence median (BP) |
|---|---|
| P1 constant-rate | ~570 |
| Starostin (accelerating) | ~1320 |
| Vasilyev (two-component) | ~1675 |

None of these is "the" answer. Which law fits basic-vocabulary decay is an
empirical question, and the honest way to choose is **held-out prediction in a
shared observable space** (`docs/validation.md` §6), not preference for the one
whose date you like. Clepsydra ports the laws so that comparison is possible; it
does not endorse one.

## References

- Swadesh, M. 1955. Towards greater accuracy in lexicostatistic dating.
- Starostin, S. 2000. Comparative-historical linguistics and lexicostatistics.
- Vasilyev's two-component glottochronological model.
- `distglot` — the distributional-glottochronology toolkit these laws are ported
  from (with the point-estimate inversion replaced by a Bayesian posterior here).
