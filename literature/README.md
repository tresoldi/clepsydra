# Required reading

The papers this implementation is built on are **not redistributed here** —
they are third-party publications under their own copyright. Each is cited
below with a DOI; `docs/paper-evaluation.md` records, per paper, what was
adopted, adapted, rejected or postponed, and why.

- **Rönchen, P., Wiklund, T. & Hammarström, H. (2024).** *Likelihood-based
  methods for lexical dating.* Language Dynamics and Change 14(1), 1–41.
  doi:10.1163/22105832-bja10032 — the primary reference; the T1 novel-state
  likelihood is an independent re-implementation of its model.
- **Pagel, M., Atkinson, Q. D. & Meade, A. (2007).** *Frequency of word-use
  predicts rates of lexical evolution throughout Indo-European history.*
  Nature 449, 717–720. doi:10.1038/nature06176 — the concept-rate prior.
- **Rama, T. & Wichmann, S. (2020).** *A test of Generalized Bayesian dating.*
  PLOS ONE. doi:10.1371/journal.pone.0147924
- **King, R. et al. (2026).** doi:10.1371/journal.pcbi.1014312 — the
  meaning-weighted partition rates used by the T2 binary likelihood.

`docs/model-spec.md` is self-contained: the mathematics, the generative models
and the pseudocode are all stated there, so none of these papers is needed to
read, build or verify the code. They are needed to judge whether the modelling
choices are the right ones.
