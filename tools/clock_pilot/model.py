"""Distances, rates, the random-effects meta-analysis, and predictives.

Notation (plan.md): b_i = ln rate per year for node i; b_i = mu + u_f + e_i,
u_f ~ N(0, sf^2), e_i ~ N(0, tau^2 + s_i^2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats

LN2 = math.log(2.0)
BOOT = 1000
SEED = 20261010


# ---------------------------------------------------------------- distances
def two_rt(p, alpha: float | None):
    """2rT from the expected match probability p (T1, or gamma with shape alpha)."""
    p = np.asarray(p, dtype=float)
    if alpha is None:
        return -np.log(p)
    return alpha * (p ** (-1.0 / alpha) - 1.0)


def pooled_p(K: float, N: float) -> tuple[float, bool]:
    if K >= N:
        return K / (N + 0.5), True
    if K <= 0:
        return 0.5 / N, True
    return K / N, False


@dataclass
class LexEstimate:
    K: float
    N: float
    p: float
    corrected: bool
    d: float                 # 2rT
    log_d: float
    s_lex: float             # SD of ln d over the concept bootstrap
    log_d_median_pairs: float


def lexical(Kc: np.ndarray, Nc: np.ndarray, pair_k: np.ndarray,
            pair_n: np.ndarray, alpha: float | None,
            seed: int = SEED) -> LexEstimate:
    K, N = float(Kc.sum()), float(Nc.sum())
    p, corr = pooled_p(K, N)
    d = float(two_rt(p, alpha))
    rng = np.random.default_rng(seed)
    C = len(Kc)
    W = rng.multinomial(C, np.full(C, 1.0 / C), size=BOOT).astype(float)
    Kb, Nb = W @ Kc, W @ Nc
    pb = np.where(Kb >= Nb, Kb / (Nb + 0.5), np.where(Kb <= 0, 0.5 / np.maximum(Nb, 1), Kb / np.maximum(Nb, 1)))
    s = float(np.std(np.log(two_rt(pb, alpha)), ddof=1))
    ok = pair_n > 0
    pp = np.array([pooled_p(k, n)[0] for k, n in zip(pair_k[ok], pair_n[ok])])
    med = float(np.median(np.log(two_rt(pp, alpha))))
    return LexEstimate(K, N, p, corr, d, math.log(d), s, med)


# ---------------------------------------------------- random-effects (grid)
def _family_arrays(b, s2, fam):
    labs = sorted(set(fam))
    idx = [np.array([i for i, f in enumerate(fam) if f == lab]) for lab in labs]
    return labs, idx


def loglik_grid(b, s2, fam, sf, tau, restricted=True):
    """Log (restricted) likelihood on a grid; sf, tau broadcastable arrays."""
    b, s2 = np.asarray(b, float), np.asarray(s2, float)
    _, idx = _family_arrays(b, s2, fam)
    sf2, tau2 = sf ** 2, tau ** 2
    A = Bc = Q = logdet = 0.0
    for ix in idx:
        D = tau2[..., None] + s2[ix]               # [..., m]
        a = (1.0 / D).sum(-1)
        c = (b[ix] / D).sum(-1)
        q = (b[ix] ** 2 / D).sum(-1)
        den = 1.0 + sf2 * a
        A = A + a / den
        Bc = Bc + c / den
        Q = Q + q - sf2 * c ** 2 / den
        logdet = logdet + np.log(D).sum(-1) + np.log(den)
    mu = Bc / A
    rss = Q - Bc ** 2 / A
    ll = -0.5 * (logdet + rss)
    if restricted:
        ll = ll - 0.5 * np.log(A)
    return ll, mu, A


@dataclass
class Fit:
    grid_sf: np.ndarray
    grid_tau: np.ndarray
    logpost: np.ndarray        # [S, T], normalised weights in .w
    w: np.ndarray
    mu: np.ndarray
    A: np.ndarray
    llr: np.ndarray            # restricted log-lik
    n: int
    F: int

    # ---- Bayesian summaries
    def _q(self, vals, qs):
        o = np.argsort(vals.ravel())
        cw = np.cumsum(self.w.ravel()[o])
        return [float(vals.ravel()[o][np.searchsorted(cw, q)]) for q in qs]

    def post(self, what: str, qs=(0.5, 0.025, 0.975)):
        S, T = np.meshgrid(self.grid_sf, self.grid_tau, indexing="ij")
        vals = {"sf": S, "tau": T, "tot": np.sqrt(S ** 2 + T ** 2),
                "mu": self.mu}[what]
        if what == "mu":
            # mixture of normals: draw-free quantiles via fine grid
            # Grid points with negligible weight are dropped, and the CDF is
            # accumulated in chunks: the dense [S, T, x] array was 2.9 GB and
            # got the process OOM-killed.
            keep = self.w > 1e-12 * self.w.max()
            w = self.w[keep] / self.w[keep].sum()
            mu = self.mu[keep]
            sd = 1.0 / np.sqrt(self.A[keep])
            xs = np.linspace(mu.min() - 6 * sd.max(), mu.max() + 6 * sd.max(), 4001)
            cdf = np.zeros_like(xs)
            for i in range(0, len(w), 2000):
                sl = slice(i, i + 2000)
                cdf += (w[sl, None] * stats.norm.cdf(
                    (xs[None, :] - mu[sl, None]) / sd[sl, None])).sum(0)
            return [float(np.interp(q, cdf, xs)) for q in qs]
        return self._q(vals, qs)

    # ---- REML
    def reml(self):
        i, j = np.unravel_index(np.argmax(self.llr), self.llr.shape)
        return (float(self.grid_sf[i]), float(self.grid_tau[j]),
                float(self.mu[i, j]), float(1.0 / self.A[i, j]))

    def profile_sf(self, drop=1.92):
        prof = self.llr.max(1)
        ok = prof >= prof.max() - drop
        g = self.grid_sf[ok]
        hit_edge = ok[-1]
        return float(g.min()), float(g.max()), bool(hit_edge)

    def profile_tau(self, drop=1.92):
        prof = self.llr.max(0)
        ok = prof >= prof.max() - drop
        g = self.grid_tau[ok]
        return float(g.min()), float(g.max())

    def predictive(self, shift: float, extra_var: float, tol=1e-10):
        """Mixture for x = shift - beta: (weights, means, sds)."""
        S, T = np.meshgrid(self.grid_sf, self.grid_tau, indexing="ij")
        keep = self.w > tol * self.w.max()
        w = self.w[keep] / self.w[keep].sum()
        m = shift - self.mu[keep]
        v = S[keep] ** 2 + T[keep] ** 2 + 1.0 / self.A[keep] + extra_var
        return w, m, np.sqrt(v)


def fit(b, s2, fam, prior_scale=1.0, gmax=3.0, ng=301) -> Fit:
    g = np.linspace(0.0, gmax, ng)
    S, T = np.meshgrid(g, g, indexing="ij")
    llr, mu, A = loglik_grid(b, s2, fam, S, T, restricted=True)
    lp = llr - 0.5 * (S / prior_scale) ** 2 - 0.5 * (T / prior_scale) ** 2
    w = np.exp(lp - lp.max())
    w /= w.sum()
    return Fit(g, g, lp, w, mu, A, llr, len(b), len(set(fam)))


# ------------------------------------------------------ censored ML (GH quad)
GH_X, GH_W = np.polynomial.hermite.hermgauss(40)


def censored_negll(theta, rows):
    """rows: list per family of (kind, value, var) with kind in obs/ge/le."""
    mu, lsf, ltau = theta
    sf, tau = math.exp(lsf), math.exp(ltau)
    u = math.sqrt(2.0) * sf * GH_X
    tot = 0.0
    for fam_rows in rows:
        logf = np.zeros_like(u)
        for kind, v, var in fam_rows:
            sd = math.sqrt(tau ** 2 + var)
            z = (v - mu - u) / sd
            if kind == "obs":
                logf += stats.norm.logpdf(z) - math.log(sd)
            elif kind == "ge":        # beta >= v
                logf += special.log_ndtr(-z)
            else:                      # beta <= v
                logf += special.log_ndtr(z)
        m = logf.max()
        tot += m + math.log((GH_W * np.exp(logf - m)).sum() / math.sqrt(math.pi))
    return -tot


def censored_fit(rows, start):
    best = optimize.minimize(censored_negll, start, args=(rows,),
                             method="Nelder-Mead",
                             options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 20000})
    return best


def censored_profile(rows, best, grid):
    out = []
    for sf in grid:
        lsf = math.log(max(sf, 1e-4))
        f = lambda th: censored_negll((th[0], lsf, th[1]), rows)
        r = optimize.minimize(f, (best.x[0], best.x[2]), method="Nelder-Mead",
                              options={"xatol": 1e-5, "fatol": 1e-7})
        out.append(-r.fun)
    out = np.array(out)
    ok = out >= -best.fun - 1.92
    return float(grid[ok].min()), float(grid[ok].max()), bool(ok[-1])


# ------------------------------------------------------------- depth slope
def depth_slope(b, logT):
    res = stats.linregress(logT, b)
    return float(res.slope), float(res.stderr)
