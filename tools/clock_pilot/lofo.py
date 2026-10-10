"""Leave-one-family-out scoring of calibration-age predictors on x = ln T."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats

from . import model

LEVELS = (0.5, 0.8, 0.95)
OLD_MU, OLD_SD = -8.8, 0.8


class Mixture:
    def __init__(self, w, m, s):
        self.w, self.m, self.s = map(np.asarray, (w, m, s))

    def cdf(self, x):
        x = np.atleast_1d(x)
        out = np.zeros(x.shape)
        for i in range(0, len(self.w), 2000):
            sl = slice(i, i + 2000)
            out += (self.w[sl, None] * stats.norm.cdf(
                (x[None, :] - self.m[sl, None]) / self.s[sl, None])).sum(0)
        return out

    def pdf(self, x):
        x = np.atleast_1d(x)
        out = np.zeros(x.shape)
        for i in range(0, len(self.w), 2000):
            sl = slice(i, i + 2000)
            out += (self.w[sl, None] * stats.norm.pdf(
                x[None, :], self.m[sl, None], self.s[sl, None])).sum(0)
        return out

    @property
    def mean(self):
        return float((self.w * self.m).sum())

    @property
    def sd(self):
        mu = self.mean
        return float(math.sqrt((self.w * (self.s ** 2 + self.m ** 2)).sum() - mu ** 2))


class T:
    def __init__(self, df, loc, scale):
        self.d = stats.t(df, loc, scale)
        self.df, self.loc, self.scale = df, loc, scale

    def cdf(self, x):
        return self.d.cdf(np.atleast_1d(x))

    def pdf(self, x):
        return self.d.pdf(np.atleast_1d(x))

    @property
    def mean(self):
        return self.loc

    @property
    def sd(self):
        return float(self.d.std()) if self.df > 2 else float("inf")


def calib_density(node, xs):
    """Calibration density q(x) on the grid xs (x = ln T)."""
    if node.use == "interval":
        a = math.log(node.lo + 50.0)
        b = math.log(node.hi + 50.0)
        q = ((xs >= a) & (xs <= b)).astype(float)
    else:
        c, s = node.log_moments()
        q = stats.norm.pdf(xs, c, s)
    q /= np.trapezoid(q, xs)
    return q


def quantile(pred, p, xs, cdf):
    return float(np.interp(p, cdf, xs))


def score(pred, node, center):
    xs = np.linspace(center - 9.0, center + 9.0, 9001)
    pdf, cdf = pred.pdf(xs), pred.cdf(xs)
    q = calib_density(node, xs)
    out = {"log_score": float(math.log(max(np.trapezoid(pdf * q, xs), 1e-300))),
           "pit": float(np.trapezoid(cdf * q, xs))}
    for L in LEVELS:
        lo = quantile(pred, (1 - L) / 2, xs, cdf)
        hi = quantile(pred, 1 - (1 - L) / 2, xs, cdf)
        inside = ((xs >= lo) & (xs <= hi)).astype(float)
        out[f"cov{int(L * 100)}"] = float(np.trapezoid(inside * q, xs))
        if L == 0.95:
            a = 0.05
            out["interval_score"] = float((hi - lo) + (2 / a) * max(lo - center, 0)
                                          + (2 / a) * max(center - hi, 0))
            out["width95_factor"] = float(math.exp((hi - lo) / 2))
            out["lo95_years"] = math.exp(lo) - 50.0
            out["hi95_years"] = math.exp(hi) - 50.0
    out["median_years"] = math.exp(quantile(pred, 0.5, xs, cdf)) - 50.0
    out["z"] = float((center - pred.mean) / pred.sd)
    return out


def censored_score(pred, node):
    x = math.log(node.bound())
    F = float(pred.cdf(x)[0])
    return math.log(max(F if node.use == "max" else 1 - F, 1e-300))




# ------------------------------------------------------------- predictors
def clock_bayes(train, node_est, prior_scale=1.0, ng=151):
    b, s2, fam = train
    f = model.fit(b, s2, fam, prior_scale=prior_scale, ng=ng)
    w, m, s = f.predictive(node_est.log_d - model.LN2, node_est.s_lex ** 2, tol=1e-8)
    return Mixture(w, m, s), f


def clock_reml(fitobj, node_est):
    sf, tau, mu, varmu = fitobj.reml()
    F = fitobj.F
    scale = math.sqrt(sf ** 2 + tau ** 2 + varmu + node_est.s_lex ** 2)
    return T(max(F - 1, 1), node_est.log_d - model.LN2 - mu, scale)


def old_default(node_est_t1):
    return Mixture([1.0], [node_est_t1.log_d - model.LN2 - OLD_MU],
                   [math.sqrt(OLD_SD ** 2 + node_est_t1.s_lex ** 2)])


def family_size(train_x, train_n, n0):
    x = np.asarray(train_x, float)
    z = np.log(np.asarray(train_n, float))
    n = len(x)
    res = stats.linregress(z, x)
    resid = x - (res.intercept + res.slope * z)
    s2 = (resid ** 2).sum() / (n - 2)
    sxx = ((z - z.mean()) ** 2).sum()
    z0 = math.log(n0)
    se = math.sqrt(s2 * (1 + 1 / n + (z0 - z.mean()) ** 2 / sxx))
    return T(n - 2, res.intercept + res.slope * z0, se)


def climatological(train_x):
    x = np.asarray(train_x, float)
    n = len(x)
    return T(n - 1, float(x.mean()), float(x.std(ddof=1) * math.sqrt(1 + 1 / n)))


@dataclass
class Row:
    family: str
    node: str
    predictor: str
    scores: dict
