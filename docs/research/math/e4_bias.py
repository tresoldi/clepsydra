import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import ln_rates
rc=ln_rates(0.15,0.8,200)
def mle_T(Ftrue, model_p, model_dp, lo=0.01, hi=60):
    # solve expected score = 0 by bisection
    def score(T):
        p=model_p(T); dp=model_dp(T); return np.sum((Ftrue-p)*dp/(p*(1-p)))
    a,b=lo,hi
    for _ in range(200):
        m=0.5*(a+b)
        if score(m)>0: a=m   # score>0 means T too small
        else: b=m
    return 0.5*(a+b)
p=lambda T: np.exp(-2*rc*T); dp=lambda T: -2*rc*np.exp(-2*rc*T)
print("(1) Unmodelled concept-x-lineage gamma heterogeneity (shape k), per-concept rates known exactly; pair MLE of T:")
for k in [1,2,4,10]:
    out=[]
    for T in [1,2,4,6,8]:
        F=(1+rc*T/k)**(-2*k)
        out.append("%g->%.2f"%(T,mle_T(F,p,dp)))
    print("  k=%-3g (cv %.2f): "%(k,1/np.sqrt(k))+"  ".join(out))
print("(2) Split-coding error: true cognate pair judged non-cognate w.p. eps (F'=F(1-eps)):")
for eps in [0.02,0.05,0.10]:
    print("  eps=%.2f: "%eps+"  ".join("%g->%.2f"%(T,mle_T(np.exp(-2*rc*T)*(1-eps),p,dp)) for T in [0.5,1,2,4,8]))
print("(3) Lump/false-positive cognacy: non-cognate pair judged cognate w.p. eta (F'=F+(1-F)eta):")
for eta in [0.01,0.03,0.05]:
    print("  eta=%.2f: "%eta+"  ".join("%g->%.2f"%(T,mle_T(np.exp(-2*rc*T)+(1-np.exp(-2*rc*T))*eta,p,dp)) for T in [1,2,4,6,8,10]))
print("(4) Undetected borrowing between the pair: fraction b of replacements copy the other lineage: F=b+(1-b)exp(-2rT):")
for b in [0.05,0.1,0.2]:
    print("  b=%.2f: "%b+"  ".join("%g->%.2f"%(T,mle_T(b+(1-b)*np.exp(-2*rc*T),p,dp)) for T in [1,2,4,6,8]))
print("(5) Split pulse kappa (extra lexical length per lineage at split, in units of mean-rate ky): true F=exp(-2r(T+kappa))")
for kap in [0.2,0.5]:
    print("  kappa=%.1f ky-equiv: "%kap+"  ".join("%g->%.2f"%(T,mle_T(np.exp(-2*rc*(T+kap)),p,dp)) for T in [0.5,1,2,4]))
