import numpy as np
rbar=0.15  # per ky per lineage, mean over meanings (retention ~0.86/ky)
def S_gamma(t,k): th=rbar/k; return (1+th*t)**(-k)
xh,wh=np.polynomial.hermite_e.hermegauss(80); wh=wh/wh.sum()
def S_logn(t,cv):
    s=np.sqrt(np.log1p(cv*cv)); r=rbar*np.exp(s*xh-s*s/2); return (wh*np.exp(-np.outer(t,r))).sum(1)
ts=np.array([0.5,1,2,3,4,6,8,10,15])
print("Pair shared-cognate fraction F(T)=S(2T) and naive date -lnF/(2 rbar) [true T in ky]")
print("T      exp    | gamma k=4 (cv.5)   | gamma k=1.56 (cv.8) | gamma k=1 (cv1) | lognormal cv.8 ")
for T in ts:
    row=["%5.1f %6.3f"%(T,np.exp(-2*rbar*T))]
    for S in [lambda t:S_gamma(t,4),lambda t:S_gamma(t,1/0.64),lambda t:S_gamma(t,1.0),lambda t:S_logn(np.atleast_1d(t),0.8)[0]]:
        F=S(2*T); row.append("%6.3f->%5.2f"%(F,-np.log(F)/(2*rbar)))
    print(" | ".join(row))
print("\nCalibrated-at-2ky exponential, extrapolated (true heterogeneity gamma cv .8):")
k=1/0.64; Fcal=S_gamma(4.0,k); rhat=-np.log(Fcal)/4.0
print("  rhat=%.4f (true mean %.3f)"%(rhat,rbar))
for T in [1,2,4,6,8,10]:
    F=S_gamma(2*T,k); print("  T=%g -> That=%.2f  (ratio %.2f)"%(T,-np.log(F)/(2*rhat),-np.log(F)/(2*rhat)/T))
print("\nRatio of node ages under naive exponential when truth gamma cv .8 (t1/t2 true=0.5):")
for t2 in [2,4,8]:
    t1=t2/2; a=-np.log(S_gamma(2*t1,k)); b=-np.log(S_gamma(2*t2,k)); print("  t2=%g: est ratio %.3f"%(t2,a/b))
print("\nAsymptotics: gamma S ~ (theta t)^-k (power law); lognormal ~ exp(-(ln t)^2/(2s^2)) slower than exp, faster than any power eventually")
