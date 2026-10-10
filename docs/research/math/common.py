import numpy as np
from statistics import NormalDist
def ln_rates(rbar,s,C):
    q=(np.arange(C)+0.5)/C
    z=np.array([NormalDist().inv_cdf(x) for x in q])
    return rbar*np.exp(s*z-s*s/2)
