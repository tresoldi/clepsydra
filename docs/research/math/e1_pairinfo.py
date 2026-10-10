import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import ln_rates


# per-meaning Fisher info about log T for a pair, known rate: f(x)=x^2/(e^x-1), x=2rT
f=lambda x: x**2/np.expm1(x)
xs=np.linspace(0.01,10,100000); i=np.argmax(f(xs)); print('max f at x=%.3f f=%.4f p=%.3f'%(xs[i],f(xs[i]),np.exp(-xs[i])))
for x in [0.2,0.5,1,1.59,2,3,4,6,8]:
    print("x=%.2f p=%.3f f=%.4f CV(C=100)=%.3f CV(C=200)=%.3f"%(x,np.exp(-x),f(x),1/np.sqrt(100*f(x)),1/np.sqrt(200*f(x))))
# heterogeneous concept rates: lognormal, mean rbar=0.15/ky (retention .86/ky), log-sd s
rbar=0.15
print("\nPair, heterogeneous lognormal concept rates (known), CV of T from data alone:")
for s in [0.0,0.5,0.8,1.0]:
    rc=ln_rates(rbar,s,200)
    row=[]
    for T in [1,2,3,4,6,8,10]:
        I=f(2*rc*T).sum()
        row.append("T=%g:%.3f/%.3f"%(T,1/np.sqrt(I),1/np.sqrt(I/2)))
    print("s=%.1f"%s," ".join(row),"  (C=200 / C=100)")
# total CV with family rate prior sd sigma (log scale) combining
print("\nTotal log-SD of T ~ sqrt(cv_data^2+sigma^2) [pair, s=0.8, C=200]")
rc=ln_rates(rbar,0.8,200)
for T in [1,2,4,6,8]:
    cv=1/np.sqrt(f(2*rc*T).sum())
    print("T=%g data=%.3f"%(T,cv)," ".join("sig=%.1f:%.3f (95%%x%.2f)"%(sg,np.hypot(cv,sg),np.exp(1.96*np.hypot(cv,sg))) for sg in [0.1,0.2,0.3,0.5]))
