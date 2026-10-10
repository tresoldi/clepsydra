"""Check E[K] (number of distinct cognate classes per meaning) on a fixed tree:
E[K] = rho(root) + sum_b (1-exp(-r t_b)) rho(child_b),  rho(v)=1-prod_c(1-exp(-r t_c) rho(c)), rho(tip)=1."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np, e2_tree_fisher as m
ages={v:m.nodes[v][1] for v in m.order}; ages.update({t:0.0 for t in m.tips})
L=sum(ages[v]-ages[c] for v in m.order for c in m.nodes[v][0])
def EK(r):
    rho={t:1.0 for t in m.tips}; s=0
    for v in m.order:
        pr=1
        for c in m.nodes[v][0]:
            e=np.exp(-r*(ages[v]-ages[c])); pr*=1-e*rho[c]; s+=(1-e)*rho[c]
        rho[v]=1-pr
    return rho['R']+s
rng=np.random.default_rng(0)
for r in [0.1,0.5,1.0,3.0]:
    tot=0;N=40000
    for _ in range(N):
        st={'R':0};nxt=1
        for v in reversed(m.order):
            for c in m.nodes[v][0]:
                if rng.random()<np.exp(-r*(ages[v]-ages[c])): st[c]=st[v]
                else: st[c]=nxt;nxt+=1
        tot+=len(set(st[t] for t in m.tips))
    print("r=%.1f  formula E[K]=%.4f  sim=%.4f  1+rL=%.4f"%(r,EK(r),tot/N,1+r*L))
# Starostin c = W(lam t^2)/(lam t^2) vs gamma
def W(y):
    w=np.log1p(y)
    for _ in range(50): w=w-(w*np.exp(w)-y)/(np.exp(w)*(1+w))
    return w
lam=0.05
print("Starostin retention c(t)=W(lam t^2)/(lam t^2), lam=0.05/ky^2:")
for t in [1,2,4,6,8,10,15,20]:
    y=lam*t*t; print("  t=%g c=%.3f"%(t,W(y)/y))
