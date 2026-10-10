"""Godambe efficiency of pairwise composite likelihood vs exact novel-state likelihood
on the 7-tip tree, strict clock, known concept rates. Params: 6 log internal ages."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools
import e2_tree_fisher as m
from common import ln_rates
C=200; rates=ln_rates(0.15,0.8,C)
tips=m.tips; order=m.order
def mrca(i,j):
    for v in order:
        if tips[i] in m.TS[v] and tips[j] in m.TS[v]: return v
pairs=list(itertools.combinations(range(7),2)); MR=[order.index(mrca(i,j)) for i,j in pairs]
for T in [1,2,4,8]:
    ages=np.array([T*m.nodes[v][1] for v in order])
    agesd={v:a for v,a in zip(order,ages)}; mult={b:1.0 for b in m.branches}
    P=6; Iex=np.zeros((P,P)); J=np.zeros((P,P)); H=np.zeros((P,P))
    h=1e-4
    # H: sum over pairs of Fisher info of each Bernoulli
    for (i,j),k in zip(pairs,MR):
        p=np.exp(-2*rates*ages[k]); dl=2*rates*ages[k]*p  # |dp/dlog a|
        H[k,k]+=np.sum(dl**2/(p*(1-p)))
    for pat in m.pats:
        p0=m.lik(pat,agesd,mult,rates)
        G=np.zeros((P,C))
        for a in range(P):
            d=dict(agesd); d[order[a]]=agesd[order[a]]*np.exp(h); up=m.lik(pat,d,mult,rates)
            d[order[a]]=agesd[order[a]]*np.exp(-h); dn=m.lik(pat,d,mult,rates)
            G[a]=(up-dn)/(2*h)
        Iex+=(G/np.where(p0>0,p0,1))@G.T
        u=np.zeros((P,C))
        for (i,j),k in zip(pairs,MR):
            p=np.exp(-2*rates*ages[k]); x=2*rates*ages[k]
            u[k]+= -x if pat[i]==pat[j] else x*p/(1-p)
        J+=(u*p0)@u.T
    Vex=np.linalg.inv(Iex); Hi=np.linalg.inv(H); Vcl=Hi@J@Hi
    naive=Hi
    print("T=%g sd(log root): exact %.4f | composite(Godambe) %.4f | composite naive(H^-1) %.4f ; mean over nodes sd ratio CL/exact %.2f"%(
        T,np.sqrt(Vex[5,5]),np.sqrt(Vcl[5,5]),np.sqrt(naive[5,5]),np.mean(np.sqrt(np.diag(Vcl)/np.diag(Vex)))))
