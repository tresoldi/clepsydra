"""Exact Fisher information for node ages on a 7-tip tree under the strict
novel-state (infinite-alleles) model, by enumerating all 877 set partitions.
Params: log internal ages (6), log branch multipliers (12), log family rate (1).
Bayesian information J = E[I] + prior precision (van Trees)."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import ln_rates
tips=list("ABCDEFG")
# internal nodes: name -> (children, relative age)
nodes={'AB':(['A','B'],0.2),'ABC':(['AB','C'],0.4),'DE':(['D','E'],0.3),
       'ABCDE':(['ABC','DE'],0.7),'FG':(['F','G'],0.5),'R':(['ABCDE','FG'],1.0)}
order=['AB','ABC','DE','ABCDE','FG','R']
branches=[c for n in order for c in nodes[n][0]]  # 12 branches by child
def partitions(n):
    if n==0: yield []; return
    for p in partitions(n-1):
        k=max(p)+1 if p else 0
        for lab in range(k+1): yield p+[lab]
pats=[tuple(p) for p in partitions(7)]
def tipset(v):
    return [v] if v in tips else [x for c in nodes[v][0] for x in tipset(c)]
TS={v:tipset(v) for v in tips+order}
def lik(pat,ages,mult,rates):
    lab={t:pat[i] for i,t in enumerate(tips)}
    allcnt={}
    for t in tips: allcnt[lab[t]]=allcnt.get(lab[t],0)+1
    age=dict(ages); 
    for t in tips: age[t]=0.0
    F={}
    R=len(rates)
    for t in tips: F[t]={lab[t]:np.ones(R),None:np.zeros(R)}
    for v in order:
        labs=set(lab[t] for t in TS[v])
        fv={s:np.ones(R) for s in list(labs)+[None]}
        for c in nodes[v][0]:
            bl=mult[c]*(age[v]-age[c])
            e=np.exp(-rates*bl)
            fc=F[c]; clabs=set(lab[t] for t in TS[c])
            g=fc[None].copy(); conf=set()
            for a in clabs:
                if sum(1 for t in TS[c] if lab[t]==a)==allcnt[a]: g=g+fc[a]; conf.add(a)
            for s in fv:
                if s in clabs:   # tips of class s below c must inherit it: no change allowed
                    fv[s]=fv[s]*(e*fc[s])
                else:
                    fv[s]=fv[s]*(e*fc[None]+(1-e)*g)
        F[v]=fv
    return sum(F['R'].values())
def run(T,rates,sb,sf,sref=0.0,nq=1):
    # rate uncertainty per concept marginalised by Gauss-Hermite on log rate
    if sref>0:
        x,w=np.polynomial.hermite_e.hermegauss(nq); w=w/w.sum()
    else: x,w=np.array([0.0]),np.array([1.0])
    rr=(rates[:,None]*np.exp(sref*x[None,:])).ravel()
    P=len(order)+len(branches)+1
    th0=np.zeros(P)
    def unpack(th):
        ages={v:T*nodes[v][1]*np.exp(th[i]) for i,v in enumerate(order)}
        mult={b:np.exp(th[len(order)+j]) for j,b in enumerate(branches)}
        return ages,mult,np.exp(th[-1])
    h=1e-4
    I=np.zeros((P,P)); tot=0
    for pat in pats:
        def L(th):
            a,m,g=unpack(th); v=lik(pat,a,m,rr*g)
            return (v.reshape(len(rates),len(x))*w).sum(1)
        p0=L(th0); tot+=p0
        G=np.zeros((P,len(rates)))
        for i in range(P):
            d=np.zeros(P); d[i]=h
            G[i]=(L(th0+d)-L(th0-d))/(2*h)
        with np.errstate(divide='ignore',invalid='ignore'):
            wgt=np.where(p0>1e-300,1/p0,0)
        I+=(G*wgt)@G.T
    return I,tot
def report(I,sb,sf,label):
    P=I.shape[0]; no=len(order)
    prior=np.zeros(P); prior[no:no+len(branches)]=1/sb**2 if sb>0 else 1e12; prior[-1]=1/sf**2 if sf>0 else 1e12
    J=I+np.diag(prior)
    S=np.linalg.inv(J)
    sd=np.sqrt(np.diag(S)[:no])
    # ratio age_AB/root etc: var(log a_i - log a_R)
    rat=[np.sqrt(S[i,i]+S[no-1,no-1]-2*S[i,no-1]) for i in range(no-1)]
    return label+" sd(log age) "+" ".join("%s=%.3f"%(v,s) for v,s in zip(order,sd))+" | sd(log ratio to root) "+" ".join("%.3f"%r for r in rat)
if __name__=="__main__":
    C=int(sys.argv[1]) if len(sys.argv)>1 else 200
    sref=float(sys.argv[2]) if len(sys.argv)>2 else 0.0
    rates=ln_rates(0.15,0.8,C)
    print("C=%d concept-rate uncertainty sref=%.2f"%(C,sref))
    for T in [1,2,4,6,8]:
        I,tot=run(T,rates,0,0,sref=sref,nq=7)
        print("T=%g ky  sum P=%.6f (should be C=%d)"%(T,tot.sum(),C))
        for sb,sf in [(0,0),(0.1,0),(0.3,0),(0,0.3),(0.3,0.3),(0.3,0.5)]:
            print("  ",report(I,sb,sf,"sb=%.1f sf=%.1f"%(sb,sf)))
