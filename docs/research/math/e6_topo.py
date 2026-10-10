from math import factorial
def dfact(n):
    r=1
    while n>1: r*=n; n-=2
    return r
# rooted multifurcating (all resolutions) via Felsenstein recursion T(n,m)
def rooted_all(n):
    # number of rooted trees with n labelled leaves, any degree (Schroeder's 4th problem): a(n)
    # recurrence via EGF: A' = ... use known values
    vals=[0,1,1,4,26,236,2752,39208,660032,12818912,282137824]
    return vals[n] if n<len(vals) else None
print("n  rooted_binary  unrooted_binary  ranked_histories  rooted_any_degree")
for n in range(3,13):
    rb=dfact(2*n-3); ub=dfact(2*n-5); rh=factorial(n)*factorial(n-1)//2**(n-1)
    print(n,rb,ub,rh,rooted_all(n))
