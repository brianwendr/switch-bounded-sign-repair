import random, time, statistics, math

def bswitch_values(signs, Bmax, start_positive=True, weights=None):
    n=len(signs)
    if weights is None: weights=[1]*n
    # signs: +1 if positive else -1
    # block b target = + if (start_positive xor b parity?)
    vals=[]
    # prefix costs per b as needed
    # F_prev for b=0
    def target(b):
        p = 1 if start_positive else -1
        return p if b%2==0 else -p
    pref=[]
    for b in range(Bmax+1):
        P=[0]*(n+1)
        t=target(b)
        for i,s in enumerate(signs):
            P[i+1]=P[i]+(weights[i] if s!=t else 0)
        pref.append(P)
    F_prev=pref[0][:]
    vals.append(F_prev[n])
    for b in range(1,Bmax+1):
        P=pref[b]
        F=[0]*(n+1)
        # running minimum of F_prev[u]-P[u]
        m=F_prev[0]-P[0]
        for t in range(0,n+1):
            if t>0:
                cand=F_prev[t]-P[t]
                if cand<m: m=cand
            F[t]=P[t]+m
        F_prev=F
        vals.append(F_prev[n])
    return vals

def bswitch_value_minrun(signs, B, ell, start_positive=True):
    # at most B switches = up to B+1 segments, each realized segment length>=ell
    # exact segments m=1..B+1, nonempty. Returns min value and m-1 switch count.
    n=len(signs)
    INF=10**9
    def cost_pref(phase):
        p=1 if start_positive else -1
        target=p if phase%2==0 else -p
        P=[0]*(n+1)
        for i,s in enumerate(signs): P[i+1]=P[i]+(1 if s!=target else 0)
        return P
    prefs=[cost_pref(j) for j in range(B+1)]
    # F for m=1 segments phase 0, valid t>=ell only; F[0] invalid unless allow empty? stable no.
    F_prev=[INF]*(n+1)
    P0=prefs[0]
    for t in range(ell, n+1): F_prev[t]=P0[t]
    best=(F_prev[n],0)
    for seg in range(2, B+2):
        phase=seg-1
        P=prefs[phase]
        F=[INF]*(n+1)
        # F[t]=P[t]+min_{u <= t-ell} F_prev[u]-P[u]
        m=INF
        for t in range(0,n+1):
            u=t-ell
            if u>=0:
                cand=F_prev[u]-P[u]
                if cand<m: m=cand
            if m<INF:
                F[t]=P[t]+m
        F_prev=F
        if F_prev[n]<best[0]: best=(F_prev[n],seg-1)
    return best

def make_planted(lengths, flip_p, rng, start=1):
    signs=[]; cur=start
    for L in lengths:
        for _ in range(L):
            s=cur
            if rng.random()<flip_p: s=-s
            signs.append(s)
        cur=-cur
    return signs

def random_signs(n,p,rng):
    return [1 if rng.random()<p else -1 for _ in range(n)]

rng=random.Random(20260605)
# Table 1 planted frontier
lengths=[80,60,100,61] # n=301, true switch count 3
for flip in [0.05,0.10,0.20]:
    vals=[]; times=[]
    for _ in range(100):
        signs=make_planted(lengths, flip, rng)
        t0=time.perf_counter()
        v=bswitch_values(signs,6,True)
        times.append((time.perf_counter()-t0)*1000)
        vals.append(v)
    means=[statistics.mean(v[b] for v in vals) for b in range(7)]
    sds=[statistics.stdev(v[b] for v in vals) for b in range(7)]
    print('frontier flip',flip)
    print('  mean_sd', ['%.1f (%.1f)'%(m,s) for m,s in zip(means,sds)])
    print('  runtime mean_sd', '%.3f (%.3f)'%(statistics.mean(times),statistics.stdev(times)), 'max', round(max(times),3))

# Table 2 penalty path example
rng2=random.Random(12345)
signs=make_planted([30,18,24,18,20],0.12,rng2) # n=110, true 4 switches
rhos=bswitch_values(signs,8,True)
print('rhos',rhos)
intervals=[]
Bmax=8
for b,rho in enumerate(rhos):
    L=0.0
    U=float('inf')
    for c,rhoc in enumerate(rhos):
        if c==b: continue
        if c>b:
            L=max(L,(rho-rhoc)/(c-b))
        else: # c<b
            U=min(U,(rhoc-rho)/(b-c))
    if L <= U + 1e-9:
        intervals.append((b,rho,L,U))
print('intervals', intervals)

# Repeated penalty-path summary for the Table 3 caption.
rng_penalty = random.Random(12345)
non_dom_counts=[]
for _ in range(100):
    signs_trial=make_planted([30,18,24,18,20],0.12,rng_penalty)
    rhos_trial=bswitch_values(signs_trial,8,True)
    count=0
    for b,rho in enumerate(rhos_trial):
        L=0.0
        U=float('inf')
        for c,rhoc in enumerate(rhos_trial):
            if c==b:
                continue
            if c>b:
                L=max(L,(rho-rhoc)/(c-b))
            else:
                U=min(U,(rhoc-rho)/(b-c))
        if L <= U + 1e-9:
            count += 1
    non_dom_counts.append(count)
print('Table 3 repeated penalty-path non-dominated count range', min(non_dom_counts), max(non_dom_counts))

# Table 3 min-run
for ell in [1,5,10,20,30]:
    vals=[]; bs=[]; times=[]
    rng3=random.Random(9876) # same set for each ell? regenerate consistently
    for _ in range(100):
        signs=make_planted([80,60,100,61],0.10,rng3)
        t0=time.perf_counter()
        if ell==1:
            v=bswitch_values(signs,8,True)
            # at most B, choose min (will be v[8] since nonincreasing)
            val=v[8]; b=8
        else:
            val,b=bswitch_value_minrun(signs,8,ell,True)
        times.append((time.perf_counter()-t0)*1000)
        vals.append(val); bs.append(b)
    print('minrun ell',ell,'meanval',round(statistics.mean(vals),1),'meanB',round(statistics.mean(bs),2),'overselect_rate',round(sum(b>3 for b in bs)/len(bs),2),'runtime',round(statistics.mean(times),3),round(max(times),3))

# Table 4 scaling Bmax on random signs maybe not use
for n in [200,500,1000,2000,5000]:
    rng4=random.Random(42+n)
    times=[]
    for _ in range(30):
        signs=random_signs(n,0.6,rng4)
        t0=time.perf_counter(); v=bswitch_values(signs,8,True); times.append((time.perf_counter()-t0)*1000)
    print('scaling n',n,'rho8',v[8],'time',round(statistics.mean(times),3),round(max(times),3))

# Exhaustive verification for Section 8.
# Compares the O(BN) dynamic program against brute-force enumeration of all
# alternating breakpoint templates for all binary strings up to n=10 and B<=5.
def brute_force_value(signs, B, start_positive=True):
    n=len(signs)
    p0=1 if start_positive else -1
    best=10**9
    # Enumerate template breakpoints using recursion. Empty blocks are allowed.
    cuts=[]
    def rec(pos, remaining, prefix):
        if remaining==0:
            yield prefix+[n]
        else:
            for u in range(pos, n+1):
                yield from rec(u, remaining-1, prefix+[u])
    for cuts in rec(0, B, [0]):
        cost=0
        for b in range(B+1):
            target=p0 if b%2==0 else -p0
            for q in range(cuts[b], cuts[b+1]):
                if signs[q]!=target:
                    cost+=1
                    if cost>=best: break
            if cost>=best: break
        if cost<best: best=cost
    return best

def exhaustive_check(max_n=10, max_B=5):
    for n in range(1,max_n+1):
        for mask in range(1<<n):
            signs=[1 if (mask>>i)&1 else -1 for i in range(n)]
            vals=bswitch_values(signs, max_B, True)
            for B in range(max_B+1):
                br=brute_force_value(signs, B, True)
                if vals[B]!=br:
                    raise AssertionError((n, mask, B, vals[B], br, signs))
    return True

print('exhaustive_check n<=10 B<=5', exhaustive_check())
