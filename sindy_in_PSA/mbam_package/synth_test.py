import sys; sys.path.insert(0,'/tmp/claude-0/mbam')
exec(open('/tmp/claude-0/mbam/load10.py').read()); exec(open('/tmp/claude-0/mbam/prep10.py').read())
from mbam_run import *
d,_=prep_gen(G[2],2)
prob=make_problem(d,2,False)
true=prob.p_base.copy()
for n,v in dict(H2=12.8,D=2.0,kf=1.8/d['vf'][0],km=0.83,c=0.1,eq0=1.0,ed0=0.4).items(): true[IDX[n]]=v
P,Q,W,Dd=prob.simulate(true)
rng=np.random.default_rng(0)
# noisy synthetic measurements (0.2 % of range)
prob.y=dict(P=P+rng.normal(0,2e-3*np.ptp(P),len(P)),Q=Q+rng.normal(0,2e-3*np.ptp(Q),len(P)),
            w=W+rng.normal(0,2e-3*np.ptp(W),len(P)),dmt=Dd-true[IDX['c']]+rng.normal(0,2e-3*np.ptp(Dd),len(P)))
prob.sig={k:max(0.005*np.ptp(v),1e-9) for k,v in prob.y.items()}
prob.w0=W[0]
# start from perturbed values (as if H, D, kf, km unknown)
start=prob.p_base.copy()
for n,v in dict(H2=20.0,D=0.5,kf=0.7*true[IDX['kf']],km=1.0,c=0.0,eq0=1.1,ed0=0.3).items(): start[IDX[n]]=v
t0=time.time(); p,cost,_=initial_fit(prob),None,None
p,cost=p
print("fit time",round(time.time()-t0,1),"s")
print(f"{'param':7s} {'true':>10s} {'fit':>10s}")
for n in NAMES: print(f"{n:7s} {true[IDX[n]]:10.4f} {p[IDX[n]]:10.4f}")
np.save('/tmp/claude-0/mbam/synth_fit.npy',p)
free=prob.free_names(); x=prob.to_x(p,free)
J,g,w,v=fisher(prob,x,free)
print("FIM eigenvalues (log-parameter metric):",np.array2string(w[::-1],precision=2,max_line_width=160))
print("sloppiest eigenvector (|component|>0.2):",[(free[i],round(v[i,0],2)) for i in range(len(free)) if abs(v[i,0])>0.2])
