import sys,warnings; sys.path.insert(0,'/tmp/claude-0/mbam'); warnings.simplefilter('ignore')
exec(open('/tmp/claude-0/mbam/load10.py').read()); exec(open('/tmp/claude-0/mbam/prep10.py').read())
from mbam_run import *
import pickle
K=int(sys.argv[1]); T0=float(sys.argv[2]) if len(sys.argv)>2 else 2.4
d,_=prep_gen(G[K],K,t0=T0)
prob=make_problem(d,K,SALIENT[K])
print(f"G{K}: window {d['t'][0]:.2f}-{d['t'][-1]:.2f} s, {len(d['t'])} samples, salient={SALIENT[K]}")
p,cost=initial_fit(prob)
r=prob.resid_p(p); print(f"full-model fit: rms {np.sqrt(np.mean(r**2)):.2f} sigma; per-signal rms (sigma):",
      [round(float(np.sqrt(np.mean(r[i*prob.N:(i+1)*prob.N]**2))),2) for i in range(4)],"(P,Q,w,load angle)")
p_full=p.copy()
t0=time.time(); p,rms,log=mbam_reduce(prob,p,max_iter=8)
print(f"\nMBAM {time.time()-t0:.0f}s; limiting parameters: {prob.fixed}")
pickle.dump(dict(p_full=p_full,p_red=p,fixed=prob.fixed,log=log,rms=rms),open(f'/tmp/claude-0/mbam/res_G{K}.pkl','wb'))
