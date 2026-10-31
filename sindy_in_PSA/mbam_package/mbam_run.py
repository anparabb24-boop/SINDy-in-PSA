import numpy as np, time
from mbam_core import *
NOMINAL={
 1:dict(H2=47.28,D=2.0,xd=0.3614,xq=0.2398,xdp=0.1505,xqp=0.2398,xdpp=0.0990,xqpp=0.1114,Td0p=8.96,Tq0p=1.0,Td0pp=0.05,Tq0pp=0.09),
 2:dict(H2=12.8,D=2.0,xd=1.7199,xq=1.6598,xdp=0.23,xqp=0.378,xdpp=0.1728,xqpp=0.192,Td0p=6.0,Tq0p=0.535,Td0pp=0.05,Tq0pp=0.04),
 3:dict(H2=6.02,D=2.0,xd=1.68,xq=1.61,xdp=0.2321,xqp=0.32,xdpp=0.128,xqpp=0.192,Td0p=5.89,Tq0p=0.6,Td0pp=0.05,Tq0pp=0.04)}
def make_problem(d,k,salient,Ob=2*np.pi*60,sig_frac=0.005,nominal_extra=None):
    nom=dict(NOMINAL[k]); nom.update(dict(kf=1.8/d['vf'][0],km=1.0,c=0.0,eq0=1.1,ed0=0.3 if not salient else 0.0))
    if nominal_extra: nom.update(nominal_extra)
    return SGProblem(d['t'],d['V'],d['th'],d['vf'],d['Pm'],d['P'],d['Q'],d['w'],d['dmt'],nom,Ob,salient=salient,sig_frac=sig_frac)
def steady_calibration(prob,p,tail=1.0):
    """Use the (settled) end of the record: algebraic steady-state relations give the angle offset c, vf scale and Pm scale."""
    m=prob.t>=prob.t[-1]-tail
    P,Q,V,dm=[float(np.mean(prob.y[k][m])) for k in ('P','Q','w','dmt')][:2]+[0,0]
    P=float(np.mean(prob.y['P'][m])); Q=float(np.mean(prob.y['Q'][m])); dm=float(np.mean(prob.y['dmt'][m]))
    V=float(np.mean(prob.Vf[-int(tail*60*prob.ns):])); xq,xd,xdp=p[IDX['xq']],p[IDX['xd']],p[IDX['xdp']]
    best=None
    for dd in np.linspace(dm-1.2,dm+1.2+0.6,2400):
        Vd,Vq=V*np.sin(dd),V*np.cos(dd); idc=(P*Vd+Q*Vq)/V**2; iqc=(P*Vq-Q*Vd)/V**2
        e=abs(Vd-xq*iqc)
        if best is None or e<best[0]: best=(e,dd,idc,iqc,Vd,Vq)
    _,dd,idc,iqc,Vd,Vq=best
    out=p.copy()
    out[IDX['c']]=dd-dm
    out[IDX['kf']]=(Vq+xd*idc)/float(np.mean(prob.vraw[-int(tail*60*prob.ns):]))
    out[IDX['km']]=P/float(np.mean(prob.Pmraw[-int(tail*60*prob.ns):]))
    out[IDX['eq0']]=Vq+xdp*idc; out[IDX['ed0']]=0.0 if prob.salient else (xq-p[IDX['xqp']])*iqc
    return out
def initial_fit(prob,verbose=True):
    p=prob.p_base.copy()
    for n,v in prob.fixed.items(): p[IDX[n]]=v
    p=steady_calibration(prob,p)
    sub=[n for n in ['c','kf','km','eq0','ed0','H2','D'] if n not in prob.fixed]
    p,cost,_=prob.fit(p,free=sub)
    if verbose: print(f"   stage1 (nuisance + H, D) cost {cost:.1f}")
    p,cost,sol=prob.fit(p)
    if verbose: print(f"   stage2 (all) cost {cost:.1f}  (rms residual {np.sqrt(2*cost/len(sol.fun)):.2f} sigma)")
    return p,cost

# ----------------------------------------------------------------------------- MBAM reduction loop (paper Fig.1, steps 3-4)
T_ALG, T_FROZEN = 1e-4, 1e5
SUPPORTED=('Td0p','Tq0p','Td0pp','Tq0pp','D')
def apply_limit(prob,name,direction):
    """direction<0: parameter -> 0 ; >0: -> infinity.  Returns text describing the limiting model or None if unsupported."""
    if name=='D' and direction<0: prob.fixed['D']=0.0; return "D -> 0: damping term removed"
    if name in ('Td0p','Tq0p','Td0pp','Tq0pp'):
        st={'Td0p':"e'q",'Tq0p':"e'd",'Td0pp':"e''q",'Tq0pp':"e''d"}[name]
        if direction<0: prob.fixed[name]=T_ALG; return f"{DESC[name].split()[0]} -> 0: {st} becomes an algebraic variable (singular-perturbation limit)"
        prob.fixed[name]=T_FROZEN; return f"{DESC[name].split()[0]} -> inf: {st} frozen (constant)"
    return None
def mbam_reduce(prob,p,max_iter=8,exc_thresh=2.5,accept_rel=0.10,accept_abs=0.05,dtau=0.25,max_steps=60,verbose=True,log=None):
    log=log if log is not None else []
    free=prob.free_names(); x=prob.to_x(p,free); r=prob.resid(x,free); rms=np.sqrt(np.mean(r**2))
    for it in range(1,max_iter+1):
        free=prob.free_names(); x=prob.to_x(p,free)
        J,g,w,v=fisher(prob,x,free)
        vmin=v[:,0]
        if verbose: print(f"\n--- MBAM iteration {it}: {len(free)} free parameters, rms residual {rms:.3f} sigma, FIM eigenvalue range {w[-1]:.1e} .. {w[0]:.1e}")
        runs=[]
        for s in (+1,-1):
            t0=time.time(); path,reason=geodesic(prob,x,free,s*vmin,dtau=dtau,max_steps=max_steps)
            xe=path[-1][1]; dx=xe-x; runs.append((s,path,reason,dx))
            if verbose: print(f"   geodesic direction {s:+d}: {len(path)-1} steps, stop: {reason}, final cost {path[-1][4]:.1f}, max|dlog p|={np.max(np.abs(dx[IS_LOG[[IDX[n] for n in free]]])):.2f}  ({time.time()-t0:.0f}s)")
        # direction that reaches the boundary soonest / largest excursion at lowest cost
        def score(run):
            s,path,reason,dx=run; lg=IS_LOG[[IDX[n] for n in free]]
            return (('boundary' in reason)*1.0, np.max(np.abs(dx[lg])) - 0.001*path[-1][4])
        s,path,reason,dx=max(runs,key=score)
        lg=IS_LOG[[IDX[n] for n in free]]
        mx=np.max(np.abs(dx[lg])) if lg.any() else 0.0
        thr=max(1.5,0.5*mx) if 'boundary' in reason else exc_thresh          # at a singular-metric boundary use the relative excursion
        moved=[(free[i],dx[i]) for i in range(len(free)) if lg[i] and abs(dx[i])>=thr]
        moved.sort(key=lambda z:-abs(z[1]))
        if verbose: print("   parameters driven to a limit along the geodesic:",[(n,round(float(d),2)) for n,d in moved] or "none")
        entry=dict(iter=it,free=len(free),rms=rms,eig_min=w[0],eig_max=w[-1],moved=moved,reason=reason,accepted=False,action=None)
        applied=False
        for name,d in moved:
            if name in SUPPORTED:
                snap_fixed=dict(prob.fixed); snap_p=p.copy()
                txt=apply_limit(prob,name,d)
                if txt is None: continue
                pn=p.copy()
                for fn,fv in prob.fixed.items(): pn[IDX[fn]]=fv
                pn,cost,sol=prob.fit(pn)
                rms_new=np.sqrt(2*cost/len(sol.fun))
                ok=rms_new<=max(rms*(1+accept_rel),rms+accept_abs)
                if verbose: print(f"   try: {txt}:  rms {rms:.3f} -> {rms_new:.3f} sigma  => {'ACCEPT' if ok else 'reject'}")
                if ok:
                    p,rms=pn,rms_new; entry.update(accepted=True,action=txt,rms_after=rms_new); applied=True; break
                else:
                    prob.fixed=snap_fixed
        log.append(entry)
        if not applied:
            if verbose: print("   no admissible reduction found -> stop")
            break
    return p,rms,log
