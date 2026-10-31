import numpy as np, pandas as pd
F0=60.0; a=np.exp(2j*np.pi/3)
SN={1:247.5e6,2:230e6,3:128e6}; VN={1:16.5e3,2:18e3,3:13.8e3}; PP={1:20,2:1,3:1}
SALIENT={1:True,2:False,3:False}
def seq(A,B,C): return (A+B+C)/3,(A+a*B+a*a*C)/3,(A+a*a*B+a*C)/3
def prep_gen(df,k,t0=2.4,n=24,sign=-1.0):
    """df: generator csv (t,Va..Ic + slow channels).  sign=-1: file current flows INTO the machine."""
    slow=[c for c in df.columns if c not in('t','Va','Vb','Vc','Ia','Ib','Ic')]
    cn=lambda key:[c for c in slow if c.lower().startswith(key)][0]
    c_th,c_w,c_vf,c_pm=cn('theta'),cn('w'),[c for c in slow if 'exc' in c.lower()][0],[c for c in slow if c.lower().startswith('pm')][0]
    sl=df.dropna(subset=[c_th]).drop_duplicates('t').sort_values('t').reset_index(drop=True)
    fa=df.drop_duplicates('t').sort_values('t')[['t','Va','Vb','Vc','Ia','Ib','Ic']]
    t=fa.t.values; fs=n*F0; tu=np.arange(0,t.max()+1e-9,1/fs)
    wv=np.exp(-2j*np.pi*np.arange(n)/n)[::-1]*np.sqrt(2)/n
    Ph={c:np.convolve(np.interp(tu,t,fa[c].values),wv)[:len(tu)]*np.exp(-2j*np.pi*F0*tu) for c in ['Va','Vb','Vc','Ia','Ib','Ic']}
    tc=tu-(n-1)/2/fs
    V1=seq(Ph['Va'],Ph['Vb'],Ph['Vc'])[1]; I1=seq(Ph['Ia'],Ph['Ib'],Ph['Ic'])[1]
    S=sign*3*V1*np.conj(I1)
    ts=sl.t.values
    ang=np.unwrap(np.angle(V1))
    I=lambda z:np.interp(ts,tc,z)
    Vpu=I(np.abs(V1))/(VN[k]/np.sqrt(3)); thV=I(ang); P=I(S.real)/SN[k]; Q=I(S.imag)/SN[k]
    w=sl[c_w].values/sl[c_w].values[0]
    dl=np.unwrap(np.deg2rad(sl[c_th].values))*PP[k]
    out=dict(t=ts,V=Vpu,th=thV,vf=sl[c_vf].values,Pm=sl[c_pm].values,P=P,Q=Q,w=w,delta=dl,dmt=dl-thV)
    m=ts>=t0
    return {kk:v[m] for kk,v in out.items()}, out
