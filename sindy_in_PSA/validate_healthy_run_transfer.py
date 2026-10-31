"""Score the previously fitted generator models on the newly added clean run.

This is an independent-run generator-input trajectory check, not a closed-loop
network DAE validation: measured generator terminal voltage remains an input.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
PKG=ROOT/"mbam_package"
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(PKG))
from prepare_healthy_data import positive_sequence_phasors
from mbam_core import IDX,NAMES,SGProblem
from mbam_run import NOMINAL

DATA=ROOT/"healthy_datasets"
OUT=PKG/"validation_results"
SN={1:247.5e6,2:230e6,3:128e6}; VN={1:16.5e3,2:18e3,3:13.8e3}
POLES={1:20,2:1,3:1}; OMEGA=2*np.pi*60

def interp_complex(tnew,t,z):
    return np.interp(tnew,t,z.real)+1j*np.interp(tnew,t,z.imag)

def load_new_run(gen):
    df=pd.read_csv(DATA/f"G{gen}_10s_H.csv")
    slow=df.dropna(subset=["theta","w","Efd","Pmpu"])[["t","theta","w","Efd","Pmpu"]].copy()
    wave=df.dropna(subset=["Va","Vb","Vc","Ia","Ib","Ic"])[["t","Va","Vb","Vc","Ia","Ib","Ic"]].copy()
    # The logger duplicates waveform rows at timestamps where slow channels
    # were joined. Remove only exact rounded-time duplicates before phasoring.
    wave["_time_key"]=wave.t.round(9)
    wave=wave.drop_duplicates("_time_key",keep="first").drop(columns="_time_key").sort_values("t")
    v=positive_sequence_phasors(wave,"voltage")
    i=positive_sequence_phasors(wave,"current")
    ph=v.merge(i,on="t_center_s",validate="one_to_one")
    z_v=ph.voltage_positive_re.to_numpy()+1j*ph.voltage_positive_im.to_numpy()
    z_i=ph.current_positive_re.to_numpy()+1j*ph.current_positive_im.to_numpy()
    t=slow.t.to_numpy(float); tph=ph.t_center_s.to_numpy(float)
    vv=interp_complex(t,tph,z_v); ii=interp_complex(t,tph,z_i)
    s=-3*vv*np.conj(ii)
    keep=(t>=5)&(t<=10)
    theta=np.unwrap(np.deg2rad(slow.theta.to_numpy(float)))*POLES[gen]
    angle=np.unwrap(np.angle(vv))
    d={"t":t[keep],"V":np.abs(vv[keep])/(VN[gen]/np.sqrt(3)),"th":angle[keep],
       "vf":slow.Efd.to_numpy(float)[keep],"Pm":slow.Pmpu.to_numpy(float)[keep],
       "P":s.real[keep]/SN[gen],"Q":s.imag[keep]/SN[gen],
       "w":slow.w.to_numpy(float)[keep]/(OMEGA/POLES[gen]),"dmt":theta[keep]-angle[keep]}
    return d

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    fitted=pd.read_csv(OUT/"fitted_parameters.csv")
    rows=[]
    for gen in (1,2,3):
        d=load_new_run(gen); p=NOMINAL[gen].copy()
        part=fitted.loc[fitted.generator==f"G{gen}"]
        for r in part.itertuples(index=False): p[r.parameter]=float(r.estimate)
        prob=SGProblem(d['t'],d['V'],d['th'],d['vf'],d['Pm'],d['P'],d['Q'],d['w'],d['dmt'],
                       p,OMEGA,salient=(gen==1),ns=1,sig_frac=.005)
        pred=list(prob.simulate(np.array([p[n] for n in NAMES])))
        pred[3]=pred[3]-p['c']
        for name,y,yhat in zip(("P_pu","Q_pu","w_pu","delta_minus_theta_rad"),
                               (d['P'],d['Q'],d['w'],d['dmt']),pred):
            rmse=float(np.sqrt(np.mean((yhat-y)**2))); span=float(np.ptp(y)); mean=float(np.mean(y))
            rows.append({"generator":f"G{gen}","output":name,"n":len(y),"rmse":rmse,
                         "mean_observed":mean,"normalization_scale_full_run_range":span,
                         "nrmse_range":rmse/max(span,1e-12),
                         "relative_rmse_mean_magnitude":rmse/max(abs(mean),1e-6),
                         "fit_parameters_reused_without_refit":True})
    result=pd.DataFrame(rows)
    result.to_csv(OUT/"independent_healthy_run_generator_transfer_metrics.csv",index=False)
    print(result.to_string(index=False,float_format=lambda x:f"{x:.6g}"))
    print("maximum_relative_rmse_mean_magnitude",float(result.relative_rmse_mean_magnitude.max()))

if __name__=="__main__": main()
