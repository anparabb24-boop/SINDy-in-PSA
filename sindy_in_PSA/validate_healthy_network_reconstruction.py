"""Reconstruct unmeasured network bus voltages over the new no-fault run.

Generator bus voltages are measured. Generator injections are evaluated both
from measured current phasors and from the previously fitted machine models.
The six unmeasured transmission bus voltages are solved from network KCL; the
three generator-bus KCL residuals are the withheld algebraic consistency check.
This is a boundary-driven DAE reconstruction, not a free-running closed-loop
simulation or direct validation against unmeasured bus voltages.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent; PKG=ROOT/"mbam_package"
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(PKG))
from prepare_healthy_data import positive_sequence_phasors
from mbam_core import NAMES,SGProblem
from mbam_run import NOMINAL
from build_ybus import build_ybus

SBASE=100.0; RATING={1:247.5,2:230.0,3:128.0}; VLL={1:16.5e3,2:18e3,3:13.8e3}
POLES={1:20,2:1,3:1}; OMEGA=2*np.pi*60
LOADS=((5,125,50),(6,90,30),(6,50,0),(8,100,35))

def interp_complex(tnew,t,z):return np.interp(tnew,t,z.real)+1j*np.interp(tnew,t,z.imag)
def load_gen(gen):
    df=pd.read_csv(ROOT/f"healthy_datasets/G{gen}_10s_H.csv")
    slow=df.dropna(subset=["theta","w","Efd","Pmpu"])[["t","theta","w","Efd","Pmpu"]].copy()
    wave=df.dropna(subset=["Va","Vb","Vc","Ia","Ib","Ic"])[["t","Va","Vb","Vc","Ia","Ib","Ic"]].copy()
    wave["key"]=wave.t.round(9);wave=wave.drop_duplicates("key").drop(columns="key").sort_values("t")
    v=positive_sequence_phasors(wave,"voltage");i=positive_sequence_phasors(wave,"current")
    ph=v.merge(i,on="t_center_s",validate="one_to_one")
    zv=ph.voltage_positive_re.to_numpy()+1j*ph.voltage_positive_im.to_numpy()
    zi=ph.current_positive_re.to_numpy()+1j*ph.current_positive_im.to_numpy()
    t=slow.t.to_numpy(float);tp=ph.t_center_s.to_numpy(float)
    vv=interp_complex(t,tp,zv);ii=interp_complex(t,tp,zi)
    keep=(t>=5)&(t<=10);t=t[keep]
    vpu=vv[keep]/(VLL[gen]/np.sqrt(3))
    # Exported generator current is positive into the machine; negate for network injection.
    ibase=SBASE*1e6/(np.sqrt(3)*VLL[gen])
    imeas=-ii[keep]/ibase
    d={"t":t,"V":np.abs(vpu),"th":np.unwrap(np.angle(vv))[keep],
       "vf":slow.Efd.to_numpy(float)[keep],"Pm":slow.Pmpu.to_numpy(float)[keep],
       "P":(-3*vv[keep]*np.conj(ii[keep])).real/(RATING[gen]*1e6),
       "Q":(-3*vv[keep]*np.conj(ii[keep])).imag/(RATING[gen]*1e6),
       "w":slow.w.to_numpy(float)[keep]/(OMEGA/POLES[gen]),
       "dmt":np.unwrap(np.deg2rad(slow.theta.to_numpy(float)))[keep]*POLES[gen]-np.unwrap(np.angle(vv))[keep]}
    return d, vpu, imeas

def main():
    OUT=PKG/"validation_results";OUT.mkdir(parents=True,exist_ok=True)
    fits=pd.read_csv(OUT/"fitted_parameters.csv"); y=build_ybus(); yaug=y.copy()
    for bus,p,q in LOADS:yaug[bus-1,bus-1]+=np.conj(complex(p,q)/SBASE)
    network_matrix=yaug[3:,3:]
    rows=[]; voltage_rows=[]; data_by_gen={}; voltage_by_gen={}; measured_by_gen={}; model_by_gen={}
    for gen in (1,2,3):
        d,v,imeas=load_gen(gen); pdata=NOMINAL[gen].copy()
        for r in fits.loc[fits.generator==f"G{gen}"].itertuples(index=False):pdata[r.parameter]=float(r.estimate)
        prob=SGProblem(d['t'],d['V'],d['th'],d['vf'],d['Pm'],d['P'],d['Q'],d['w'],d['dmt'],pdata,OMEGA,
                       salient=(gen==1),ns=1,sig_frac=.005)
        pp=prob.simulate(np.array([pdata[n] for n in NAMES])); pmach=pp[0]+1j*pp[1]
        # Undo the fitted angle output offset only for the scored angle metric; it does not affect P/Q.
        imodel=np.conj((pmach*(RATING[gen]/SBASE))/v)
        data_by_gen[gen]=d;voltage_by_gen[gen]=v;measured_by_gen[gen]=imeas;model_by_gen[gen]=imodel
    t=data_by_gen[1]['t']; n=len(t)
    measured=np.column_stack([measured_by_gen[g] for g in (1,2,3)])
    model=np.column_stack([model_by_gen[g] for g in (1,2,3)])
    vg=np.column_stack([voltage_by_gen[g] for g in (1,2,3)])
    for key,igen in (("measured_current",measured),("generator_model",model)):
        residual=[]; vnet=[]
        for k in range(n):
            rhs=-yaug[3:,:3]@vg[k]
            vn=np.linalg.solve(network_matrix,rhs)
            fullv=np.r_[vg[k],vn]
            igen9=np.r_[igen[k],np.zeros(6,complex)]
            residual.append(yaug[:3,:]@fullv-igen9[:3]);vnet.append(vn)
        residual=np.asarray(residual);vnet=np.asarray(vnet)
        for b in (1,2,3):
            r=residual[:,b-1]
            rows.append({"injection_source":key,"generator_bus":b,"n":n,
                         "rms_kcl_residual_system_pu":float(np.sqrt(np.mean(abs(r)**2))),
                         "max_kcl_residual_system_pu":float(np.max(abs(r))),
                         "rms_kcl_residual_mva_equivalent":float(np.sqrt(np.mean(abs(r)**2))*SBASE)})
        for j,b in enumerate(range(4,10)):
            vm=np.abs(vnet[:,j])
            voltage_rows.append({"injection_source":key,"bus":b,"n":n,"vm_min_pu":float(vm.min()),
                                 "vm_mean_pu":float(vm.mean()),"vm_max_pu":float(vm.max())})
    pd.DataFrame(rows).to_csv(OUT/"healthy_network_reconstruction_metrics.csv",index=False)
    pd.DataFrame(voltage_rows).to_csv(OUT/"healthy_network_reconstructed_bus_voltage_ranges.csv",index=False)
    print(pd.DataFrame(rows).to_string(index=False,float_format=lambda x:f"{x:.7g}"))
    print(pd.DataFrame(voltage_rows).to_string(index=False,float_format=lambda x:f"{x:.6g}"))

if __name__=="__main__":main()
