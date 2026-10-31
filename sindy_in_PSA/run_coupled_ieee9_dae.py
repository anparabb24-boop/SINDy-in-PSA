"""Integrate fixed fitted generator models with the assembled IEEE9 DAE.

Measured mechanical-power and field-voltage channels drive the machines. Bus
voltages are solved from algebraic KCL at every RHS evaluation. No parameters
are refitted. The output scores this coupled simulation against the separate
healthy candidate's generator terminal channels.
"""
from __future__ import annotations
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent; PKG=ROOT/"mbam_package"
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(PKG))
from prepare_healthy_data import positive_sequence_phasors
from mbam_core import NAMES,SGProblem
from mbam_run import NOMINAL
from ieee9_dae import GeneratorState,generator_rhs,constant_impedance_load_current
from build_ybus import build_ybus

SBASE=100.0; RATING={1:247.5,2:230.0,3:128.0}; VLL={1:16.5e3,2:18e3,3:13.8e3}
POLES={1:20,2:1,3:1}; OMEGA=2*np.pi*60
LOADS=((5,125,50),(6,90,30),(6,50,0),(8,100,35))

def _interp_complex(tnew,t,z):return np.interp(tnew,t,z.real)+1j*np.interp(tnew,t,z.imag)
def load_generator(gen):
    df=pd.read_csv(ROOT/f"healthy_datasets/G{gen}_10s_H.csv")
    slow=df.dropna(subset=["theta","w","Efd","Pmpu"])[["t","theta","w","Efd","Pmpu"]].copy()
    wave=df.dropna(subset=["Va","Vb","Vc","Ia","Ib","Ic"])[["t","Va","Vb","Vc","Ia","Ib","Ic"]].copy()
    wave["_key"]=wave.t.round(9);wave=wave.drop_duplicates("_key").drop(columns="_key").sort_values("t")
    vph=positive_sequence_phasors(wave,"voltage");iph=positive_sequence_phasors(wave,"current")
    ph=vph.merge(iph,on="t_center_s",validate="one_to_one")
    zv=ph.voltage_positive_re.to_numpy()+1j*ph.voltage_positive_im.to_numpy()
    zi=ph.current_positive_re.to_numpy()+1j*ph.current_positive_im.to_numpy()
    tslow=slow.t.to_numpy(float);tp=ph.t_center_s.to_numpy(float)
    vv=_interp_complex(tslow,tp,zv);ii=_interp_complex(tslow,tp,zi)
    keep=(tslow>=5)&(tslow<=10);t=tslow[keep]
    vpu=vv[keep]/(VLL[gen]/np.sqrt(3))
    ibase=SBASE*1e6/(np.sqrt(3)*VLL[gen])
    imeas=-ii[keep]/ibase # exported current; source sensors point into generator
    theta=np.unwrap(np.deg2rad(slow.theta.to_numpy(float)))*POLES[gen]
    thv=np.unwrap(np.angle(vv))
    s=-3*vv[keep]*np.conj(ii[keep])
    return {"t":t,"v_meas":vpu,"i_meas_sys":imeas,"V":np.abs(vpu),"th":thv[keep],
      "vf_raw":slow.Efd.to_numpy(float)[keep],"pm_raw":slow.Pmpu.to_numpy(float)[keep],
      "P_meas_machine":s.real/(RATING[gen]*1e6),"Q_meas_machine":s.imag/(RATING[gen]*1e6),
      "w_meas":slow.w.to_numpy(float)[keep]/(OMEGA/POLES[gen]),
      "dmt_meas":theta[keep]-thv[keep]}

def load_parameters():
    table=pd.read_csv(PKG/"validation_results"/"fitted_parameters.csv")
    params={}
    for gen in (1,2,3):
        params[gen]=NOMINAL[gen].copy()
        for r in table.loc[table.generator==f"G{gen}"].itertuples(index=False):params[gen][r.parameter]=float(r.estimate)
    return params

def initial_state(gen,d,p):
    dd0=float(d["dmt_meas"][0]+p["c"]);delta=float(d["th"][0]+dd0)
    w=float(d["w_meas"][0]);eqp=float(p["eq0"]);edp=float(p["ed0"])
    a3=(p["xdp"]-p["xdpp"])/p["xdpp"]
    a4=(p["xqp"]-p["xqpp"])/p["xqpp"]
    vd=d["V"][0]*np.sin(dd0);vq=d["V"][0]*np.cos(dd0)
    eqpp=(eqp+a3*vq)/(1+a3);edpp=(edp+a4*vd)/(1+a4)
    return GeneratorState(delta,w,eqp,edp,eqpp,edpp)

def pack_state(states):
    return np.array([getattr(states[g],n) for g in (1,2,3)
                     for n in ("delta","omega","eqp","edp","eqpp","edpp")],float)
def unpack_state(x):
    names=("delta","omega","eqp","edp","eqpp","edpp")
    return {g:GeneratorState(*x[(g-1)*6:g*6]) for g in (1,2,3)}
def pack_voltage(v):return np.r_[v.real,v.imag]
def unpack_voltage(x):return x[:9]+1j*x[9:]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rtol',type=float,default=2e-6)
    parser.add_argument('--atol',type=float,default=2e-8)
    parser.add_argument('--max-step',type=float,default=1/600)
    parser.add_argument('--tag',default='')
    args=parser.parse_args()
    outdir=PKG/"validation_results";outdir.mkdir(parents=True,exist_ok=True)
    data={g:load_generator(g) for g in (1,2,3)};t_eval=data[1]["t"]
    for g in (2,3):
        if not np.allclose(data[g]["t"],t_eval):raise ValueError("generator slow-channel grids differ")
    params=load_parameters(); ybus=build_ybus();yaug=ybus.copy()
    for b,p,q in LOADS:yaug[b-1,b-1]+=np.conj(complex(p,q)/SBASE)
    states0={g:initial_state(g,data[g],params[g]) for g in (1,2,3)}
    x0=pack_state(states0)
    # Seed transmission buses from the measured-current KCL reconstruction.
    vg0=np.array([data[g]["v_meas"][0] for g in (1,2,3)])
    vnet0=np.linalg.solve(yaug[3:,3:],-yaug[3:,:3]@vg0)
    vguess=np.r_[vg0,vnet0]
    counters={"algebraic_calls":0,"algebraic_failures":0,"max_root_residual":0.0}
    def solve_voltage(states,t,guess):
        def residual(x):
            v=unpack_voltage(x);igen=np.zeros(9,complex)
            for g in (1,2,3):
                _,ig,_,_=generator_rhs(states[g],v[g-1],
                    params[g]["kf"]*np.interp(t,data[g]["t"],data[g]["vf_raw"]),
                    params[g]["km"]*np.interp(t,data[g]["t"],data[g]["pm_raw"]),
                    params[g],OMEGA)
                igen[g-1]=ig*(RATING[g]/SBASE)
            rg=yaug@v-igen
            return pack_voltage(rg)
        sol=root(residual,pack_voltage(guess),method="hybr",options={"xtol":1e-10})
        rn=float(np.linalg.norm(residual(sol.x),ord=np.inf))
        counters["algebraic_calls"]+=1;counters["max_root_residual"]=max(counters["max_root_residual"],rn)
        if not sol.success and rn>1e-8:counters["algebraic_failures"]+=1
        return unpack_voltage(sol.x),rn
    vcache=[vguess.copy()]
    def rhs(t,x):
        states=unpack_state(x);v,rn=solve_voltage(states,t,vcache[0]);vcache[0]=v
        dx=[]
        for g in (1,2,3):
            dstate,_,_,_=generator_rhs(states[g],v[g-1],
                params[g]["kf"]*np.interp(t,data[g]["t"],data[g]["vf_raw"]),
                params[g]["km"]*np.interp(t,data[g]["t"],data[g]["pm_raw"]),
                params[g],OMEGA)
            dx.extend(getattr(dstate,n) for n in ("delta","omega","eqp","edp","eqpp","edpp"))
        return np.asarray(dx)
    sol=solve_ivp(rhs,(t_eval[0],t_eval[-1]),x0,method="RK45",t_eval=t_eval,
                  rtol=args.rtol,atol=args.atol,max_step=args.max_step)
    if not sol.success:raise RuntimeError(sol.message)
    rows=[]; voltage_guess=vcache[0]
    for k,t in enumerate(sol.t):
        states=unpack_state(sol.y[:,k]);v,rn=solve_voltage(states,t,voltage_guess);voltage_guess=v
        igen=np.zeros(9,complex);powers={};dstate_map={}
        for g in (1,2,3):
            vf=params[g]["kf"]*np.interp(t,data[g]["t"],data[g]["vf_raw"])
            pm=params[g]["km"]*np.interp(t,data[g]["t"],data[g]["pm_raw"])
            dstate,ig,pg,qg=generator_rhs(states[g],v[g-1],vf,pm,params[g],OMEGA)
            igen[g-1]=ig*(RATING[g]/SBASE);powers[g]=(pg,qg)
            dstate_map[g]=dstate
        residual=yaug@v-igen
        row={"t_s":float(t),"algebraic_solver_residual_inf_pu":rn}
        for b in range(1,10):
            row[f"V{b}_re_pu"]=v[b-1].real;row[f"V{b}_im_pu"]=v[b-1].imag;row[f"V{b}_mag_pu"]=abs(v[b-1])
            row[f"KCL{b}_re_pu"]=residual[b-1].real;row[f"KCL{b}_im_pu"]=residual[b-1].imag;row[f"KCL{b}_mag_pu"]=abs(residual[b-1])
        for g in (1,2,3):
            d=data[g];state=states[g];ig=igen[g-1];pg,qg=powers[g]
            idx=int(np.argmin(np.abs(d["t"]-t)))
            dmt_pred=state.delta-np.angle(v[g-1])-params[g]["c"]
            row.update({f"G{g}_Imeas_re_sys_pu":d["i_meas_sys"][idx].real,
                        f"G{g}_Imeas_im_sys_pu":d["i_meas_sys"][idx].imag,
                        f"G{g}_Ipred_re_sys_pu":ig.real,f"G{g}_Ipred_im_sys_pu":ig.imag,
                        f"G{g}_Ierr_mag_sys_pu":abs(ig-d["i_meas_sys"][idx]),
                        f"G{g}_Vmeas_re_pu":d["v_meas"][idx].real,f"G{g}_Vmeas_im_pu":d["v_meas"][idx].imag,
                        f"G{g}_Vpred_re_pu":v[g-1].real,f"G{g}_Vpred_im_pu":v[g-1].imag,
                        f"G{g}_Pmeas_machine_pu":d["P_meas_machine"][idx],f"G{g}_Ppred_machine_pu":pg,
                        f"G{g}_Qmeas_machine_pu":d["Q_meas_machine"][idx],f"G{g}_Qpred_machine_pu":qg,
                        f"G{g}_Vmeas_mag_pu":abs(d["v_meas"][idx]),f"G{g}_Vpred_mag_pu":abs(v[g-1]),
                        f"G{g}_Verr_mag_pu":abs(v[g-1])-abs(d["v_meas"][idx]),
                        f"G{g}_delta_state_rad":state.delta,
                        f"G{g}_w_meas_pu":d["w_meas"][idx],f"G{g}_w_pred_pu":state.omega,
                        f"G{g}_delta_minus_theta_meas_rad":d["dmt_meas"][idx],
                        f"G{g}_delta_minus_theta_pred_rad":dmt_pred})
        rows.append(row)
    frame=pd.DataFrame(rows)
    metrics=[]
    for b in range(1,10):
        e=frame[f"KCL{b}_mag_pu"].to_numpy()
        metrics.append({"quantity":"bus_kcl","asset":f"Bus{b}","n":len(e),"rmse":float(np.sqrt(np.mean(e**2))),"max_abs":float(e.max()),
                        "measurement_status":"generator buses voltage measured; transmission buses voltage reconstructed" if b<=3 else "transmission voltage reconstructed; no independent healthy bus channel"})
    for g in (1,2,3):
        phase=np.unwrap(np.angle(frame[f"G{g}_Vpred_re_pu"].to_numpy()+1j*frame[f"G{g}_Vpred_im_pu"].to_numpy()))
        measured_phase=data[g]["th"]
        phase+=2*np.pi*np.round((measured_phase[0]-phase[0])/(2*np.pi))
        frame[f"G{g}_Vpred_phase_aligned_rad"]=phase
        frame[f"G{g}_delta_minus_theta_pred_rad"]=frame[f"G{g}_delta_state_rad"].to_numpy()-phase-params[g]["c"]
        vp=frame[f"G{g}_Vpred_re_pu"].to_numpy()+1j*frame[f"G{g}_Vpred_im_pu"].to_numpy()
        vm=frame[f"G{g}_Vmeas_re_pu"].to_numpy()+1j*frame[f"G{g}_Vmeas_im_pu"].to_numpy()
        ip=frame[f"G{g}_Ipred_re_sys_pu"].to_numpy()+1j*frame[f"G{g}_Ipred_im_sys_pu"].to_numpy()
        im=frame[f"G{g}_Imeas_re_sys_pu"].to_numpy()+1j*frame[f"G{g}_Imeas_im_sys_pu"].to_numpy()
        for label,col in (("current_complex_error_system_pu",np.abs(ip-im)),
                          ("current_magnitude_error_system_pu",np.abs(ip)-np.abs(im)),
                          ("current_phase_error_rad",np.angle(ip*np.conj(im))),
                          ("terminal_voltage_complex_error_pu",np.abs(vp-vm)),
                          ("terminal_voltage_magnitude_error_pu",np.abs(vp)-np.abs(vm)),
                          ("terminal_voltage_phase_error_rad",np.angle(vp*np.conj(vm))),
                          ("active_power_error_machine_pu",None),("reactive_power_error_machine_pu",None),
                          ("rotor_speed_error_pu",None),("rotor_relative_angle_error_rad",None)):
            if isinstance(col,str):err=frame[col].to_numpy()
            elif isinstance(col,np.ndarray):err=col
            elif label.startswith("active"):err=frame[f"G{g}_Ppred_machine_pu"].to_numpy()-frame[f"G{g}_Pmeas_machine_pu"].to_numpy()
            elif label.startswith("reactive"):err=frame[f"G{g}_Qpred_machine_pu"].to_numpy()-frame[f"G{g}_Qmeas_machine_pu"].to_numpy()
            elif label.startswith("rotor_speed"):err=frame[f"G{g}_w_pred_pu"].to_numpy()-frame[f"G{g}_w_meas_pu"].to_numpy()
            else:err=frame[f"G{g}_delta_minus_theta_pred_rad"].to_numpy()-frame[f"G{g}_delta_minus_theta_meas_rad"].to_numpy()
            metrics.append({"quantity":label,"asset":f"G{g}","n":len(err),"rmse":float(np.sqrt(np.mean(err**2))),"mean_error":float(np.mean(err)),"max_abs":float(np.max(np.abs(err))),
                            "measurement_status":"direct healthy generator channel comparison"})
    tag=f"_{args.tag}" if args.tag else ""
    frame.to_csv(outdir/f"coupled_dae_healthy_candidate_5s_10s{tag}.csv",index=False)
    metrics_df=pd.DataFrame(metrics);metrics_df.to_csv(outdir/f"coupled_dae_healthy_candidate_metrics{tag}.csv",index=False)
    (outdir/f"coupled_dae_solver_diagnostics{tag}.json").write_text(__import__("json").dumps({"success":bool(sol.success),"message":sol.message,
      "ode_nfev":int(sol.nfev),"algebraic_calls":counters["algebraic_calls"],"algebraic_failures":counters["algebraic_failures"],
      "max_algebraic_root_residual_inf_pu":counters["max_root_residual"],"system_base_mva":SBASE,
      "ode_rtol":args.rtol,"ode_atol":args.atol,"ode_max_step_s":args.max_step,
      "parameter_refit":False,"controls":"measured healthy-run Efd and Pm time series with previously fitted kf/km scaling"},indent=2),encoding="utf-8")
    print(metrics_df.to_string(index=False,float_format=lambda x:f"{x:.7g}"))
    print("solver",sol.success,sol.message,"nfev",sol.nfev,"algebraic",counters)

if __name__=="__main__":main()
