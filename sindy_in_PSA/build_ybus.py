"""Build the positive-sequence Ybus from the source-derived network data.

The common 100 MVA base is an explicit analyst-selected reference. Transformer
phase_shift_deg is the angle of V_HV/V_LV at no load. Delta 1 clock means the
LV winding lags the HV winding by 30 degrees, so this tap angle is +30 degrees.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BUS_ORDER = tuple(range(1, 10))

def _network_data():
    params = json.loads((HERE / "network_parameters.json").read_text())
    sbase = params["system_base_mva"]
    if not isinstance(sbase, (int, float)) or sbase <= 0:
        raise ValueError("network_parameters.json must declare a positive DAE reference base")
    vb = np.array([params["voltage_bases_kv_ll"][str(b)] for b in BUS_ORDER], float)
    return float(sbase), vb

def line_positive_sequence(row):
    """Return full-line (R1 ohm, X1 ohm, B1 siemens), retaining mutual coupling.

    For the documented three-phase line reduction, R1=R-Rm, L1=L-M,
    C1=Cg+3Cl. Total shunt susceptance is split across both pi ends.
    """
    length = float(row.length_km)
    resistance = (float(row.R_ohm_per_km)-float(row.Rm_ohm_per_km))*length
    inductance = (float(row.Lself_mH_per_km)-float(row.M_mH_per_km))*1e-3*length
    capacitance = (float(row.Cg_uF_per_km)+3*float(row.Cl_uF_per_km))*1e-6*length
    f = float(row.frequency_hz)
    return resistance, 2*np.pi*f*inductance, 2*np.pi*f*capacitance

def build_ybus(base_mva=None, line_csv=None, transformer_csv=None, units="pu"):
    """Return Ybus in bus order 1..9.

    `units` is `pu` or `siemens`. Per-unit conversion uses each bus's LL voltage
    base and the declared common system MVA base. The transformer tap stamp uses
    the standard HV-side complex tap convention V_HV/a=V_LV.
    """
    declared_base, vb = _network_data()
    sbase = declared_base if base_mva is None else float(base_mva)
    if sbase <= 0: raise ValueError("base_mva must be positive")
    if units not in ("pu", "siemens"): raise ValueError("units must be 'pu' or 'siemens'")
    lines = pd.read_csv(line_csv or HERE/"branch_parameters.csv")
    transformers = pd.read_csv(transformer_csv or HERE/"transformer_parameters.csv")
    y = np.zeros((9,9),complex); idx={b:b-1 for b in BUS_ORDER}
    for row in lines.itertuples(index=False):
        i,j=idx[int(row.from_bus)],idx[int(row.to_bus)]
        r,x,b=line_positive_sequence(row)
        zbase=vb[i]**2/sbase
        ys=zbase/complex(r,x)
        bpu=b*zbase
        y[i,i]+=ys+0.5j*bpu; y[j,j]+=ys+0.5j*bpu
        y[i,j]-=ys; y[j,i]-=ys
    for row in transformers.itertuples(index=False):
        if pd.isna(row.phase_shift_deg):
            raise ValueError(f"{row.name}: transformer phase shift is unresolved")
        i,j=idx[int(row.hv_bus)],idx[int(row.lv_bus)]
        own_z=complex(float(row.rw1_pu_own_base+row.rw2_pu_own_base),
                      float(row.xl1_pu_own_base+row.xl2_pu_own_base))
        # Convert transformer pu impedance from its 100 MVA rating and rated
        # winding voltage to the selected common base and bus voltage base.
        own_vbase=float(row.hv_kv_ll)
        zpu=own_z*(sbase/float(row.rating_mva))*(own_vbase/vb[i])**2
        ys=1/zpu
        tap=float(row.tap_pu)*np.exp(1j*np.deg2rad(float(row.phase_shift_deg)))
        y[i,i]+=ys/(abs(tap)**2)
        y[i,j]-=ys/np.conj(tap)
        y[j,i]-=ys/tap
        y[j,j]+=ys
        # Linear magnetizing Rm and Xm are connected on winding 1 (HV).
        y[i,i]+=(complex(1/float(row.rm_pu_own_base),-1/float(row.xm_pu_own_base))
                 *sbase/float(row.rating_mva)*(own_vbase/vb[i])**2)
    if units=="pu": return y
    return y*sbase/np.outer(vb,vb)

def branch_stamp(from_bus,to_bus,y_series,b_shunt_pu=0j,tap=1+0j):
    """Independent 2x2 pi-branch stamp, useful for verification and reuse."""
    sh=0.5*b_shunt_pu
    return np.array([[(y_series+sh)/abs(tap)**2,-y_series/np.conj(tap)],
                     [-y_series/tap,y_series+sh]],complex)

if __name__=="__main__":
    print(build_ybus())
