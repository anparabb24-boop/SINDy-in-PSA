"""Check the independently tabulated operating point against the built Ybus."""
from __future__ import annotations
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
import numpy as np
import pandas as pd
from build_ybus import build_ybus, BUS_ORDER

ROOT=Path(__file__).resolve().parent
def main():
    table=pd.read_csv(ROOT/"IEEE9BUS_loadflow_results.csv")
    buses=table.loc[table.BlockType=="Busbar"]
    voltage=[]
    for bus in BUS_ORDER:
        row=buses.loc[buses.Row.str.startswith(f"Bus{bus} ")].iloc[0]
        voltage.append(float(row.ActualVoltageMagnitude_Pu)*np.exp(1j*np.deg2rad(float(row.VoltageAngle_Deg))))
    voltage=np.asarray(voltage)
    base=100.0
    loads=np.zeros(9,complex)
    # Fixed RLC load powers scale with |V|^2 from rated-voltage settings.
    for bus,p,q in ((5,125,50),(6,90,30),(6,50,0),(8,100,35)):
        loads[bus-1]+=(p+1j*q)*abs(voltage[bus-1])**2/base
    generation=np.zeros(9,complex)
    for row in table.loc[table.BlockType=="Synchronous Machine"].itertuples(index=False):
        bus=int(row.Row.split("@Bus")[1].split()[0])
        generation[bus-1]=(float(row.ActualGenerationP_MW)+1j*float(row.ActualGenerationQ_Mvar))/base
    y=build_ybus(base_mva=base)
    mismatch=(voltage*np.conj(y@voltage)+loads-generation)*base
    transformers=pd.read_csv(ROOT/"transformer_parameters.csv")
    transformers["phase_shift_deg"]=-transformers["phase_shift_deg"]
    with NamedTemporaryFile(suffix=".csv") as alternate_file:
        transformers.to_csv(alternate_file.name,index=False)
        y_alt=build_ybus(base_mva=base,transformer_csv=alternate_file.name)
    mismatch_alt=(voltage*np.conj(y_alt@voltage)+loads-generation)*base
    report={"bus_order":list(BUS_ORDER),"system_base_mva":base,
            "validation_reference":"IEEE9BUS_loadflow_results.csv bus voltage, generation and demand snapshot",
            "load_model":"fixed impedance; P,Q rated at 1 pu voltage; 50 MW block is resistive (Q=0)",
            "max_abs_bus_mismatch_mva":float(np.max(np.abs(mismatch))),
            "rms_complex_bus_mismatch_mva":float(np.sqrt(np.mean(np.abs(mismatch)**2))),
            "opposite_phase_max_abs_bus_mismatch_mva":float(np.max(np.abs(mismatch_alt))),
            "opposite_phase_rms_complex_bus_mismatch_mva":float(np.sqrt(np.mean(np.abs(mismatch_alt)**2))),
            "bus_mismatch_mva_mvar":[{"bus":b,"p_mw":float(m.real),"q_mvar":float(m.imag)} for b,m in zip(BUS_ORDER,mismatch)],
            "ybus_pu_rows":[[f"{v.real:.6f}{v.imag:+.6f}j" for v in row] for row in y]}
    out=ROOT/"mbam_package"/"validation_results"/"ybus_reference_check.json"
    out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("ybus_pu_rows","bus_mismatch_mva_mvar")},indent=2))
    print(json.dumps(report["bus_mismatch_mva_mvar"],indent=2))
    print(np.array2string(y,precision=5,suppress_small=True))
if __name__=="__main__":main()
