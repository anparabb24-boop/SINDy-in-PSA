"""Modular positive-sequence IEEE 9-bus network DAE.

The selected 100 MVA base is an analyst modelling reference. Generators retain
their machine-base differential equations; terminal currents are converted to
the system base for network KCL. Wye RLC loads are represented as constant-Z
at nominal frequency, consistent with the fixed passive components.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from build_ybus import BUS_ORDER

@dataclass
class GeneratorState:
    delta: float
    omega: float
    eqp: float
    edp: float
    eqpp: float
    edpp: float

def constant_impedance_load_current(v, p_pu, q_pu, v_nominal_pu=1.0):
    """Current drawn by rated-power fixed RLC load at rated voltage/frequency."""
    if v_nominal_pu<=0: raise ValueError('v_nominal_pu must be positive')
    y_load=np.conj(complex(p_pu,q_pu))/(v_nominal_pu**2)
    return y_load*v

def loads_pu_from_mw(load_rows, system_base_mva=100.0):
    """Aggregate SLX load ratings into constant-impedance Spu by bus."""
    if system_base_mva<=0: raise ValueError('system_base_mva must be positive')
    out={}
    for row in load_rows:
        bus=int(row['bus'])
        out[bus]=out.get(bus,0j)+complex(float(row['p_mw']),float(row['q_mvar']))/system_base_mva
    return out

def machine_impedance_to_system_base(x_pu, machine_mva, system_mva=100.0,
                                    machine_kv=None, system_kv=None):
    """Convert a machine-base pu impedance to the network base."""
    if machine_mva<=0 or system_mva<=0: raise ValueError('MVA bases must be positive')
    voltage_factor=1.0
    if machine_kv is not None or system_kv is not None:
        if machine_kv is None or system_kv is None or machine_kv<=0 or system_kv<=0:
            raise ValueError('supply both positive voltage bases or neither')
        voltage_factor=(machine_kv/system_kv)**2
    return float(x_pu)*(system_mva/machine_mva)*voltage_factor

def generator_terminal(state: GeneratorState, v, p):
    """Subtransient stator relation consistent with mbam_core's dq convention.

    Returns injected current and electrical P/Q on the machine base. The caller
    converts terminal current to the common network base.
    """
    dd=state.delta-np.angle(v); vm=abs(v)
    vd=vm*np.sin(dd); vq=vm*np.cos(dd)
    id_=(state.eqpp-vq)/p['xdpp']; iq=(vd-state.edpp)/p['xqpp']
    pg=vd*id_+vq*iq; qg=vq*id_-vd*iq
    # Network positive sequence complex convention consistent with S=V*conj(I).
    s=complex(pg,qg); iinj=np.conj(s/v)
    return iinj,pg,qg

def generator_rhs(state: GeneratorState, v, vf, pm, p, omega_base=2*np.pi*60):
    iinj,pg,qg=generator_terminal(state,v,p)
    inertia_2h=p.get('H2',p.get('H'))
    if inertia_2h is None or inertia_2h<=0:
        raise ValueError("generator parameters must include positive H2 (2H inertia)")
    dd=state.delta-np.angle(v); vd=abs(v)*np.sin(dd); vq=abs(v)*np.cos(dd)
    id_=(state.eqpp-vq)/p['xdpp']; iq=(vd-state.edpp)/p['xqpp']
    d=GeneratorState(
      omega_base*(state.omega-1),
      (pm-pg-p.get('D',0)*(state.omega-1))/inertia_2h,
      (vf-state.eqp-(p['xd']-p['xdp'])*id_)/p['Td0p'],
      (-state.edp+(p['xq']-p['xqp'])*iq)/p['Tq0p'],
      (-state.eqpp+state.eqp-(p['xdp']-p['xdpp'])*id_)/p['Td0pp'],
      (-state.edpp+state.edp+(p['xqp']-p['xqpp'])*iq)/p['Tq0pp'])
    return d, iinj, pg, qg

class IEEE9DAE:
    """Algebraic current-balance assembly: YV + Iload(V) - Igen(x,V)=0."""
    def __init__(self, ybus, loads, generator_buses=(1,2,3), system_base_mva=100.0,
                 generator_ratings_mva=None):
        self.ybus=np.asarray(ybus,complex)
        if self.ybus.shape!=(9,9): raise ValueError('Ybus must be 9x9 in bus order 1..9')
        self.loads={int(k):complex(v) for k,v in loads.items()}
        self.generator_buses=tuple(generator_buses)
        self.system_base_mva=float(system_base_mva)
        if self.system_base_mva<=0: raise ValueError('system_base_mva must be positive')
        self.generator_ratings_mva={int(k):float(v) for k,v in (generator_ratings_mva or {}).items()}
    def algebraic_residual(self, voltage, generator_injections):
        voltage=np.asarray(voltage,complex); gen=np.asarray(generator_injections,complex)
        if voltage.shape!=(9,) or gen.shape!=(9,): raise ValueError('voltage and current vectors must have length 9')
        load=np.zeros(9,complex)
        for bus,s in self.loads.items(): load[bus-1]=constant_impedance_load_current(voltage[bus-1],s.real,s.imag)
        return self.ybus@voltage+load-gen
    def differential_rhs(self, states_by_bus, voltage, controls, parameters_by_bus):
        out={}; current=np.zeros(9,complex); power={}
        for bus,state in states_by_bus.items():
            if bus not in self.generator_ratings_mva:
                raise ValueError(f'missing MVA rating for generator bus {bus}')
            rhs,i,p,q=generator_rhs(state,voltage[bus-1],**controls[bus],p=parameters_by_bus[bus])
            ratio=self.generator_ratings_mva[bus]/self.system_base_mva
            out[bus]=rhs; current[bus-1]=i*ratio; power[bus]=(p*ratio,q*ratio)
        return out,current,power

def dae_residual(xdot, f, voltage, network_residual):
    """Return stacked dx/dt-f and real/imag algebraic residuals."""
    return np.r_[np.asarray(xdot)-np.asarray(f), np.real(network_residual), np.imag(network_residual)]
