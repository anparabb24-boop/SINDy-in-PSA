import unittest
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ieee9_dae import (constant_impedance_load_current, loads_pu_from_mw,
                       machine_impedance_to_system_base, GeneratorState,
                       generator_rhs, IEEE9DAE)

class DAETests(unittest.TestCase):
    def test_constant_impedance_current_scales_linearly(self):
        v=0.98*np.exp(1j*0.2); s=0.7+0.25j
        i=constant_impedance_load_current(v,s.real,s.imag)
        self.assertAlmostEqual((v*np.conj(i)).real,s.real*abs(v)**2)
        self.assertAlmostEqual((v*np.conj(i)).imag,s.imag*abs(v)**2)
    def test_load_ratings_aggregate_and_convert_to_system_base(self):
        rows=[{'bus':6,'p_mw':90,'q_mvar':30},{'bus':6,'p_mw':50,'q_mvar':0}]
        self.assertEqual(loads_pu_from_mw(rows)[6],1.4+0.3j)
    def test_machine_impedance_base_conversion(self):
        self.assertAlmostEqual(machine_impedance_to_system_base(.2,200,100),.1)
    def test_generator_rhs_uses_validated_h2_parameter(self):
        p={'xdpp':.25,'xqpp':.25,'H2':6.0,'D':0.0,'xd':1.8,'xdp':.3,
           'xq':1.7,'xqp':.3,'Td0p':5.0,'Tq0p':.5,'Td0pp':.05,'Tq0pp':.05}
        s=GeneratorState(0.0,1.0,1.0,0.0,1.0,0.0)
        rhs,_,_,_=generator_rhs(s,1+0j,1.0,.8,p)
        self.assertEqual(len((rhs.delta,rhs.omega,rhs.eqp,rhs.edp,rhs.eqpp,rhs.edpp)),6)
    def test_algebraic_residual_shape(self):
        dae=IEEE9DAE(np.eye(9,dtype=complex),{5:0.1+0.02j})
        v=np.ones(9,dtype=complex); gen=np.zeros(9,dtype=complex)
        r=dae.algebraic_residual(v,gen)
        self.assertEqual(r.shape,(9,)); self.assertTrue(np.isfinite(r).all())

if __name__=='__main__': unittest.main()
