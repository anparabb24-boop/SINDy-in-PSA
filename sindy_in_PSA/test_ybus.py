import unittest
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_ybus import line_positive_sequence, build_ybus, branch_stamp

class YbusTests(unittest.TestCase):
    def test_positive_sequence_reduction_and_pi_capacitance(self):
        row=SimpleNamespace(length_km=10,R_ohm_per_km=0.2,Rm_ohm_per_km=0.05,
          Lself_mH_per_km=3.0,M_mH_per_km=0.5,Cl_uF_per_km=0.002,Cg_uF_per_km=0.006,frequency_hz=60)
        r,x,b=line_positive_sequence(row)
        self.assertAlmostEqual(r,1.5)
        self.assertAlmostEqual(x,2*np.pi*60*0.025)
        self.assertAlmostEqual(b,2*np.pi*60*1.2e-7)
    def test_network_ybus_is_finite_and_expected_shape(self):
        y=build_ybus()
        self.assertEqual(y.shape,(9,9))
        self.assertTrue(np.isfinite(y).all())
        self.assertTrue(np.all(np.abs(np.diag(y))>0))
    def test_branch_stamp_matches_pi_model(self):
        ys=2-3j; b=0.4j; tap=np.exp(1j*np.deg2rad(30))
        m=branch_stamp(1,2,ys,b,tap)
        self.assertAlmostEqual(m[0,0],(ys+b/2)/abs(tap)**2)
        self.assertAlmostEqual(m[0,1],-ys/np.conj(tap))
        self.assertAlmostEqual(m[1,0],-ys/tap)

if __name__=='__main__': unittest.main()
