"""Unit tests for the healthy-only generator fault experiment utilities."""
import unittest
import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))

import generator_fault_experiment as exp


class ExperimentTests(unittest.TestCase):
    def test_threshold_is_calibration_quantile(self):
        scores = np.array([0., 1., 2., 3., 4.])
        self.assertEqual(exp.calibrate_threshold(scores), float(np.quantile(scores, .99, method="higher")))

    def test_residual_metric_is_same_two_target_metric(self):
        residual = np.array([[3., 4.], [0., 2.]])
        scales = np.array([3., 4.])
        got = exp.residual_metric(residual, scales)
        np.testing.assert_allclose(got, [np.sqrt(.5*(1+1)), np.sqrt(.5*(0+.25))])

    def test_split_windows_are_disjoint(self):
        t = np.array([5., 7.999, 8., 9.999, 10.])
        fit = (t >= exp.FIT[0]) & (t < exp.FIT[1])
        cal = (t >= exp.CAL[0]) & (t < exp.CAL[1])
        self.assertEqual(t[fit].tolist(), [5., 7.999])
        self.assertEqual(t[cal].tolist(), [8., 9.999])
        self.assertFalse(np.any(fit & cal))

    def test_stream_split_preserves_native_streams(self):
        path = exp.ROOT / "healthy_datasets" / "G1_10s_H.csv"
        _, wave, slow = exp.split_data(path, "Pmpu", "Efd")
        self.assertEqual(len(wave), 14401)
        self.assertEqual(len(slow), 601)
        self.assertAlmostEqual(float(np.median(np.diff(wave.t))), 1/1440, places=8)
        self.assertAlmostEqual(float(np.median(np.diff(slow.t))), 1/60, places=8)
        self.assertNotIn("theta", wave.columns)
        self.assertIn("theta", slow.columns)


if __name__ == "__main__":
    unittest.main()
