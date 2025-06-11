import os
import unittest
from datetime import datetime
import siren.core
import numpy as np


class TestGetAcqDateTime(unittest.TestCase):

    def test_acq_datetime_8_3V_1(self):
        dcm_path = os.path.join(
            'test', '8_3V',
            'Patient_test_Study_10_Scan_10_Bed_1_Dyn_1.dcm')
        dt = siren.core.get_acq_datetime(dcm_path)
        self.assertEqual(dt, datetime(2023, 12, 1, 13, 30, 28, 0))

    def test_acq_datetime_8_3V_5(self):
        dcm_path = os.path.join(
            'test', '8_3V',
            'Patient_test_Study_10_Scan_10_Bed_1_Dyn_5.dcm')
        dt = siren.core.get_acq_datetime(dcm_path)
        self.assertEqual(dt, datetime(2023, 12, 1, 13, 30, 40, 800000))


class TestGetTAC(unittest.TestCase):

    def test_series_roi_8_3V_tac(self):
        dcm_path = os.path.join('test', '8_3V')
        roi_path = os.path.join(
            'test', '8_3V_seg', 'Segmentation.nrrd')
        dyn = siren.core.get_tac_from_paths(dcm_path, [roi_path])

        tacq_exp = np.array([0, 3.0, 6.3, 9.5, 12.8, 16.0, 19.3, 22.5, 25.8])
        self.assertFalse(np.any(dyn['tacq'] - tacq_exp))

        r1 = dyn['Segmentation.nrrd']
        r1_exp = np.array([0.0, 188.081845, 301254.45, 2944728.5,
                           2954430.5, 312867.45, 3303.139, 183.26686, 0.0])
        test_arr = np.nan_to_num(abs(r1 - r1_exp) / r1_exp)
        self.assertTrue(np.all(test_arr < 0.00001))
