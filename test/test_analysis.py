import unittest
import siren.analysis
import numpy as np
import os


class TestFindPeak(unittest.TestCase):

    def test_find_peak_no_start(self):
        t = np.array([0.0, 1.0, 2.0,  3.0,   4.0,   5.0,   6.0])
        a = np.array([0.0, 0.0, 10.5, 4.5e2, 3.7e3, 3.9e3, 1.2e2])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak, apeak = siren.analysis.find_peak(tac, 'aorta')
        self.assertEqual(tpeak, 5.0)
        self.assertEqual(apeak, 3.9e3)

    def test_find_peak2_no_start(self):
        t = np.array([0.0, 1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0])
        a = np.array([0.0, 1.0e1, 1.5e2, 4.5e3, 3.7e3, 3.9e3, 1.2e2, 9.2e1])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak, apeak = siren.analysis.find_peak(tac, 'aorta')
        self.assertEqual(tpeak, 3.0)
        self.assertEqual(apeak, 4.5e3)

    def test_find_peak_start(self):
        t = np.array([0.0, 1.0,   2.0,   3.0,   4.0,   5.0,   6.0])
        a = np.array([0.0, 1.0e4, 1.0e1, 4.5e2, 3.7e3, 3.9e3, 1.2e2])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak, apeak = siren.analysis.find_peak(tac, 'aorta', start=1.5)
        self.assertEqual(tpeak, 5.0)
        self.assertEqual(apeak, 3.9e3)


class TestFindCrossUnder(unittest.TestCase):

    def test_find_cross_under_no_start(self):
        t = np.array([0.0, 1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0])
        a = np.array([4.0e3, 1.0e4, 1.5e4, 4.5e3, 3.7e3, 2.9e3, 2.1e3, 9.2e2])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak = siren.analysis.find_first_under(tac,
                                                value=2.5e3,
                                                label='aorta')
        self.assertEqual(tpeak, 6.0)

    def test_find_cross_under_start(self):
        t = np.array([0.0, 1.0,   2.0,   3.0,   4.0,   5.0,   6.0])
        a = np.array([0.0, 1.0e4, 1.0e1, 4.5e3, 3.0e3, 2.0e3, 1.2e2])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak = siren.analysis.find_first_under(tac,
                                                label='aorta',
                                                value=2.9e3,
                                                start=2.1)
        self.assertEqual(tpeak, 5.0)

    def test_find_cross_under_no_half(self):
        t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        a = np.array([0.0, 1.0e1, 1.5e2, 4.5e3, 3.7e3, 3.9e3, 3.1e3, 3.2e3])
        tac = {
            'tacq': t,
            'aorta': a
        }
        tpeak = siren.analysis.find_first_under(tac,
                                                label='aorta',
                                                value=3.0e3,
                                                start=2.5)
        self.assertEqual(tpeak, -1.0)


class TestIntegrate(unittest.TestCase):

    def test_integrate_at_smpls_rect(self):
        t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        k = np.array([0.0, 1.0, 2.0, 2.0, 2.0, 2.0, 1.0])
        tac = {
            'tacq': t,
            'kidney': k

        }
        i = siren.analysis.integrate(tac, 'kidney', 2.0, 5.0)
        self.assertEqual(6.0, i)

    def test_integrate_at_smpls_trap(self):
        t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        k = np.array([0.0, 1.0, 2.0, 3.0, 3.0, 2.0, 1.0])
        tac = {
            'tacq': t,
            'kidney': k

        }
        i = siren.analysis.integrate(tac, 'kidney', 2.0, 5.0)
        self.assertEqual(8.0, i)

    def test_integrate_out_smpls_rect(self):
        t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        k = np.array([0.0, 1.0, 2.0, 2.0, 2.0, 2.0, 1.0])
        tac = {
            'tacq': t,
            'kidney': k

        }
        i = siren.analysis.integrate(tac, 'kidney', 1.1, 5.9)
        self.assertEqual(6.0, i)


class TestAvg(unittest.TestCase):

    def test_avg(self):
        t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        k = np.array([0.0, 100.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        tac = {
            'tacq': t,
            'kidney': k

        }
        i = siren.analysis.avg(tac,
                               label='kidney',
                               start=1.1)
        self.assertEqual(3.0, i)


class TestVolume(unittest.TestCase):

    def test_volume(self):
        v = siren.get_volume(
            os.path.join('test', '8_3V_seg', 'Segmentation.nrrd'))
        self.assertAlmostEqual(v, 29.1784, places=4)
