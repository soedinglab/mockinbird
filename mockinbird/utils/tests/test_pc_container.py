import unittest
from os import path
import os

import numpy.testing as npt

from stammp.utils import ParclipSiteContainer
from stammp.obj.parclipsites import ParclipSites

DATA_DIR = path.join(path.dirname(__file__), 'data')
TABLE_DIR = path.join(DATA_DIR, 'test.table')


class TestClass(unittest.TestCase):

    def _get_container(self):
        pc_table = ParclipSiteContainer.from_file(TABLE_DIR)
        return pc_table

    def test_overlap(self):
        pc = self._get_container()
        ovl_pos = []
        for rec in pc.get_overlaps('chrI', 5, 9, '+'):
            ovl_pos.append(rec.position)
        self.assertEqual(ovl_pos, [5, 9])

        ovl_pos = []
        for rec in pc.get_overlaps('chrI', 5, 9, '-'):
            ovl_pos.append(rec.position)
        self.assertEqual(ovl_pos, [6])

        ovl_pos = []
        for rec in pc.get_overlaps('chrI', 10, 10, '-'):
            ovl_pos.append(rec.position)
        self.assertEqual(ovl_pos, [10])

        ovl_pos = []
        for rec in pc.get_overlaps('chrII', 8, 13, '+'):
            ovl_pos.append(rec.position)
        self.assertEqual(ovl_pos, [12])

    def test_remove_ovl(self):
        pc = self._get_container()
        pc.remove_ovl_sites('chrI', 5, 12, '-')
        ovl_pos = []
        for rec in pc:
            ovl_pos.append(rec.position)
        self.assertEqual(ovl_pos, [5, 9, 12])

    def test_occ_profile(self):
        pc = self._get_container()
        profile = pc.get_occ_profile('chrI', 4, 10, '+')
        npt.assert_allclose(profile, [0, 0.3, 0, 0, 0, 0.5, 0])

        profile = pc.get_occ_profile('chrI', 4, 10, '-')
        npt.assert_allclose(profile, [0, 0, 0.2, 0, 0, 0, 0.1])

    def test_slicing(self):
        pc = self._get_container()
        obs_pos = []
        for rec in pc[:3]:
            obs_pos.append(rec.position)
        self.assertEqual(obs_pos, [5, 6, 9])

    def test_slicing_sort(self):
        pc = self._get_container()
        pc.sort(by='occupancy', ascending=False)
        obs_pos = []
        for rec in pc[:3]:
            obs_pos.append(rec.position)
        self.assertEqual(obs_pos, [9, 5, 6])

        pc.sort(by='score', ascending=False)
        obs_pos = []
        for rec in pc[:3]:
            obs_pos.append(rec.position)
        self.assertEqual(obs_pos, [12, 6, 9])

        pc.sort()
        print(pc._data_list)
        obs_pos = []
        for rec in pc[:3]:
            obs_pos.append(rec.position)
        self.assertEqual(obs_pos, [5, 6, 9])
