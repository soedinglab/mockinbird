import unittest

from mockinbird.ivtree import GenomicInterval


class IntervalTestCase(unittest.TestCase):

    def test_iv_one_based(self):
        GenomicInterval(1, 4, one_based=True)
        GenomicInterval(1, 1, one_based=True)
        with self.assertRaises(ValueError):
            GenomicInterval(2, 1, one_based=True)

    def test_iv_zero_based(self):
        GenomicInterval(1, 4, one_based=False)
        GenomicInterval(1, 2, one_based=False)
        with self.assertRaises(ValueError):
            GenomicInterval(1, 1, one_based=False)

    def test_iv_overlap(self):
        iv1 = GenomicInterval(1, 3)
        iv2 = GenomicInterval(4, 5)
        iv3 = GenomicInterval(5, 6)
        iv4 = GenomicInterval(1, 2)
        iv5 = GenomicInterval(1, 1)
        iv6 = GenomicInterval(3, 3)
        self.assertFalse(iv1.overlaps(iv2))
        self.assertTrue(iv2.overlaps(iv3))
        self.assertFalse(iv4.overlaps(iv2))
        self.assertTrue(iv2.overlaps(iv3))
        self.assertTrue(iv5.overlaps(iv1))
        self.assertTrue(iv6.overlaps(iv1))
