import tempfile
import unittest

from stammp.utils import native_wordcount


class TestClass(unittest.TestCase):

    def test_native_wc(self):
        with tempfile.NamedTemporaryFile('wt') as tmp_file:
            for word in 'this is a test'.split():
                print(word, file=tmp_file)
            tmp_file.flush()
            line_count = native_wordcount(tmp_file.name)
            self.assertEqual(line_count, 4)
