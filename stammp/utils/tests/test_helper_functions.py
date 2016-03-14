import tempfile
import unittest
import os

from stammp.utils import native_wordcount, prepare_output_dir


class TestClass(unittest.TestCase):

    def test_native_wc(self):
        with tempfile.NamedTemporaryFile('wt') as tmp_file:
            for word in 'this is a test'.split():
                print(word, file=tmp_file)
            tmp_file.flush()
            line_count = native_wordcount(tmp_file.name)
            self.assertEqual(line_count, 4)

    def test_prepare_outdir(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            prepare_output_dir(tmp_dir)

        with tempfile.TemporaryDirectory() as tmp_dir:
            final_dir = os.path.join(tmp_dir, 'out_dir')
            prepare_output_dir(final_dir)
            self.assertTrue(os.path.isdir(final_dir))

        with tempfile.TemporaryDirectory() as tmp_dir:
            final_path = os.path.join(tmp_dir, 'foo.tmp')
            with open(final_path, 'w') as tmp_file:
                print('Hello world!', file=tmp_file)
            with self.assertRaisesRegex(ValueError, '.*not .* a directory'):
                prepare_output_dir(final_path)

        with tempfile.TemporaryDirectory() as tmp_dir:
            # make it non-writable
            os.chmod(tmp_dir, 0o555)
            with self.assertRaisesRegex(ValueError, '.*not writable'):
                prepare_output_dir(tmp_dir)

        with tempfile.TemporaryDirectory() as tmp_dir:
            os.chmod(tmp_dir, 0o555)
            final_dir = os.path.join(tmp_dir, 'out_dir')
            with self.assertRaisesRegex(ValueError, '.*cannot be created'):
                prepare_output_dir(final_dir)
