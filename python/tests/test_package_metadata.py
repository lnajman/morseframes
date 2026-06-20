import morseframes as mf
import unittest


class PackageMetadataTest(unittest.TestCase):
    def test_package_exposes_version_and_backend_status(self):
        self.assertEqual(mf.__version__, "0.1.0a1")
        self.assertIsInstance(mf.cpp_backend_available(), bool)

        complex_ = mf.FilteredComplex.from_simplices([([0], 0.0)])
        self.assertIsInstance(mf.cpp_backend_active(complex_), bool)


if __name__ == "__main__":
    unittest.main()
