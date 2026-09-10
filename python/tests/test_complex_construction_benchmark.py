import copy
from contextlib import redirect_stderr
from io import StringIO
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import benchmark_complex_construction as bench


class ComplexConstructionBenchmarkTest(unittest.TestCase):
    def samples(self):
        def sample(scale):
            return dict(construction_seconds=scale * 10,
                        enumeration_and_insertion_seconds=scale * 3,
                        finalize_seconds=scale * 7, fmax_seconds=scale * 2,
                        rk_seconds=scale)
        return [dict(order=["baseline", "candidate"] if i % 2 else ["candidate", "baseline"],
                     baseline=[sample(i + 1)], candidate=[sample((i + 1) / 2)])
                for i in range(4)]

    def test_paired_summary_and_partition(self):
        samples = self.samples()
        original = copy.deepcopy(samples)
        summary = bench.summarize_samples(samples)
        self.assertEqual(summary["construction_seconds"]["median_paired_ratio"], 0.5)
        self.assertEqual(summary["construction_seconds"]["median_seconds"]["baseline"], 25)
        self.assertEqual(samples, original)

    def test_reject_incomplete_nonfinite_and_overlapping_samples(self):
        for mutate in (
            lambda s: s.pop(),
            lambda s: s[0]["baseline"].clear(),
            lambda s: s[0]["baseline"][0].update(construction_seconds=11),
            lambda s: s[0]["baseline"][0].update(rk_seconds=math.nan),
            lambda s: s[0]["baseline"][0].update(rk_seconds=0),
            lambda s: s[0]["baseline"][0].update(diagnostic_seconds=100),
            lambda s: s[0].update(order=["baseline", "candidate"]),
        ):
            samples = self.samples()
            mutate(samples)
            with self.assertRaises(ValueError):
                bench.summarize_samples(samples)

    def test_argument_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.txt"
            source.touch()
            args = ["--baseline", "HEAD", "--inputs", str(source),
                    "--output", str(Path(directory) / "new.json")]
            self.assertEqual(bench.parse_args(args).workers, [1, 8])
            for extra in (["--blocks", "3"], ["--workers", "1", "1"],
                          ["--diagnostics", "0"], ["--output", str(source)]):
                with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                    bench.parse_args(args + extra)

    @unittest.skipUnless(os.environ.get("MORSEFRAMES_CONSTRUCTION_BENCHMARK"),
                         "optional native construction driver")
    def test_native_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.txt"
            path.write_text("morseframes-ttk-v1 2 3 1\n0 1 2\n0 1 2\n")
            binary = os.environ["MORSEFRAMES_CONSTRUCTION_BENCHMARK"]
            result = subprocess.run([binary, str(path), str(Path(directory) / "dump")],
                                    input="run 2 2\nprofile\nquit\n", text=True,
                                    capture_output=True, check=True)
            metadata, runs, diagnostic = map(json.loads, result.stdout.splitlines())
            self.assertEqual(metadata["simplices"], 7)
            self.assertEqual(metadata["rk_critical_counts"], [1, 0, 0])
            self.assertEqual(metadata["fmax_critical_counts"], [1, 0, 0])
            self.assertEqual(len(runs), 2)
            self.assertEqual(diagnostic["complex_fingerprint"], metadata["complex_fingerprint"])
            if "generated_faces" in diagnostic: # Archived headers use the legacy adapter.
                self.assertEqual(diagnostic["generated_faces"], 7)
                self.assertEqual(diagnostic["unique_faces_submitted"], 7)
                self.assertGreaterEqual(diagnostic["bulk_sort_and_dedup_seconds"], 0)
            for run in runs:
                self.assertAlmostEqual(run["construction_seconds"],
                                       run["enumeration_and_insertion_seconds"] + run["finalize_seconds"])
            memory = json.loads(subprocess.check_output([binary, str(path), "--memory"], text=True))
            self.assertGreaterEqual(memory["construction_peak_bytes"], memory["input_peak_bytes"])
            self.assertEqual(memory["simplices"], 7)


if __name__ == "__main__":
    unittest.main()
