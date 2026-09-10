"""Binary-isolation design, accounting and native ordinary-only contract."""
from collections import Counter
import copy
import itertools
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import benchmark_rk_binary_isolation as isolation
import benchmark_rk_same_binary as same_binary
import render_rk_binary_isolation as render


class BinaryIsolationTests(unittest.TestCase):
    def test_independent_summary_and_corruption(self):
        blocks = self.blocks()
        data = dict(cases=[dict(blocks=blocks, summary=isolation.summaries(blocks))])
        self.assertEqual(render.independently_check(data), 36)
        data['cases'][0]['summary']['candidate_ordinary/baseline_ordinary']['f_max']['algorithm_seconds']['median_paired_ratio'] += .01
        with self.assertRaises(AssertionError):
            render.independently_check(data)

    def test_same_binary_schedule_and_ratios(self):
        for reverse in (0, 1):
            orders = same_binary.orders(reverse)
            self.assertEqual(Counter(tuple(r) for r in orders), {('a', 'b'): 2, ('b', 'a'): 2})
        source = self.blocks()
        blocks = [dict(order=['a', 'b'], runs={'a': b['runs']['baseline_ordinary'], 'b': b['runs']['candidate_ordinary']}) for b in source[:4]]
        for phases in same_binary.summarize(blocks).values():
            self.assertEqual(phases['algorithm_seconds']['median_paired_ratio'], 2)

    def test_trace_symbol_detection(self):
        self.assertEqual(isolation.trace_symbols('xbuild_reduction_kernel_implILb0E'), 0)
        self.assertEqual(isolation.trace_symbols('xbuild_reduction_kernel_implILb1E\nxrun_levels'), 2)

    def test_crossover(self):
        for count in (8, 16):
            for reverse in (False, True):
                orders = isolation.schedule(count, reverse)
                for pos in range(4):
                    self.assertEqual(Counter(r[pos] for r in orders), {v: count // 4 for v in isolation.VARIANTS})
                for a, b in itertools.combinations(isolation.VARIANTS, 2):
                    self.assertEqual(sum(r.index(a) < r.index(b) for r in orders), count // 2)
            self.assertEqual(isolation.schedule(count, True), [r[::-1] for r in isolation.schedule(count)])
        for count in (0, 4, 10, -8):
            with self.assertRaises(ValueError):
                isolation.schedule(count)

    @staticmethod
    def blocks():
        blocks = []
        factors = dict(baseline_ordinary=1, candidate_ordinary=2, baseline_full=4, candidate_full=8)
        for i, order in enumerate(isolation.schedule(8)):
            runs = {}
            for v in order:
                runs[v] = []
                for methods in itertools.permutations(isolation.common.ALGORITHMS):
                    t = (i + 1) * factors[v]
                    runs[v].append({m: dict(builder_seconds=.1*t, kernel_seconds=.9*t, algorithm_seconds=t) for m in methods})
            blocks.append(dict(order=order, runs=runs))
        return blocks

    def test_ratio_direction(self):
        rows = isolation.summaries(self.blocks())
        for key, methods in rows.items():
            expected = 2 if key.startswith('candidate_') and '/baseline_' in key else 4
            for phases in methods.values():
                self.assertEqual(phases['algorithm_seconds']['median_paired_ratio'], expected)
                self.assertEqual(phases['algorithm_seconds']['paired_ratio_bootstrap_95_interval'], [expected, expected])

    def test_coverage_and_timing_rejections(self):
        with tempfile.TemporaryDirectory() as temp:
            dump = Path(temp) / 'reference.dump'; dump.write_text('reference')
            identity = dict(simplices=3, lower_star_sizes=[1, 2], euler=1,
                            algorithms={m: dict(steps=2, critical_counts=[1, 0]) for m in isolation.common.ALGORITHMS})
            blocks = self.blocks()
            case = dict(blocks=blocks, launch_order=list(isolation.VARIANTS),
                        metadata={v: dict(identity=identity) for v in isolation.VARIANTS},
                        reference_dumps={v: str(dump) for v in isolation.VARIANTS},
                        reference_sha256=isolation.common.digest(dump), summary=isolation.summaries(blocks))
            settings = dict(blocks=8, repeats=6, reverse=False)
            isolation.validate_case(case, settings)
            mutations = [lambda c: c['blocks'].pop(),
                         lambda c: c['blocks'][0]['order'].reverse(),
                         lambda c: c['blocks'][0]['runs'].pop('candidate_full'),
                         lambda c: c['blocks'][0]['runs']['candidate_ordinary'][0]['f_max'].update(algorithm_seconds=0),
                         lambda c: c.update(reference_sha256='bad')]
            for mutate in mutations:
                bad = copy.deepcopy(case); mutate(bad)
                with self.assertRaises((ValueError, AssertionError)):
                    isolation.validate_case(bad, settings)

    @unittest.skipUnless(os.getenv('MORSEFRAMES_ORDINARY_BENCHMARK') and os.getenv('MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK'), 'optional paired native drivers')
    def test_ordinary_native_contract(self):
        binaries = [Path(os.environ[k]).resolve() for k in ('MORSEFRAMES_ORDINARY_BENCHMARK', 'MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK')]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); path = root / 'input.txt'
            path.write_text(isolation.common.grid_input(4, 2, 1))
            identities, dumps = [], []
            for i, binary in enumerate(binaries):
                dump = root / f'{i}.dump'; dumps.append(dump)
                worker = isolation.common.Worker(binary, path, dump)
                try:
                    identities.append(worker.metadata['identity'])
                    for count in (1, 8):
                        isolation.common.check_runs(worker.run(count, 6), 6)
                finally:
                    worker.close()
            self.assertEqual(*identities)
            self.assertEqual(dumps[0].read_bytes(), dumps[1].read_bytes())
            for command in ('profile 1', 'rk_detailed 1', 'rk_levels_coarse 1'):
                result = subprocess.run([str(binaries[0]), str(path), str(root / 'rejected.dump')],
                                        input=command + '\n', capture_output=True, text=True)
                self.assertEqual(result.returncode, 1)
                self.assertIn('Unknown worker command', result.stderr)
            self.assertEqual(isolation.trace_symbols(isolation.read_symbols(binaries[0])), 0)
            self.assertGreater(isolation.trace_symbols(isolation.read_symbols(binaries[1])), 0)


if __name__ == '__main__':
    unittest.main()
