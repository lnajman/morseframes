"""Level-trace accounting and crossover design; no performance thresholds."""
from collections import Counter
import copy
import itertools
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import profile_reduction_kernel_levels as levels
from benchmark_reduction_kernel_ab import Worker
from benchmark_simplicial_gradients import grid_input


class LevelProfileTests(unittest.TestCase):
    def test_parallel_closure_elapsed_and_cumulative(self):
        from profile_reduction_kernel_simplicial import validate_parallel_closure, PARALLEL_CLOSURE_FIELDS
        row = {k: 0 for k in PARALLEL_CLOSURE_FIELDS}
        row.update({'closure_' + k + '_seconds': 0 for k in ('initial','traversal','sort','materialize')})
        row.update(closure_seconds=1., closure_parallel_seconds=.8,
                   closure_parallel_merge_seconds=.1, closure_parallel_traversal_seconds=2.,
                   closure_parallel_batches=1, closure_parallel_tasks=4)
        validate_parallel_closure(row, 4) # Cumulative worker time may exceed elapsed time.
        for changes in (dict(closure_parallel_seconds=2), dict(closure_parallel_tasks=5),
                        dict(closure_parallel_batches=0), dict(closure_parallel_merge_seconds=.9)):
            with self.assertRaises(ValueError):
                validate_parallel_closure(dict(row, **changes), 4)

    def test_schedule(self):
        for reverse in (False, True):
            orders = levels.mode_orders(reverse)
            for pos in range(5):
                self.assertEqual(Counter(row[pos] for row in orders), {m: 2 for m in levels.MODES})
            for a, b in itertools.combinations(levels.MODES, 2):
                self.assertEqual(sum(row.index(a) < row.index(b) for row in orders), 5)
        self.assertEqual(levels.mode_orders(True), levels.mode_orders()[::-1])

    def test_owner_rank_mapping(self):
        self.assertEqual(levels.level_owners([4, 1, 8, 2]), [1, 3, 0, 2])
        for values in ([], [1, 1], [float('nan')], [float('inf')]):
            with self.assertRaises(ValueError):
                levels.level_owners(values)

    @staticmethod
    def fixture(detailed=True):
        identity = {'simplices': 3, 'lower_star_sizes': [2, 1],
                    'algorithms': {'reduction_kernel': {'steps': 2}}}
        rows = []
        for i in range(2):
            row = dict(level=i, task=0, simplices=i + 1, events=1, completed=True,
                       start_seconds=.1 + .2 * i, duration_seconds=.1)
            row.update({k + '_seconds': 0 for k in levels.TIMES})
            row.update({k: 0 for k in levels.COUNTS})
            if detailed:
                row.update(reductions=i, perforations=1-i, kernel_rounds=1)
            rows.append(row)
        return identity, dict(builder_seconds=.1, kernel_seconds=.9, algorithm_seconds=1.,
                              setup_seconds=.1, level_wall_seconds=.6, replay_seconds=.1,
                              level_trace=dict(completed=True, detailed=detailed, executor_workers=2,
                                               level_tasks=2, level_wall_seconds=.5, levels=rows))

    def test_valid_and_corrupt_traces(self):
        for detailed in (False, True):
            identity, row = self.fixture(detailed)
            mode = 'rk_levels_detailed' if detailed else 'rk_levels_coarse'
            levels.validate(row, mode, identity, [1, 0], 2)
        identity, row = self.fixture()
        mutations = [lambda r: r['level_trace'].update(completed=False),
                     lambda r: r['level_trace']['levels'].pop(),
                     lambda r: r['level_trace']['levels'][0].update(simplices=2),
                     lambda r: r['level_trace']['levels'][1].update(start_seconds=.15),
                     lambda r: r['level_trace']['levels'][0].update(task=2),
                     lambda r: r['level_trace']['levels'][0].update(reductions=1),
                     lambda r: r['level_trace']['levels'][0].update(closure_seconds=.2),
                     lambda r: r['level_trace']['levels'][0].update(duration_seconds=float('nan')),
                     lambda r: r.update(rounds=100),
                     lambda r: r.update(algorithm_seconds=2)]
        for mutate in mutations:
            bad = copy.deepcopy(row); mutate(bad)
            with self.assertRaises(ValueError):
                levels.validate(bad, 'rk_levels_detailed', identity, [1, 0], 2)

    @unittest.skipUnless(os.getenv('MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK'), 'optional native worker')
    def test_native_roundtrip(self):
        binary = Path(os.environ['MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK']).resolve()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.txt'
            path.write_text(grid_input(4, 2, 3))
            tokens = path.read_text().split()
            values = list(map(float, tokens[4:4 + int(tokens[2])]))
            worker = Worker(binary, path, Path(directory) / 'reference.dump')
            try:
                for count in (1, 4, 8):
                    for mode in levels.MODES:
                        worker.process.stdin.write(f'{mode} {count}\n'); worker.process.stdin.flush()
                        levels.validate(worker.read(), mode, worker.metadata['identity'], values, count)
            finally:
                worker.close()


if __name__ == '__main__':
    unittest.main()
