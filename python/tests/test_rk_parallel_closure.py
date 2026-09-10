"""Parallel closure study: balanced modes and optional native toggle contract."""
from collections import Counter
import itertools
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import benchmark_rk_parallel_closure as study


class ParallelClosureTests(unittest.TestCase):
    def test_balanced_modes(self):
        for reverse in (False, True):
            orders = study.orders(reverse)
            self.assertEqual({tuple(p) for p in orders}, set(itertools.permutations(study.MODES)))
            for i in range(3):
                self.assertEqual(Counter(p[i] for p in orders), {m: 2 for m in study.MODES})
        self.assertEqual(study.orders(True), [p[::-1] for p in study.orders()])

    def test_summary_direction_independent_check(self):
        blocks = []
        for i, order in enumerate(study.orders()):
            runs = {}
            for mode, factor in zip(study.MODES, (4., 2., 1.)):
                runs[mode] = [{m: dict(builder_seconds=.1*factor, kernel_seconds=.9*factor,
                                     algorithm_seconds=factor) for m in permutation}
                              for permutation in itertools.permutations(study.common.ALGORITHMS)]
            blocks.append(dict(order=order, runs=runs))
        summary = study.summarize(blocks)
        for n, d in study.PAIRS:
            for m in study.common.ALGORITHMS:
                for phase in study.common.PHASES:
                    study.check_equal(summary[n+'/'+d][m][phase], study.recompute(blocks,n,d,m,phase))
        self.assertEqual(summary['candidate_on/baseline']['reduction_kernel']['algorithm_seconds']['median_paired_ratio'], .25)

    @unittest.skipUnless(os.getenv('MORSEFRAMES_PARALLEL_CLOSURE_BENCHMARK'), 'optional closure-toggle native binary')
    def test_native_toggle(self):
        binary = Path(os.environ['MORSEFRAMES_PARALLEL_CLOSURE_BENCHMARK']).resolve()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.txt'; path.write_text(study.common.grid_input(6,2,4))
            worker = study.common.Worker(binary,path,Path(tmp)/'reference.dump')
            try:
                for workers in (1,8):
                    for mode in ('candidate_on','candidate_off','candidate_on'):
                        study.common.check_runs(study.run_mode({'candidate':worker},mode,workers,6),6)
            finally: worker.close()


if __name__ == '__main__': unittest.main()
