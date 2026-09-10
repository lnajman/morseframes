"""Parallel closure study: balanced modes and optional native toggle contract."""
from collections import Counter
import itertools
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import benchmark_rk_parallel_closure as study


class ParallelClosureTests(unittest.TestCase):
    def test_baseline_switch_is_explicit_and_backward_compatible(self):
        for baseline_serial in (False, True):
            for mode in study.MODES:
                stream = io.StringIO()
                worker = SimpleNamespace(process=SimpleNamespace(stdin=stream), read=lambda: [])
                study.run_mode({'baseline':worker,'candidate':worker},mode,8,6,baseline_serial)
                serial = mode == 'candidate_off' or (mode == 'baseline' and baseline_serial)
                self.assertEqual(stream.getvalue(), ('run_closure_serial' if serial else 'run')+' 8 6\n')

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

    @unittest.skipUnless(os.getenv('MORSEFRAMES_SELECTIVE_CLOSURE_BENCHMARK'), 'optional selective-closure full native binary')
    def test_native_selective_profiles(self):
        binary = Path(os.environ['MORSEFRAMES_SELECTIVE_CLOSURE_BENCHMARK']).resolve()
        with tempfile.TemporaryDirectory() as tmp:
            for dimension,seed,eligible in ((6,4,False),(7,3,True)):
                path = Path(tmp)/'input.txt'; path.write_text(study.common.grid_input(dimension,2,seed))
                worker = study.common.Worker(binary,path,Path(tmp)/'reference.dump')
                try:
                    for command in ('rk_detailed','rk_levels_detailed'):
                        worker.process.stdin.write(command+' 8\n'); worker.process.stdin.flush()
                        row = worker.read()
                        self.assertEqual(row['closure_parallel_tasks']>0,eligible)
                        if 'level_trace' in row:
                            self.assertEqual(sum(l['closure_parallel_tasks'] for l in row['level_trace']['levels']),row['closure_parallel_tasks'])
                            for level in row['level_trace']['levels']:
                                if level['simplices']<32768:
                                    self.assertEqual(level['closure_parallel_tasks'],0)
                finally: worker.close()


if __name__ == '__main__': unittest.main()
