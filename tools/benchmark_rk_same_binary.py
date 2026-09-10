#!/usr/bin/env python3
"""Post-hoc A/A control: one frozen ordinary binary, independent processes.

Three fresh process pairs at each worker count. Keep their results separate;
within-process block intervals do not measure between-process variability.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import benchmark_simplicial_gradients as common
import benchmark_rk_binary_isolation as isolation
from validate_simplicial_gradient_ab import equivalent


def orders(reverse):
    return [['b', 'a'] if (i + reverse) % 2 else ['a', 'b'] for i in range(4)]


def summarize(blocks):
    return {m: {p: common.summarize([{'baseline': [r[m][p] for r in b['runs']['a']],
                                     'candidate': [r[m][p] for r in b['runs']['b']]} for b in blocks])
                for p in common.PHASES} for m in common.ALGORITHMS}


def audit(path):
    data = json.loads(path.read_text())
    assert data['completed'] and data['schema'] == 'rk-same-binary-v1' and not data['source_status']
    provenance = isolation.audit(Path(data['binary_study']))
    build = provenance['builds']['candidate_ordinary']
    assert data['binary'] == build['binary'] and data['binary_sha256'] == build['binary_sha256']
    assert common.digest(Path(data['binary_study'])) == data['binary_study_sha256']
    assert common.digest(Path(data['binary'])) == data['binary_sha256']
    assert common.digest(Path(data['input'])) == data['input_sha256']
    for source, sha in data['source_sha256'].items():
        import hashlib
        raw = subprocess.check_output(['git', 'show', data['source_revision'] + ':' + source], cwd=common.ROOT)
        assert hashlib.sha256(raw).hexdigest() == sha
    assert len(data['pairs']) == 6
    assert [(p['workers'], p['pair']) for p in data['pairs']] == [(w, p) for w in (1, 8) for p in range(3)]
    identities = []
    for pair in data['pairs']:
        expected = orders(pair['reverse'])
        assert pair['launch_order'] == expected[0] and len(pair['blocks']) == 4
        for block, order in zip(pair['blocks'], expected):
            assert block['order'] == order and set(block['runs']) == {'a', 'b'}
            for rows in block['runs'].values():
                common.check_runs(rows, 6)
                assert len({tuple(r) for r in rows}) == 6
        assert pair['metadata']['a']['identity'] == pair['metadata']['b']['identity']
        identities.append(pair['metadata']['a']['identity'])
        assert all(common.digest(Path(d)) == pair['reference_sha256'] for d in pair['dumps'].values())
        equivalent(pair['summary'], summarize(pair['blocks']))
    assert all(i == identities[0] for i in identities)
    print('VALIDATED', path.name, 'six independent A/A process pairs', common.digest(path))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary-study', type=Path)
    p.add_argument('--input', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--audit', type=Path)
    a = p.parse_args()
    if a.audit:
        audit(a.audit); return
    if not a.binary_study or not a.input or not a.output or a.output.exists():
        p.error('A frozen binary study, input and fresh output are required')
    status = common.command_output('git', 'status', '--porcelain')
    if status:
        p.error('Freeze sources first')
    provenance = isolation.audit(a.binary_study)
    build = provenance['builds']['candidate_ordinary']
    binary = Path(build['binary'])
    artifacts = a.output.resolve().with_suffix('.artifacts'); artifacts.mkdir()
    sources = tuple(dict.fromkeys(('tools/benchmark_rk_same_binary.py',) + isolation.SOURCES))
    data = dict(schema='rk-same-binary-v1', completed=False, source_status=status,
                source_revision=common.command_output('git', 'rev-parse', 'HEAD'),
                source_sha256={s: common.digest(common.ROOT / s) for s in sources},
                binary_study=str(a.binary_study.resolve()), binary_study_sha256=common.digest(a.binary_study),
                binary=str(binary), binary_sha256=build['binary_sha256'],
                input=str(a.input.resolve()), input_sha256=common.digest(a.input),
                started_utc=datetime.now(timezone.utc).isoformat(), pairs=[])
    with a.output.open('x') as output:
        try:
            for wi, workers in enumerate((1, 8)):
                for pair in range(3):
                    reverse = (wi + pair) % 2
                    launch = orders(reverse)[0]
                    with ExitStack() as stack:
                        processes, dumps = {}, {}
                        for label in launch:
                            dump = artifacts / f'w{workers}-pair{pair}-{label}.dump'
                            dumps[label] = str(dump)
                            processes[label] = common.Worker(binary, a.input.resolve(), dump)
                            stack.callback(processes[label].close)
                        assert common.digest(Path(dumps['a'])) == common.digest(Path(dumps['b']))
                        for label in launch:
                            common.check_runs(processes[label].run(workers, 2), 2)
                        blocks = []
                        for order in orders(reverse):
                            blocks.append(dict(order=order, runs={v: processes[v].run(workers, 6) for v in order}))
                        data['pairs'].append(dict(workers=workers, pair=pair, reverse=reverse,
                            launch_order=launch, dumps=dumps, reference_sha256=common.digest(Path(dumps['a'])),
                            metadata={v: w.metadata for v, w in processes.items()}, blocks=blocks))
                        print('Measured identical binary:', workers, 'workers, process pair', pair, flush=True)
            data['finished_utc'] = datetime.now(timezone.utc).isoformat()
            for pair in data['pairs']:
                pair['summary'] = summarize(pair['blocks'])
            data['completed'] = True
        finally:
            json.dump(data, output, indent=2); output.write('\n')
    audit(a.output)


if __name__ == '__main__':
    main()
