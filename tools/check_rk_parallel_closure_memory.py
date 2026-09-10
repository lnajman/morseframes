#!/usr/bin/env python3
"""Supplementary fresh-process peak RSS observations using the frozen binaries.

This is not a gradient-only allocation measurement: construction contributes to
the process high-water mark. Do not run alongside performance measurements.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import benchmark_simplicial_gradients as common


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists(): p.error('Fresh output required')
    study=json.loads(args.study.read_text())
    assert study['completed'] and study['schema']=='rk-parallel-closure-v1'
    data=dict(schema='rk-parallel-closure-memory-v1',study=str(args.study.resolve()),
              study_sha256=common.digest(args.study),source_revision=common.command_output('git','rev-parse','HEAD'),
              source_status=common.command_output('git','status','--porcelain'),
              script_source=Path(__file__).read_text(),
              script_sha256=common.digest(Path(__file__)),
              started_utc=datetime.now(timezone.utc).isoformat(),cases=[])
    seen=set()
    for c in study['cases']:
        name=Path(c['input']).name
        if name not in ('grid-d7-n2-seed3.txt','grid-d6-n2-seed4.txt','volume-n32-seed0.txt'): continue
        key=(c['input'],c['workers'])
        if key in seen: continue
        seen.add(key)
        assert common.digest(Path(c['input']))==c['input_sha256']
        identity=c['metadata']['baseline']['identity']
        row=dict(input=c['input'],input_sha256=c['input_sha256'],workers=c['workers'],samples=[])
        for repeat in range(3):
            order=['candidate','baseline'] if repeat%2 else ['baseline','candidate']
            values={}
            for v in order:
                binary=Path(study['builds'][v]['binary'])
                assert common.digest(binary)==study['builds'][v]['binary_sha256']
                result=json.loads(subprocess.check_output([str(binary),c['input'],'--memory','2',str(c['workers'])],text=True))
                assert result['simplices']==identity['simplices'] and result['steps']==identity['algorithms']['reduction_kernel']['steps']
                assert 0<result['complex_peak_bytes']<=result['gradient_peak_bytes']
                values[v]=result
            row['samples'].append(dict(order=order,values=values))
        data['cases'].append(row)
        print('Memory observations:',name,c['workers'],'workers',flush=True)
    assert len(data['cases'])==6
    data['finished_utc']=datetime.now(timezone.utc).isoformat()
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print('Retained 36 fresh-process observations:',common.digest(args.output))


if __name__=='__main__': main()
