#!/usr/bin/env python3
"""Pack completed process proofs without duplicating thousands of module maps.

Run after all audit commands finish. Raw per-process evidence remains available
as compressed JSON Lines; the result files retain a digest and process summary.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re

from audit import ROOT, archive_evidence, read, render, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'audits/downstream')
    parser.add_argument('--cache', type=Path, default=Path('/tmp/symox-featured'))
    args = parser.parse_args()
    rows = [(p, read(p)) for p in sorted((args.out/'results').glob('*.json'))]
    if any(row['status'].startswith('running_') for _, row in rows):
        parser.error('Wait for all audit runs to finish before packing their evidence')
    for path, row in rows:
        name = row['project']['id']
        # A separately recorded reference confirmation can add new raw logs.
        if 'reference_confirmation' in row:
            row.setdefault('evidence', {}).update(archive_evidence(args, name,
                (args.cache/'logs'/name).glob('reference-confirm*')))
        for stage in ('reference', 'native', 'reference_confirmation'):
            execution = row.get(stage, {})
            proofs = execution.pop('process_provenance', None)
            if stage == 'reference_confirmation' and execution and not proofs:
                proofs = [read(p) for p in sorted((args.cache/'logs'/name).glob('reference-confirm-provenance.*.json'))]
            if not isinstance(proofs, list) or not proofs:
                continue
            target = args.out/'evidence'/name/(stage+'-processes.jsonl.gz')
            target.parent.mkdir(parents=True, exist_ok=True)
            data = b'\n'.join(json.dumps(p, sort_keys=True, separators=(',', ':')).encode() for p in proofs)+b'\n'
            target.write_bytes(gzip.compress(data, mtime=0))
            relative = str(target.relative_to(args.out))
            execution['process_provenance_artifact'] = relative
            execution['process_provenance_summary'] = dict(processes=len(proofs),
                backends=dict(Counter(p['backend'] for p in proofs)),
                blocked_reference_execution=sum(len(p['blocked_reference_execution']) for p in proofs))
            row.setdefault('evidence', {})[relative] = hashlib.sha256(target.read_bytes()).hexdigest()
        row['evidence'] = {p:digest for p,digest in row.get('evidence', {}).items()
                           if not re.search(r'provenance\.\d+\.json\.gz$', p)}
        write(path, row)
    referenced = {p for _,row in rows for p in row.get('evidence', {})}
    validation = args.out/'validation.json'
    if validation.exists(): referenced.update(read(validation).get('evidence', {}))
    for path in (args.out/'evidence').rglob('*.gz'):
        if str(path.relative_to(args.out)) not in referenced:
            path.unlink()
    render(args)


if __name__ == '__main__':
    main()
