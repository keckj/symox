#!/usr/bin/env python3
"""Attest that the recorded tests match their pinned upstream sources."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import tarfile
import urllib.request

from audit import ROOT, SUPPORT, read, write


def sha(data):
    return hashlib.sha256(data).hexdigest()


def verify(result, project, cache):
    selected = result.get('selected_test_sha256', result.get('selected_tests', {}))
    if not selected:
        return dict(status='no_selected_test_hashes', selected_files=0)
    mismatches = [p for p, digest in selected.items() if not Path(p).is_file() or sha(Path(p).read_bytes()) != digest]
    if project.get('source_git'):
        git = project['source_git']
        if project['id'] == 'torchinductor':
            for filename, digest in selected.items():
                relative = Path(filename).relative_to(cache/'git/torch').as_posix()
                url = 'https://raw.githubusercontent.com/pytorch/pytorch/'+git['commit']+'/'+relative
                if sha(urllib.request.urlopen(url, timeout=60).read()) != digest:
                    mismatches.append(filename)
        else:
            source = Path(result.get('source', cache/'git'/project['id']))
            head = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
            clean = subprocess.run(['git', '-C', str(source), 'diff', '--quiet', git['commit'], '--']).returncode == 0
            if head != git['commit'] or not clean:
                mismatches.append('Git revision or tracked source differs from the pin')
            for filename, digest in selected.items():
                path = Path(filename)
                if source in path.parents:
                    continue
                # Compiled wheels can supply the same tests at a different path.
                module_path = path.as_posix().split('/site-packages/', 1)[-1]
                counterpart = source/module_path
                if not counterpart.is_file() or sha(counterpart.read_bytes()) != digest:
                    mismatches.append(filename)
    elif project.get('source_artifact'):
        artifact = project['source_artifact']
        archive = cache/'archives'/artifact['filename']
        if sha(archive.read_bytes()) != artifact['digests']['sha256']:
            mismatches.append('Source archive digest differs from the pin')
        else:
            source = Path(result['source'])
            with tarfile.open(archive) as package:
                contents = {Path(m.name).parts[1:]: m for m in package.getmembers() if m.isfile()}
                for filename, digest in selected.items():
                    path = Path(filename)
                    relative = (path.relative_to(source) if source in path.parents else
                                Path(path.as_posix().split('/site-packages/', 1)[-1]))
                    member = contents.get(relative.parts)
                    if member is None or sha(package.extractfile(member).read()) != digest:
                        mismatches.append(filename)
    else:
        # Sage ships its tests as part of the explicitly locked Conda runtime.
        return dict(status='installed_runtime_hashes_match' if not mismatches else 'mismatch',
                    selected_files=len(selected), mismatches=mismatches,
                    note='Runtime source is identified by the recorded Conda explicit package list')
    return dict(status='pinned_upstream_tests_match' if not mismatches else 'mismatch',
                selected_files=len(selected), mismatches=mismatches)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('/tmp/symox-featured'))
    parser.add_argument('--out', type=Path, default=ROOT/'audits/downstream')
    args = parser.parse_args()
    projects = {p['id']:p for p in read(SUPPORT/'projects.json')['projects']}
    results = {}
    for path in sorted((args.out/'results').glob('*.json')):
        result = read(path); name = result['project']['id']
        results[name] = verify(result, projects[name], args.cache.resolve())
        print(name+': '+results[name]['status'], flush=True)
    write(args.out/'source-integrity.json', results)
    if any(row['status'] == 'mismatch' for row in results.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
