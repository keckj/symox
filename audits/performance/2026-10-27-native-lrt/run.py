import pathlib,subprocess,os,sys,json,time,statistics
root=pathlib.Path(__file__).resolve().parent
backends={'cb5':('/home/keckj/symox-ci/wt/native-derivative-r1b-20261027/python','cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e'),'r2':('/home/keckj/symox-ci/wt/native-lrt-limits-20261027/python','daef424363f37e580e8c55dc0ba5ef22631c5ec1c093e282d8774fee7fa15007')}
started=time.monotonic();reports={k:[] for k in backends}
for round in range(3):
 for name in (('cb5','r2') if round%2==0 else ('r2','cb5')):
  package,sha=backends[name];out=root/f'{round}-{name}.json';assert not out.exists()
  subprocess.run([sys.executable,str(root/'worker.py'),'--sha',sha,'--out',str(out),'--cpu','24'],env=dict(os.environ,PYTHONPATH=package),check=True,timeout=120)
  reports[name].append(json.loads(out.read_text()));print(round,name,flush=True)
summary={'fresh_processes_per_backend':3,'cpu':24,'reference_runtime':'blocked in every worker','iterations':50,'rounds_per_process':5,'certificates':'independent SymPy1.14 exact coefficient fixtures; before/after timers; domains and signs preserved','elapsed_seconds':time.monotonic()-started,'cases':{}}
for case in reports['cb5'][0]['cases']:
 digests={p['cases'][case]['digest'] for rows in reports.values() for p in rows};assert len(digests)==1
 old=[p['cases'][case]['median_ns_per_pair'] for p in reports['cb5']];new=[p['cases'][case]['median_ns_per_pair'] for p in reports['r2']]
 summary['cases'][case]={'cb5_ns':statistics.median(old),'r2_ns':statistics.median(new),'observed_ratio':statistics.median(old)/statistics.median(new),'cb5_process_medians':old,'r2_process_medians':new,'digest':digests.pop()}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
