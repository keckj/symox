from pathlib import Path
import json,os,statistics,subprocess,sys,time
root=Path(__file__).resolve().parent;started=time.monotonic()
backends={'cb5':('/home/keckj/symox-ci/wt/native-derivative-r1b-20261027/python','cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e'),'r2':('/home/keckj/symox-ci/wt/native-lrt-limits-20261027/python','daef424363f37e580e8c55dc0ba5ef22631c5ec1c093e282d8774fee7fa15007')}
reports={name:[] for name in backends}
for i in range(3):
 for name in (['cb5','r2'] if i%2==0 else ['r2','cb5']):
  package,sha=backends[name];out=root/f'{i}-{name}.json';assert not out.exists()
  subprocess.run([sys.executable,str(root/'worker.py'),'--sha',sha,'--out',str(out)],env=dict(os.environ,PYTHONPATH=package),check=True,timeout=180)
  reports[name].append(json.loads(out.read_text()));print(i,name,flush=True)
summary={'fresh_processes_per_backend':3,'cpu':24,'reference_runtime':'blocked in all six workers','iterations':3,'rounds_per_process':5,'families':5,'distinct_matrices':15,'full_coefficients_checked_before_after':180,'independent_Cayley_Hamilton_certificates':15,'elapsed_seconds':time.monotonic()-started,'cases':{}}
for name in reports['cb5'][0]['cases']:
 digests={r['cases'][name]['digest'] for rows in reports.values() for r in rows};assert len(digests)==1
 before=[r['cases'][name]['median_ns'] for r in reports['cb5']];after=[r['cases'][name]['median_ns'] for r in reports['r2']]
 summary['cases'][name]={'cb5_ns':statistics.median(before),'r2_ns':statistics.median(after),'ratio':statistics.median(before)/statistics.median(after),'cb5_process_medians':before,'r2_process_medians':after,'digest':digests.pop()}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
