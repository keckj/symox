import argparse,pathlib,subprocess,os,sys,json,time,statistics
p=argparse.ArgumentParser();p.add_argument('--after-package',required=True);p.add_argument('--after-source',required=True);a=p.parse_args();root=pathlib.Path(__file__).resolve().parent
backends={'r2':('/home/keckj/symox-ci/wt/native-lrt-limits-20261027/python','daef424363f37e580e8c55dc0ba5ef22631c5ec1c093e282d8774fee7fa15007'),'r3b':(a.after_package,'8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a')}
started=time.monotonic();reports={}
for variant in ('plain','vector'):
 reports[variant]={key:[] for key in backends}
 for round in range(3):
  for name in (('r2','r3b') if round%2==0 else ('r3b','r2')):
   package,sha=backends[name];out=root/f'{variant}-{round}-{name}.json';assert not out.exists()
   subprocess.run([sys.executable,str(root/'worker.py'),'--sha',sha,'--out',str(out),'--variant',variant,'--cpu','24'],env=dict(os.environ,PYTHONPATH=package),check=True,timeout=180)
   reports[variant][name].append(json.loads(out.read_text()));print(variant,round,name,flush=True)
summary={'fresh_processes_per_variant_backend':3,'cpu':24,'reference_runtime':'blocked in every worker','after_source':a.after_source,'elapsed_seconds':time.monotonic()-started,'variants':{}}
for variant,backends in reports.items():
 summary['variants'][variant]={}
 for case in backends['r2'][0]['cases']:
  digests={p['cases'][case]['digest'] for rows in backends.values() for p in rows};assert len(digests)==1
  old=[p['cases'][case]['median_ns_per_pair'] for p in backends['r2']];new=[p['cases'][case]['median_ns_per_pair'] for p in backends['r3b']]
  assert all(p['cases'][case]['native_selected_outside_timers'] for p in backends['r3b'])
  summary['variants'][variant][case]={'r2_ns':statistics.median(old),'r3b_ns':statistics.median(new),'ratio_r2_over_r3b':statistics.median(old)/statistics.median(new),'r2_process_medians':old,'r3b_process_medians':new,'digest':digests.pop(),'r2_native_selection':[p['cases'][case]['native_selected_outside_timers'] for p in backends['r2']]}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
