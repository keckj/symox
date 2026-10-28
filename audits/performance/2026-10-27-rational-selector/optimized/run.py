import argparse,pathlib,subprocess,os,sys,json,time,statistics
p=argparse.ArgumentParser();p.add_argument('--optimized-adapter-sha',required=True);a=p.parse_args();root=pathlib.Path(__file__).resolve().parent
backends={
 'r2':('/home/keckj/symox-ci/wt/native-lrt-limits-20261027/python','daef424363f37e580e8c55dc0ba5ef22631c5ec1c093e282d8774fee7fa15007','add84089dc8297efa1d1dfbf25125629db2df6b498183a809a66fa1bd2bb4c77'),
 'r3b':('/home/keckj/symox-ci/wt/native-inexact-r1-20261027/python','8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a','525ea4127a4f216750625575296c32ad28672cd91aa51f6d3157b69eaf0e0b44'),
 'optimized':('/home/keckj/symox-ci/wt/selector-proof-release-r3b/python','8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a',a.optimized_adapter_sha)}
started=time.monotonic();reports={}
for variant in ('plain','vector'):
 reports[variant]={key:[] for key in backends}
 for round,order in enumerate((('r2','r3b','optimized'),('optimized','r2','r3b'),('r3b','optimized','r2'))):
  for name in order:
   package,sha,adapter=backends[name];out=root/f'{variant}-{round}-{name}.json';assert not out.exists()
   subprocess.run([sys.executable,str(root/'worker.py'),'--sha',sha,'--adapter-sha',adapter,'--out',str(out),'--variant',variant,'--cpu','24'],env=dict(os.environ,PYTHONPATH=package),check=True,timeout=180)
   reports[variant][name].append(json.loads(out.read_text()));print(variant,round,name,flush=True)
summary={'fresh_processes_per_variant_backend':3,'cpu':24,'reference_runtime':'blocked in every worker','elapsed_seconds':time.monotonic()-started,'backend_identities':backends,'variants':{}}
for variant,backends in reports.items():
 summary['variants'][variant]={}
 for case in backends['r2'][0]['cases']:
  digests={p['cases'][case]['digest'] for rows in backends.values() for p in rows};assert len(digests)==1
  values={name:[p['cases'][case]['median_ns_per_pair'] for p in rows] for name,rows in backends.items()}
  medians={name:statistics.median(v) for name,v in values.items()}
  assert all(p['cases'][case]['native_selected_outside_timers'] for p in backends['optimized'])
  summary['variants'][variant][case]={'median_ns':medians,'ratio_r3b_over_optimized':medians['r3b']/medians['optimized'],'ratio_r2_over_optimized':medians['r2']/medians['optimized'],'process_medians_ns':values,'digest':digests.pop()}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
