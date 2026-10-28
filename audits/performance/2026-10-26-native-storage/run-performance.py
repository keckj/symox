from pathlib import Path
import hashlib,json,subprocess,time
cache=Path('/home/keckj/.cache/symox-validation/native-storage-20261026')
after='/home/keckj/symox-ci/wt/native-storage-20261026/python'
checks="""import os,subprocess,json
rows=[]
for line in subprocess.check_output(['ps','-eo','pid,comm,args'],text=True).splitlines()[1:]:
 p=line.split(None,2)
 if len(p)<3 or int(p[0])==os.getpid():continue
 if p[1] in ('rustc','cargo','cargo-clippy','clippy-driver') or (p[1].startswith('python') and 'symox-ci' in p[2] and any(k in p[2] for k in ('pytest','run_suite.py','rust-checks.py','benchmark-'))):rows.append(line)
print(json.dumps(rows))
"""
# Pass code through stdin; the process command line cannot match its own patterns.
while True:
 process=subprocess.run(['ssh','-o','BatchMode=yes','demeter','python3 -'],input=checks,capture_output=True,text=True,check=True)
 active=json.loads(process.stdout)
 if not active:break
 print('Waiting for quiet host:',len(active),flush=True);time.sleep(15)
identities={'linprog':{'before':'/home/keckj/symox-ci/wt/native-corrections-20261026-r1/python','after':after},'normal':{'before':'/home/keckj/symox-ci/wt/native-algebra-20261026/python','after':after}}
commands=[]
for trial in range(3):
 for tag in ['before','after']:
  for kind in ['linprog','normal']:
   package=identities[kind][tag]
   script='/home/keckj/symox-ci/simplex-algebra-benchmark.py' if kind=='linprog' else '/home/keckj/symox-ci/benchmark-normal-raw-common.py'
   command=f'SYMOX_CI_PY={package} PYTHONPATH={package} taskset -c 24 /home/keckj/symox-ci/venv/bin/python {script}'
   if kind=='normal':command+=f' /home/keckj/symox-ci/normal-common-{tag}-{trial}.json'
   output=cache/f'{kind}-{tag}-{trial}.json';assert not output.exists()
   with output.open('w') as stdout,(cache/f'{kind}-{tag}-{trial}.stderr').open('w') as stderr:
    subprocess.run(['ssh','-o','BatchMode=yes','demeter',command],stdout=stdout,stderr=stderr,check=True)
   record=json.loads(output.read_text());expected='570f1832' if tag=='after' else ('d9b14020' if kind=='linprog' else '9405dccb')
   assert record['native_sha256'].startswith(expected)
   commands.append(command);print(kind,tag,trial,flush=True)
(cache/'performance-commands.json').write_text(json.dumps(commands,indent=2))
