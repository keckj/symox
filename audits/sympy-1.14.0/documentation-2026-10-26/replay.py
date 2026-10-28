import concurrent.futures,collections,gzip,hashlib,json,os,resource,subprocess,sys,time
from pathlib import Path
root=Path(__file__).parent;out=root/'results-576a';out.mkdir(exist_ok=True)
source=Path('/home/keckj/symox-ci/wt/native-products-20261026/python/symox')
expected_native='576a9f76ab3be1cc0d3896d326e4b6615d1a1ffb1d5f0351f23c522bb37d8ad8'
assert hashlib.sha256((source/'_core.abi3.so').read_bytes()).hexdigest()==expected_native
package=root/'package';package.mkdir(exist_ok=True)
if not (package/'symox').exists():(package/'symox').symlink_to(source,target_is_directory=True)
ledger=json.loads(gzip.decompress((root/'catalog.json.gz').read_bytes()))
seen=set();records=[]
for row in sorted(ledger['records'],key=lambda r:r['kind']=='method'):
 if row.get('kind')=='data' or not row.get('examples') or row['canonical'] in seen:continue
 seen.add(row['canonical']);records.append(row)
expected=collections.Counter(r['id'] for r in records);assert len(records)==2916 and all(n==1 for n in expected.values())
for shard in range(32):
 part=dict(ledger,records=records[shard::32]);(out/f'catalog-{shard:02d}.json.gz').write_bytes(gzip.compress(json.dumps(part).encode(),mtime=0))
start=time.time()
def run(job):
 mode,shard=job;name=f'{mode}-{shard:02d}'
 cmd=[sys.executable,str(root/'tools/sympy_parity.py'),'examples','--catalog',str(out/f'catalog-{shard:02d}.json.gz'),'--out',str(out/(name+'.json.gz')),'--snapshot',str(package),'--timeout','60','--fresh']
 if mode=='keywords':cmd.append('--keywords')
 # Inherited by disposable workers, as well as the controller for this shard.
 def limits():resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
 env=dict(os.environ,SYMOX_CI_PY=str(source.parent),PYTHONPATH=str(source.parent),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
 with (out/(name+'.log')).open('w') as log:
  p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
 print(name,p.returncode,flush=True);return name,p.returncode
jobs=[(mode,shard) for mode in ['ordinary','keywords'] for shard in range(32)]
with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:completed=list(pool.map(run,jobs))
summary={'native_sha256':expected_native,'reference_version':'1.14.0','limits':{'seconds_per_object_per_backend':60,'address_space_gib':8},'runners':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'tools').glob('*.py')},'jobs':completed,'coverage':{},'counts':{},'elapsed_seconds':time.time()-start}
for mode in ['ordinary','keywords']:
 results={};actual=collections.Counter()
 for shard in range(32):
  p=out/f'{mode}-{shard:02d}.json.gz'
  if not p.exists():continue
  data=json.loads(gzip.decompress(p.read_bytes()));results.update(data['results'])
 # Count ids rather than report values.
 actual=collections.Counter(key for shard in range(32) if (out/f'{mode}-{shard:02d}.json.gz').exists() for key in json.loads(gzip.decompress((out/f'{mode}-{shard:02d}.json.gz').read_bytes()))['results'])
 summary['coverage'][mode]={'expected':len(records),'finished':sum(actual.values()),'exactly_once':actual==expected,'missing':list((expected-actual).elements()),'unexpected':list((actual-expected).elements())}
 summary['counts'][mode]=dict(collections.Counter(row['status'] for row in results.values()))
assert hashlib.sha256((source/'_core.abi3.so').read_bytes()).hexdigest()==expected_native
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
assert all(code==0 for _,code in completed) and all(row['exactly_once'] for row in summary['coverage'].values())
