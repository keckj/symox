"""Inject certified packed LU into unchanged public Phase2, no engine edits."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,importlib,json
from derive_correction import corrected
from certificates import certify,encoded
parser=argparse.ArgumentParser();parser.add_argument('--library',choices=['sympy','symox'],required=True);args=parser.parse_args()
s=importlib.import_module(args.library);DDM=importlib.import_module(args.library+'.polys.matrices.ddm').DDM
if args.library=='sympy':assert s.__version__=='1.14.0'
class Phase2Input(DDM):
 def _fflu(self):return self._prepared
folder=Path(__file__).resolve().parent;cases=[c for name in ['fixtures-neighbours.json','fixtures-rational.json'] for c in json.loads((folder/name).read_text())['cases']];records=[]
for case in cases:
 m,n=case['shape'];a=[[Fraction(v) for v in row] for row in case['rows']]
 for domain in case.get('domains',['ZZ','QQ']):
  packed,perm,expected,certificate,trace,divisions=corrected(a,m,n,domain=='ZZ');k=getattr(s,domain)
  def cells(rows):return [[k(Fraction(v).numerator,Fraction(v).denominator) if domain=='QQ' else k(int(v)) for v in row] for row in rows]
  source=Phase2Input(cells(a),(m,n),k);source._prepared=DDM(cells(packed),(m,n),k),perm[:]
  output=source.fflu();actual=[[[Fraction(str(v)) for v in row] for row in factor] for factor in output]
  assert actual==list(expected)
  proof=certify(a,m,n,*actual);assert proof['complete_certificate'] and all(factor.domain==k for factor in output)
  records.append({'case':case['id'],'domain':domain,'shape':[m,n],'all_factor_entries_equal_math_prototype':True,'certificate':proof,'domains_preserved':True})
root=Path(s.__file__).resolve().parent;source=root/'polys/matrices/ddm.py'
result={'status':'unchanged original Phase2 succeeds for every corrected packed Phase1 input','library':args.library,'package':str(root),'version':getattr(s,'__version__',None),'DDM_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'native_sha256':hashlib.sha256(Path(s._core.__file__).read_bytes()).hexdigest() if args.library=='symox' else None,'domain_cases':len(records),'records':records}
(folder/('phase2-'+args.library+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}))
