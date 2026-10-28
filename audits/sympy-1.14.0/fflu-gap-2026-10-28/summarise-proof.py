from pathlib import Path
from fractions import Fraction
from collections import Counter
import hashlib,json
folder=Path(__file__).resolve().parent
cases=[c for name in ['fixtures-neighbours.json','fixtures-rational.json'] for c in json.loads((folder/name).read_text())['cases']]
assert len(cases)==203 and len({c['id'] for c in cases})==203
assert sum(len(c.get('domains',['ZZ','QQ'])) for c in cases)==404
unique={json.dumps([c['shape'],[[str(Fraction(v)) for v in row] for row in c['rows']]],sort_keys=True) for c in cases}
corrected=json.loads((folder/'correction-math.json').read_text());fixed={(r['case'],r['domain']):r for r in corrected['records']}
summary={'fixture_ids':203,'numeric_matrices_distinct_by_shape_and_exact_entries':len(unique),'integer_fixture_ids_with_both_domains':201,'rational_fixture_ids_QQ_only':2,'domain_case_count':404,'group_counts':dict(Counter(c['group'] for c in cases)),'all_small_shape_pairs':sorted({tuple(c['shape']) for c in cases if c['group']=='all-small-shapes-ranks'}),'legacy_calls_per_engine':1616,'legacy':{},'prototype':{k:v for k,v in corrected.items() if k!='records'}}
for lib in ['sympy','symox']:
 records=[r for prefix in ['neighbours','rational'] for r in json.loads((folder/(prefix+'-'+lib+'.json')).read_text())['records']]
 assert len(records)==1616 and len({(r['case'],r['domain'],r['entry']) for r in records})==1616
 assert all(r['input_unchanged'] and r.get('input_rows_cells_unchanged',True) for r in records)
 assert all(r.get('domain_preserved',True) and r.get('private_shapes_domain',True) for r in records)
 passed=[];failures=[]
 for r in records:
  expected=fixed[(r['case'],r['domain'])]
  if r['status']=='returned' and r['certificate']['complete_certificate'] and r['entry']!='private':
   assert r['factors']==expected['factors'],r
   passed.append((r['case'],r['domain'],r['entry']))
  if r['status']=='raised' or not r['certificate']['complete_certificate']:failures.append({'case':r['case'],'domain':r['domain'],'entry':r['entry'],'status':r['status'],'input_rank':expected['certificate']['rank_A'],'error':r.get('error'),'certificate':r.get('certificate')})
 phase2=json.loads((folder/('phase2-'+lib+'.json')).read_text());assert phase2['domain_cases']==404 and all(r['certificate']['complete_certificate'] for r in phase2['records'])
 summary['legacy'][lib]={'returned_certified':sum(r.get('certificate',{}).get('complete_certificate',False) for r in records),'returned_wrong_identity_or_shape':sum(r['status']=='returned' and not r['certificate']['complete_certificate'] for r in records),'raised':sum(r['status']=='raised' for r in records),'unchanged_inputs_and_domains':True,'original_correct_public_factors_unchanged':len(passed),'original_phase2_corrected_inputs_certified':404,'failures':failures}
summary['reference_intentional_correction_required']=True
(folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('legacy','prototype')}))
