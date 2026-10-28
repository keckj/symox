"""Review every archived observation; normalize only literal set element order.

No execution/evaluation of displayed expressions. Numeric and string lexemes
are retained verbatim, including precision/display digits. All other AST field
ordering and every exception/status field remain unchanged.
"""
from pathlib import Path
import ast,gzip,hashlib,json,sys
root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[3];source=root/'audits/sympy-1.14.0/documentation-2026-10-26';out=Path(__file__).parent

def canonical_output(text):
    try:tree=ast.parse(text)
    except (SyntaxError,ValueError,RecursionError):return None
    if not any(isinstance(node,ast.Set) for node in ast.walk(tree)):return None
    def visit(node):
        if isinstance(node,ast.Constant):
            # Do not silently equate differently printed Float precisions.
            return ['Constant',ast.get_source_segment(text,node)]
        if isinstance(node,ast.AST):
            values=[]
            for field,value in ast.iter_fields(node):
                child=visit(value)
                if isinstance(node,ast.Set) and field=='elts':child=sorted(child,key=lambda item:json.dumps(item,sort_keys=True))
                values.append([field,child])
            return [type(node).__name__,values]
        if isinstance(node,list):return [visit(value) for value in node]
        return node
    return visit(tree)

review=json.loads((source/'review-required.json').read_text());inputs={};records={};checks=[]
for path in sorted((source/'raw').glob('*.json.gz')):
    inputs[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
    mode=path.name.split('-')[0]
    data=json.loads(gzip.decompress(path.read_bytes()))
    for identifier,record in data['results'].items():
        assert (identifier,mode) not in records
        records[identifier,mode]=(path,record)
assert len(records)==2916*2
classified=[];pending=[]
for object_review in review:
    identifier=object_review['id'];all_differences=[];remaining=[]
    for mode in ['ordinary','keywords']:
        path,record=records[identifier,mode];ref=record['reference'];native=record['native']
        if ref.get('status')!='completed' or native.get('status')!='completed':
            remaining.append({'mode':mode,'reason':'backend incomplete'});continue
        a,b=ref.get('examples',[]),native.get('examples',[])
        if len(a)!=len(b):remaining.append({'mode':mode,'reason':'different example counts'});continue
        found=[]
        for index,(left,right) in enumerate(zip(a,b)):
            if left==right:continue
            found.append(index)
            c,d=dict(left),dict(right);lo=c.pop('output',None);ro=d.pop('output',None)
            lc=canonical_output(lo) if type(lo) is str else None
            rc=canonical_output(ro) if type(ro) is str else None
            if c==d and left.get('status')=='ok' and lc is not None and lc==rc:
                all_differences.append({'mode':mode,'index':index,'source':record['sources'][index],'reference_output':lo,'native_output':ro,'archive':str(path.relative_to(root))})
            else:remaining.append({'mode':mode,'index':index,'reason':'difference beyond literal set permutation'})
        if found!=record.get('differences',[]):remaining.append({'mode':mode,'reason':'auditor differences not identical to full observation differences','observed':found,'reported':record.get('differences')})
    if all_differences and not remaining:classified.append({'id':identifier,'verdict':'archived output differs only by literal set element order','observations':all_differences})
    else:pending.append({'id':identifier,'remaining':remaining,'set_order_observations':len(all_differences)})
result={'date':'2026-10-27','source_native_sha256':'576a9f76ab3be1cc0d3896d326e4b6615d1a1ffb1d5f0351f23c522bb37d8ad8','scope':'Historical 26 October displayed observations only; no new runtime/domain/precision/completeness certificate and no current release replay. Every ordinary/keyword example inspected; numeric/string lexemes preserved; exceptions/statuses compared in full.','objects_with_differences':len(review),'objects_only_set_order':len(classified),'observations_only_set_order':sum(len(r['observations']) for r in classified),'objects_requiring_further_review':len(pending),'classified':classified,'remaining':pending,'input_sha256':inputs}
(out/'set-order-review.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ['classified','remaining','input_sha256']})
