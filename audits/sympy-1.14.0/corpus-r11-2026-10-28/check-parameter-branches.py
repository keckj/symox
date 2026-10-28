"""Targeted tree/substitution controls, separate from the exactly-once replay."""
import ast,hashlib,importlib,json,pathlib,resource,sys
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
class NoReference:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden in native branch controls')
sys.meta_path.insert(0,NoReference())
import symox as S
root=pathlib.Path(__file__).resolve().parent
out=root/'integration-r11'
rows={r['id']:r for r in json.loads((out/'classified.json').read_text())}
cases={r['id']:r for r in json.loads((out/'cases.json').read_text())}
plan=json.loads((root/'parameter-gap-classification-plan.json').read_text())
ns0={};exec('from symox import *',ns0)
ns0.update(Complexes=S.S.Complexes,Reals=S.S.Reals,ExprCondPair=lambda value,condition:(value,condition))
checks=[];failures=[]
def present(value):return S.srepr(value,order='none')
def apply(value,a,b,method):
    if method=='subs':return value.subs(a,b)
    if method=='xreplace':return value.xreplace({a:b})
    return value.subs({a:b},simultaneous=True)
for item in plan['cases']:
    ident=item['id'];row=rows[ident];case=cases[ident];result={'id':ident,'expected':item['expected_classification'],'controls':[]}
    try:
        assert row['symox']['ok'],row['symox']
        ns=dict(ns0);exec(case['base'],ns);exec(case['setup'],ns)
        a,x=ns['a'],ns['x']
        # Reconstruct the recorded result; do not evaluate the corpus call again.
        value=eval(row['symox']['srepr'],ns)
        assert present(value)==row['symox']['srepr'],'native srepr roundtrip changed'
        call=ast.parse(case['call'],mode='eval').body
        assert isinstance(call,ast.Call) and not call.keywords
        arguments=[eval(compile(ast.Expression(arg),'<recorded-integrand>','eval'),ns) for arg in call.args]
        original=S.Integral(*arguments)
        family=item['family']
        if item['expected_classification']=='retained_affine_nonzero_capability_gap':
            assert value==original
            result['classification']='retained_affine_nonzero_capability_gap'
        else:
            assert value.is_Piecewise and value.args[-1][0]==original and value.args[-1][1]==S.true
            if family in ('1/(x-a)','1/(x+a)'):
                c=-a if family=='1/(x-a)' else a
                expected=S.log(1+c)-S.log(c)
                condition=S.And(S.Contains(c,S.S.Complexes),S.Not(S.Contains(-c,S.Interval(0,1))))
                assert value.args[0][0]==expected and value.args[0][1]==condition
            elif family=='1/(1+a*x)':
                expected=(S.log(1+1/a)-S.log(1/a))/a
                condition=S.And(S.Contains(a,S.S.Complexes),S.Ne(a,0),S.Not(S.Contains(-1/a,S.Interval(0,1))))
                assert value.args[0][0]==1 and value.args[0][1]==S.Eq(a,0)
                assert value.args[1][0]==expected
                assert S.And(value.args[1][1],S.Ne(a,0))==condition
            else:
                expected=-S.log(a) if family=='1/x' else S.pi/(2*S.sqrt(a))
                condition=S.And(S.Contains(a,S.S.Complexes),S.Not(S.Contains(a,S.Interval(-S.oo,0))))
                assert value.args[0][0]==expected and value.args[0][1]==condition
            result.update(classification='certified_conditional_value_with_original_integral_fallback',formula=present(expected),condition=present(condition))
            for iteration in range(2):
                again=value.doit();assert again==value;value=again
            result['controls'].append({'kind':'doit_fixed_point','iterations':2})
        for method in ('subs','xreplace','simultaneous'):
            b=S.Symbol('audit_imaginary',imaginary=True)
            for replacement in (b,b+1):
                transformed=apply(value,a,replacement,method)
                established=apply(original,a,replacement,method)
                assert transformed==established
                result['controls'].append({'kind':'symbolic_chain_restores_original','method':method,'replacement':present(replacement),'actual':present(transformed)})
            for point in (S.Integer(-2),S.Integer(-1),S.Rational(-1,2),S.Integer(0),S.Integer(2),S.I,S.oo,-S.oo,S.zoo,S.nan):
                transformed=apply(value,a,point,method)
                established=apply(original,a,point,method)
                if transformed.has(S.Integral):assert transformed==established
                elif family=='1/(1+a*x)' and point==0:assert transformed==1
                elif item['expected_classification']!='retained_affine_nonzero_capability_gap':
                    assert transformed==apply(expected,a,point,method)
                result['controls'].append({'kind':'numeric_specialization_control','method':method,'point':present(point),'actual':present(transformed),'original':present(established)})
        assert value.subs(x,S.Symbol('unused_bound'))==value
        checks.append(result)
    except Exception as error:
        result['failure']=type(error).__name__+': '+str(error);failures.append(result)
assert 'sympy' not in sys.modules
report={'scope':'targeted controls on reconstructed recorded outputs; no replay of original corpus calls','native_module':S.__file__,'native_sha256':hashlib.sha256((pathlib.Path(S.__file__).parent/'_core.abi3.so').read_bytes()).hexdigest(),'reference_imported':False,'expected_cases':21,'passed_cases':len(checks),'failed_cases':len(failures),'checks':checks,'failures':failures,'interpretation':'Analytic proof in parameter-branch-proof.md establishes domains/branches for all admitted parameters. These exact controls exercise reconstruction, specialization and fixed points; they are not sampled convergence/equivalence proofs.'}
(root/'parameter-branch-controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('native_sha256','reference_imported','expected_cases','passed_cases','failed_cases')},indent=2))
if failures:raise SystemExit(1)
