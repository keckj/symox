"""The parity auditor must not turn permissive signatures into certificates."""
import inspect
import csv
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from sympy_parity import binding_probes, compare_signatures, function_ledger, keyword_binding_plan, worker


def test_concurrent_snapshots_are_complete_and_detect_later_corruption(tmp_path,monkeypatch):
    import sympy_parity as audit
    root=tmp_path/'workspace';package=root/'python'/'symox';package.mkdir(parents=True)
    (package/'__init__.py').write_text('from .module import value\n')
    (package/'module.py').write_text('value=42\n')
    monkeypatch.setattr(audit,'ROOT',root)
    revision=audit.native_fingerprint();snapshot=tmp_path/'snapshot'
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(audit.prepare_snapshot,snapshot,revision) for _ in range(2)]
        for future in futures:future.result()
    assert audit.native_fingerprint(snapshot/'symox')==revision
    (snapshot/'symox'/'module.py').write_text('value=43\n')
    with pytest.raises(ValueError,match='snapshot does not match'):
        audit.prepare_snapshot(snapshot,revision)


def test_keyword_names_and_positional_options_are_checked():
    def reference(expr, ratio=1.7, *, doit=True): pass
    def native(e, **kwargs): pass
    result = compare_signatures(reference, native)
    assert result['status'] == 'binding_gap'
    assert any(v['positional'] == 2 for v in result['rejected'])
    assert any('expr' in v['keywords'] for v in result['rejected'])


def test_permissive_kwargs_are_not_runtime_certification():
    def reference(x, *, method=None): pass
    def native(x, **kwargs): pass
    result = compare_signatures(reference, native)
    assert result['status'] == 'binding_probes_pass'
    assert result['requires_runtime_option_audit'] is True


def test_positional_only_and_required_keyword_binding():
    def reference(a, /, b=2, *, c): pass
    probes = list(binding_probes(inspect.signature(reference)))
    assert (1, ('c',)) in probes
    assert (2, ('c',)) in probes
    assert all('a' not in keys for _, keys in probes)


def test_keyword_conversion_preserves_positional_only_and_variadic_binding():
    def function(a,/,b,*rest,flag=True):pass
    assert keyword_binding_plan(function,(1,2),{})=={'retained':1,'names':['b']}
    assert keyword_binding_plan(function,(1,2,3),{}) is None
    assert keyword_binding_plan(function,(1,),{'b':2}) is None


def test_native_worker_rewrites_dotted_import_and_retains_state():
    sources = ['import sympy.utilities.iterables\n',
               'sympy.utilities.iterables.flatten([[1], [2]])\n',
               'from sympy import Symbol as Variable\n',
               "x = Variable('x')\n", '(x + 1)**2\n',
               'raise ValueError("expected")\n']
    payload = {'examples': [{'source': s, 'skip': False} for s in sources]}
    reference = worker(dict(payload, library='sympy'), 10)
    native = worker(dict(payload, library='symox'), 10)
    assert reference['status'] == native['status'] == 'completed'
    assert reference['examples'] == native['examples']
    assert native['examples'][-1]['error']['type'] == 'ValueError'


def test_literal_exec_imports_are_rewritten_without_reference_imports():
    sources=['ns={}\n',"exec('from sympy.abc import x', ns)\n",
             "ns['x'] + 1\n","exec('import sympy as s', ns)\n",
             "ns['s'].Rational(1,3)\n"]
    payload={'examples':[{'source':source,'skip':False} for source in sources]}
    reference=worker(dict(payload,library='sympy'),10)
    native=worker(dict(payload,library='symox'),10)
    assert reference['status']==native['status']=='completed'
    assert reference['examples']==native['examples']
    assert all(example['status']=='ok' for example in native['examples'])


def test_real_keyword_calls_use_reference_contracts_and_preserve_callable_identity():
    sources=['from sympy import Symbol, Rational, Integer, Float, Integral, sin, pi\n',
             "x = Symbol('x', positive=True)\n",'Rational(1, 3)\n',
             'for constructor in [Integer, Float]:\n    print(constructor(1))\n',
             'Integral(function=x).doit()\n','pi.n(20)\n','sin(x)\n']
    payload={'examples':[{'source':source,'skip':False} for source in sources],
             'keyword_calls':True}
    reference=worker(dict(payload,library='sympy'),10)
    assert reference['status']=='completed'
    native=worker(dict(payload,library='symox',keyword_plans=reference['keyword_plans']),10)
    assert native['status']=='completed'
    assert reference['examples']==native['examples']
    assert sum(reference['keyword_usage'].values())>=5
    assert reference['keyword_usage']==native['keyword_usage']
    assert reference['successful_keyword_calls']==native['successful_keyword_calls']


@pytest.mark.parametrize('review_sha',['current','outdated'])
def test_scoped_reviews_retain_scope_and_fingerprint_without_universal_certification(tmp_path,review_sha):
    identifier='sympy.combinatorics.Permutation.rank'
    catalog={'native_source_sha256':'current','records':[{
        'id':identifier,'canonical':identifier,'kind':'method',
        'signature':{'status':'binding_probes_pass','reference':'(self)','native':'(self)'}}]}
    review={'native_source_sha256':review_sha,'records':[{
        'id':identifier,'status':'closed_integer_domain_review',
        'precision':'exact integers','complexity':'bounded index operations',
        'scope':'validated finite permutations','universal_function_parity_certified':True}]}
    catalog_path=tmp_path/'catalog.json';catalog_path.write_text(json.dumps(catalog))
    review_path=tmp_path/'review.json';review_path.write_text(json.dumps(review))
    out=tmp_path/'functions.csv'
    function_ledger(SimpleNamespace(catalog=catalog_path,examples=None,certificates=None,
                                   reviews=[review_path],out=out))
    row=next(csv.DictReader(out.open()))
    assert row['review_stale']==str(review_sha!='current')
    assert row['review_scope']=='validated finite permutations'
    assert row['universal_function_parity_certified']=='False'
    assert row['algorithm_review'].startswith('stale_')==(review_sha!='current')


def test_algorithm_review_is_stale_when_rust_changes_without_a_rebuild(tmp_path,monkeypatch):
    import hashlib
    import sympy_parity as audit
    monkeypatch.setattr(audit,'ROOT',tmp_path)
    source=tmp_path/'algorithm.rs';source.write_text('exact_old_algorithm')
    identifier='sympy.example'
    catalog={'native_source_sha256':'unchanged_binary','records':[{
        'id':identifier,'canonical':identifier,'kind':'callable',
        'signature':{'status':'binding_probes_pass','reference':'(x)','native':'(x)'}}]}
    review={'native_source_sha256':'unchanged_binary','records':[{
        'id':identifier,'status':'scoped_invariant','precision':'exact',
        'complexity':'pending','scope':'closed family',
        'native_sources_sha256':{'algorithm.rs':hashlib.sha256(source.read_bytes()).hexdigest()}}]}
    catalog_path=tmp_path/'catalog.json';catalog_path.write_text(json.dumps(catalog))
    review_path=tmp_path/'review.json';review_path.write_text(json.dumps(review))
    out=tmp_path/'functions.csv'
    args=SimpleNamespace(catalog=catalog_path,examples=None,certificates=None,reviews=[review_path],out=out)
    audit.function_ledger(args)
    assert next(csv.DictReader(out.open()))['review_source_stale']=='False'
    source.write_text('changed_algorithm_not_compiled_yet')
    audit.function_ledger(args)
    row=next(csv.DictReader(out.open()))
    assert row['review_source_stale']==row['review_stale']=='True'
    assert row['algorithm_review']=='stale_scoped_invariant'
    assert row['universal_function_parity_certified']=='False'


def test_targeted_replays_override_aliases_without_refreshing_other_evidence(tmp_path):
    from sympy_parity import observation_evidence
    canonical='sympy.function';alias='sympy.module.function';other='sympy.other'
    catalog={'native_source_sha256':'new','reference_version':'1.14.0','records':[
        dict(id=name,canonical=canonical if name!=other else other,kind='callable',
             signature={'status':'binding_probes_pass'}) for name in [canonical,alias,other]]}
    old=dict(status='observation_difference',native_source_sha256='old')
    new=dict(status='examples_match',native_source_sha256='new')
    paths=[]
    for name,records in [('base',{canonical:old,alias:old,other:old}),('target',{canonical:new})]:
        path=tmp_path/(name+'.json')
        path.write_text(json.dumps(dict(reference_version='1.14.0',results=records)))
        paths.append(path)
    catalog_path=tmp_path/'catalog.json';catalog_path.write_text(json.dumps(catalog))
    out=tmp_path/'ledger.csv'
    args=SimpleNamespace(catalog=catalog_path,examples=paths[0],extra_examples=[paths[1]],
                         keyword_examples=paths[0],extra_keyword_examples=[paths[1]],
                         certificates=None,reviews=[],out=out)
    merged=observation_evidence(args,catalog)
    assert merged['native_source_sha256']=='mixed; consult each record'
    assert merged['results']=={canonical:new,other:old}
    function_ledger(args)
    rows={r['reference_callable']:r for r in csv.DictReader(out.open())}
    assert rows[canonical]['example_stale']==rows[canonical]['keyword_stale']=='False'
    assert rows[other]['example_stale']==rows[other]['keyword_stale']=='True'
    assert all(row['universal_function_parity_certified']=='False' for row in rows.values())


def test_targeted_replays_reject_a_different_reference_version(tmp_path):
    from sympy_parity import observation_evidence
    evidence=tmp_path/'evidence.json'
    evidence.write_text(json.dumps(dict(reference_version='different',results={})))
    with pytest.raises(ValueError,match='reference version'):
        observation_evidence(SimpleNamespace(examples=evidence),dict(reference_version='1.14.0',records=[]))
