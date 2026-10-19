"""Versioned, resumable behavioral audit of the documented SymPy API.

``catalog`` examines every exported callable and its public class methods.
``examples`` differentially replays reference docstrings in isolated processes.
Neither binding compatibility nor matching examples is a universal certificate.
No reference implementation is used by the native execution worker.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import contextlib
import csv
import doctest
import gzip
import hashlib
import importlib
from importlib import metadata
import base64
import inspect
import io
import json
import os
from pathlib import Path
import random
import re
import signal
import shutil
import shlex
import subprocess
import sys
import tempfile
import time
import warnings

from sympy_coverage import MODULES, ROOT, exports, load_inventory

AUDITOR_SOURCE=Path(__file__).read_bytes()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(path.suffix + '.pending')
    data = (json.dumps(value, indent=2, ensure_ascii=False) + '\n').encode()
    pending.write_bytes(gzip.compress(data, mtime=0) if path.suffix == '.gz' else data)
    pending.replace(path)


def read_json(path):
    path = Path(path)
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def reference_environment():
    """Pin actual installed sources and numeric backends, beyond a version label."""
    import mpmath
    from sympy.external.gmpy import GROUND_TYPES
    from mpmath.libmp import BACKEND
    distribution=metadata.distribution('sympy')
    record=distribution.read_text('RECORD')
    digest=hashlib.sha256();matched=0;mismatched=[];unhashed=[]
    for name,expected,_ in sorted(csv.reader(record.splitlines())):
        if not name.endswith('.py'):continue
        data=distribution.locate_file(name).read_bytes()
        digest.update(name.encode()+b'\0'+hashlib.sha256(data).digest())
        if not expected:unhashed.append(name);continue
        algorithm,value=expected.split('=',1)
        actual=base64.urlsafe_b64encode(hashlib.new(algorithm,data).digest()).decode().rstrip('=')
        if actual==value:matched+=1
        else:mismatched.append(name)
    result=dict(reference_version=distribution.version,reference_python_source_sha256=digest.hexdigest(),
        installed_record_sha256=hashlib.sha256(record.encode()).hexdigest(),
        manifest_matching_python_files=matched,manifest_mismatches=mismatched,unhashed_python_files=unhashed,
        python_version=sys.version,mpmath_version=mpmath.__version__,
        sympy_ground_types=GROUND_TYPES,mpmath_backend=BACKEND,
        scope='Installed Python sources verified against their distribution RECORD; not a signature over an external wheel or a universal parity certificate.')
    attestation=ROOT/'audits'/('sympy-'+distribution.version)/'reference-wheel-integrity.json'
    if attestation.exists():
        wheel=json.loads(attestation.read_text())
        if wheel.get('reference_python_source_sha256')==result['reference_python_source_sha256'] and not wheel.get('mismatches'):
            result['published_wheel_integrity']=wheel
    return result


def native_fingerprint(package=None):
    package=Path(package) if package is not None else ROOT/'python/symox'
    digest=hashlib.sha256()
    for source in sorted(package.rglob('*')):
        if source.suffix in ('.py','.so'):
            digest.update(str(Path('python/symox')/source.relative_to(package)).encode())
            digest.update(source.read_bytes())
    return digest.hexdigest()


def prepare_snapshot(snapshot,revision):
    """Publish a complete package atomically, including concurrent audit runs."""
    snapshot=Path(snapshot)
    target=snapshot/'symox'
    snapshot.mkdir(parents=True,exist_ok=True)
    if not target.exists():
        with tempfile.TemporaryDirectory(prefix='.snapshot-',dir=snapshot) as staging:
            pending=Path(staging)/'symox'
            shutil.copytree(ROOT/'python/symox',pending,
                            ignore=shutil.ignore_patterns('__pycache__','*.pending'))
            if native_fingerprint(pending)!=revision:
                raise ValueError('Native sources changed while creating the audit snapshot')
            try:pending.rename(target)
            except OSError:
                if not target.exists():raise
    if native_fingerprint(target)!=revision:
        raise ValueError('Native snapshot does not match its recorded source fingerprint')


def atomic_bytes(path,data):
    path=Path(path)
    with tempfile.NamedTemporaryFile(prefix='.'+path.name,dir=path.parent,delete=False) as pending:
        pending.write(data);source=Path(pending.name)
    try:source.replace(path)
    finally:source.unlink(missing_ok=True)


def module_names(inventory):
    names = set(inventory['modules']) | {'sympy.' + m for _, m, _ in MODULES}
    for name in list(names):
        while '.' in name:
            name = name.rsplit('.', 1)[0]
            names.add(name)
    return sorted(names)


def signature(obj):
    try:
        return inspect.signature(obj)
    except (TypeError, ValueError):
        return None


def binding_probes(sig):
    """Witness positional/keyword forms without executing dummy arguments.

    Variadics are bounded probes, not a proof for all argument counts/names.
    """
    P = inspect.Parameter
    params = list(sig.parameters.values())
    pos = [p for p in params if p.kind in (P.POSITIONAL_ONLY, P.POSITIONAL_OR_KEYWORD)]
    required_pos = max((i + 1 for i, p in enumerate(pos) if p.default is P.empty), default=0)
    kwonly = [p for p in params if p.kind is P.KEYWORD_ONLY]
    required_kw = {p.name for p in kwonly if p.default is P.empty}
    trials = [(n, required_kw) for n in range(required_pos, len(pos) + 1)]
    # Every positional-or-keyword name matters, including required parameters.
    for index, p in enumerate(pos):
        if p.kind is P.POSITIONAL_OR_KEYWORD:
            remaining = {q.name for q in pos[index:] if q.default is P.empty}
            trials.append((index, required_kw | remaining | {p.name}))
    for p in kwonly:
        trials.append((required_pos, required_kw | {p.name}))
    trials.append((0, required_kw | {p.name for p in pos if p.kind is P.POSITIONAL_OR_KEYWORD}))
    trials.append((len(pos), {p.name for p in kwonly}))
    if any(p.kind is P.VAR_POSITIONAL for p in params):
        trials.extend((len(pos) + extra, required_kw) for extra in (1, 3))
    if any(p.kind is P.VAR_KEYWORD for p in params):
        trials.append((required_pos, required_kw | {'_parity_extra_keyword'}))
    seen = set()
    for n, keys in trials:
        key = n, tuple(sorted(keys))
        if key in seen:
            continue
        seen.add(key)
        try:
            sig.bind(*([None] * n), **dict.fromkeys(keys))
        except TypeError:
            continue
        yield n, key[1]


def compare_signatures(reference, native):
    a, b = signature(reference), signature(native)
    result = {'reference': str(a) if a is not None else None,
              'native': str(b) if b is not None else None}
    if a is None or b is None:
        result['status'] = 'opaque_signature'
        return result
    rejected = []
    for n, names in binding_probes(a):
        try:
            b.bind(*([None] * n), **dict.fromkeys(names))
        except TypeError as error:
            rejected.append({'positional': n, 'keywords': list(names), 'error': str(error)})
    defaults = []
    for name, p in a.parameters.items():
        q = b.parameters.get(name)
        if q is not None and p.default is not inspect.Parameter.empty:
            if repr(p.default) != repr(q.default):
                defaults.append({'parameter': name, 'reference': repr(p.default),
                                 'native': repr(q.default)})
    result['status'] = 'binding_gap' if rejected else 'binding_probes_pass'
    if rejected:
        result['rejected'] = rejected
    if defaults:
        result['defaults_to_review'] = defaults
    # **kwargs may reject or silently ignore an option inside the function.
    result['requires_runtime_option_audit'] = True
    return result


def catalog(path):
    import sympy
    inventory = load_inventory(sympy.__version__)
    records, modules, identities, keepalive = [], [], {}, []
    implementations, class_surfaces = {}, {}
    def entry(key, reference, native, kind, owner=None, static_unresolved=False):
        canonical = getattr(reference, '__func__', reference)
        identity = id(canonical)
        # Hold references: temporary bound methods must not reuse identities.
        keepalive.extend((reference, canonical, native))
        record = {'id': key, 'kind': kind, 'canonical': identities.setdefault(identity, key),
                  'certified': False}
        if owner:
            record['owner'] = owner
        native_identity = id(getattr(native, '__func__', native))
        pair = identity, native_identity
        if pair in implementations:
            record['evidence_from'] = implementations[pair]
        else:
            implementations[pair] = key
        if native is None:
            record['signature'] = {'status': 'instance_probe_required' if static_unresolved else 'missing'}
        else:
            record['signature'] = compare_signatures(reference, native)
        doc = inspect.getdoc(reference) or ''
        record['doc_sha256'] = hashlib.sha256(doc.encode()).hexdigest()
        try:
            record['examples'] = len(doctest.DocTestParser().get_examples(doc))
        except ValueError as error:
            record['examples'] = 0
            record['doc_error'] = str(error)
        records.append(record)
    for name in module_names(inventory):
        status = {'id': name}
        try:
            reference = importlib.import_module(name)
        except Exception as error:
            status.update(status='reference_import_error', error=f'{type(error).__name__}: {error}')
            modules.append(status)
            continue
        try:
            native = importlib.import_module('symox' + name[5:])
        except Exception as error:
            native = None
            status.update(status='native_import_error', error=f'{type(error).__name__}: {error}')
        else:
            status['status'] = 'imported'
        names = exports(reference, inventory['documented_exports'].get(name, ()))
        status['exports'] = len(names)
        modules.append(status)
        for public in names:
            ref = getattr(reference, public)
            actual = getattr(native, public, None) if native is not None else None
            if not callable(ref):
                records.append({'id': name + '.' + public, 'kind': 'data',
                                'present': actual is not None, 'certified': False})
                continue
            key = name + '.' + public
            entry(key, ref, actual, 'class' if inspect.isclass(ref) else 'callable')
            if inspect.isclass(ref):
                pair = id(ref), id(actual)
                if pair in class_surfaces:
                    records[-1]['methods_from'] = class_surfaces[pair]
                    continue
                class_surfaces[pair] = key
                for method in sorted(n for n in dir(ref) if not n.startswith('_')):
                    try:
                        value = getattr(ref, method)
                    except Exception:
                        continue
                    if not (inspect.isroutine(value) or inspect.ismethoddescriptor(value)):
                        continue
                    try:
                        native_method = getattr(actual, method) if actual is not None else None
                    except Exception:
                        native_method = None
                    entry(key + '.' + method, value, native_method, 'method', key,
                          static_unresolved=actual is not None)
        print(f'{name}: {len(records)} bindings', file=sys.stderr, flush=True)
    counts = Counter(r.get('signature', {}).get('status', 'data') for r in records)
    report = dict(schema=1, reference_version=sympy.__version__,
                  reference_environment=reference_environment(),
                  native_source_sha256=native_fingerprint(),
                  inventory_sha256=hashlib.sha256((ROOT/'tools/sympy_public_modules.json').read_bytes()).hexdigest(),
                  scope='module exports and public callable class methods; no universal certificates',
                  modules=modules, counts=dict(counts), records=records)
    atomic_json(path, report)
    print(json.dumps({'modules': len(modules), 'bindings': len(records),
                      'canonical_callables': len(identities), 'counts': counts}))


class NativeImports(ast.NodeTransformer):
    def visit_Call(self,node):
        node=self.generic_visit(node)
        # A few reference docs demonstrate imports inside literal exec/eval
        # snippets. Rewrite their import statements too; the native import
        # guard must remain strict and must never load the reference library.
        if (isinstance(node.func,ast.Name) and node.func.id in ('exec','eval')
                and node.args and isinstance(node.args[0],ast.Constant)
                and isinstance(node.args[0].value,str)):
            try:
                source=ast.parse(node.args[0].value,mode='exec' if node.func.id=='exec' else 'eval')
            except SyntaxError:return node
            source=ast.fix_missing_locations(NativeImports().visit(source))
            node.args[0]=ast.copy_location(ast.Constant(ast.unparse(source)),node.args[0])
        return node

    def visit_Import(self, node):
        # Preserve the local binding of ``import sympy`` and dotted imports.
        # Explicit aliases work directly; unaliased imports need sympy bound
        # to the native root so later ``sympy.foo`` uses native code.
        replacements = []
        for alias in node.names:
            if alias.name == 'sympy' or alias.name.startswith('sympy.'):
                rewritten = 'symox' + alias.name[5:]
                if alias.asname:
                    replacements.append(ast.Import(names=[ast.alias(name=rewritten, asname=alias.asname)]))
                else:
                    replacements.append(ast.Import(names=[ast.alias(name=rewritten)]))
                    replacements.append(ast.Import(names=[ast.alias(name='symox', asname='sympy')]))
            else:
                replacements.append(ast.Import(names=[alias]))
        return [ast.copy_location(n, node) for n in replacements]

    def visit_ImportFrom(self, node):
        if node.module and (node.module == 'sympy' or node.module.startswith('sympy.')):
            node.module = 'symox' + node.module[5:]
        return node


def resolve(module, name):
    obj = importlib.import_module(module)
    for part in name.split('.'):
        obj = getattr(obj, part)
    return obj


class KeywordCalls(ast.NodeTransformer):
    """Route each doc call through a runtime argument-binding probe."""
    def __init__(self,example):self.example=example

    def visit_Call(self,node):
        node=self.generic_visit(node)
        spelling=(node.func.id if isinstance(node.func,ast.Name) else
                  node.func.attr if isinstance(node.func,ast.Attribute) else '')
        site=f'{self.example}:{node.lineno}:{node.col_offset}'
        return ast.copy_location(ast.Call(
            func=ast.Name(id='__symox_audit_keyword_call__',ctx=ast.Load()),
            args=[ast.Constant(site),ast.Constant(spelling),node.func,*node.args],
            keywords=node.keywords),node)


def keyword_binding_plan(function,args,kwargs):
    sig=signature(function)
    if sig is None:return None
    try:sig.bind(*args,**kwargs)
    except TypeError:return None
    P=inspect.Parameter
    positional=[p for p in sig.parameters.values() if p.kind in (P.POSITIONAL_ONLY,P.POSITIONAL_OR_KEYWORD)]
    if len(args)>len(positional):return None
    retained=sum(p.kind is P.POSITIONAL_ONLY for p in positional[:len(args)])
    names=[p.name for p in positional[retained:len(args)]]
    return {'retained':retained,'names':names} if names else None


def execute_document(payload):
    """One library in one disposable worker, with no oracle imports."""
    # Documentation preview examples can launch desktop viewers. Keep the
    # actual TeX/code-generation pipeline available, but never open a GUI.
    original_popen = subprocess.Popen
    compiler_tools = {'latex', 'pdflatex', 'xelatex', 'lualatex', 'dvipng',
                      'dvips', 'dvisvgm', 'convert', 'magick', 'gs', 'gcc',
                      'cc', 'g++', 'c++', 'clang', 'clang++', 'gfortran',
                      'f95', 'f77', 'ninja', 'make', 'cmake', 'pkg-config',
                      'python', 'python3', Path(sys.executable).name}
    class ExternalProcessBlocked(RuntimeError):
        pass
    def guarded_popen(command, *args, **kwargs):
        words = shlex.split(command) if isinstance(command, str) else command
        program = Path(os.fsdecode(words[0])).name if words else ''
        if program not in compiler_tools:
            raise ExternalProcessBlocked('desktop/external process blocked by audit: ' + program)
        return original_popen(command, *args, **kwargs)
    subprocess.Popen = guarded_popen
    def guarded_system(command):
        return subprocess.call(shlex.split(command))
    os.system = guarded_system
    import webbrowser
    def blocked_browser(*args, **kwargs):
        raise ExternalProcessBlocked('browser blocked by audit')
    webbrowser.open = webbrowser.open_new = webbrowser.open_new_tab = blocked_browser
    native = payload['library'] == 'symox'
    if native:
        class NoReference:
            def find_spec(self, fullname, path=None, target=None):
                if fullname == 'sympy' or fullname.startswith('sympy.'):
                    raise ImportError('reference library forbidden in the native audit worker')
        sys.meta_path.insert(0, NoReference())
    namespace = {'__name__': '__parity_doc__'}
    exec('from ' + payload['library'] + ' import *', namespace)
    # SymPy's documentation examples conventionally have these symbols.
    exec('from ' + payload['library'] + '.abc import x, y, z, t, k, n', namespace)
    random.seed(1729)
    importlib.import_module(payload['library']+'.core.random').seed(1729)
    keyword_mode=payload.get('keyword_calls',False)
    plans=dict(payload.get('keyword_plans',{}));usage=Counter();successful=Counter()
    if keyword_mode:
        symbolic_function=namespace['Function']
        def keyword_call(site,spelling,function,/,*args,**kwargs):
            name=getattr(function,'__name__',type(function).__name__)
            if spelling=='n' and name in ('n','evalf'):name='n'
            key=json.dumps([site,spelling,name,len(args),sorted(kwargs)])
            if native:plan=plans.get(key)
            else:
                module=getattr(function,'__module__','') or ''
                eligible=(module=='sympy' or module.startswith('sympy.')) and (
                    inspect.isroutine(function) or inspect.isclass(function))
                # Function heads display named signatures but their runtime
                # call interface interprets those keywords as options instead.
                if inspect.isclass(function) and issubclass(function,symbolic_function):eligible=False
                plan=keyword_binding_plan(function,args,kwargs) if eligible else None
                if plan is not None:plans[key]=plan
            if plan is None:return function(*args,**kwargs)
            usage[key]+=1
            retained=plan['retained'];options=dict(kwargs)
            options.update(zip(plan['names'],args[retained:]))
            result=function(*args[:retained],**options)
            successful[key]+=1
            return result
        namespace['__symox_audit_keyword_call__']=keyword_call
    outcomes = []
    for sample_index,sample in enumerate(payload['examples']):
        if sample['skip']:
            outcomes.append({'status': 'skipped'})
            continue
        output = io.StringIO()
        old_display = sys.displayhook
        def display(value):
            if value is not None:
                namespace['_'] = value
                print(repr(value))
        sys.displayhook = display
        error = None
        try:
            tree = ast.parse(sample['source'])
            if native:
                tree = NativeImports().visit(tree)
            if keyword_mode:tree=KeywordCalls(sample_index).visit(tree)
            tree = ast.fix_missing_locations(ast.Interactive(body=tree.body))
            with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
                exec(compile(tree, '<doc-example>', 'single'), namespace)
        except BaseException as exc:
            error = {'type': type(exc).__name__, 'message': str(exc)}
        finally:
            sys.displayhook = old_display
        # Addresses and library prefixes carry no mathematical meaning.
        value = re.sub(r'0x[0-9a-fA-F]+', '0xADDRESS', output.getvalue())
        value = value.replace('symox.', 'sympy.')
        result = {'status': 'exception' if error else 'ok', 'output': value}
        if error:
            result['error'] = error
        outcomes.append(result)
    if keyword_mode:
        return {'examples':outcomes,'keyword_plans':plans,
                'keyword_usage':dict(usage),'successful_keyword_calls':dict(successful)}
    return outcomes


def worker(payload, timeout, package_root=None, audit_script=None):
    env = dict(os.environ, MPLBACKEND='Agg', PYTHONHASHSEED='0',
               DISPLAY='', WAYLAND_DISPLAY='')
    if package_root:
        env['PYTHONPATH'] = str(package_root)
    with tempfile.TemporaryDirectory(prefix='symox-parity-') as cwd:
        try:
            result = subprocess.run([sys.executable, str(audit_script or Path(__file__).resolve()), '_worker'],
                input=json.dumps(payload), text=True, capture_output=True, timeout=timeout,
                cwd=cwd, env=env)
        except subprocess.TimeoutExpired:
            return {'status': 'timeout', 'seconds': timeout}
        if result.returncode:
            return {'status': 'worker_error', 'code': result.returncode,
                    'stderr': result.stderr[-4000:]}
        try:
            value=json.loads(result.stdout)
            return dict(status='completed',**value) if isinstance(value,dict) else {'status':'completed','examples':value}
        except ValueError:
            return {'status': 'invalid_worker_output', 'stdout': result.stdout[-2000:],
                    'stderr': result.stderr[-2000:]}


def examples(args):
    ledger = read_json(args.catalog)
    import sympy
    if sympy.__version__ != ledger['reference_version']:
        raise ValueError('Reference version differs from the catalog')
    path = Path(args.out)
    report = read_json(path) if path.exists() and not args.fresh else {
        'schema': 1, 'reference_version': sympy.__version__, 'results': {}}
    environment=reference_environment()
    if environment['manifest_mismatches']:
        raise ValueError('Installed reference sources differ from their distribution RECORD')
    report['reference_environment']=environment
    # Fingerprint native Python sources and backend, so stale evidence is visible.
    existing_snapshot=Path(args.snapshot)/'symox' if args.snapshot else None
    revision = native_fingerprint(existing_snapshot if existing_snapshot and existing_snapshot.exists() else None)
    report['native_source_sha256'] = revision
    # Audit an immutable native package snapshot. Edits made while a long audit
    # runs cannot acquire a misleading fingerprint from its start time.
    snapshot = Path(args.snapshot) if args.snapshot else path.parent / ('native-' + revision[:16])
    prepare_snapshot(snapshot,revision)
    report['snapshot'] = str(snapshot.resolve())
    auditor_sha=hashlib.sha256(AUDITOR_SOURCE).hexdigest()
    report['auditor_sha256'] = auditor_sha
    report['mode']='keyword_calls' if args.keywords else 'ordinary_examples'
    audit_tools=snapshot/'_audit_tools'/auditor_sha[:16]
    audit_tools.mkdir(parents=True,exist_ok=True)
    audit_script=audit_tools/'sympy_parity.py'
    atomic_bytes(audit_script,AUDITOR_SOURCE)
    atomic_bytes(audit_tools/'sympy_coverage.py',Path(__file__).with_name('sympy_coverage.py').read_bytes())
    records = sorted(ledger['records'], key=lambda r: r['kind'] == 'method')
    seen, count = set(), 0
    for record in records:
        key = record['id']
        if record.get('kind') == 'data' or not record.get('examples'):
            continue
        if args.module and not any(key == m or key.startswith(m + '.') for m in args.module):
            continue
        # Deduplicate only after selecting the requested public surfaces.
        canonical = record['canonical']
        if canonical in seen:
            continue
        seen.add(canonical)
        previous = report['results'].get(key)
        if args.pending_only and previous:
            continue
        if (previous and previous.get('native_source_sha256') == revision
                and previous.get('auditor_sha256') == auditor_sha
                and previous.get('mode','ordinary_examples') == report['mode']
                and previous.get('doc_sha256') == record['doc_sha256'] and not args.retry):
            continue
        # Find the longest actual module prefix; nested classes are allowed.
        module = max((m['id'] for m in ledger['modules']
                      if m['status']!='reference_import_error' and key.startswith(m['id'] + '.')), key=len)
        obj = resolve(module, key[len(module) + 1:])
        doc = inspect.getdoc(obj) or ''
        samples = doctest.DocTestParser().get_examples(doc)
        payload = {'keyword_calls':args.keywords,
                   'examples': [{'source': s.source, 'skip': bool(s.options.get(doctest.SKIP))}
                                for s in samples]}
        start = time.monotonic()
        reference = worker(dict(payload, library='sympy'), args.timeout, audit_script=audit_script)
        if args.keywords:payload['keyword_plans']=reference.get('keyword_plans',{})
        native = worker(dict(payload, library='symox'), args.timeout, snapshot, audit_script)
        result = {'native_source_sha256': revision, 'doc_sha256': record['doc_sha256'],
                  'reference_python_source_sha256':environment['reference_python_source_sha256'],
                  'auditor_sha256':auditor_sha,
                  'mode':report['mode'],
                  'seconds': round(time.monotonic() - start, 3), 'certified': False,
                  'reference': reference, 'native': native,
                  'sources': [s.source for s in samples]}
        if reference['status'] != 'completed' or native['status'] != 'completed':
            result['status'] = 'incomplete'
        else:
            differences = [i for i, (a, b) in enumerate(zip(reference['examples'], native['examples'])) if a != b]
            result['differences'] = differences
            result['status'] = 'observation_difference' if differences else 'examples_match'
            # Reference failures are explicit: equal failures do not establish support.
            result['reference_exceptions'] = [i for i, value in enumerate(reference['examples'])
                                               if value['status'] == 'exception']
            categories = Counter()
            for a, b in zip(reference['examples'], native['examples']):
                if a['status'] == b['status'] == 'skipped':
                    categories['skipped'] += 1
                elif a == b:
                    categories['matched_exception' if a['status'] == 'exception' else 'matched_success'] += 1
                elif a['status'] == 'ok' and b['status'] == 'exception':
                    categories['native_exception'] += 1
                elif a['status'] == 'exception':
                    categories['reference_exception_difference'] += 1
                else:
                    categories['output_difference'] += 1
            result['case_counts'] = dict(categories)
            if args.keywords:
                a,b=reference.get('keyword_usage',{}),native.get('keyword_usage',{})
                result['keyword_counts']={
                    'reference_calls':sum(a.values()),'native_calls':sum(b.values()),
                    'reference_successful_calls':sum(reference.get('successful_keyword_calls',{}).values()),
                    'native_successful_calls':sum(native.get('successful_keyword_calls',{}).values()),
                    'unmatched_reference_calls':sum(max(0,count-b.get(key,0)) for key,count in a.items())}
            if not categories['matched_success'] and not differences:
                result['status'] = 'no_successful_examples'
        report['results'][key] = result
        report['counts'] = dict(Counter(r['status'] for r in report['results'].values()))
        atomic_json(path, report)
        print(f"{key}: {result['status']} ({result['seconds']}s)", flush=True)
        count += 1
        if args.limit and count >= args.limit:
            break
    print(json.dumps(report['counts']))


def observation_evidence(args, catalog, keyword=False):
    """Overlay targeted runs by reference identity, retaining record provenance."""
    attribute='keyword_examples' if keyword else 'examples'
    base=getattr(args,attribute,None)
    paths=([base] if base else [])+list(getattr(args,'extra_'+attribute,None) or [])
    identities={row['id']:row.get('canonical',row['id']) for row in catalog['records']}
    selected={};metadata={}
    for path in paths:
        package=read_json(path)
        if (package.get('reference_version') and catalog.get('reference_version') and
                package['reference_version']!=catalog['reference_version']):
            raise ValueError('Example evidence uses a different reference version')
        metadata=package
        for key,result in package['results'].items():
            selected[identities.get(key,key)]=(key,result)
    result=dict(metadata,results=dict(selected.values()),evidence_files=[str(path) for path in paths])
    if len(paths)>1:
        fingerprints={row.get('native_source_sha256','unrecorded') for _,row in selected.values()}
        result['native_source_sha256']=(next(iter(fingerprints)) if len(fingerprints)==1
                                      else 'mixed; consult each record')
    return result


def render_report(args):
    ledger=read_json(args.catalog)
    evidence=observation_evidence(args,ledger)
    counts=Counter(row['status'] for row in evidence['results'].values())
    keyword=observation_evidence(args,ledger,keyword=True)
    keyword_counts=Counter(row['status'] for row in keyword['results'].values())
    calls=Counter()
    for row in keyword['results'].values():calls.update(row.get('keyword_counts',{}))
    reference=ledger.get('reference_environment',{})
    lines=['# Versioned SymPy parity audit','',
        f"Reference: **SymPy {ledger['reference_version']}**. Global algorithmic parity is **not certified**.",
        '', 'The catalog examines every documented module export and public callable class method.',
        'A binding probe checks Python argument binding, including keyword names and positional options.',
        'A passing binding probe does not establish that an option is implemented. A rejected probe',
        'can also expose an overly permissive reference signature whose implementation rejects that call.',
        '', 'Documentation examples execute the reference and native libraries in separate processes,',
        'with timeouts and an immutable native package snapshot. Reference imports and desktop viewers',
        'are blocked in native workers. Output differences require review: different printing is not',
        'automatically a mathematical error. Equal exceptions are recorded separately from successful examples.',
        '', '## Evidence','',
        f"- Module paths: {len(ledger['modules'])}.",
        f"- Export/method bindings, including aliases: {len(ledger['records'])}.",
        f"- Distinct reference callables: {len({r['canonical'] for r in ledger['records'] if 'canonical' in r})}.",
        f"- Documented objects replayed: {len(evidence['results'])}.",
        f"- Matching example observations: {counts['examples_match']} objects.",
        f"- Different observations to review: {counts['observation_difference']} objects.",
        f"- Incomplete replays: {counts['incomplete']} objects.",
        f"- Named-call objects replayed: {len(keyword['results'])}.",
        f"- Matching named-call observations: {keyword_counts['examples_match']} objects.",
        f"- Different named-call observations to review: {keyword_counts['observation_difference']} objects.",
        f"- Incomplete named-call replays: {keyword_counts['incomplete']} objects.",
        '', f"Reference Python source fingerprint: `{reference.get('reference_python_source_sha256','not recorded')}`.",
        f"Installed source files matching distribution RECORD: {reference.get('manifest_matching_python_files','not recorded')}; mismatches: {len(reference.get('manifest_mismatches',[]))}.",
        'The separate [published-wheel integrity record](../audits/sympy-1.14.0/reference-wheel-integrity.json)',
        'compares these sources with the [SymPy 1.14.0 PyPI release](https://pypi.org/project/sympy/1.14.0/#files), without modifying the reference installation.',
        '', f"Catalog native fingerprint: `{ledger.get('native_source_sha256','not recorded in initial snapshot')}`.",
        '', f"Example snapshot fingerprint: `{evidence.get('native_source_sha256','not run')}`.",
        '', f"Named-call snapshot fingerprint: `{keyword.get('native_source_sha256','not run')}`.",
        '', 'Named-call replays use actual reference argument binding to move eligible positional inputs',
        'to their public keyword names, preserving positional-only and variadic calls. They record',
        f"{calls['reference_calls']} reference calls and {calls['native_calls']} native calls, with",
        f"{calls['reference_successful_calls']} and {calls['native_successful_calls']} successful calls respectively.",
        'Cascading exceptions and repeated calls affect these counts; they are not a compliance percentage.',
        '', 'These fingerprints may differ while fixes are being made. Each example record keeps its',
        'own fingerprint; stale evidence must be rerun before approving a changed implementation.',
        'Targeted supplementary runs replace older observations by reference callable identity,',
        'including aliases. Other records retain their original snapshots and stale status.',
        'Ordinary evidence files: '+', '.join('`'+p+'`' for p in evidence['evidence_files'])+'.',
        'Named-call evidence files: '+', '.join('`'+p+'`' for p in keyword['evidence_files'])+'.',
        'Scoped review source hashes are also checked against the working tree, so edited Rust',
        'cannot retain a fresh review merely because an older compiled extension is still installed.',
        '', 'Scoped source-equivalence checks for container algorithms are recorded in',
        '[source-certificates.json](../audits/sympy-1.14.0/source-certificates.json).',
        'These compare signatures, normalized control flow and closed standard-library dependencies,',
        'including error propagation and time/space behavior. Symbolic input-object protocols are a separate obligation.',
        '', 'Scoped invariant, exact-arithmetic and complexity reviews are recorded in',
        '[permutation-review.json](../audits/sympy-1.14.0/permutation-review.json) and',
        '[group-review.json](../audits/sympy-1.14.0/group-review.json).',
        'Complete elementary solution-lattice invariants are recorded separately in',
        '[solver-review.json](../audits/sympy-1.14.0/solver-review.json).',
        'Radical controls and principal complex-branch invariants are recorded in',
        '[radical-review.json](../audits/sympy-1.14.0/radical-review.json).',
        'Scoped gamma-moment convergence and public Piecewise branch protocols are recorded in',
        '[integration-review.json](../audits/sympy-1.14.0/integration-review.json).',
        'Their finite-domain assumptions and unresolved complexity obligations are retained in each record;',
        'they do not certify an entire constructor, inherited protocol or unsupported optional interface.',
        '', '## Module traversal','',
        '“Instance” means that a method could not be inspected on its native constructor; it needs',
        'an instance probe. This is not a claim that the method is absent on returned objects.',
        '', '| Module | Exports | Method bindings | Binding gaps | Instance probes |',
        '| --- | ---: | ---: | ---: | ---: |']
    for module in ledger['modules']:
        name=module['id']
        # Methods carry their exported owner; longest module matching excludes
        # child-module rows from the parent package's counts.
        rows=[r for r in ledger['records'] if r['id'].rsplit('.',2 if r['kind']=='method' else 1)[0]==name]
        statuses=Counter(r.get('signature',{}).get('status') for r in rows)
        if module['status']!='imported':
            lines.append(f"| `{name}` | {module['status']} | — | — | — |")
        else:
            lines.append(f"| `{name}` | {module.get('exports',0)} | {sum(r['kind']=='method' for r in rows)} | {statuses['binding_gap']} | {statuses['instance_probe_required']} |")
    lines.extend(['','## Function and constructor binding gaps','',
        'Aliases share their reference identity below. The compressed catalog retains every public binding,',
        'all method signatures and concrete rejected-call witnesses. These rows are an open work queue.',
        '', '| Reference callable | Reference signature | Native signature |','| --- | --- | --- |'])
    seen=set()
    for row in ledger['records']:
        if row['kind']=='method' or row.get('signature',{}).get('status')!='binding_gap':continue
        if row['canonical'] in seen:continue
        seen.add(row['canonical'])
        sig=row['signature']
        clean=lambda text:re.sub(r'0x[0-9a-fA-F]+','0x…',str(text)).replace('|','\\|')
        lines.append(f"| `{row['id']}` | `{clean(sig['reference'])}` | `{clean(sig['native'])}` |")
    lines.extend(['','## Reproduce and resume','', '```sh',
        'python tools/sympy_parity.py catalog --out audits/sympy-1.14.0/catalog.json.gz',
        'python tools/sympy_parity.py examples --catalog audits/sympy-1.14.0/catalog.json.gz \\',
        '  --out /tmp/symox-parity/nightly.json --timeout 15',
        'python tools/sympy_parity.py report --catalog audits/sympy-1.14.0/catalog.json.gz \\',
        '  --examples /tmp/symox-parity/nightly.json --out docs/sympy-parity-audit.md',
        'python tools/sympy_stdlib_certificates.py --out audits/sympy-1.14.0/source-certificates.json',
        'python tools/sympy_parity.py examples --catalog audits/sympy-1.14.0/catalog.json.gz --keywords \\',
        '  --out /tmp/symox-parity/keywords.json --timeout 20',
        '```','', 'Use `--module sympy.utilities.iterables` to target a surface, `--limit N` to bound a batch,',
        '`--retry` to rerun existing evidence, and `--pending-only` to continue traversal while retaining',
        'explicitly stale records. Normally a changed native fingerprint causes previous examples to be replayed.',
        '', 'See [the capability contract and remaining algorithmic gaps](algorithmic-parity.md).',''])
    Path(args.out).write_text('\n'.join(lines))


def function_ledger(args):
    catalog=read_json(args.catalog)
    observations=observation_evidence(args,catalog)['results']
    keyword_observations=observation_evidence(args,catalog,keyword=True)['results']
    by_id={r['id']:r for r in catalog['records']}
    groups={}
    for row in catalog['records']:
        if 'canonical' in row:groups.setdefault(row['canonical'],[]).append(row)
    evidence={}
    for key,result in observations.items():
        if key in by_id:evidence[by_id[key]['canonical']]=(key,result)
    keyword_evidence={}
    for key,result in keyword_observations.items():
        if key in by_id:keyword_evidence[by_id[key]['canonical']]=(key,result)
    certificate_package=read_json(args.certificates) if args.certificates else {}
    certificates=certificate_package.get('records',[])
    certified={by_id[r['id']]['canonical']:r for r in certificates if r['id'] in by_id}
    reviews={}
    source_digests={}
    def changed_review_sources(record):
        sources=record.get('native_sources_sha256',{})
        if not isinstance(sources,dict):return True
        for name,expected in sources.items():
            if name not in source_digests:
                try:source_digests[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                except OSError:source_digests[name]=None
            if source_digests[name]!=expected:return True
        return False
    for filename in args.reviews or []:
        review=read_json(filename)
        for record in review['records']:
            if record['id'] in by_id:
                reviews[by_id[record['id']]['canonical']]=(record,review,filename)
    columns=['reference_callable','kind','aliases','reference_signature','native_signature',
        'binding_states','documentation_examples','example_observations','example_evidence_id',
        'native_snapshot_sha256','example_stale','first_difference_source','first_reference_observation',
        'keyword_observations','keyword_evidence_id','keyword_snapshot_sha256','keyword_stale',
        'keyword_reference_calls','keyword_native_calls','keyword_reference_successful_calls',
        'keyword_native_successful_calls','keyword_unmatched_reference_calls','keyword_first_difference_source',
        'first_native_observation','algorithm_review','domain_review','precision_review',
        'complexity_review','algorithm_evidence','review_native_sha256','review_stale','review_source_stale','review_scope',
        'universal_function_parity_certified']
    output=io.StringIO();writer=csv.DictWriter(output,fieldnames=columns);writer.writeheader()
    for canonical,rows in sorted(groups.items()):
        row=by_id[canonical]
        key,result=evidence.get(canonical,('',{}))
        difference=(result.get('differences') or [None])[0]
        sig=row.get('signature',{})
        item=dict(reference_callable=canonical,kind=row['kind'],
            aliases=json.dumps([r['id'] for r in rows if r['id']!=canonical]),
            reference_signature=sig.get('reference'),native_signature=sig.get('native'),
            binding_states='|'.join(sorted({r['signature']['status'] for r in rows})),
            documentation_examples=row.get('examples',0),
            example_observations=result.get('status','not_replayed'),example_evidence_id=key,
            native_snapshot_sha256=result.get('native_source_sha256',''),
            example_stale=bool(result and result.get('native_source_sha256')!=catalog.get('native_source_sha256')),
            algorithm_review='source_equivalent' if canonical in certified else 'pending',
            domain_review='standard_input_protocols' if canonical in certified else 'pending',
            precision_review='opaque_objects_preserved' if canonical in certified else 'pending',
            complexity_review='same_standard_library_operations' if canonical in certified else 'pending',
            universal_function_parity_certified=False)
        if canonical in certified:
            proof=certified[canonical]
            stale=certificate_package.get('native_source_sha256')!=catalog.get('native_source_sha256')
            item.update(algorithm_review=('stale_' if stale else '')+'source_equivalent',
                algorithm_evidence=str(args.certificates),
                review_native_sha256=certificate_package.get('native_source_sha256',''),
                review_stale=stale,review_scope=proof.get('scope',''))
        if canonical in reviews:
            review,package,filename=reviews[canonical]
            source_stale=changed_review_sources(review)
            stale=source_stale or package.get('native_source_sha256')!=catalog.get('native_source_sha256')
            item.update(algorithm_review=('stale_' if stale else '')+review['status'],
                domain_review='closed_domain_only',precision_review=review['precision'],
                complexity_review=review['complexity'],algorithm_evidence=filename,
                review_native_sha256=package.get('native_source_sha256'),review_stale=stale,review_source_stale=source_stale,
                review_scope=review['scope'])
        keyword_key,keyword=keyword_evidence.get(canonical,('',{}))
        keyword_counts=keyword.get('keyword_counts',{})
        item.update(keyword_observations=keyword.get('status','not_replayed'),
            keyword_evidence_id=keyword_key,
            keyword_snapshot_sha256=keyword.get('native_source_sha256',''),
            keyword_stale=bool(keyword and keyword.get('native_source_sha256')!=catalog.get('native_source_sha256')))
        for field in ('reference_calls','native_calls','reference_successful_calls',
                      'native_successful_calls','unmatched_reference_calls'):
            item['keyword_'+field]=keyword_counts.get(field,'')
        if keyword.get('differences'):
            item['keyword_first_difference_source']=keyword['sources'][keyword['differences'][0]]
        if difference is not None:
            item['first_difference_source']=result['sources'][difference]
            item['first_reference_observation']=json.dumps(result['reference']['examples'][difference])
            item['first_native_observation']=json.dumps(result['native']['examples'][difference])
        writer.writerow(item)
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    data=output.getvalue().encode()
    path.write_bytes(gzip.compress(data,mtime=0) if path.suffix=='.gz' else data)
    print(f'{len(groups)} function-level review rows written')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('catalog')
    p.add_argument('--out', required=True)
    p = sub.add_parser('examples')
    p.add_argument('--catalog', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--module', action='append')
    p.add_argument('--limit', type=int)
    p.add_argument('--timeout', type=float, default=30)
    p.add_argument('--fresh', action='store_true')
    p.add_argument('--retry', action='store_true')
    p.add_argument('--pending-only', action='store_true')
    p.add_argument('--snapshot', help='Directory for the immutable native package snapshot')
    p.add_argument('--keywords',action='store_true',help='Replay real doc inputs with reference positional-or-keyword names')
    p=sub.add_parser('report')
    p.add_argument('--catalog',required=True)
    p.add_argument('--examples')
    p.add_argument('--keyword-examples')
    p.add_argument('--extra-examples',action='append',help='Overlay a targeted replay, keeping per-record fingerprints')
    p.add_argument('--extra-keyword-examples',action='append',help='Overlay targeted named-call observations')
    p.add_argument('--out',required=True)
    p=sub.add_parser('ledger')
    p.add_argument('--catalog',required=True)
    p.add_argument('--examples')
    p.add_argument('--keyword-examples',help='Separate runtime named-call observations with independent fingerprints')
    p.add_argument('--extra-examples',action='append',help='Overlay a targeted replay, keeping per-record fingerprints')
    p.add_argument('--extra-keyword-examples',action='append',help='Overlay targeted named-call observations')
    p.add_argument('--certificates')
    p.add_argument('--reviews',action='append',help='Scoped algorithm review JSON, with its native fingerprint')
    p.add_argument('--out',required=True)
    sub.add_parser('_worker')
    args = parser.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if args.command == 'catalog':
            catalog(args.out)
        elif args.command == '_worker':
            print(json.dumps(execute_document(json.load(sys.stdin))))
        elif args.command == 'examples':
            examples(args)
        elif args.command == 'report':
            render_report(args)
        else:
            function_ledger(args)


if __name__ == '__main__':
    main()
