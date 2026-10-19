"""Measure API presence and a curated differential corpus, not global parity.

Run with a development installation and SymPy installed::

    python tools/sympy_coverage.py --corpus --json /tmp/symox-coverage.json

API names say nothing about signatures, domains or behavior. The expression
corpus comes from symox's tests and is intentionally biased toward implemented
features. Neither percentage estimates the fraction of SymPy implemented.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import inspect
import importlib
import json
from pathlib import Path
import runpy
import subprocess
import warnings


ROOT = Path(__file__).resolve().parents[1]

# Public import surfaces, including modules that symox has not started yet.
# A missing module has zero available names; it is never skipped.
MODULES = (
    ('Geometry', 'geometry', 'geometry'),
    ('Ordinary differential equations', 'solvers.ode', 'solvers.ode'),
    ('Partial differential equations', 'solvers.pde', 'solvers.pde'),
    ('Equation solvers', 'solvers.solvers', 'solvers.solvers'),
    ('Inequalities', 'solvers.inequalities', 'solvers.inequalities'),
    ('Arrays', 'tensor.array', 'tensor.array'),
    ('Indexed tensor algebra', 'tensor.tensor', 'tensor.tensor'),
    ('Series and limits', 'series', 'series'),
    ('Sequences', 'series.sequences', 'series.sequences'),
    ('Real analysis', 'calculus.util', 'calculus.util'),
    ('Singularities and monotonicity', 'calculus.singularities', 'calculus.singularities'),
    ('Coefficient domains', 'polys.domains', 'polys.domains'),
    ('Polynomial ring constructors', 'polys.rings', 'polys.rings'),
    ('Fraction field constructors', 'polys.fields', 'polys.fields'),
    ('Number fields', 'polys.numberfields', 'polys.numberfields'),
    ('Combinatorics', 'combinatorics', 'combinatorics'),
    ('Symbolic sets', 'sets', None),
    ('Simplification', 'simplify', None),
    ('Integral transforms', 'integrals.transforms', '_transforms'),
    ('Special functions: gamma', 'functions.special.gamma_functions', None),
    ('Special functions: hypergeometric', 'functions.special.hyper', None),
    ('Special functions: Bessel', 'functions.special.bessel', None),
    ('Special functions: Mathieu', 'functions.special.mathieu_functions', None),
    ('Complex functions and branches', 'functions.elementary.complexes', None),
    ('Number theory', 'ntheory', 'ntheory'),
    ('Iterable utilities', 'utilities.iterables', 'utilities.iterables'),
    ('Printing and code generation', 'printing', 'printing'),
    ('Plotting', 'plotting', 'plotting'),
    ('Physics: units', 'physics.units', 'physics.units'),
    ('Physics: vectors', 'physics.vector', 'physics.vector'),
    ('Physics: mechanics', 'physics.mechanics', 'physics.mechanics'),
    ('Physics: quantum mechanics', 'physics.quantum', 'physics.quantum'),
    ('Physics: optics', 'physics.optics', 'physics.optics'),
    ('Physics: Pauli algebra', 'physics.paulialgebra', 'physics.paulialgebra'),
    ('Physics: control', 'physics.control', 'physics.control'),
    ('Physics: continuum beams', 'physics.continuum_mechanics.beam', 'physics.continuum_mechanics.beam'),
    ('Physics: continuum trusses', 'physics.continuum_mechanics.truss', 'physics.continuum_mechanics.truss'),
    ('Physics: continuum cables', 'physics.continuum_mechanics.cable', 'physics.continuum_mechanics.cable'),
    ('Physics: continuum arches', 'physics.continuum_mechanics.arch', 'physics.continuum_mechanics.arch'),
    ('Physics: biomechanics', 'physics.biomechanics', 'physics.biomechanics'),
    ('Physics: high energy', 'physics.hep.gamma_matrices', 'physics.hep.gamma_matrices'),
    ('Physics: Wigner symbols', 'physics.wigner', 'physics.wigner'),
    ('Physics: second quantization', 'physics.secondquant', 'physics.secondquant'),
    ('Statistics', 'stats', 'stats'),
    ('Cryptography', 'crypto.crypto', 'crypto.crypto'),
    ('Holonomic functions', 'holonomic', 'holonomic'),
    ('Differential geometry', 'diffgeom', 'diffgeom'),
    ('Lie algebras', 'liealgebras', 'liealgebras'),
    ('Categories', 'categories', 'categories'),
    ('Discrete transforms and convolutions', 'discrete', 'discrete'),
    ('Coordinate vector calculus', 'vector', 'vector'),
    ('Expression parsing', 'parsing.sympy_parser', 'parsing.sympy_parser'),
)


def exports(module, documented=()):
    if hasattr(module, '__all__'):
        names = list(module.__all__)
    # Do not count imported helpers such as Symbol, exp or simplify as
    # APIs of a submodule with no explicit export list.
    else:
        names = [name for name in dir(module) if not name.startswith('_')
                 and callable(getattr(module, name))
                 and (getattr(getattr(module, name), '__module__', '') or '').startswith(module.__name__)]
    # The official documentation also exposes APIs absent from __all__, such
    # as evaluate_pauli_product. Only count bindings present in this version.
    return sorted(set(names).union(name for name in documented if hasattr(module, name)))


def load_inventory(version):
    inventory = json.loads((ROOT/'tools/sympy_public_modules.json').read_text())
    if inventory['sympy_version'] != version:
        raise ValueError('The documented module inventory targets SymPy '
                         f"{inventory['sympy_version']}; installed version is {version}")
    return inventory


def native_module(qualified):
    try:
        return importlib.import_module(qualified), 'available', None
    except ModuleNotFoundError as error:
        if error.name == qualified or qualified.startswith((error.name or '')+'.'):
            return None, 'missing_module', None
        return None, 'native_import_error', f'{type(error).__name__}: {error}'
    except Exception as error:
        return None, 'native_import_error', f'{type(error).__name__}: {error}'


def documented_surfaces(sx, inventory):
    """Inventory every documented module, parent package and curated surface.

    Import failures remain visible with an explicit status. In particular an
    unavailable reference module must never masquerade as 0% or 100% coverage.
    """
    modules = set(inventory['modules'])
    modules.update('sympy.'+reference for _, reference, _ in MODULES)
    for name in list(modules):
        while '.' in name:
            name = name.rsplit('.', 1)[0]
            modules.add(name)
    results = {}
    for name in sorted(modules):
        qualified = 'symox'+name[len('sympy'):]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                reference = importlib.import_module(name)
        except Exception as error:
            results[name] = dict(reference_module=name, actual_module=qualified,
                                 status='reference_import_error', present=None,
                                 total=None, presence_percent=None, missing=[],
                                 error=f'{type(error).__name__}: {error}')
            continue
        actual, status, error = native_module(qualified)
        names = exports(reference, inventory['documented_exports'].get(name, ()))
        result = surface(reference, actual, names)
        result.update(reference_module=name, actual_module=qualified, status=status)
        if error:
            result['error'] = error
        if not names:
            result['status'] = 'no_exports' if actual is not None else status
        results[name] = result
    return results


def module_surfaces(sx, inventory=None):
    results = {}
    for label, reference_name, actual_name in MODULES:
        reference_path = 'sympy.'+reference_name
        actual_path = 'symox' if actual_name is None else 'symox.'+actual_name
        try:
            reference = importlib.import_module(reference_path)
        except Exception as error:
            results[label+' exports'] = dict(reference_module=reference_path,
                actual_module=actual_path, status='reference_import_error',
                present=None, total=None, presence_percent=None, missing=[],
                error=f'{type(error).__name__}: {error}')
            continue
        actual = sx
        status, error = 'available', None
        if actual_name is not None:
            actual, status, error = native_module(actual_path)
        documented = ()
        if inventory is not None and (reference_name.startswith('physics.') or
                                      reference_name in ('stats', 'crypto.crypto')):
            documented = inventory['documented_exports'].get(reference.__name__, ())
        result = surface(reference, actual, exports(reference, documented))
        result.update(reference_module=reference.__name__,
                      actual_module=actual_path, status=status)
        if error:
            result['error'] = error
        results[label+' exports'] = result
    return results


def percentage(result):
    value = result['presence_percent']
    return 'N/A' if value is None else f'{value:.2f}%'


def surface_family(label):
    name = label.removesuffix(' exports')
    if name.startswith('Physics:'):
        return 'Physics'
    if name.startswith('Special functions:') or name == 'Complex functions and branches':
        return 'Special functions and complex branches'
    groups = {
        'Core and package exports': ('SymPy top-level', 'Expr methods'),
        'Polynomials, algebra and number theory': (
            'Poly methods', 'AlgebraicNumber methods', 'Rational domain methods',
            'Coefficient domains', 'Polynomial ring constructors',
            'Fraction field constructors', 'Number fields', 'Number theory', 'Simplification'),
        'Matrices, arrays and tensors': (
            'Matrix methods', 'MatrixExpr methods', 'Arrays', 'Indexed tensor algebra'),
        'Calculus and equation solvers': (
            'Ordinary differential equations', 'Partial differential equations',
            'Equation solvers', 'Inequalities', 'Series and limits', 'Sequences',
            'Real analysis', 'Singularities and monotonicity', 'Integral transforms',
            'Holonomic functions', 'Coordinate vector calculus'),
        'Geometry and discrete mathematics': (
            'Geometry', 'Differential geometry', 'Combinatorics', 'Symbolic sets',
            'Lie algebras', 'Categories', 'Discrete transforms and convolutions'),
        'Statistics': ('Statistics',),
        'Cryptography': ('Cryptography',),
        'Utilities, parsing, printing and plotting': (
            'Iterable utilities', 'Expression parsing',
            'Printing and code generation', 'Plotting'),
    }
    return next((family for family, names in groups.items() if name in names), 'Other APIs')


def surfaces_markdown(surfaces):
    groups = {}
    for label, result in surfaces.items():
        groups.setdefault(surface_family(label), []).append((label, result))
    families = ('Core and package exports', 'Polynomials, algebra and number theory',
                'Matrices, arrays and tensors', 'Calculus and equation solvers',
                'Geometry and discrete mathematics', 'Special functions and complex branches',
                'Physics', 'Statistics', 'Cryptography',
                'Utilities, parsing, printing and plotting', 'Other APIs')
    labels = ['SymPy top-level exports', 'Expr methods', 'Poly methods',
              'Matrix methods', 'MatrixExpr methods', 'AlgebraicNumber methods',
              'Rational domain methods']+[label+' exports' for label, _, _ in MODULES]
    order = {label: index for index, label in enumerate(labels)}
    lines = ['<!-- sympy-api-surfaces:start -->']
    for family in families:
        entries = groups.get(family)
        if not entries:
            continue
        lines += ['', '<details>', f'<summary>{family}</summary>', '',
                  '| Surface | Present | Missing | API presence |', '|---|---:|---:|---:|']
        for label, result in sorted(entries, key=lambda entry: (order.get(entry[0], len(order)), entry[0])):
            fraction = 'N/A' if result['total'] is None else f"{result['present']} / {result['total']}"
            missing = 'N/A' if result['total'] is None else str(len(result['missing']))
            lines.append(f"| {label} | {fraction} | {missing} | {percentage(result)} |")
        lines += ['', '</details>']
    return '\n'.join(lines+['', '<!-- sympy-api-surfaces:end -->', ''])


def inventory_markdown(report):
    inventory = report['module_inventory']
    modules = report['documented_modules']
    counts = Counter(result['status'] for result in modules.values())
    lines = [
        '# SymPy module API inventory', '', '[README](../README.md)', '',
        f"Measured on {report['measured_at'][:10]} against SymPy {report['sympy_version']}.", '',
        f"Source: [official Sphinx module inventory]({inventory['source_url']}) "
        f"({inventory['documented_count']} documented modules), extended with parent "
        'packages and the README import surfaces.', '',
        f"All {len(modules)} resulting module paths are listed below. "
        f"{counts['missing_module']} native module paths are missing; "
        f"{counts['native_import_error']} native imports failed; "
        f"{counts['reference_import_error']} reference imports could not be measured.", '',
        'Names come from `__all__`, or public callables defined in the module when '
        '`__all__` is absent, plus documented functions, classes and data actually '
        'bound in the installed reference module. Imported helper names are excluded '
        'unless explicitly exported or documented. Constants absent from these '
        'sources and individual class methods are outside this module export measure.', '',
        'The native column measures the corresponding `symox.*` module path. '
        'An API available only at the package root is still absent at its compatible '
        'module path. Existing README rows with an explicit root mapping retain '
        'their original meaning; they can differ from these module rows.', '',
        '`missing_module` means the module or a parent package is absent. '
        '`native_import_error` means an import failed for another reason. '
        '`reference_import_error` means the reference requires an unavailable '
        'dependency or the upstream documentation names an unimportable module. '
        '`no_exports` has no countable names under the stated rule. '
        'N/A is never treated as 100%.', '',
        'Known documentation inventory corrections are recorded in '
        '`tools/sympy_public_modules.json`: `Arch` is a class exported by '
        '`sympy.physics.continuum_mechanics.arch`, not a module named `Arch`.', '',
        'These are name-presence measurements. They do not check signatures, '
        'domains, algorithms or behavior. Modules overlap; rows must not be added '
        'together. This inventory is bounded by the documented modules and listed '
        'parents, not a claim to enumerate every internal SymPy module.', '',
        'Reproduce offline after installing the pinned reference version:', '',
        '```bash', 'python tools/sympy_coverage.py --corpus --json /tmp/symox-coverage.json --update-readme',
        '```', '', 'The JSON includes missing names and import errors for every module.', '',
    ]
    family = None
    for name, result in modules.items():
        current = name.split('.')[1] if '.' in name else 'root'
        if current != family:
            if family is not None:
                lines += ['', '</details>', '']
            family = current
            lines += ['<details>', f'<summary>{family}</summary>', '',
                      '| Reference module | Present | Missing | API presence | Status |',
                      '|---|---:|---:|---:|---|']
        url = inventory['documentation_urls'].get(name)
        label = f'[{name}]({url})' if url else f'`{name}`'
        total, present = result['total'], result['present']
        fraction = 'N/A' if total is None else f'{present} / {total}'
        missing = 'N/A' if total is None else str(len(result['missing']))
        lines.append(f"| {label} | {fraction} | {missing} | {percentage(result)} | {result['status']} |")
    if family is not None:
        lines += ['', '</details>', '']
    errors = [(name, result['error']) for name, result in modules.items() if result.get('error')]
    if errors:
        lines += ['', '## Import errors', '']
        lines += [f'- `{name}`: `{error.replace(chr(10), " ").replace("`", chr(39))}`' for name, error in errors]
    return '\n'.join(lines).rstrip()+'\n'


def update_readme(report):
    path = ROOT/'README.md'
    contents = path.read_text()
    marker = '<!-- sympy-api-surfaces:start -->'
    if marker in contents:
        start = contents.index(marker)
        end_marker = '<!-- sympy-api-surfaces:end -->'
        end = contents.index(end_marker, start)+len(end_marker)
        end += int(contents[end:end+1] == '\n')
    else:
        start = contents.index('| Surface |')
        end = contents.index('\nThese percentages', start)
    path.write_text(contents[:start]+surfaces_markdown(report['surfaces'])+contents[end:])
    (ROOT/'docs/sympy-api.md').write_text(inventory_markdown(report))


def surface(reference, actual, names):
    names = sorted(set(names))
    missing = [name for name in names if not hasattr(actual, name)]
    return {
        "present": len(names) - len(missing),
        "total": len(names),
        "presence_percent": round(100 * (len(names) - len(missing)) / len(names), 2) if names else None,
        "missing": missing,
    }


def methods(value):
    for name in dir(type(value)):
        if not name.startswith('_'):
            member = inspect.getattr_static(value, name)
            if callable(member) or isinstance(member, (classmethod, staticmethod)):
                yield name


def converted(value):
    if hasattr(value, 'to_sympy'):
        return value.to_sympy()
    if isinstance(value, (list, tuple)):
        return type(value)(converted(item) for item in value)
    if isinstance(value, dict):
        return {converted(key): converted(item) for key, item in value.items()}
    return value


def corpus_report(sp, sx):
    suite = runpy.run_path(str(ROOT / 'tests' / 'compare_sympy.py'))
    reference, native = suite['namespace'](sp), suite['namespace'](sx)
    results, details = Counter(), []
    for case in suite['CASES']:
        try:
            expected = eval(case, reference)
        except Exception as error:
            status, detail = 'reference_error', f'{type(error).__name__}: {error}'
        else:
            try:
                actual = eval(case, native)
                back = converted(actual)
                if expected == back:
                    status, detail = 'same_tree', None
                elif isinstance(expected, sp.Expr) and isinstance(back, sp.Expr) and sp.simplify(expected-back) == 0:
                    status, detail = 'proved_equivalent', None
                else:
                    # Do not count numerical agreement or matching text as proof.
                    status, detail = 'unverified_difference', {'symox': str(actual), 'sympy': str(expected)}
            except Exception as error:
                status, detail = 'native_error', f'{type(error).__name__}: {error}'
        results[status] += 1
        if detail is not None:
            details.append({'case': case, 'status': status, 'detail': detail,
                            'known_deviation': case in suite['KNOWN_DEVIATIONS']})
    total = sum(results.values()) - results['reference_error']
    confirmed = results['same_tree'] + results['proved_equivalent']
    return {'total': total, 'confirmed': confirmed,
            'confirmed_percent': round(100*confirmed/total, 2) if total else None,
            'counts': dict(results), 'details': details}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', action='store_true', help='also run exact differential expression checks')
    parser.add_argument('--json', type=Path, help='write the complete report, including missing API names')
    parser.add_argument('--update-readme', action='store_true', help='refresh the README table and full documented module inventory')
    options = parser.parse_args()
    import sympy as sp
    import symox as sx
    inventory = load_inventory(sp.__version__)

    x, X = sx.Symbol('x'), sp.Symbol('x')
    samples = {
        'Expr methods': (X, x),
        'Poly methods': (sp.Poly(X*X+1, X), sx.Poly(x*x+1, x)),
        'Matrix methods': (sp.Matrix([[1, 2], [3, 4]]), sx.Matrix([[1, 2], [3, 4]])),
        'MatrixExpr methods': (sp.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('A', 2, 2)),
        'AlgebraicNumber methods': (sp.AlgebraicNumber(sp.sqrt(2)), sx.AlgebraicNumber(sx.sqrt(2))),
        'Rational domain methods': (sp.QQ, sx.QQ),
    }
    report = {
        'measured_at': datetime.now(timezone.utc).isoformat(),
        'sympy_version': sp.__version__, 'symox_version': sx.__version__,
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'source_dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()),
        'limitations': [
            'API presence does not verify signatures, domains or behavior.',
            'The expression corpus is curated for symox; its percentage is not global SymPy coverage.',
            'There is no combined global parity percentage.',
            'The full module inventory follows documented modules and parent packages, not every internal module.',
        ],
        'surfaces': {'SymPy top-level exports': surface(sp, sx, sp.__all__)},
    }
    for name, (reference, actual) in samples.items():
        report['surfaces'][name] = surface(reference, actual, methods(reference))
    report['surfaces'].update(module_surfaces(sx, inventory))
    report['module_inventory'] = {
        'source_url': inventory['source_url'], 'source_sha256': inventory['source_sha256'],
        'documented_count': len(inventory['modules']),
        'documentation_urls': inventory['documentation_urls'],
    }
    report['documented_modules'] = documented_surfaces(sx, inventory)
    if options.corpus:
        report['corpus'] = corpus_report(sp, sx)
    for name, values in report['surfaces'].items():
        fraction = 'N/A' if values['total'] is None else f"{values['present']}/{values['total']}"
        print(f"{name}: {fraction} names present ({percentage(values)})")
    print('Module inventory:', len(report['documented_modules']), 'paths;',
          json.dumps(dict(Counter(result['status'] for result in report['documented_modules'].values())), sort_keys=True))
    if options.corpus:
        values = report['corpus']
        print(f"Curated corpus: {values['confirmed']}/{values['total']} confirmed ({values['confirmed_percent']}%)")
        print('Corpus counts:', json.dumps(values['counts'], sort_keys=True))
    print('These are separate indicators, not a percentage of global SymPy parity.')
    if options.json:
        options.json.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    if options.update_readme:
        update_readme(report)


if __name__ == '__main__':
    main()
