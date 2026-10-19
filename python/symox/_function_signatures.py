"""SymPy 1.14.0 named-argument contracts for native function heads.

The same Signature objects both expose the interface and bind runtime keyword
calls. Omitted optional arguments stay omitted, including in held expressions.
No reference library is imported by this adapter.
"""
import inspect
from . import _core

_PARAMETERS = {'Abs': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Chi': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Ci': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'DiracDelta': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['k', 'POSITIONAL_OR_KEYWORD', {'default': 0}]],
 'E1': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Ei': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'FallingFactorial': [['x', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                      ['k', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Heaviside': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}],
               ['H0', 'POSITIONAL_OR_KEYWORD', {'rational': [1, 2]}]],
 'LambertW': [['x', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['k', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'Li': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Max': [['args', 'VAR_POSITIONAL', {'required': True}]],
 'Min': [['args', 'VAR_POSITIONAL', {'required': True}]],
 'RisingFactorial': [['x', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                     ['k', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Shi': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Si': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'Ynm': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
         ['m', 'POSITIONAL_OR_KEYWORD', {'required': True}],
         ['theta', 'POSITIONAL_OR_KEYWORD', {'required': True}],
         ['phi', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'acos': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'acosh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'acot': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'airyai': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'airyaiprime': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'airybi': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'airybiprime': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'andre': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'arg': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'asin': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'asinh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'assoc_laguerre': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                    ['alpha', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                    ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'assoc_legendre': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                    ['m', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                    ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'atan': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'atan2': [['y', 'POSITIONAL_OR_KEYWORD', {'required': True}],
           ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'atanh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'bernoulli': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
               ['x', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'besseli': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'besselj': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'besselk': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'bessely': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'beta': [['x', 'POSITIONAL_OR_KEYWORD', {'required': True}],
          ['y', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'binomial': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['k', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'carmichael': [['args', 'VAR_POSITIONAL', {'required': True}]],
 'catalan': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'ceiling': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'chebyshevt': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'chebyshevu': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'conjugate': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'cos': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'cosh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'cot': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'coth': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'csc': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'digamma': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'dirichlet_eta': [['s', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                   ['a', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'divisor_sigma': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                   ['k', 'POSITIONAL_OR_KEYWORD', {'default': 1}]],
 'elliptic_e': [['m', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['z', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'elliptic_f': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['m', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'elliptic_k': [['m', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'elliptic_pi': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                 ['m', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                 ['z', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'erf': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'erf2': [['x', 'POSITIONAL_OR_KEYWORD', {'required': True}],
          ['y', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'erfc': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'erfcinv': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'erfi': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'erfinv': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'euler': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
           ['x', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'exp': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'expint': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
            ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'factorial': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'factorial2': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'fibonacci': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
               ['sym', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'floor': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'fresnelc': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'fresnels': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'gamma': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'gegenbauer': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'genocchi': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['x', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'hankel1': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'hankel2': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'harmonic': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['m', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'hermite': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'im': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'jacobi': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
            ['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
            ['b', 'POSITIONAL_OR_KEYWORD', {'required': True}],
            ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'jacobi_symbol': [['m', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                   ['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'jn': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
        ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'kronecker_symbol': [['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                      ['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'laguerre': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'legendre': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
              ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'legendre_symbol': [['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                     ['p', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'lerchphi': [['args', 'VAR_POSITIONAL', {'required': True}]],
 'li': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'loggamma': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'lowergamma': [['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['x', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'lucas': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'mobius': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'motzkin': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'partition': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'polygamma': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
               ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'polylog': [['s', 'POSITIONAL_OR_KEYWORD', {'required': True}],
             ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'primenu': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'primeomega': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'primepi': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 're': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'reduced_totient': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'sec': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'sign': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'sin': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'sinh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'subfactorial': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'tan': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'tanh': [['arg', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'totient': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'tribonacci': [['n', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['sym', 'POSITIONAL_OR_KEYWORD', {'default': None}]],
 'trigamma': [['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'uppergamma': [['a', 'POSITIONAL_OR_KEYWORD', {'required': True}],
                ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'yn': [['nu', 'POSITIONAL_OR_KEYWORD', {'required': True}],
        ['z', 'POSITIONAL_OR_KEYWORD', {'required': True}]],
 'zeta': [['s', 'POSITIONAL_OR_KEYWORD', {'required': True}],
          ['a', 'POSITIONAL_OR_KEYWORD', {'default': None}]]}


def signature(name):
    fields=_PARAMETERS.get(name)
    if fields is None:return None
    parameters=[]
    for name,kind,options in fields:
        if 'rational' in options:default=_core.Rational(*options['rational'])
        else:default=options.get('default',inspect.Parameter.empty)
        parameters.append(inspect.Parameter(name,getattr(inspect.Parameter,kind),default=default))
    return inspect.Signature(parameters)


def bind(contract,args,keywords):
    bound=contract.bind(*args,**keywords)
    positional=[]
    for parameter in contract.parameters.values():
        if parameter.kind==parameter.VAR_POSITIONAL:
            positional.extend(bound.arguments.get(parameter.name,()))
        elif parameter.kind in (parameter.POSITIONAL_ONLY,parameter.POSITIONAL_OR_KEYWORD):
            positional.append(bound.arguments.get(parameter.name,parameter.default))
        else:
            raise TypeError('unsupported native function argument contract')
    # Native optional-None conventions mean the omitted function overload.
    while positional:
        parameter=list(contract.parameters.values())[len(positional)-1] if len(positional)<=len(contract.parameters) else None
        if parameter is not None and parameter.default is None and positional[-1] is None:
            positional.pop()
        elif parameter is not None and parameter.name not in bound.arguments and parameter.default is not inspect.Parameter.empty:
            positional.pop()
        else:break
    return tuple(positional)
