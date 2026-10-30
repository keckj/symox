"""Differential tests of the native finite Markov chain kernels
(symox::stats::markov) against SymPy 1.14.

The expected values are SymPy 1.14 results (srepr) over the Markov cases of
sympy/stats/tests/test_stochastic_process.py and 20 random chains (sizes 2-8,
rational entries and one-symbol entries). Probabilities conditioned on
several observation times are checked by the Markov property instead: SymPy
1.14 ignores the later conditions there (see docs/algorithmic-parity.md,
"Wave 21 spec7" and "Wave 33 markovrust")."""
import ast
import warnings

import pytest

import symox as S
from symox.stats import DiscreteMarkovChain, ContinuousMarkovChain, P, E
from symox.stats import _markov_native

NS = {name: getattr(S, name) for name in dir(S) if not name.startswith('_')}

CASES = {'TO2': {'E2': 'Integer(0)',
         'P2': 'Integer(0)',
         'Pne': 'Integer(0)',
         'abschain': 'true',
         'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)]])',
         'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Integer(1), "
                      'Integer(0), Integer(0)], [Rational(1, 2), Rational(1, 4), Rational(1, 4)], '
                      "[Rational(1, 3), Rational(1, 3), Rational(1, 3)]])']",
         'classes': "[[['Integer(0)'], 'true', 'Integer(1)'], [['Integer(1)', 'Integer(2)'], 'false', "
                    "'Integer(1)']]",
         'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                      "'ImmutableDenseMatrix([[Rational(1, 2)], [Rational(1, 3)]])', "
                      "'ImmutableDenseMatrix([[Rational(1, 4), Rational(1, 4)], [Rational(1, 3), Rational(1, "
                      "3)]])']",
         'ergodic': 'false',
         'fundamental': 'ImmutableDenseMatrix([[Rational(8, 5), Rational(3, 5)], [Rational(4, 5), '
                        'Rational(9, 5)]])',
         'regular': 'false',
         'rows': [['1', '0', '0'], ['1/2', '1/4', '1/4'], ['1/3', '1/3', '1/3']],
         'stationary': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)]])'},
 'TO3': {'E2': 'Rational(15, 16)',
         'P2': 'Rational(1, 4)',
         'Pne': 'Rational(11, 16)',
         'abschain': 'false',
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Rational(1, 4), "
                      'Rational(3, 4), Integer(0)], [Rational(1, 3), Rational(1, 3), Rational(1, 3)], '
                      "[Integer(0), Rational(1, 4), Rational(3, 4)]])']",
         'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'true', 'Integer(1)']]",
         'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Rational(1, 4), "
                      'Rational(3, 4), Integer(0)], [Rational(1, 3), Rational(1, 3), Rational(1, 3)], '
                      "[Integer(0), Rational(1, 4), Rational(3, 4)]])', 'ImmutableDenseMatrix(0, 3, [])', "
                      "'ImmutableDenseMatrix([])']",
         'ergodic': 'true',
         'fundamental': 'ImmutableDenseMatrix([[Rational(176, 125), Rational(81, 125), Rational(-132, 125)], '
                        '[Rational(36, 125), Rational(141, 125), Rational(-52, 125)], [Rational(-44, 125), '
                        'Rational(-39, 125), Rational(208, 125)]])',
         'regular': 'true',
         'rows': [['1/4', '3/4', '0'], ['1/3', '1/3', '1/3'], ['0', '1/4', '3/4']],
         'stationary': 'ImmutableDenseMatrix([[Rational(4, 25), Rational(9, 25), Rational(12, 25)]])'},
 'TO6': {'E2': 'Integer(0)',
         'P2': 'Integer(0)',
         'Pne': 'Integer(0)',
         'abschain': 'true',
         'absorbing': 'ImmutableDenseMatrix([[Rational(3, 4), Rational(1, 4)], [Rational(1, 2), Rational(1, '
                      '2)], [Rational(1, 4), Rational(3, 4)]])',
         'canonical': "[['Integer(0)', 'Integer(4)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                      "'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], "
                      '[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], [Rational(1, 2), '
                      'Integer(0), Integer(0), Rational(1, 2), Integer(0)], [Integer(0), Integer(0), '
                      'Rational(1, 2), Integer(0), Rational(1, 2)], [Integer(0), Rational(1, 2), Integer(0), '
                      "Rational(1, 2), Integer(0)]])']",
         'classes': "[[['Integer(0)'], 'true', 'Integer(1)'], [['Integer(4)'], 'true', 'Integer(1)'], "
                    "[['Integer(1)', 'Integer(2)', 'Integer(3)'], 'false', 'Integer(2)']]",
         'decompose': "[['Integer(0)', 'Integer(4)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                      "'ImmutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(1)]])', "
                      "'ImmutableDenseMatrix([[Rational(1, 2), Integer(0)], [Integer(0), Integer(0)], "
                      "[Integer(0), Rational(1, 2)]])', 'ImmutableDenseMatrix([[Integer(0), Rational(1, 2), "
                      'Integer(0)], [Rational(1, 2), Integer(0), Rational(1, 2)], [Integer(0), Rational(1, '
                      "2), Integer(0)]])']",
         'ergodic': 'false',
         'fundamental': 'ImmutableDenseMatrix([[Rational(3, 2), Integer(1), Rational(1, 2)], [Integer(1), '
                        'Integer(2), Integer(1)], [Rational(1, 2), Integer(1), Rational(3, 2)]])',
         'regular': 'false',
         'rows': [['1', '0', '0', '0', '0'],
                  ['1/2', '0', '1/2', '0', '0'],
                  ['0', '1/2', '0', '1/2', '0'],
                  ['0', '0', '1/2', '0', '1/2'],
                  ['0', '0', '0', '0', '1']],
         'stationary': "ImmutableDenseMatrix([[Add(Integer(1), Mul(Integer(-1), Symbol('tau0'))), "
                       "Integer(0), Integer(0), Integer(0), Symbol('tau0')]])"},
 'TO7': {'E2': 'Rational(15, 16)',
         'P2': 'Rational(3, 8)',
         'Pne': 'Rational(9, 16)',
         'abschain': 'false',
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Rational(1, 2), "
                      'Rational(1, 4), Rational(1, 4)], [Rational(1, 2), Integer(0), Rational(1, 2)], '
                      "[Rational(1, 4), Rational(1, 4), Rational(1, 2)]])']",
         'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'true', 'Integer(1)']]",
         'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'ImmutableDenseMatrix([[Rational(1, 2), "
                      'Rational(1, 4), Rational(1, 4)], [Rational(1, 2), Integer(0), Rational(1, 2)], '
                      "[Rational(1, 4), Rational(1, 4), Rational(1, 2)]])', 'ImmutableDenseMatrix(0, 3, "
                      "[])', 'ImmutableDenseMatrix([])']",
         'ergodic': 'true',
         'fundamental': 'ImmutableDenseMatrix([[Rational(86, 75), Rational(1, 25), Rational(-14, 75)], '
                        '[Rational(2, 25), Rational(21, 25), Rational(2, 25)], [Rational(-14, 75), '
                        'Rational(1, 25), Rational(86, 75)]])',
         'regular': 'true',
         'rows': [['1/2', '1/4', '1/4'], ['1/2', '0', '1/2'], ['1/4', '1/4', '1/2']],
         'stationary': 'ImmutableDenseMatrix([[Rational(2, 5), Rational(1, 5), Rational(2, 5)]])'},
 'Y7': {'E2': 'Integer(2)',
        'P2': 'Integer(0)',
        'Pne': 'Integer(1)',
        'abschain': 'false',
        'absorbing': "ERR:ValueError:The fundamental matrix doesn't exist.",
        'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)', 'Integer(4)', 'Integer(5)', "
                     "'Integer(6)'], 'ImmutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0), "
                     'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(1), Integer(0), '
                     'Integer(0), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(1), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(1), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(0), Integer(1)], [Integer(0), Integer(0), Integer(0), Integer(1), '
                     "Integer(0), Integer(0), Integer(0)]])']",
        'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)'], 'true', 'Integer(3)'], [['Integer(3)', "
                   "'Integer(4)', 'Integer(5)', 'Integer(6)'], 'true', 'Integer(4)']]",
        'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)', 'Integer(4)', 'Integer(5)', "
                     "'Integer(6)'], 'ImmutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0), "
                     'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(1), Integer(0), '
                     'Integer(0), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(1), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(1), Integer(0)], [Integer(0), Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(0), Integer(1)], [Integer(0), Integer(0), Integer(0), Integer(1), '
                     "Integer(0), Integer(0), Integer(0)]])', 'ImmutableDenseMatrix(0, 7, [])', "
                     "'ImmutableDenseMatrix([])']",
        'ergodic': 'false',
        'fundamental': "ERR:ValueError:The fundamental matrix doesn't exist.",
        'regular': 'false',
        'rows': [['0', '1', '0', '0', '0', '0', '0'],
                 ['0', '0', '1', '0', '0', '0', '0'],
                 ['1', '0', '0', '0', '0', '0', '0'],
                 ['0', '0', '0', '0', '1', '0', '0'],
                 ['0', '0', '0', '0', '0', '1', '0'],
                 ['0', '0', '0', '0', '0', '0', '1'],
                 ['0', '0', '0', '1', '0', '0', '0']],
        'stationary': 'ImmutableDenseMatrix([[Add(Rational(1, 3), Mul(Integer(-1), Rational(4, 3), '
                      "Symbol('tau0'))), Add(Rational(1, 3), Mul(Integer(-1), Rational(4, 3), "
                      "Symbol('tau0'))), Add(Rational(1, 3), Mul(Integer(-1), Rational(4, 3), "
                      "Symbol('tau0'))), Symbol('tau0'), Symbol('tau0'), Symbol('tau0'), Symbol('tau0')]])"},
 'Y9': {'E2': 'Integer(1)',
        'P2': 'Integer(0)',
        'Pne': 'Rational(1, 2)',
        'abschain': 'false',
        'absorbing': 'ImmutableDenseMatrix([[Integer(1), Integer(0)], [Integer(1), Integer(0)], [Integer(1), '
                     'Integer(0)]])',
        'canonical': "[['Integer(3)', 'Integer(4)', 'Integer(2)', 'Integer(0)', 'Integer(1)'], "
                     "'ImmutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], "
                     '[Rational(1, 2), Rational(1, 2), Integer(0), Integer(0), Integer(0)], [Integer(1), '
                     'Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(0), '
                     'Integer(0), Integer(1)], [Integer(0), Integer(0), Rational(1, 2), Rational(1, 2), '
                     "Integer(0)]])']",
        'classes': "[[['Integer(3)', 'Integer(4)'], 'true', 'Integer(1)'], [['Integer(2)'], 'false', "
                   "'Integer(1)'], [['Integer(0)', 'Integer(1)'], 'false', 'Integer(2)']]",
        'decompose': "[['Integer(3)', 'Integer(4)', 'Integer(2)', 'Integer(0)', 'Integer(1)'], "
                     "'ImmutableDenseMatrix([[Integer(0), Integer(1)], [Rational(1, 2), Rational(1, 2)]])', "
                     "'ImmutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(0)], "
                     "[Integer(0), Integer(0)]])', 'ImmutableDenseMatrix([[Integer(0), Integer(0), "
                     'Integer(0)], [Integer(0), Integer(0), Integer(1)], [Rational(1, 2), Rational(1, 2), '
                     "Integer(0)]])']",
        'ergodic': 'false',
        'fundamental': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)], [Integer(1), Integer(2), '
                       'Integer(2)], [Integer(1), Integer(1), Integer(2)]])',
        'regular': 'false',
        'rows': [['0', '1', '0', '0', '0'],
                 ['1/2', '0', '1/2', '0', '0'],
                 ['0', '0', '0', '1', '0'],
                 ['0', '0', '0', '0', '1'],
                 ['0', '0', '0', '1/2', '1/2']],
        'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(0), Rational(1, 3), '
                      'Rational(2, 3)]])'},
 'float': {'E2': "Float('0.75', precision=53)",
           'P2': "Float('0.75', precision=53)",
           'Pne': "Float('0.75', precision=53)",
           'abschain': 'false',
           'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
           'canonical': '[[\'Integer(0)\', \'Integer(1)\'], "ImmutableDenseMatrix([[Float(\'0.5\', '
                        "precision=53), Float('0.5', precision=53)], [Float('0.0', precision=53), "
                        'Float(\'1.0\', precision=53)]])"]',
           'classes': "[[['Integer(0)', 'Integer(1)'], 'true', 'Integer(1)']]",
           'decompose': '[[\'Integer(0)\', \'Integer(1)\'], "ImmutableDenseMatrix([[Float(\'0.5\', '
                        "precision=53), Float('0.5', precision=53)], [Float('0.0', precision=53), "
                        'Float(\'1.0\', precision=53)]])", \'ImmutableDenseMatrix(0, 2, [])\', '
                        "'ImmutableDenseMatrix([])']",
           'ergodic': 'true',
           'fundamental': "ImmutableDenseMatrix([[Float('2.0', precision=53), Float('-1.0', precision=53)], "
                          "[Integer(0), Float('1.0', precision=53)]])",
           'regular': 'true',
           'rows': [['0.5', '0.5'], ['0.0', '1.0']],
           'stationary': "ImmutableDenseMatrix([[Integer(0), Float('1.0', precision=53)]])"},
 'per': {'E2': 'Integer(2)',
         'P2': 'Integer(0)',
         'Pne': 'Integer(1)',
         'abschain': 'false',
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
         'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                      "'ImmutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(0), "
                      'Integer(0), Integer(1), Integer(0)], [Integer(0), Integer(0), Integer(0), '
                      "Integer(1)], [Rational(1, 2), Integer(0), Rational(1, 2), Integer(0)]])']",
         'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], 'true', 'Integer(2)']]",
         'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                      "'ImmutableDenseMatrix([[Integer(0), Integer(1), Integer(0), Integer(0)], [Integer(0), "
                      'Integer(0), Integer(1), Integer(0)], [Integer(0), Integer(0), Integer(0), '
                      "Integer(1)], [Rational(1, 2), Integer(0), Rational(1, 2), Integer(0)]])', "
                      "'ImmutableDenseMatrix(0, 4, [])', 'ImmutableDenseMatrix([])']",
         'ergodic': 'true',
         'fundamental': 'ImmutableDenseMatrix([[Rational(25, 36), Rational(19, 36), Rational(1, 18), '
                        'Rational(-5, 18)], [Rational(-5, 36), Rational(25, 36), Rational(7, 18), '
                        'Rational(1, 18)], [Rational(1, 36), Rational(-5, 36), Rational(13, 18), Rational(7, '
                        '18)], [Rational(7, 36), Rational(1, 36), Rational(1, 18), Rational(13, 18)]])',
         'regular': 'false',
         'rows': [['0', '1', '0', '0'], ['0', '0', '1', '0'], ['0', '0', '0', '1'], ['1/2', '0', '1/2', '0']],
         'stationary': 'ImmutableDenseMatrix([[Rational(1, 6), Rational(1, 6), Rational(1, 3), Rational(1, '
                       '3)]])'},
 'rand0': {'E2': 'Rational(26, 49)',
           'P2': 'Rational(26, 49)',
           'Pne': 'Rational(26, 49)',
           'abschain': 'false',
           'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
           'canonical': "[['Integer(0)', 'Integer(1)'], 'ImmutableDenseMatrix([[Rational(3, 7), Rational(4, "
                        "7)], [Rational(1, 2), Rational(1, 2)]])']",
           'classes': "[[['Integer(0)', 'Integer(1)'], 'true', 'Integer(1)']]",
           'decompose': "[['Integer(0)', 'Integer(1)'], 'ImmutableDenseMatrix([[Rational(3, 7), Rational(4, "
                        "7)], [Rational(1, 2), Rational(1, 2)]])', 'ImmutableDenseMatrix(0, 2, [])', "
                        "'ImmutableDenseMatrix([])']",
           'ergodic': 'true',
           'fundamental': 'ImmutableDenseMatrix([[Rational(217, 225), Rational(8, 225)], [Rational(7, 225), '
                          'Rational(218, 225)]])',
           'regular': 'true',
           'rows': [['3/7', '4/7'], ['1/2', '1/2']],
           'stationary': 'ImmutableDenseMatrix([[Rational(7, 15), Rational(8, 15)]])'},
 'rand1': {'E2': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
           'P2': 'Integer(0)',
           'Pne': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
           'abschain': 'true',
           'absorbing': 'ImmutableDenseMatrix([[Mul(Integer(-1), Add(Integer(1), Mul(Integer(-1), '
                        "Symbol('p'))), Pow(Add(Symbol('p'), Integer(-1)), Integer(-1))), Integer(0)]])",
           'canonical': "[['Integer(1)', 'Integer(2)', 'Integer(0)'], "
                        '"ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)], [Integer(0), '
                        "Integer(1), Integer(0)], [Add(Integer(1), Mul(Integer(-1), Symbol('p'))), "
                        'Integer(0), Symbol(\'p\')]])"]',
           'classes': "[[['Integer(1)'], 'true', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)'], "
                      "[['Integer(2)'], 'true', 'Integer(1)']]",
           'decompose': "[['Integer(1)', 'Integer(2)', 'Integer(0)'], 'ImmutableDenseMatrix([[Integer(1), "
                        'Integer(0)], [Integer(0), Integer(1)]])\', "ImmutableDenseMatrix([[Add(Integer(1), '
                        'Mul(Integer(-1), Symbol(\'p\'))), Integer(0)]])", '
                        '"ImmutableDenseMatrix([[Symbol(\'p\')]])"]',
           'ergodic': 'false',
           'fundamental': "ImmutableDenseMatrix([[Mul(Integer(-1), Pow(Add(Symbol('p'), Integer(-1)), "
                          'Integer(-1)))]])',
           'regular': 'false',
           'rows': [['p', '1 - p', '0'], ['0', '1', '0'], ['0', '0', '1']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Add(Integer(1), Mul(Integer(-1), '
                         "Symbol('tau0'))), Symbol('tau0')]])"},
 'rand10': {'E2': 'Integer(3)',
            'P2': 'Integer(0)',
            'Pne': 'Integer(1)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)], [Integer(1)], [Integer(1)]])',
            'canonical': "[['Integer(3)', 'Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(4)'], "
                         "'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), "
                         'Integer(0)], [Integer(1), Integer(0), Integer(0), Integer(0), Integer(0)], '
                         '[Integer(0), Rational(1, 3), Integer(0), Rational(2, 3), Integer(0)], [Rational(1, '
                         '3), Integer(0), Rational(4, 9), Rational(2, 9), Integer(0)], [Rational(1, 7), '
                         "Rational(4, 7), Rational(2, 7), Integer(0), Integer(0)]])']",
            'classes': "[[['Integer(3)'], 'true', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)'], "
                       "[['Integer(1)', 'Integer(2)'], 'false', 'Integer(1)'], [['Integer(4)'], 'false', "
                       "'Integer(1)']]",
            'decompose': "[['Integer(3)', 'Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(4)'], "
                         "'ImmutableDenseMatrix([[Integer(1)]])', 'ImmutableDenseMatrix([[Integer(1)], "
                         "[Integer(0)], [Rational(1, 3)], [Rational(1, 7)]])', "
                         "'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(0), Integer(0)], "
                         '[Rational(1, 3), Integer(0), Rational(2, 3), Integer(0)], [Integer(0), Rational(4, '
                         '9), Rational(2, 9), Integer(0)], [Rational(4, 7), Rational(2, 7), Integer(0), '
                         "Integer(0)]])']",
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], '
                           '[Rational(7, 13), Rational(21, 13), Rational(18, 13), Integer(0)], [Rational(4, '
                           '13), Rational(12, 13), Rational(27, 13), Integer(0)], [Rational(66, 91), '
                           'Rational(6, 13), Rational(36, 91), Integer(1)]])',
            'regular': 'false',
            'rows': [['0', '0', '0', '1', '0'],
                     ['1/3', '0', '2/3', '0', '0'],
                     ['0', '4/9', '2/9', '1/3', '0'],
                     ['0', '0', '0', '1', '0'],
                     ['4/7', '2/7', '0', '1/7', '0']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(0), Integer(1), '
                          'Integer(0)]])'},
 'rand11': {'E2': "Add(Rational(82, 33), Mul(Integer(-1), Rational(7, 108), Symbol('p')))",
            'P2': 'Rational(79, 594)',
            'Pne': 'Rational(77, 81)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Mul(Integer(-1), Integer(12), Add(Rational(7, 12), '
                         "Mul(Integer(-1), Rational(7, 12), Symbol('p'))), Pow(Add(Mul(Integer(7), "
                         "Symbol('p')), Integer(-7)), Integer(-1)))], [Mul(Integer(-1), Integer(12), "
                         "Add(Rational(7, 12), Mul(Integer(-1), Rational(7, 12), Symbol('p'))), "
                         "Pow(Add(Mul(Integer(7), Symbol('p')), Integer(-7)), Integer(-1)))], "
                         '[Mul(Integer(-1), Integer(12), Add(Rational(7, 12), Mul(Integer(-1), Rational(7, '
                         "12), Symbol('p'))), Pow(Add(Mul(Integer(7), Symbol('p')), Integer(-7)), "
                         'Integer(-1)))], [Mul(Integer(-1), Integer(12), Add(Rational(7, 12), '
                         "Mul(Integer(-1), Rational(7, 12), Symbol('p'))), Pow(Add(Mul(Integer(7), "
                         "Symbol('p')), Integer(-7)), Integer(-1)))], [Add(Mul(Integer(-1), Integer(48), "
                         "Add(Rational(7, 12), Mul(Integer(-1), Rational(7, 12), Symbol('p'))), "
                         "Pow(Add(Mul(Integer(49), Symbol('p')), Integer(-49)), Integer(-1))), Rational(3, "
                         '7))]])',
            'canonical': "[['Integer(2)', 'Integer(1)', 'Integer(3)', 'Integer(5)', 'Integer(4)', "
                         '\'Integer(0)\'], "ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), '
                         'Integer(0), Integer(0), Integer(0)], [Integer(0), Rational(3, 11), Rational(2, '
                         '11), Rational(2, 11), Rational(4, 11), Integer(0)], [Integer(0), Rational(1, 4), '
                         'Rational(1, 4), Rational(1, 2), Integer(0), Integer(0)], [Integer(0), Rational(1, '
                         '2), Rational(1, 2), Integer(0), Integer(0), Integer(0)], [Add(Rational(7, 12), '
                         "Mul(Integer(-1), Rational(7, 12), Symbol('p'))), Mul(Rational(7, 12), "
                         "Symbol('p')), Integer(0), Rational(1, 3), Rational(1, 12), Integer(0)], "
                         '[Rational(1, 3), Rational(2, 9), Rational(1, 9), Integer(0), Rational(1, 9), '
                         'Rational(2, 9)]])"]',
            'classes': "[[['Integer(2)'], 'true', 'Integer(1)'], [['Integer(1)', 'Integer(3)', 'Integer(5)', "
                       "'Integer(4)'], 'false', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(2)', 'Integer(1)', 'Integer(3)', 'Integer(5)', 'Integer(4)', "
                         "'Integer(0)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                         '"ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)], [Add(Rational(7, '
                         '12), Mul(Integer(-1), Rational(7, 12), Symbol(\'p\')))], [Rational(1, 3)]])", '
                         '"ImmutableDenseMatrix([[Rational(3, 11), Rational(2, 11), Rational(2, 11), '
                         'Rational(4, 11), Integer(0)], [Rational(1, 4), Rational(1, 4), Rational(1, 2), '
                         'Integer(0), Integer(0)], [Rational(1, 2), Rational(1, 2), Integer(0), Integer(0), '
                         "Integer(0)], [Mul(Rational(7, 12), Symbol('p')), Integer(0), Rational(1, 3), "
                         'Rational(1, 12), Integer(0)], [Rational(2, 9), Rational(1, 9), Integer(0), '
                         'Rational(1, 9), Rational(2, 9)]])"]',
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Mul(Integer(-1), Integer(121), Pow(Add(Mul(Integer(28), '
                           "Symbol('p')), Integer(-28)), Integer(-1))), Mul(Integer(-1), Integer(41), "
                           "Pow(Add(Mul(Integer(14), Symbol('p')), Integer(-14)), Integer(-1))), "
                           "Mul(Integer(-1), Integer(79), Pow(Add(Mul(Integer(28), Symbol('p')), "
                           'Integer(-28)), Integer(-1))), Mul(Integer(-1), Integer(12), '
                           "Pow(Add(Mul(Integer(7), Symbol('p')), Integer(-7)), Integer(-1))), Integer(0)], "
                           "[Mul(Integer(-1), Integer(121), Pow(Add(Mul(Integer(28), Symbol('p')), "
                           "Integer(-28)), Integer(-1))), Mul(Pow(Add(Mul(Integer(14), Symbol('p')), "
                           "Integer(-14)), Integer(-1)), Add(Mul(Integer(28), Symbol('p')), Integer(-69))), "
                           "Mul(Add(Mul(Integer(28), Symbol('p')), Integer(-107)), Pow(Add(Mul(Integer(28), "
                           "Symbol('p')), Integer(-28)), Integer(-1))), Mul(Integer(-1), Integer(12), "
                           "Pow(Add(Mul(Integer(7), Symbol('p')), Integer(-7)), Integer(-1))), Integer(0)], "
                           "[Mul(Integer(-1), Integer(121), Pow(Add(Mul(Integer(28), Symbol('p')), "
                           "Integer(-28)), Integer(-1))), Mul(Add(Mul(Integer(14), Symbol('p')), "
                           "Integer(-55)), Pow(Add(Mul(Integer(14), Symbol('p')), Integer(-14)), "
                           "Integer(-1))), Mul(Pow(Add(Mul(Integer(28), Symbol('p')), Integer(-28)), "
                           "Integer(-1)), Add(Mul(Integer(42), Symbol('p')), Integer(-121))), "
                           "Mul(Integer(-1), Integer(12), Pow(Add(Mul(Integer(7), Symbol('p')), "
                           'Integer(-7)), Integer(-1))), Integer(0)], [Mul(Add(Mul(Integer(-1), Integer(77), '
                           "Symbol('p')), Integer(-44)), Pow(Add(Mul(Integer(28), Symbol('p')), "
                           'Integer(-28)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(21), '
                           "Symbol('p')), Integer(-20)), Pow(Add(Mul(Integer(14), Symbol('p')), "
                           'Integer(-14)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(35), '
                           "Symbol('p')), Integer(-44)), Pow(Add(Mul(Integer(28), Symbol('p')), "
                           'Integer(-28)), Integer(-1))), Mul(Integer(-1), Integer(12), '
                           "Pow(Add(Mul(Integer(7), Symbol('p')), Integer(-7)), Integer(-1))), Integer(0)], "
                           "[Mul(Add(Mul(Integer(-1), Integer(77), Symbol('p')), Integer(-407)), "
                           "Pow(Add(Mul(Integer(196), Symbol('p')), Integer(-196)), Integer(-1))), "
                           "Mul(Add(Mul(Integer(7), Symbol('p')), Integer(-171)), Pow(Add(Mul(Integer(98), "
                           "Symbol('p')), Integer(-98)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(7), "
                           "Symbol('p')), Integer(-309)), Pow(Add(Mul(Integer(196), Symbol('p')), "
                           'Integer(-196)), Integer(-1))), Mul(Integer(-1), Integer(48), '
                           "Pow(Add(Mul(Integer(49), Symbol('p')), Integer(-49)), Integer(-1))), Rational(9, "
                           '7)]])',
            'regular': 'false',
            'rows': [['2/9', '2/9', '1/3', '1/9', '1/9', '0'],
                     ['0', '3/11', '0', '2/11', '4/11', '2/11'],
                     ['0', '0', '1', '0', '0', '0'],
                     ['0', '1/4', '0', '1/4', '0', '1/2'],
                     ['0', '7*p/12', '7/12 - 7*p/12', '0', '1/12', '1/3'],
                     ['0', '1/2', '0', '1/2', '0', '0']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(1), Integer(0), '
                          'Integer(0), Integer(0)]])'},
 'rand12': {'E2': 'Integer(0)',
            'P2': 'Integer(0)',
            'Pne': 'Integer(0)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)], [Integer(1)], [Integer(1)], '
                         '[Integer(1)], [Integer(1)]])',
            'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(3)', 'Integer(4)', 'Integer(2)', "
                         "'Integer(6)', 'Integer(5)'], 'ImmutableDenseMatrix([[Integer(1), Integer(0), "
                         'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
                         'Rational(1, 8), Rational(1, 2), Integer(0), Integer(0), Rational(1, 4), '
                         'Rational(1, 8)], [Integer(0), Rational(2, 7), Rational(2, 7), Rational(1, 7), '
                         'Integer(0), Rational(3, 14), Rational(1, 14)], [Integer(0), Rational(1, 7), '
                         'Rational(3, 7), Integer(0), Rational(2, 7), Rational(1, 7), Integer(0)], '
                         '[Rational(2, 3), Rational(1, 3), Integer(0), Integer(0), Integer(0), Integer(0), '
                         'Integer(0)], [Integer(0), Rational(3, 13), Rational(4, 13), Rational(3, 13), '
                         'Integer(0), Rational(2, 13), Rational(1, 13)], [Rational(1, 2), Rational(1, 8), '
                         "Rational(1, 4), Integer(0), Integer(0), Integer(0), Rational(1, 8)]])']",
            'classes': "[[['Integer(0)'], 'true', 'Integer(1)'], [['Integer(1)', 'Integer(3)', 'Integer(4)', "
                       "'Integer(2)', 'Integer(6)', 'Integer(5)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(3)', 'Integer(4)', 'Integer(2)', "
                         "'Integer(6)', 'Integer(5)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                         "'ImmutableDenseMatrix([[Integer(0)], [Integer(0)], [Integer(0)], [Rational(2, 3)], "
                         "[Integer(0)], [Rational(1, 2)]])', 'ImmutableDenseMatrix([[Rational(1, 8), "
                         'Rational(1, 2), Integer(0), Integer(0), Rational(1, 4), Rational(1, 8)], '
                         '[Rational(2, 7), Rational(2, 7), Rational(1, 7), Integer(0), Rational(3, 14), '
                         'Rational(1, 14)], [Rational(1, 7), Rational(3, 7), Integer(0), Rational(2, 7), '
                         'Rational(1, 7), Integer(0)], [Rational(1, 3), Integer(0), Integer(0), Integer(0), '
                         'Integer(0), Integer(0)], [Rational(3, 13), Rational(4, 13), Rational(3, 13), '
                         'Integer(0), Rational(2, 13), Rational(1, 13)], [Rational(1, 8), Rational(1, 4), '
                         "Integer(0), Integer(0), Integer(0), Rational(1, 8)]])']",
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Rational(7382, 1631), Rational(1383, 233), Rational(729, '
                           '466), Rational(729, 1631), Rational(10127, 3262), Rational(2290, 1631)], '
                           '[Rational(18082, 4893), Rational(1587, 233), Rational(787, 466), Rational(787, '
                           '1631), Rational(30355, 9786), Rational(6638, 4893)], [Rational(15538, 4893), '
                           'Rational(1201, 233), Rational(1093, 466), Rational(1093, 1631), Rational(25831, '
                           '9786), Rational(5414, 4893)], [Rational(7382, 4893), Rational(461, 233), '
                           'Rational(243, 466), Rational(1874, 1631), Rational(10127, 9786), Rational(2290, '
                           '4893)], [Rational(5870, 1631), Rational(1341, 233), Rational(813, 466), '
                           'Rational(813, 1631), Rational(13039, 3262), Rational(2178, 1631)], '
                           '[Rational(1190, 699), Rational(651, 233), Rational(329, 466), Rational(47, 233), '
                           'Rational(1859, 1398), Rational(1210, 699)]])',
            'regular': 'false',
            'rows': [['1', '0', '0', '0', '0', '0', '0'],
                     ['0', '1/8', '0', '1/2', '0', '1/8', '1/4'],
                     ['2/3', '1/3', '0', '0', '0', '0', '0'],
                     ['0', '2/7', '0', '2/7', '1/7', '1/14', '3/14'],
                     ['0', '1/7', '2/7', '3/7', '0', '0', '1/7'],
                     ['1/2', '1/8', '0', '1/4', '0', '1/8', '0'],
                     ['0', '3/13', '0', '4/13', '3/13', '1/13', '2/13']],
            'stationary': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), '
                          'Integer(0), Integer(0), Integer(0)]])'},
 'rand13': {'E2': "Add(Mul(Integer(-1), Rational(9, 49), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), "
                  "Rational(600, 833), Symbol('p')), Rational(1129, 357))",
            'P2': "Add(Rational(9, 119), Mul(Integer(-1), Rational(9, 119), Symbol('p')))",
            'Pne': "Add(Integer(1), Mul(Integer(-1), Rational(9, 49), Pow(Symbol('p'), Integer(2))))",
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Rational(189, 391), Rational(202, 391)], [Rational(189, '
                         '391), Rational(202, 391)], [Rational(189, 391), Rational(202, 391)], '
                         '[Rational(201, 391), Rational(190, 391)], [Rational(189, 391), Rational(202, '
                         '391)], [Rational(105, 391), Rational(286, 391)]])',
            'canonical': "[['Integer(2)', 'Integer(6)', 'Integer(0)', 'Integer(1)', 'Integer(4)', "
                         "'Integer(5)', 'Integer(3)', 'Integer(7)'], "
                         '"ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), '
                         'Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(1), '
                         'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
                         "[Integer(0), Integer(0), Mul(Rational(3, 7), Symbol('p')), Add(Rational(3, 7), "
                         "Mul(Integer(-1), Rational(3, 7), Symbol('p'))), Integer(0), Integer(0), "
                         'Rational(4, 7), Integer(0)], [Rational(3, 17), Rational(2, 17), Integer(0), '
                         'Rational(4, 17), Rational(3, 17), Rational(2, 17), Integer(0), Rational(3, 17)], '
                         '[Integer(0), Integer(0), Rational(2, 3), Integer(0), Rational(1, 3), Integer(0), '
                         'Integer(0), Integer(0)], [Rational(1, 4), Rational(1, 6), Rational(1, 12), '
                         'Integer(0), Rational(1, 12), Integer(0), Rational(1, 3), Rational(1, 12)], '
                         '[Integer(0), Integer(0), Integer(0), Rational(1, 3), Rational(1, 3), Integer(0), '
                         'Rational(1, 3), Integer(0)], [Integer(0), Rational(4, 9), Rational(1, 3), '
                         'Integer(0), Integer(0), Integer(0), Rational(2, 9), Integer(0)]])"]',
            'classes': "[[['Integer(2)'], 'true', 'Integer(1)'], [['Integer(6)'], 'true', 'Integer(1)'], "
                       "[['Integer(0)', 'Integer(1)', 'Integer(4)', 'Integer(5)', 'Integer(3)', "
                       "'Integer(7)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(2)', 'Integer(6)', 'Integer(0)', 'Integer(1)', 'Integer(4)', "
                         "'Integer(5)', 'Integer(3)', 'Integer(7)'], 'ImmutableDenseMatrix([[Integer(1), "
                         "Integer(0)], [Integer(0), Integer(1)]])', 'ImmutableDenseMatrix([[Integer(0), "
                         'Integer(0)], [Rational(3, 17), Rational(2, 17)], [Integer(0), Integer(0)], '
                         '[Rational(1, 4), Rational(1, 6)], [Integer(0), Integer(0)], [Integer(0), '
                         'Rational(4, 9)]])\', "ImmutableDenseMatrix([[Mul(Rational(3, 7), Symbol(\'p\')), '
                         "Add(Rational(3, 7), Mul(Integer(-1), Rational(3, 7), Symbol('p'))), Integer(0), "
                         'Integer(0), Rational(4, 7), Integer(0)], [Integer(0), Rational(4, 17), Rational(3, '
                         '17), Rational(2, 17), Integer(0), Rational(3, 17)], [Rational(2, 3), Integer(0), '
                         'Rational(1, 3), Integer(0), Integer(0), Integer(0)], [Rational(1, 12), Integer(0), '
                         'Rational(1, 12), Integer(0), Rational(1, 3), Rational(1, 12)], [Integer(0), '
                         'Rational(1, 3), Rational(1, 3), Integer(0), Rational(1, 3), Integer(0)], '
                         '[Rational(1, 3), Integer(0), Integer(0), Integer(0), Rational(2, 9), '
                         'Integer(0)]])"]',
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Mul(Integer(-1), Integer(4655), '
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           "Rational(54, 23), Mul(Add(Mul(Integer(936), Symbol('p')), Integer(-3555)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           "Rational(108, 391), Mul(Add(Mul(Integer(333), Symbol('p')), Integer(-4545)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           'Rational(171, 391)], [Mul(Integer(-1), Integer(1918), Pow(Add(Mul(Integer(1173), '
                           "Symbol('p')), Integer(-1955)), Integer(-1))), Rational(54, 23), "
                           "Mul(Add(Mul(Integer(936), Symbol('p')), Integer(-2382)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           "Rational(108, 391), Mul(Add(Mul(Integer(333), Symbol('p')), Integer(-2199)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           'Rational(171, 391)], [Mul(Integer(-1), Integer(4655), Pow(Add(Mul(Integer(1173), '
                           "Symbol('p')), Integer(-1955)), Integer(-1))), Rational(54, 23), "
                           "Mul(Pow(Add(Mul(Integer(2346), Symbol('p')), Integer(-3910)), Integer(-1)), "
                           "Add(Mul(Integer(5391), Symbol('p')), Integer(-12975))), Rational(108, 391), "
                           "Mul(Add(Mul(Integer(333), Symbol('p')), Integer(-4545)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           'Rational(171, 391)], [Mul(Integer(-1), Integer(4123), Pow(Add(Mul(Integer(2346), '
                           "Symbol('p')), Integer(-3910)), Integer(-1))), Rational(59, 46), "
                           "Mul(Add(Mul(Integer(1935), Symbol('p')), Integer(-4992)), "
                           "Pow(Add(Mul(Integer(2346), Symbol('p')), Integer(-3910)), Integer(-1))), "
                           "Rational(450, 391), Mul(Add(Mul(Integer(801), Symbol('p')), Integer(-3102)), "
                           "Pow(Add(Mul(Integer(1173), Symbol('p')), Integer(-1955)), Integer(-1))), "
                           'Rational(126, 391)], [Mul(Integer(-1), Integer(6573), Pow(Add(Mul(Integer(2346), '
                           "Symbol('p')), Integer(-3910)), Integer(-1))), Rational(54, 23), "
                           "Mul(Pow(Add(Mul(Integer(4692), Symbol('p')), Integer(-7820)), Integer(-1)), "
                           "Add(Mul(Integer(7263), Symbol('p')), Integer(-17739))), Rational(108, 391), "
                           "Mul(Pow(Add(Mul(Integer(2346), Symbol('p')), Integer(-3910)), Integer(-1)), "
                           "Add(Mul(Integer(4185), Symbol('p')), Integer(-12609))), Rational(171, 391)], "
                           "[Mul(Integer(-1), Integer(2282), Pow(Add(Mul(Integer(1173), Symbol('p')), "
                           'Integer(-1955)), Integer(-1))), Rational(30, 23), Mul(Add(Mul(Integer(1431), '
                           "Symbol('p')), Integer(-4341)), Pow(Add(Mul(Integer(2346), Symbol('p')), "
                           'Integer(-3910)), Integer(-1))), Rational(60, 391), Mul(Add(Mul(Integer(576), '
                           "Symbol('p')), Integer(-2916)), Pow(Add(Mul(Integer(1173), Symbol('p')), "
                           'Integer(-1955)), Integer(-1))), Rational(486, 391)]])',
            'regular': 'false',
            'rows': [['3*p/7', '3/7 - 3*p/7', '0', '4/7', '0', '0', '0', '0'],
                     ['0', '4/17', '3/17', '0', '3/17', '2/17', '2/17', '3/17'],
                     ['0', '0', '1', '0', '0', '0', '0', '0'],
                     ['0', '1/3', '0', '1/3', '1/3', '0', '0', '0'],
                     ['2/3', '0', '0', '0', '1/3', '0', '0', '0'],
                     ['1/12', '0', '1/4', '1/3', '1/12', '0', '1/6', '1/12'],
                     ['0', '0', '0', '0', '0', '0', '1', '0'],
                     ['1/3', '0', '0', '2/9', '0', '0', '4/9', '0']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Add(Integer(1), Mul(Integer(-1), '
                          "Symbol('tau0'))), Integer(0), Integer(0), Integer(0), Symbol('tau0'), "
                          'Integer(0)]])'},
 'rand14': {'E2': 'Rational(4, 25)',
            'P2': 'Rational(4, 25)',
            'Pne': 'Rational(4, 25)',
            'abschain': 'false',
            'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
            'canonical': "[['Integer(0)', 'Integer(1)'], 'ImmutableDenseMatrix([[Rational(4, 5), Rational(1, "
                         "5)], [Integer(1), Integer(0)]])']",
            'classes': "[[['Integer(0)', 'Integer(1)'], 'true', 'Integer(1)']]",
            'decompose': "[['Integer(0)', 'Integer(1)'], 'ImmutableDenseMatrix([[Rational(4, 5), Rational(1, "
                         "5)], [Integer(1), Integer(0)]])', 'ImmutableDenseMatrix(0, 2, [])', "
                         "'ImmutableDenseMatrix([])']",
            'ergodic': 'true',
            'fundamental': 'ImmutableDenseMatrix([[Rational(35, 36), Rational(1, 36)], [Rational(5, 36), '
                           'Rational(31, 36)]])',
            'regular': 'true',
            'rows': [['4/5', '1/5'], ['1', '0']],
            'stationary': 'ImmutableDenseMatrix([[Rational(5, 6), Rational(1, 6)]])'},
 'rand15': {'E2': 'Rational(10, 9)',
            'P2': 'Rational(5, 9)',
            'Pne': 'Rational(5, 9)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)]])',
            'canonical': "[['Integer(2)', 'Integer(0)', 'Integer(1)'], "
                         '"ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)], [Rational(1, 3), '
                         "Rational(2, 3), Integer(0)], [Integer(0), Symbol('p'), Add(Integer(1), "
                         'Mul(Integer(-1), Symbol(\'p\')))]])"]',
            'classes': "[[['Integer(2)'], 'true', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)'], "
                       "[['Integer(1)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(2)', 'Integer(0)', 'Integer(1)'], "
                         "'ImmutableDenseMatrix([[Integer(1)]])', 'ImmutableDenseMatrix([[Rational(1, 3)], "
                         '[Integer(0)]])\', "ImmutableDenseMatrix([[Rational(2, 3), Integer(0)], '
                         '[Symbol(\'p\'), Add(Integer(1), Mul(Integer(-1), Symbol(\'p\')))]])"]',
            'ergodic': 'false',
            'fundamental': "ImmutableDenseMatrix([[Integer(3), Integer(0)], [Integer(3), Pow(Symbol('p'), "
                           'Integer(-1))]])',
            'regular': 'false',
            'rows': [['2/3', '0', '1/3'], ['p', '1 - p', '0'], ['0', '0', '1']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(1)]])'},
 'rand16': {'E2': 'Rational(62, 135)',
            'P2': 'Rational(2, 27)',
            'Pne': 'Rational(14, 45)',
            'abschain': 'false',
            'absorbing': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)]])',
            'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(3)', 'Integer(2)'], "
                         "'ImmutableDenseMatrix([[Rational(1, 3), Rational(4, 9), Rational(2, 9), "
                         'Integer(0)], [Integer(1), Integer(0), Integer(0), Integer(0)], [Rational(3, 5), '
                         'Rational(2, 5), Integer(0), Integer(0)], [Integer(1), Integer(0), Integer(0), '
                         "Integer(0)]])']",
            'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(3)'], 'true', 'Integer(1)'], "
                       "[['Integer(2)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(3)', 'Integer(2)'], "
                         "'ImmutableDenseMatrix([[Rational(1, 3), Rational(4, 9), Rational(2, 9)], "
                         '[Integer(1), Integer(0), Integer(0)], [Rational(3, 5), Rational(2, 5), '
                         "Integer(0)]])', 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0)]])', "
                         "'ImmutableDenseMatrix([[Integer(0)]])']",
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Integer(1)]])',
            'regular': 'false',
            'rows': [['1/3', '4/9', '0', '2/9'],
                     ['1', '0', '0', '0'],
                     ['1', '0', '0', '0'],
                     ['3/5', '2/5', '0', '0']],
            'stationary': 'ImmutableDenseMatrix([[Rational(45, 79), Rational(24, 79), Integer(0), '
                          'Rational(10, 79)]])'},
 'rand17': {'E2': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Rational(7, 4), "
                  "Symbol('p')), Rational(11, 4))",
            'P2': "Add(Rational(1, 4), Mul(Integer(-1), Rational(1, 4), Symbol('p')))",
            'Pne': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Rational(2, 3), Rational(1, 3)], [Rational(2, 3), '
                         'Rational(1, 3)], [Rational(2, 3), Rational(1, 3)]])',
            'canonical': "[['Integer(2)', 'Integer(4)', 'Integer(0)', 'Integer(1)', 'Integer(3)'], "
                         '"ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), '
                         'Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), Integer(0)], '
                         "[Integer(0), Integer(0), Symbol('p'), Add(Integer(1), Mul(Integer(-1), "
                         "Symbol('p'))), Integer(0)], [Rational(1, 2), Rational(1, 4), Integer(0), "
                         'Integer(0), Rational(1, 4)], [Integer(0), Integer(0), Integer(1), Integer(0), '
                         'Integer(0)]])"]',
            'classes': "[[['Integer(2)'], 'true', 'Integer(1)'], [['Integer(4)'], 'true', 'Integer(1)'], "
                       "[['Integer(0)', 'Integer(1)', 'Integer(3)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(2)', 'Integer(4)', 'Integer(0)', 'Integer(1)', 'Integer(3)'], "
                         "'ImmutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(1)]])', "
                         "'ImmutableDenseMatrix([[Integer(0), Integer(0)], [Rational(1, 2), Rational(1, 4)], "
                         '[Integer(0), Integer(0)]])\', "ImmutableDenseMatrix([[Symbol(\'p\'), '
                         "Add(Integer(1), Mul(Integer(-1), Symbol('p'))), Integer(0)], [Integer(0), "
                         'Integer(0), Rational(1, 4)], [Integer(1), Integer(0), Integer(0)]])"]',
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Mul(Integer(-1), Integer(4), Pow(Add(Mul(Integer(3), '
                           "Symbol('p')), Integer(-3)), Integer(-1))), Rational(4, 3), Rational(1, 3)], "
                           "[Mul(Integer(-1), Pow(Add(Mul(Integer(3), Symbol('p')), Integer(-3)), "
                           'Integer(-1))), Rational(4, 3), Rational(1, 3)], [Mul(Integer(-1), Integer(4), '
                           "Pow(Add(Mul(Integer(3), Symbol('p')), Integer(-3)), Integer(-1))), Rational(4, "
                           '3), Rational(4, 3)]])',
            'regular': 'false',
            'rows': [['p', '1 - p', '0', '0', '0'],
                     ['0', '0', '1/2', '1/4', '1/4'],
                     ['0', '0', '1', '0', '0'],
                     ['1', '0', '0', '0', '0'],
                     ['0', '0', '0', '0', '1']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Add(Integer(1), Mul(Integer(-1), '
                          "Symbol('tau0'))), Integer(0), Symbol('tau0')]])"},
 'rand18': {'E2': 'Integer(0)',
            'P2': 'Integer(0)',
            'Pne': 'Integer(0)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)], [Integer(1)], [Integer(1)], '
                         '[Integer(1)]])',
            'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(5)', 'Integer(4)', "
                         "'Integer(3)'], 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), "
                         'Integer(0), Integer(0), Integer(0)], [Integer(0), Rational(1, 8), Rational(1, 2), '
                         'Rational(3, 8), Integer(0), Integer(0)], [Rational(1, 2), Integer(0), Integer(0), '
                         'Rational(1, 2), Integer(0), Integer(0)], [Integer(0), Integer(0), Rational(1, 3), '
                         'Rational(1, 3), Rational(1, 3), Integer(0)], [Integer(0), Integer(0), Rational(4, '
                         '11), Rational(1, 11), Rational(3, 11), Rational(3, 11)], [Integer(0), Rational(1, '
                         "5), Integer(0), Rational(4, 5), Integer(0), Integer(0)]])']",
            'classes': "[[['Integer(0)'], 'true', 'Integer(1)'], [['Integer(1)', 'Integer(2)', 'Integer(5)', "
                       "'Integer(4)', 'Integer(3)'], 'false', 'Integer(1)']]",
            'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(5)', 'Integer(4)', "
                         "'Integer(3)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                         "'ImmutableDenseMatrix([[Integer(0)], [Rational(1, 2)], [Integer(0)], [Integer(0)], "
                         "[Integer(0)]])', 'ImmutableDenseMatrix([[Rational(1, 8), Rational(1, 2), "
                         'Rational(3, 8), Integer(0), Integer(0)], [Integer(0), Integer(0), Rational(1, 2), '
                         'Integer(0), Integer(0)], [Integer(0), Rational(1, 3), Rational(1, 3), Rational(1, '
                         '3), Integer(0)], [Integer(0), Rational(4, 11), Rational(1, 11), Rational(3, 11), '
                         'Rational(3, 11)], [Rational(1, 5), Integer(0), Rational(4, 5), Integer(0), '
                         "Integer(0)]])']",
            'ergodic': 'false',
            'fundamental': 'ImmutableDenseMatrix([[Rational(11, 9), Integer(2), Rational(25, 9), '
                           'Rational(275, 216), Rational(25, 72)], [Rational(1, 18), Integer(2), '
                           'Rational(35, 18), Rational(385, 432), Rational(35, 144)], [Rational(1, 9), '
                           'Integer(2), Rational(35, 9), Rational(385, 216), Rational(35, 72)], [Rational(1, '
                           '6), Integer(2), Rational(17, 6), Rational(385, 144), Rational(35, 48)], '
                           '[Rational(1, 3), Integer(2), Rational(11, 3), Rational(121, 72), Rational(35, '
                           '24)]])',
            'regular': 'false',
            'rows': [['1', '0', '0', '0', '0', '0'],
                     ['0', '1/8', '1/2', '0', '0', '3/8'],
                     ['1/2', '0', '0', '0', '0', '1/2'],
                     ['0', '1/5', '0', '0', '0', '4/5'],
                     ['0', '0', '4/11', '3/11', '3/11', '1/11'],
                     ['0', '0', '1/3', '0', '1/3', '1/3']],
            'stationary': 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0), '
                          'Integer(0), Integer(0)]])'},
 'rand19': {'E2': 'Rational(352, 245)',
            'P2': 'Rational(4, 35)',
            'Pne': 'Rational(158, 245)',
            'abschain': 'true',
            'absorbing': 'ImmutableDenseMatrix([[Add(Mul(Rational(1, 7), Pow(Add(Mul(Integer(30), '
                         "Symbol('p')), Integer(2430)), Integer(-1)), Add(Mul(Integer(210), Symbol('p')), "
                         "Integer(6503))), Mul(Integer(6), Pow(Add(Mul(Integer(5), Symbol('p')), "
                         'Integer(405)), Integer(-1))), Mul(Rational(616, 5), Pow(Add(Mul(Integer(3), '
                         "Symbol('p')), Integer(243)), Integer(-1)))), Mul(Rational(233, 10), "
                         "Pow(Add(Mul(Integer(3), Symbol('p')), Integer(243)), Integer(-1)))], "
                         "[Add(Mul(Rational(1, 7), Pow(Add(Mul(Integer(20), Symbol('p')), Integer(1620)), "
                         "Integer(-1)), Add(Mul(Integer(140), Symbol('p')), Integer(833))), Mul(Integer(9), "
                         "Pow(Add(Mul(Integer(5), Symbol('p')), Integer(405)), Integer(-1))), "
                         "Mul(Rational(466, 5), Pow(Add(Mul(Integer(2), Symbol('p')), Integer(162)), "
                         "Integer(-1))), Mul(Integer(15), Pow(Add(Symbol('p'), Integer(81)), Integer(-1)))), "
                         "Mul(Rational(233, 10), Pow(Add(Mul(Integer(2), Symbol('p')), Integer(162)), "
                         'Integer(-1)))], [Add(Mul(Rational(2, 5), Add(Integer(416), Mul(Integer(-1), '
                         "Integer(15), Symbol('p'))), Pow(Add(Mul(Integer(6), Symbol('p')), Integer(486)), "
                         "Integer(-1))), Mul(Rational(1, 7), Pow(Add(Mul(Integer(60), Symbol('p')), "
                         "Integer(4860)), Integer(-1)), Add(Mul(Integer(945), Symbol('p')), Integer(2996))), "
                         "Mul(Integer(21), Pow(Add(Mul(Integer(5), Symbol('p')), Integer(405)), "
                         "Integer(-1))), Mul(Integer(105), Pow(Add(Mul(Integer(3), Symbol('p')), "
                         'Integer(243)), Integer(-1)))), Mul(Rational(1, 10), Add(Integer(416), '
                         "Mul(Integer(-1), Integer(15), Symbol('p'))), Pow(Add(Mul(Integer(6), Symbol('p')), "
                         'Integer(486)), Integer(-1)))], [Add(Mul(Rational(3, 11), Add(Integer(264), '
                         "Mul(Integer(-1), Integer(11), Symbol('p'))), Pow(Add(Mul(Integer(3), Symbol('p')), "
                         'Integer(243)), Integer(-1))), Mul(Rational(2, 5), Add(Integer(438), '
                         "Mul(Integer(-1), Integer(55), Symbol('p'))), Pow(Add(Mul(Integer(6), Symbol('p')), "
                         'Integer(486)), Integer(-1))), Mul(Rational(1, 7), Pow(Add(Mul(Integer(60), '
                         "Symbol('p')), Integer(4860)), Integer(-1)), Add(Mul(Integer(2765), Symbol('p')), "
                         "Integer(3318))), Mul(Integer(63), Pow(Add(Mul(Integer(5), Symbol('p')), "
                         'Integer(405)), Integer(-1)))), Mul(Rational(1, 10), Add(Integer(438), '
                         "Mul(Integer(-1), Integer(55), Symbol('p'))), Pow(Add(Mul(Integer(6), Symbol('p')), "
                         'Integer(486)), Integer(-1)))], [Add(Mul(Rational(2, 5), Add(Integer(659), '
                         "Mul(Integer(-1), Integer(12), Symbol('p'))), Pow(Add(Mul(Integer(12), "
                         "Symbol('p')), Integer(972)), Integer(-1))), Mul(Rational(1, 7), "
                         "Pow(Add(Mul(Integer(120), Symbol('p')), Integer(9720)), Integer(-1)), "
                         "Add(Mul(Integer(1008), Symbol('p')), Integer(8099))), Mul(Rational(3, 11), "
                         "Pow(Add(Mul(Integer(10), Symbol('p')), Integer(810)), Integer(-1)), "
                         "Add(Mul(Integer(11), Symbol('p')), Integer(968))), Mul(Integer(105), "
                         "Pow(Add(Mul(Integer(6), Symbol('p')), Integer(486)), Integer(-1)))), "
                         'Mul(Rational(1, 10), Add(Integer(659), Mul(Integer(-1), Integer(12), '
                         "Symbol('p'))), Pow(Add(Mul(Integer(12), Symbol('p')), Integer(972)), "
                         'Integer(-1)))]])',
            'canonical': "[['Integer(1)', 'Integer(3)', 'Integer(0)', 'Integer(2)', 'Integer(6)', "
                         '\'Integer(4)\', \'Integer(5)\'], "ImmutableDenseMatrix([[Integer(1), Integer(0), '
                         'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
                         'Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
                         '[Rational(1, 7), Integer(0), Rational(4, 7), Rational(2, 7), Integer(0), '
                         'Integer(0), Integer(0)], [Rational(2, 5), Rational(1, 10), Rational(1, 10), '
                         'Integer(0), Rational(2, 5), Integer(0), Integer(0)], [Rational(3, 11), Integer(0), '
                         'Rational(1, 11), Rational(4, 11), Integer(0), Rational(3, 11), Integer(0)], '
                         "[Integer(0), Integer(0), Mul(Rational(1, 5), Symbol('p')), Add(Rational(1, 5), "
                         "Mul(Integer(-1), Rational(1, 5), Symbol('p'))), Rational(2, 5), Integer(0), "
                         'Rational(2, 5)], [Rational(3, 11), Integer(0), Rational(2, 11), Rational(1, 11), '
                         'Rational(4, 11), Integer(0), Rational(1, 11)]])"]',
            'classes': "[[['Integer(1)'], 'true', 'Integer(1)'], [['Integer(3)'], 'true', 'Integer(1)'], "
                       "[['Integer(0)', 'Integer(2)', 'Integer(6)', 'Integer(4)', 'Integer(5)'], 'false', "
                       "'Integer(1)']]",
            'decompose': "[['Integer(1)', 'Integer(3)', 'Integer(0)', 'Integer(2)', 'Integer(6)', "
                         "'Integer(4)', 'Integer(5)'], 'ImmutableDenseMatrix([[Integer(1), Integer(0)], "
                         "[Integer(0), Integer(1)]])', 'ImmutableDenseMatrix([[Rational(1, 7), Integer(0)], "
                         '[Rational(2, 5), Rational(1, 10)], [Rational(3, 11), Integer(0)], [Integer(0), '
                         "Integer(0)], [Rational(3, 11), Integer(0)]])', "
                         '"ImmutableDenseMatrix([[Rational(4, 7), Rational(2, 7), Integer(0), Integer(0), '
                         'Integer(0)], [Rational(1, 10), Integer(0), Rational(2, 5), Integer(0), '
                         'Integer(0)], [Rational(1, 11), Rational(4, 11), Integer(0), Rational(3, 11), '
                         "Integer(0)], [Mul(Rational(1, 5), Symbol('p')), Add(Rational(1, 5), "
                         "Mul(Integer(-1), Rational(1, 5), Symbol('p'))), Rational(2, 5), Integer(0), "
                         'Rational(2, 5)], [Rational(2, 11), Rational(1, 11), Rational(4, 11), Integer(0), '
                         'Rational(1, 11)]])"]',
            'ergodic': 'false',
            'fundamental': "ImmutableDenseMatrix([[Mul(Pow(Add(Mul(Integer(30), Symbol('p')), "
                           "Integer(2430)), Integer(-1)), Add(Mul(Integer(210), Symbol('p')), "
                           "Integer(6503))), Mul(Integer(233), Pow(Add(Mul(Integer(3), Symbol('p')), "
                           'Integer(243)), Integer(-1))), Mul(Integer(110), Pow(Add(Mul(Integer(3), '
                           "Symbol('p')), Integer(243)), Integer(-1))), Mul(Integer(10), "
                           "Pow(Add(Symbol('p'), Integer(81)), Integer(-1))), Mul(Integer(22), "
                           "Pow(Add(Mul(Integer(5), Symbol('p')), Integer(405)), Integer(-1)))], "
                           "[Mul(Pow(Add(Mul(Integer(20), Symbol('p')), Integer(1620)), Integer(-1)), "
                           "Add(Mul(Integer(140), Symbol('p')), Integer(833))), Mul(Integer(233), "
                           "Pow(Add(Mul(Integer(2), Symbol('p')), Integer(162)), Integer(-1))), "
                           "Mul(Integer(55), Pow(Add(Symbol('p'), Integer(81)), Integer(-1))), "
                           "Mul(Integer(15), Pow(Add(Symbol('p'), Integer(81)), Integer(-1))), "
                           "Mul(Integer(33), Pow(Add(Mul(Integer(5), Symbol('p')), Integer(405)), "
                           "Integer(-1)))], [Mul(Pow(Add(Mul(Integer(60), Symbol('p')), Integer(4860)), "
                           "Integer(-1)), Add(Mul(Integer(945), Symbol('p')), Integer(2996))), "
                           "Mul(Add(Integer(416), Mul(Integer(-1), Integer(15), Symbol('p'))), "
                           "Pow(Add(Mul(Integer(6), Symbol('p')), Integer(486)), Integer(-1))), "
                           "Mul(Integer(385), Pow(Add(Mul(Integer(3), Symbol('p')), Integer(243)), "
                           "Integer(-1))), Mul(Integer(35), Pow(Add(Symbol('p'), Integer(81)), "
                           "Integer(-1))), Mul(Integer(77), Pow(Add(Mul(Integer(5), Symbol('p')), "
                           "Integer(405)), Integer(-1)))], [Mul(Pow(Add(Mul(Integer(60), Symbol('p')), "
                           "Integer(4860)), Integer(-1)), Add(Mul(Integer(2765), Symbol('p')), "
                           'Integer(3318))), Mul(Add(Integer(438), Mul(Integer(-1), Integer(55), '
                           "Symbol('p'))), Pow(Add(Mul(Integer(6), Symbol('p')), Integer(486)), "
                           'Integer(-1))), Mul(Add(Integer(264), Mul(Integer(-1), Integer(11), '
                           "Symbol('p'))), Pow(Add(Mul(Integer(3), Symbol('p')), Integer(243)), "
                           "Integer(-1))), Mul(Integer(105), Pow(Add(Symbol('p'), Integer(81)), "
                           "Integer(-1))), Mul(Integer(231), Pow(Add(Mul(Integer(5), Symbol('p')), "
                           "Integer(405)), Integer(-1)))], [Mul(Pow(Add(Mul(Integer(120), Symbol('p')), "
                           "Integer(9720)), Integer(-1)), Add(Mul(Integer(1008), Symbol('p')), "
                           'Integer(8099))), Mul(Add(Integer(659), Mul(Integer(-1), Integer(12), '
                           "Symbol('p'))), Pow(Add(Mul(Integer(12), Symbol('p')), Integer(972)), "
                           "Integer(-1))), Mul(Integer(385), Pow(Add(Mul(Integer(6), Symbol('p')), "
                           'Integer(486)), Integer(-1))), Mul(Integer(35), Pow(Add(Mul(Integer(2), '
                           "Symbol('p')), Integer(162)), Integer(-1))), Mul(Pow(Add(Mul(Integer(10), "
                           "Symbol('p')), Integer(810)), Integer(-1)), Add(Mul(Integer(11), Symbol('p')), "
                           'Integer(968)))]])',
            'regular': 'false',
            'rows': [['4/7', '1/7', '2/7', '0', '0', '0', '0'],
                     ['0', '1', '0', '0', '0', '0', '0'],
                     ['1/10', '2/5', '0', '1/10', '0', '0', '2/5'],
                     ['0', '0', '0', '1', '0', '0', '0'],
                     ['p/5', '0', '1/5 - p/5', '0', '0', '2/5', '2/5'],
                     ['2/11', '3/11', '1/11', '0', '0', '1/11', '4/11'],
                     ['1/11', '3/11', '4/11', '0', '3/11', '0', '0']],
            'stationary': 'ImmutableDenseMatrix([[Integer(0), Add(Integer(1), Mul(Integer(-1), '
                          "Symbol('tau0'))), Integer(0), Symbol('tau0'), Integer(0), Integer(0), "
                          'Integer(0)]])'},
 'rand2': {'E2': 'Rational(2477, 1848)',
           'P2': 'Rational(4, 77)',
           'Pne': 'Rational(531, 616)',
           'abschain': 'false',
           'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
           'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                        "'ImmutableDenseMatrix([[Integer(0), Rational(1, 2), Rational(3, 8), Rational(1, "
                        '8)], [Integer(0), Rational(2, 3), Rational(1, 3), Integer(0)], [Rational(3, 11), '
                        'Rational(3, 11), Rational(4, 11), Rational(1, 11)], [Rational(2, 7), Integer(0), '
                        "Rational(4, 7), Rational(1, 7)]])']",
           'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], 'true', 'Integer(1)']]",
           'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)'], "
                        "'ImmutableDenseMatrix([[Integer(0), Rational(1, 2), Rational(3, 8), Rational(1, "
                        '8)], [Integer(0), Rational(2, 3), Rational(1, 3), Integer(0)], [Rational(3, 11), '
                        'Rational(3, 11), Rational(4, 11), Rational(1, 11)], [Rational(2, 7), Integer(0), '
                        "Rational(4, 7), Rational(1, 7)]])', 'ImmutableDenseMatrix(0, 4, [])', "
                        "'ImmutableDenseMatrix([])']",
           'ergodic': 'true',
           'fundamental': 'ImmutableDenseMatrix([[Rational(162048, 177419), Rational(-2169, 177419), '
                          'Rational(437, 16129), Rational(12733, 177419)], [Rational(-27944, 177419), '
                          'Rational(232527, 177419), Rational(-960, 16129), Rational(-16604, 177419)], '
                          '[Rational(33016, 177419), Rational(-4596, 16129), Rational(16566, 16129), '
                          'Rational(12733, 177419)], [Rational(52320, 177419), Rational(-131328, 177419), '
                          'Rational(4374, 16129), Rational(208313, 177419)]])',
           'regular': 'true',
           'rows': [['0', '1/2', '3/8', '1/8'],
                    ['0', '2/3', '1/3', '0'],
                    ['3/11', '3/11', '4/11', '1/11'],
                    ['2/7', '0', '4/7', '1/7']],
           'stationary': 'ImmutableDenseMatrix([[Rational(160, 1397), Rational(654, 1397), Rational(46, '
                         '127), Rational(7, 127)]])'},
 'rand3': {'E2': 'Rational(25, 16)',
           'P2': 'Rational(1, 6)',
           'Pne': 'Rational(31, 48)',
           'abschain': 'false',
           'absorbing': 'ERR:ShapeError:Matrix size mismatch: (5, 5) * (0, 5).',
           'canonical': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)', 'Integer(4)'], "
                        '"ImmutableDenseMatrix([[Rational(1, 4), Rational(1, 4), Integer(0), Integer(0), '
                        'Rational(1, 2)], [Rational(1, 6), Integer(0), Rational(2, 3), Integer(0), '
                        "Rational(1, 6)], [Symbol('p'), Integer(0), Integer(0), Add(Integer(1), "
                        "Mul(Integer(-1), Symbol('p'))), Integer(0)], [Rational(1, 3), Rational(1, 6), "
                        'Rational(1, 6), Integer(0), Rational(1, 3)], [Rational(1, 2), Integer(0), '
                        'Rational(1, 2), Integer(0), Integer(0)]])"]',
           'classes': "[[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)', 'Integer(4)'], 'true', "
                      "'Integer(1)']]",
           'decompose': "[['Integer(0)', 'Integer(1)', 'Integer(2)', 'Integer(3)', 'Integer(4)'], "
                        '"ImmutableDenseMatrix([[Rational(1, 4), Rational(1, 4), Integer(0), Integer(0), '
                        'Rational(1, 2)], [Rational(1, 6), Integer(0), Rational(2, 3), Integer(0), '
                        "Rational(1, 6)], [Symbol('p'), Integer(0), Integer(0), Add(Integer(1), "
                        "Mul(Integer(-1), Symbol('p'))), Integer(0)], [Rational(1, 3), Rational(1, 6), "
                        'Rational(1, 6), Integer(0), Rational(1, 3)], [Rational(1, 2), Integer(0), '
                        'Rational(1, 2), Integer(0), Integer(0)]])", \'ImmutableDenseMatrix(0, 5, [])\', '
                        "'ImmutableDenseMatrix([])']",
           'ergodic': 'true',
           'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(352), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(-1), Integer(10364), Symbol('p')), Integer(92092)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                          "Integer(12), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(822), "
                          "Symbol('p')), Integer(9798)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(1194), Symbol('p')), Integer(-6693)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(192), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(25338), Symbol('p')), "
                          "Integer(-25530)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(-1), Integer(48), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(-1), Integer(2190), Symbol('p')), Integer(19734)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1)))], [Mul(Pow(Add(Mul(Integer(484), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), "
                          "Integer(-1)), Add(Mul(Integer(2288), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(7940), Symbol('p')), Integer(-16744))), Mul(Pow(Add(Mul(Integer(484), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), "
                          "Integer(-1)), Add(Mul(Integer(516), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(11370), Symbol('p')), Integer(77970))), Mul(Add(Mul(Integer(3834), "
                          "Symbol('p')), Integer(29187)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(-1), Integer(2448), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(-1), Integer(7902), Symbol('p')), Integer(10350)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(128), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(2086), Symbol('p')), "
                          "Integer(-11362)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1)))], "
                          "[Mul(Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1)), Add(Mul(Integer(5192), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(47408), Symbol('p')), "
                          "Integer(-16744))), Mul(Add(Mul(Integer(252), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(1710), Symbol('p')), Integer(-4554)), Pow(Add(Mul(Integer(484), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), "
                          "Integer(-1))), Mul(Add(Mul(Integer(6870), Symbol('p')), Integer(70449)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                          "Integer(5484), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(46128), "
                          "Symbol('p')), Integer(51612)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1)), Add(Mul(Integer(524), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(3296), Symbol('p')), "
                          "Integer(-11362)))], [Mul(Add(Mul(Integer(352), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(-1), Integer(16656), Symbol('p')), Integer(6578)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                          "Integer(12), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(1218), "
                          "Symbol('p')), Integer(4416)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(1678), Symbol('p')), Integer(-115)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(192), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(32400), Symbol('p')), Integer(70449)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                          "Integer(48), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(3048), "
                          "Symbol('p')), Integer(8073)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1)))], "
                          "[Mul(Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1)), Add(Mul(Integer(1320), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(2928), Symbol('p')), "
                          "Integer(14352))), Mul(Add(Mul(Integer(-1), Integer(12), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(-1), Integer(2010), Symbol('p')), Integer(-6348)), "
                          "Pow(Add(Mul(Integer(484), Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), "
                          "Symbol('p')), Integer(89401)), Integer(-1))), Mul(Add(Mul(Integer(2646), "
                          "Symbol('p')), Integer(13041)), Pow(Add(Mul(Integer(484), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(-1), Integer(1260), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(7056), Symbol('p')), Integer(-5796)), Pow(Add(Mul(Integer(484), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), "
                          "Integer(-1))), Mul(Add(Mul(Integer(436), Pow(Symbol('p'), Integer(2))), "
                          "Mul(Integer(8392), Symbol('p')), Integer(74152)), Pow(Add(Mul(Integer(484), "
                          "Pow(Symbol('p'), Integer(2))), Mul(Integer(13156), Symbol('p')), Integer(89401)), "
                          'Integer(-1)))]])',
           'regular': 'true',
           'rows': [['1/4', '1/4', '0', '0', '1/2'],
                    ['1/6', '0', '2/3', '0', '1/6'],
                    ['p', '0', '0', '1 - p', '0'],
                    ['1/3', '1/6', '1/6', '0', '1/3'],
                    ['1/2', '0', '1/2', '0', '0']],
           'stationary': "ImmutableDenseMatrix([[Mul(Pow(Add(Mul(Integer(22), Symbol('p')), Integer(299)), "
                         "Integer(-1)), Add(Mul(Integer(66), Symbol('p')), Integer(78))), "
                         "Mul(Add(Mul(Integer(6), Symbol('p')), Integer(30)), Pow(Add(Mul(Integer(22), "
                         "Symbol('p')), Integer(299)), Integer(-1))), Mul(Integer(63), "
                         "Pow(Add(Mul(Integer(22), Symbol('p')), Integer(299)), Integer(-1))), "
                         "Mul(Add(Integer(63), Mul(Integer(-1), Integer(63), Symbol('p'))), "
                         "Pow(Add(Mul(Integer(22), Symbol('p')), Integer(299)), Integer(-1))), "
                         "Mul(Add(Mul(Integer(13), Symbol('p')), Integer(65)), Pow(Add(Mul(Integer(22), "
                         "Symbol('p')), Integer(299)), Integer(-1)))]])"},
 'rand4': {'E2': 'Rational(93, 25)',
           'P2': 'Rational(4, 15)',
           'Pne': 'Integer(1)',
           'abschain': 'true',
           'absorbing': 'ImmutableDenseMatrix([[Integer(1)], [Integer(1)], [Integer(1)], [Integer(1)], '
                        '[Integer(1)]])',
           'canonical': "[['Integer(4)', 'Integer(3)', 'Integer(1)', 'Integer(5)', 'Integer(2)', "
                        "'Integer(0)'], 'ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), "
                        'Integer(0), Integer(0), Integer(0)], [Rational(4, 5), Rational(1, 5), Integer(0), '
                        'Integer(0), Integer(0), Integer(0)], [Integer(0), Rational(1, 6), Rational(1, 6), '
                        'Rational(2, 3), Integer(0), Integer(0)], [Integer(0), Rational(2, 5), Integer(0), '
                        'Integer(0), Rational(3, 5), Integer(0)], [Rational(2, 5), Rational(2, 5), '
                        'Rational(1, 5), Integer(0), Integer(0), Integer(0)], [Integer(0), Rational(2, 5), '
                        "Rational(2, 5), Integer(0), Rational(1, 5), Integer(0)]])']",
           'classes': "[[['Integer(4)'], 'true', 'Integer(1)'], [['Integer(3)'], 'false', 'Integer(1)'], "
                      "[['Integer(1)', 'Integer(5)', 'Integer(2)'], 'false', 'Integer(1)'], [['Integer(0)'], "
                      "'false', 'Integer(1)']]",
           'decompose': "[['Integer(4)', 'Integer(3)', 'Integer(1)', 'Integer(5)', 'Integer(2)', "
                        "'Integer(0)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                        "'ImmutableDenseMatrix([[Rational(4, 5)], [Integer(0)], [Integer(0)], [Rational(2, "
                        "5)], [Integer(0)]])', 'ImmutableDenseMatrix([[Rational(1, 5), Integer(0), "
                        'Integer(0), Integer(0), Integer(0)], [Rational(1, 6), Rational(1, 6), Rational(2, '
                        '3), Integer(0), Integer(0)], [Rational(2, 5), Integer(0), Integer(0), Rational(3, '
                        '5), Integer(0)], [Rational(2, 5), Rational(1, 5), Integer(0), Integer(0), '
                        'Integer(0)], [Rational(2, 5), Rational(2, 5), Integer(0), Rational(1, 5), '
                        "Integer(0)]])']",
           'ergodic': 'false',
           'fundamental': 'ImmutableDenseMatrix([[Rational(5, 4), Integer(0), Integer(0), Integer(0), '
                          'Integer(0)], [Rational(445, 452), Rational(150, 113), Rational(100, 113), '
                          'Rational(60, 113), Integer(0)], [Rational(415, 452), Rational(18, 113), '
                          'Rational(125, 113), Rational(75, 113), Integer(0)], [Rational(315, 452), '
                          'Rational(30, 113), Rational(20, 113), Rational(125, 113), Integer(0)], '
                          '[Rational(467, 452), Rational(66, 113), Rational(44, 113), Rational(49, 113), '
                          'Integer(1)]])',
           'regular': 'false',
           'rows': [['0', '2/5', '1/5', '2/5', '0', '0'],
                    ['0', '1/6', '0', '1/6', '0', '2/3'],
                    ['0', '1/5', '0', '2/5', '2/5', '0'],
                    ['0', '0', '0', '1/5', '4/5', '0'],
                    ['0', '0', '0', '0', '1', '0'],
                    ['0', '0', '3/5', '2/5', '0', '0']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), '
                         'Integer(0)]])'},
 'rand5': {'E2': 'Rational(15, 16)',
           'P2': 'Integer(0)',
           'Pne': 'Rational(13, 16)',
           'abschain': 'true',
           'absorbing': 'ImmutableDenseMatrix([[Rational(11, 13), Rational(2, 13)], [Rational(12, 13), '
                        "Rational(1, 13)], [Add(Mul(Rational(1, 2), Add(Mul(Integer(630), Symbol('p')), "
                        "Integer(360)), Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), "
                        "Integer(-1))), Mul(Rational(1, 2), Pow(Add(Mul(Integer(1365), Symbol('p')), "
                        "Integer(884)), Integer(-1)), Add(Mul(Integer(1680), Symbol('p')), Integer(544))), "
                        "Mul(Integer(16), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                        "Integer(-1)))), Add(Mul(Rational(1, 8), Pow(Add(Mul(Integer(1365), Symbol('p')), "
                        "Integer(884)), Integer(-1)), Add(Mul(Integer(1680), Symbol('p')), Integer(544))), "
                        "Mul(Integer(12), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                        "Integer(-1))))], [Add(Mul(Rational(2, 7), Add(Mul(Integer(35), Symbol('p')), "
                        "Integer(28)), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), Integer(-1))), "
                        "Mul(Rational(1, 2), Add(Mul(Integer(600), Symbol('p')), Integer(384)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))), "
                        "Mul(Rational(1, 2), Add(Mul(Integer(1340), Symbol('p')), Integer(816)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1)))), "
                        "Add(Mul(Rational(3, 14), Add(Mul(Integer(35), Symbol('p')), Integer(28)), "
                        "Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), Integer(-1))), "
                        "Mul(Rational(1, 8), Add(Mul(Integer(1340), Symbol('p')), Integer(816)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))))], "
                        "[Add(Mul(Rational(1, 2), Add(Mul(Integer(540), Symbol('p')), Integer(336)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))), "
                        "Mul(Rational(1, 2), Add(Mul(Integer(660), Symbol('p')), Integer(272)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))), "
                        "Mul(Rational(2, 7), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                        "Integer(-1)), Add(Mul(Integer(105), Symbol('p')), Integer(84)))), "
                        "Add(Mul(Rational(1, 8), Add(Mul(Integer(660), Symbol('p')), Integer(272)), "
                        "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))), "
                        "Mul(Rational(3, 14), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                        "Integer(-1)), Add(Mul(Integer(105), Symbol('p')), Integer(84))))]])",
           'canonical': "[['Integer(1)', 'Integer(2)', 'Integer(0)', 'Integer(3)', 'Integer(4)', "
                        '\'Integer(5)\', \'Integer(6)\'], "ImmutableDenseMatrix([[Integer(1), Integer(0), '
                        'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
                        'Integer(1), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], '
                        '[Rational(1, 2), Rational(1, 8), Integer(0), Rational(3, 8), Integer(0), '
                        'Integer(0), Integer(0)], [Rational(1, 2), Integer(0), Rational(1, 2), Integer(0), '
                        'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Mul(Rational(5, 9), '
                        "Symbol('p')), Integer(0), Add(Rational(5, 9), Mul(Integer(-1), Rational(5, 9), "
                        "Symbol('p'))), Rational(2, 9), Rational(2, 9)], [Integer(0), Integer(0), "
                        'Rational(2, 3), Integer(0), Integer(0), Integer(0), Rational(1, 3)], [Rational(2, '
                        '7), Rational(3, 14), Integer(0), Rational(3, 14), Rational(2, 7), Integer(0), '
                        'Integer(0)]])"]',
           'classes': "[[['Integer(1)'], 'true', 'Integer(1)'], [['Integer(2)'], 'true', 'Integer(1)'], "
                      "[['Integer(0)', 'Integer(3)'], 'false', 'Integer(2)'], [['Integer(4)', 'Integer(5)', "
                      "'Integer(6)'], 'false', 'Integer(1)']]",
           'decompose': "[['Integer(1)', 'Integer(2)', 'Integer(0)', 'Integer(3)', 'Integer(4)', "
                        "'Integer(5)', 'Integer(6)'], 'ImmutableDenseMatrix([[Integer(1), Integer(0)], "
                        "[Integer(0), Integer(1)]])', 'ImmutableDenseMatrix([[Rational(1, 2), Rational(1, "
                        '8)], [Rational(1, 2), Integer(0)], [Integer(0), Integer(0)], [Integer(0), '
                        "Integer(0)], [Rational(2, 7), Rational(3, 14)]])', "
                        '"ImmutableDenseMatrix([[Integer(0), Rational(3, 8), Integer(0), Integer(0), '
                        'Integer(0)], [Rational(1, 2), Integer(0), Integer(0), Integer(0), Integer(0)], '
                        "[Mul(Rational(5, 9), Symbol('p')), Integer(0), Add(Rational(5, 9), Mul(Integer(-1), "
                        "Rational(5, 9), Symbol('p'))), Rational(2, 9), Rational(2, 9)], [Rational(2, 3), "
                        'Integer(0), Integer(0), Integer(0), Rational(1, 3)], [Integer(0), Rational(3, 14), '
                        'Rational(2, 7), Integer(0), Integer(0)]])"]',
           'ergodic': 'false',
           'fundamental': 'ImmutableDenseMatrix([[Rational(16, 13), Rational(6, 13), Integer(0), Integer(0), '
                          'Integer(0)], [Rational(8, 13), Rational(16, 13), Integer(0), Integer(0), '
                          "Integer(0)], [Mul(Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), "
                          "Integer(-1)), Add(Mul(Integer(1680), Symbol('p')), Integer(544))), "
                          "Mul(Add(Mul(Integer(630), Symbol('p')), Integer(360)), Pow(Add(Mul(Integer(1365), "
                          "Symbol('p')), Integer(884)), Integer(-1))), Mul(Integer(189), "
                          "Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), Integer(-1))), "
                          "Mul(Integer(42), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                          "Integer(-1))), Mul(Integer(56), Pow(Add(Mul(Integer(105), Symbol('p')), "
                          "Integer(68)), Integer(-1)))], [Mul(Add(Mul(Integer(1340), Symbol('p')), "
                          "Integer(816)), Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), "
                          "Integer(-1))), Mul(Add(Mul(Integer(600), Symbol('p')), Integer(384)), "
                          "Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), Integer(-1))), "
                          "Mul(Integer(18), Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                          "Integer(-1))), Mul(Pow(Add(Mul(Integer(105), Symbol('p')), Integer(68)), "
                          "Integer(-1)), Add(Mul(Integer(105), Symbol('p')), Integer(72))), "
                          "Mul(Add(Mul(Integer(35), Symbol('p')), Integer(28)), Pow(Add(Mul(Integer(105), "
                          "Symbol('p')), Integer(68)), Integer(-1)))], [Mul(Add(Mul(Integer(660), "
                          "Symbol('p')), Integer(272)), Pow(Add(Mul(Integer(1365), Symbol('p')), "
                          "Integer(884)), Integer(-1))), Mul(Add(Mul(Integer(540), Symbol('p')), "
                          "Integer(336)), Pow(Add(Mul(Integer(1365), Symbol('p')), Integer(884)), "
                          "Integer(-1))), Mul(Integer(54), Pow(Add(Mul(Integer(105), Symbol('p')), "
                          'Integer(68)), Integer(-1))), Mul(Integer(12), Pow(Add(Mul(Integer(105), '
                          "Symbol('p')), Integer(68)), Integer(-1))), Mul(Pow(Add(Mul(Integer(105), "
                          "Symbol('p')), Integer(68)), Integer(-1)), Add(Mul(Integer(105), Symbol('p')), "
                          'Integer(84)))]])',
           'regular': 'false',
           'rows': [['0', '1/2', '1/8', '3/8', '0', '0', '0'],
                    ['0', '1', '0', '0', '0', '0', '0'],
                    ['0', '0', '1', '0', '0', '0', '0'],
                    ['1/2', '1/2', '0', '0', '0', '0', '0'],
                    ['5*p/9', '0', '0', '0', '5/9 - 5*p/9', '2/9', '2/9'],
                    ['2/3', '0', '0', '0', '0', '0', '1/3'],
                    ['0', '2/7', '3/14', '3/14', '2/7', '0', '0']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Add(Integer(1), Mul(Integer(-1), '
                         "Symbol('tau0'))), Symbol('tau0'), Integer(0), Integer(0), Integer(0), "
                         'Integer(0)]])'},
 'rand6': {'E2': 'Rational(101, 20)',
           'P2': 'Rational(1, 2)',
           'Pne': 'Rational(127, 140)',
           'abschain': 'true',
           'absorbing': 'ImmutableDenseMatrix([[Rational(1855, 11321), Rational(7724, 11321), Rational(1742, '
                        '11321)], [Rational(4906, 11321), Rational(4780, 11321), Rational(1635, 11321)], '
                        '[Rational(5809, 11321), Rational(4866, 11321), Rational(646, 11321)], '
                        '[Rational(1916, 11321), Rational(8808, 11321), Rational(597, 11321)], '
                        '[Rational(1060, 11321), Rational(6031, 11321), Rational(4230, 11321)]])',
           'canonical': "[['Integer(3)', 'Integer(7)', 'Integer(6)', 'Integer(0)', 'Integer(2)', "
                        "'Integer(1)', 'Integer(4)', 'Integer(5)'], 'ImmutableDenseMatrix([[Integer(1), "
                        'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), Integer(0), '
                        'Integer(0)], [Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), '
                        'Integer(0), Integer(0), Integer(0)], [Integer(0), Integer(0), Integer(1), '
                        'Integer(0), Integer(0), Integer(0), Integer(0), Integer(0)], [Integer(0), '
                        'Rational(2, 5), Rational(1, 10), Integer(0), Rational(3, 10), Integer(0), '
                        'Rational(1, 5), Integer(0)], [Rational(2, 7), Integer(0), Integer(0), Rational(1, '
                        '7), Integer(0), Rational(1, 7), Rational(1, 7), Rational(2, 7)], [Rational(1, 4), '
                        'Integer(0), Integer(0), Rational(1, 12), Rational(1, 12), Rational(1, 3), '
                        'Rational(1, 4), Integer(0)], [Integer(0), Rational(1, 2), Integer(0), Rational(1, '
                        '4), Integer(0), Rational(1, 4), Integer(0), Integer(0)], [Integer(0), Rational(1, '
                        '7), Rational(2, 7), Rational(4, 7), Integer(0), Integer(0), Integer(0), '
                        "Integer(0)]])']",
           'classes': "[[['Integer(3)'], 'true', 'Integer(1)'], [['Integer(7)'], 'true', 'Integer(1)'], "
                      "[['Integer(6)'], 'true', 'Integer(1)'], [['Integer(0)', 'Integer(2)', 'Integer(1)', "
                      "'Integer(4)', 'Integer(5)'], 'false', 'Integer(1)']]",
           'decompose': "[['Integer(3)', 'Integer(7)', 'Integer(6)', 'Integer(0)', 'Integer(2)', "
                        "'Integer(1)', 'Integer(4)', 'Integer(5)'], 'ImmutableDenseMatrix([[Integer(1), "
                        'Integer(0), Integer(0)], [Integer(0), Integer(1), Integer(0)], [Integer(0), '
                        "Integer(0), Integer(1)]])', 'ImmutableDenseMatrix([[Integer(0), Rational(2, 5), "
                        'Rational(1, 10)], [Rational(2, 7), Integer(0), Integer(0)], [Rational(1, 4), '
                        'Integer(0), Integer(0)], [Integer(0), Rational(1, 2), Integer(0)], [Integer(0), '
                        "Rational(1, 7), Rational(2, 7)]])', 'ImmutableDenseMatrix([[Integer(0), Rational(3, "
                        '10), Integer(0), Rational(1, 5), Integer(0)], [Rational(1, 7), Integer(0), '
                        'Rational(1, 7), Rational(1, 7), Rational(2, 7)], [Rational(1, 12), Rational(1, 12), '
                        'Rational(1, 3), Rational(1, 4), Integer(0)], [Rational(1, 4), Integer(0), '
                        'Rational(1, 4), Integer(0), Integer(0)], [Rational(4, 7), Integer(0), Integer(0), '
                        "Integer(0), Integer(0)]])']",
           'ergodic': 'false',
           'fundamental': 'ImmutableDenseMatrix([[Rational(13860, 11321), Rational(4361, 11321), '
                          'Rational(2436, 11321), Rational(4004, 11321), Rational(1246, 11321)], '
                          '[Rational(5470, 11321), Rational(13328, 11321), Rational(4392, 11321), '
                          'Rational(4096, 11321), Rational(3808, 11321)], [Rational(4100, 11321), '
                          'Rational(2891, 11321), Rational(19932, 11321), Rational(6216, 11321), '
                          'Rational(826, 11321)], [Rational(4490, 11321), Rational(1813, 11321), '
                          'Rational(5592, 11321), Rational(13876, 11321), Rational(518, 11321)], '
                          '[Rational(7920, 11321), Rational(2492, 11321), Rational(1392, 11321), '
                          'Rational(2288, 11321), Rational(12033, 11321)]])',
           'regular': 'false',
           'rows': [['0', '0', '3/10', '0', '1/5', '0', '1/10', '2/5'],
                    ['1/12', '1/3', '1/12', '1/4', '1/4', '0', '0', '0'],
                    ['1/7', '1/7', '0', '2/7', '1/7', '2/7', '0', '0'],
                    ['0', '0', '0', '1', '0', '0', '0', '0'],
                    ['1/4', '1/4', '0', '0', '0', '0', '0', '1/2'],
                    ['4/7', '0', '0', '0', '0', '0', '2/7', '1/7'],
                    ['0', '0', '0', '0', '0', '0', '1', '0'],
                    ['0', '0', '0', '0', '0', '0', '0', '1']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(0), Integer(0), Add(Mul(Integer(-1), '
                         "Symbol('tau0')), Mul(Integer(-1), Symbol('tau1')), Integer(1)), Integer(0), "
                         "Integer(0), Symbol('tau0'), Symbol('tau1')]])"},
 'rand7': {'E2': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
           'P2': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
           'Pne': "Add(Integer(1), Mul(Integer(-1), Pow(Symbol('p'), Integer(2))))",
           'abschain': 'true',
           'absorbing': 'ImmutableDenseMatrix([[Mul(Integer(-1), Add(Integer(1), Mul(Integer(-1), '
                        "Symbol('p'))), Pow(Add(Symbol('p'), Integer(-1)), Integer(-1)))]])",
           'canonical': '[[\'Integer(1)\', \'Integer(0)\'], "ImmutableDenseMatrix([[Integer(1), Integer(0)], '
                        '[Add(Integer(1), Mul(Integer(-1), Symbol(\'p\'))), Symbol(\'p\')]])"]',
           'classes': "[[['Integer(1)'], 'true', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)']]",
           'decompose': "[['Integer(1)', 'Integer(0)'], 'ImmutableDenseMatrix([[Integer(1)]])', "
                        '"ImmutableDenseMatrix([[Add(Integer(1), Mul(Integer(-1), Symbol(\'p\')))]])", '
                        '"ImmutableDenseMatrix([[Symbol(\'p\')]])"]',
           'ergodic': 'false',
           'fundamental': "ImmutableDenseMatrix([[Mul(Integer(-1), Pow(Add(Symbol('p'), Integer(-1)), "
                          'Integer(-1)))]])',
           'regular': 'false',
           'rows': [['p', '1 - p'], ['0', '1']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Integer(1)]])'},
 'rand8': {'E2': 'Integer(1)',
           'P2': 'Rational(1, 4)',
           'Pne': 'Rational(3, 4)',
           'abschain': 'false',
           'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
           'canonical': "[['Integer(0)', 'Integer(2)', 'Integer(1)'], 'ImmutableDenseMatrix([[Integer(0), "
                        'Integer(1), Integer(0)], [Rational(1, 4), Rational(1, 4), Rational(1, 2)], '
                        "[Rational(1, 5), Integer(0), Rational(4, 5)]])']",
           'classes': "[[['Integer(0)', 'Integer(2)', 'Integer(1)'], 'true', 'Integer(1)']]",
           'decompose': "[['Integer(0)', 'Integer(2)', 'Integer(1)'], 'ImmutableDenseMatrix([[Integer(0), "
                        'Integer(1), Integer(0)], [Rational(1, 4), Rational(1, 4), Rational(1, 2)], '
                        "[Rational(1, 5), Integer(0), Rational(4, 5)]])', 'ImmutableDenseMatrix(0, 3, [])', "
                        "'ImmutableDenseMatrix([])']",
           'ergodic': 'true',
           'fundamental': 'ImmutableDenseMatrix([[Rational(257, 289), Rational(-220, 289), Rational(252, '
                          '289)], [Rational(2, 289), Rational(375, 289), Rational(-88, 289)], [Rational(19, '
                          '289), Rational(-50, 289), Rational(320, 289)]])',
           'regular': 'true',
           'rows': [['0', '0', '1'], ['1/5', '4/5', '0'], ['1/4', '1/2', '1/4']],
           'stationary': 'ImmutableDenseMatrix([[Rational(3, 17), Rational(10, 17), Rational(4, 17)]])'},
 'rand9': {'E2': 'Rational(8, 9)',
           'P2': 'Integer(0)',
           'Pne': 'Rational(8, 9)',
           'abschain': 'true',
           'absorbing': "ImmutableDenseMatrix([[Integer(1), Integer(0)], [Symbol('p'), Add(Integer(1), "
                        "Mul(Integer(-1), Symbol('p')))]])",
           'canonical': "[['Integer(1)', 'Integer(3)', 'Integer(0)', 'Integer(2)'], "
                        '"ImmutableDenseMatrix([[Integer(1), Integer(0), Integer(0), Integer(0)], '
                        '[Integer(0), Integer(1), Integer(0), Integer(0)], [Rational(2, 3), Integer(0), '
                        'Rational(1, 3), Integer(0)], [Integer(0), Add(Integer(1), Mul(Integer(-1), '
                        'Symbol(\'p\'))), Symbol(\'p\'), Integer(0)]])"]',
           'classes': "[[['Integer(1)'], 'true', 'Integer(1)'], [['Integer(0)'], 'false', 'Integer(1)'], "
                      "[['Integer(3)'], 'true', 'Integer(1)'], [['Integer(2)'], 'false', 'Integer(1)']]",
           'decompose': "[['Integer(1)', 'Integer(3)', 'Integer(0)', 'Integer(2)'], "
                        "'ImmutableDenseMatrix([[Integer(1), Integer(0)], [Integer(0), Integer(1)]])', "
                        '"ImmutableDenseMatrix([[Rational(2, 3), Integer(0)], [Integer(0), Add(Integer(1), '
                        'Mul(Integer(-1), Symbol(\'p\')))]])", "ImmutableDenseMatrix([[Rational(1, 3), '
                        'Integer(0)], [Symbol(\'p\'), Integer(0)]])"]',
           'ergodic': 'false',
           'fundamental': 'ImmutableDenseMatrix([[Rational(3, 2), Integer(0)], [Mul(Rational(3, 2), '
                          "Symbol('p')), Integer(1)]])",
           'regular': 'false',
           'rows': [['1/3', '2/3', '0', '0'],
                    ['0', '1', '0', '0'],
                    ['p', '0', '0', '1 - p'],
                    ['0', '0', '0', '1']],
           'stationary': 'ImmutableDenseMatrix([[Integer(0), Add(Integer(1), Mul(Integer(-1), '
                         "Symbol('tau0'))), Integer(0), Symbol('tau0')]])"}}


def _value(text):
    if text.startswith('['):
        return _tree(ast.literal_eval(text))
    return eval(text, NS)


def _tree(v):
    if isinstance(v, list):
        return [_tree(x) for x in v]
    if isinstance(v, bool):
        return v
    return eval(v, NS)


def _chain(rows):
    T = S.Matrix([[S.sympify(x) for x in r] for r in rows])
    return DiscreteMarkovChain('X', list(range(T.rows)), T), T


def _same(got, expected):
    if got == expected:
        return True
    if isinstance(got, (list, tuple)) and isinstance(expected, (list, tuple)):
        return len(got) == len(expected) and all(_same(a, b) for a, b in zip(got, expected))
    if isinstance(expected, S.MatrixBase):
        return (isinstance(got, S.MatrixBase) and got.shape == expected.shape
                and all(_same(a, b) for a, b in zip(got, expected)))
    # symbolic entries: equal up to simplification (the reference form is
    # produced by linsolve/inv over the polynomial domain)
    return expected.free_symbols and S.simplify(got - expected) == 0


def _call(f):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return f()
    except Exception as e:  # noqa
        return 'ERR:%s:%s' % (type(e).__name__, str(e)[:120])


METHODS = {
    'classes': lambda X: X.communication_classes(),
    'decompose': lambda X: X.decompose(),
    'canonical': lambda X: X.canonical_form(),
    'stationary': lambda X: X.stationary_distribution(),
    'fundamental': lambda X: X.fundamental_matrix(),
    'absorbing': lambda X: X.absorbing_probabilities(),
    'regular': lambda X: X.is_regular(),
    'ergodic': lambda X: X.is_ergodic(),
    'abschain': lambda X: X.is_absorbing_chain(),
    'P2': lambda X: P(S.Eq(X[2], X.number_of_states - 1), S.Eq(X[0], 0)),
    'E2': lambda X: E(X[2], S.Eq(X[0], 0)),
    'Pne': lambda X: P(S.Ne(X[2], 0), S.Eq(X[0], 0)),
}


@pytest.mark.parametrize('name', sorted(CASES))
def test_reference_values(name):
    case = CASES[name]
    X, T = _chain(case['rows'])
    for method, f in METHODS.items():
        expected = case[method]
        got = _call(lambda: f(X))
        if expected.startswith('ERR:'):
            assert isinstance(got, str) and got.split(':')[:2] == expected.split(':')[:2], (method, got, expected)
            if 'fundamental matrix' in expected:
                assert got == expected
            continue
        assert not isinstance(got, str), (method, got)
        assert _same(got, _value(expected)), (method, got, expected)


@pytest.mark.parametrize('name', sorted(CASES))
def test_two_time_conditions_follow_the_markov_property(name):
    X, T = _chain(CASES[name]['rows'])
    n = T.rows
    got = _call(lambda: P(S.Eq(X[3], 0), S.Eq(X[0], n - 1) & S.Eq(X[1], 0)))
    if T[n - 1, 0].is_zero:  # also Float(0.0)
        assert got == 'ERR:ValueError:conditioning event has zero probability'
    else:
        assert _same(got, (T**2)[0, 0])


@pytest.mark.parametrize('name', sorted(CASES))
def test_native_kernels_answer_rational_chains(name):
    # Rational, Float (RR partial pivoting) and symbolic (fraction field)
    # chains are solved natively when the stationary vector is unique; the
    # non-unique ones (tau parameters) decline to the linsolve fallback.
    X, T = _chain(CASES[name]['rows'])
    T = X.transition_probabilities
    assert _markov_native.communication_classes(T) is not None
    assert _markov_native.decompose(T) is not None
    expected = CASES[name]['stationary']
    unique = 'tau' not in expected
    got = _call(lambda: _markov_native.stationary_distribution(T))
    if expected.startswith('ERR:'):
        assert isinstance(got, str) and got.split(':')[:2] == expected.split(':')[:2], (got, expected)
        return
    assert (got is not None) == unique
    if unique:
        assert _same(got, _value(expected)), (got, expected)


def test_reference_test_file_cases():
    R = S.Rational
    TO2 = S.Matrix([[1, 0, 0], [R(1, 2), R(1, 4), R(1, 4)], [R(1, 3), R(1, 3), R(1, 3)]])
    TO3 = S.Matrix([[R(1, 4), R(3, 4), 0], [R(1, 3), R(1, 3), R(1, 3)], [0, R(1, 4), R(3, 4)]])
    Y2 = DiscreteMarkovChain('Y', trans_probs=TO2)
    Y3 = DiscreteMarkovChain('Y', trans_probs=TO3)
    assert Y3.fundamental_matrix() == S.ImmutableMatrix([[176, 81, -132], [36, 141, -52], [-44, -39, 208]])/125
    assert Y2.is_absorbing_chain() == True
    assert Y3.is_absorbing_chain() == False
    assert Y2.canonical_form() == ([0, 1, 2], TO2)
    assert Y3.canonical_form() == ([0, 1, 2], TO3)
    assert Y2.decompose() == ([0, 1, 2], TO2[0:1, 0:1], TO2[1:3, 0:1], TO2[1:3, 1:3])
    assert Y3.decompose() == ([0, 1, 2], TO3, S.Matrix(0, 3, []), S.Matrix(0, 0, []))
    assert type(Y3.decompose()[1]) is S.ImmutableDenseMatrix
    assert type(Y2.canonical_form()[1]) is S.ImmutableDenseMatrix
    assert type(Y3.fundamental_matrix()) is S.ImmutableDenseMatrix
    TO4 = S.Matrix([[R(1, 5), R(2, 5), R(2, 5)], [R(1, 10), R(1, 2), R(2, 5)], [R(3, 5), R(3, 10), R(1, 10)]])
    Y4 = DiscreteMarkovChain('Y', trans_probs=TO4)
    w = S.ImmutableMatrix([[R(11, 39), R(16, 39), R(4, 13)]])
    assert Y4.limiting_distribution == w
    assert Y4.is_regular() == True
    assert Y4.is_ergodic() == True
    h = S.S.Half
    TO6 = S.Matrix([[1, 0, 0, 0, 0], [h, 0, h, 0, 0], [0, h, 0, h, 0], [0, 0, h, 0, h], [0, 0, 0, 0, 1]])
    Y6 = DiscreteMarkovChain('Y', trans_probs=TO6)
    assert Y6.absorbing_probabilities() == S.ImmutableMatrix([[R(3, 4), R(1, 4)], [h, h], [R(1, 4), R(3, 4)]])
    with pytest.warns(S.SymPyDeprecationWarning if hasattr(S, 'SymPyDeprecationWarning') else Warning):
        Y6.absorbing_probabilites()
    assert Y6.is_absorbing_state(0) is True
    assert Y6.is_absorbing_state(1) is False
    X = DiscreteMarkovChain('X', trans_probs=S.Matrix([]))
    assert X.communication_classes() == []
    assert X.decompose() == ([], S.Matrix([]), S.Matrix([]), S.Matrix([]))
    assert X.is_regular() == False
    assert X.is_ergodic() == False
    T = S.Matrix([[0, R(1, 2), R(1, 2), 0, 0], [0, 0, 0, 1, 0], [R(1, 2), 0, R(1, 2), 0, 0],
                  [0, 1, 0, 0, 0], [0, R(3, 10), 0, R(3, 10), R(2, 5)]])
    Z = DiscreteMarkovChain('Z', ['a', 'b', 'c', 'd', 'e'], T)
    assert Z.communication_classes() == [([S.Symbol('b'), S.Symbol('d')], S.true, 2),
                                         ([S.Symbol('a'), S.Symbol('c')], S.false, 1),
                                         ([S.Symbol('e')], S.false, 1)]


def test_continuous_chain():
    R = S.Rational
    T1 = S.Matrix([[-2, 2, 0], [0, -1, 1], [R(3, 2), R(3, 2), -3]])
    C1 = ContinuousMarkovChain('C', [0, 1, 2], T1)
    assert C1.limiting_distribution() == S.ImmutableMatrix([[R(3, 19), R(12, 19), R(4, 19)]])
    T2 = S.Matrix([[-1, 1, 0], [1, -1, 0], [0, 1, -1]])
    C2 = ContinuousMarkovChain('C', [0, 1, 2], T2)
    A, t = C2.generator_matrix, S.symbols('t', positive=True)
    e, h = S.exp, S.S.Half
    assert C2.transition_probabilities(A)(t) == S.Matrix([[h + e(-2*t)/2, h - e(-2*t)/2, 0],
                                                          [h - e(-2*t)/2, h + e(-2*t)/2, 0],
                                                          [h - e(-t) + e(-2*t)/2, h - e(-2*t)/2, e(-t)]])
    G = S.MatrixSymbol('G', 3, 3)
    CS1 = ContinuousMarkovChain('C', [0, 1, 2], G)
    assert CS1.transition_probabilities(CS1.generator_matrix)(t) == S.exp(t*G)
    assert C2.transition_probabilities() is None


def test_custom_matrix_subclass_declines():
    class M(S.ImmutableDenseMatrix):
        pass
    T = M([[S.Rational(1, 2), S.Rational(1, 2)], [0, 1]])
    assert _markov_native.communication_classes(T) is None
    assert _markov_native.stationary_distribution(T) is None


def test_period_search_that_never_terminates_in_sympy():
    # SymPy 1.14 loops forever on this class: its breadth-first search keeps
    # the successors of the last frontier state only and never reaches state
    # 4. The cycles 0-1-3-4-0 and 0-2-0 give the period gcd(4, 2) = 2.
    h = S.S.Half
    T = S.Matrix([[0, h, h, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0]])
    X = DiscreteMarkovChain('X', list(range(5)), T)
    # the class order is the one of SymPy's strongly connected components,
    # computed before the period search
    assert X.communication_classes() == [([0, 1, 3, 4, 2], S.true, 2)]
    assert X.is_regular() == False
    assert X.is_ergodic() == True


def test_queries_use_the_path_kernel(monkeypatch):
    calls = []
    original = _markov_native.path_sums

    def spy(*args):
        result = original(*args)
        calls.append(result)
        return result
    monkeypatch.setattr(_markov_native, 'path_sums', spy)
    X, T = _chain(CASES['TO6']['rows'])
    assert P(S.Eq(X[3], 4), S.Eq(X[0], 1)) == S.Rational(1, 8)
    assert E(X[2], S.Eq(X[0], 1)) == 1
    assert len(calls) == 2 and all(c is not None for c in calls)
