"""Round-2 differential tests of the native Markov chain kernels against
SymPy 1.14 (expected values: SymPy srepr or exception type): float (RR)
and symbolic (fraction field) stationary and limiting distributions, the
fundamental and absorbing matrices built on them, and the query values
assembled as SymPy's unsimplified sums."""
import warnings

import pytest

import symox as S
from symox.stats import DiscreteMarkovChain, ContinuousMarkovChain, P, E

p, q, r = S.symbols('p q r')
a, b, c = S.symbols('a b c', positive=True)
LOCALS = {'a': a, 'b': b, 'c': c, 'p': p, 'q': q, 'r': r}

CASES = {'ctmcq': {'q2': 'Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), exp(Integer(-4))))',
           'q612': 'Add(Mul(Rational(1, 2), exp(Integer(-2))), Rational(1, 2))',
           'qE': 'Add(Mul(Rational(1, 2), exp(Integer(-2))), Rational(1, 2))'},
 'ctmcq2': {'q2': 'Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), exp(Integer(-4))))',
            'q612': 'Add(Mul(Rational(1, 2), exp(Integer(-2))), Rational(1, 2))',
            'qE': 'Add(Mul(Rational(1, 2), exp(Integer(-2))), Rational(1, 2))',
            'qgt': 'Add(Mul(Integer(-1), Rational(1, 2), exp(Integer(-2))), exp(Integer(-1)), Rational(1, '
                   '2))',
            'qne': 'Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), exp(Integer(-2))))'},
 'dec2': {'E2': "Float('0.27999999999999997', precision=53)",
          'E2sq': "Float('0.10000000000000001', precision=53)",
          'P2': "Float('0.27999999999999997', precision=53)",
          'P3': "Float('0.36000000000000004', precision=53)",
          'Pand': "Float('0.20999999999999999', precision=53)",
          'Pgt': "Float('0.27999999999999997', precision=53)",
          'Pne': "Float('0.27999999999999997', precision=53)",
          'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
          'fundamental': "ImmutableDenseMatrix([[Float('0.8359375', precision=53), "
                         "Float('0.16406249999999997', precision=53)], [Float('0.2109375', precision=53), "
                         "Float('0.7890625', precision=53)]])",
          'rows': [['0.3', '0.7'], ['0.9', '0.1']],
          'stationary': "ImmutableDenseMatrix([[Float('0.5625', precision=53), Float('0.4375', "
                        'precision=53)]])'},
 'dec3': {'E2': "Float('0.9375', precision=53)",
          'E2sq': "Float('2.375', precision=53)",
          'P2': "Float('0.375', precision=53)",
          'P3': "Float('0.28125', precision=53)",
          'Pand': 'Integer(0)',
          'Pgt': "Float('0.5625', precision=53)",
          'Pne': "Float('0.5625', precision=53)",
          'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
          'fundamental': "ImmutableDenseMatrix([[Float('1.0830449826989623', precision=53), "
                         "Float('0.31141868512110732', precision=53), Float('-0.39446366782006936', "
                         "precision=53)], [Float('0.14186851211072671', precision=53), "
                         "Float('0.7820069204152249', precision=53), Float('0.076124567474048374', "
                         "precision=53)], [Float('-0.21107266435986155', precision=53), "
                         "Float('-0.041522491349480939', precision=53), Float('1.2525951557093424', "
                         'precision=53)]])',
          'rows': [['0.25', '0.75', '0'], ['0.5', '0', '0.5'], ['0.125', '0.375', '0.5']],
          'stationary': "ImmutableDenseMatrix([[Float('0.29411764705882348', precision=53), "
                        "Float('0.3529411764705882', precision=53), Float('0.35294117647058826', "
                        'precision=53)]])'},
 'fl2': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
         'classes': '[([Integer(0), Integer(1)], true, Integer(1))]',
         'fundamental': "ImmutableDenseMatrix([[Float('0.87573964497041423', precision=53), "
                        "Float('0.12426035502958571', precision=53)], [Float('0.10650887573964503', "
                        "precision=53), Float('0.89349112426035504', precision=53)]])",
         'rows': [['0.3', '0.7'], ['0.6', '0.4']],
         'stationary': "ImmutableDenseMatrix([[Float('0.46153846153846145', precision=53), "
                       "Float('0.53846153846153855', precision=53)]])"},
 'fl3': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'classes': '[([Integer(0), Integer(1), Integer(2)], true, Integer(1))]',
         'fundamental': "ImmutableDenseMatrix([[Float('1.0', precision=53), Integer(0), Integer(0)], "
                        "[Integer(0), Float('1.0', precision=53), Integer(0)], [Integer(0), Integer(0), "
                        "Float('1.0', precision=53)]])",
         'rows': [['0.5', '0.25', '0.25'], ['0.5', '0.25', '0.25'], ['0.5', '0.25', '0.25']],
         'stationary': "ImmutableDenseMatrix([[Float('0.5', precision=53), Float('0.25', precision=53), "
                       "Float('0.25', precision=53)]])"},
 'fl4': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
         'classes': '[([Integer(0), Integer(1), Integer(2), Integer(3)], true, Integer(1))]',
         'fundamental': "ImmutableDenseMatrix([[Float('1.1992337164750959', precision=53), "
                        "Float('0.095785440613026865', precision=53), Float('-0.21647509578544064', "
                        "precision=53), Float('-0.078544061302682003', precision=53)], "
                        "[Float('-0.088122605363984682', precision=53), Float('1.0153256704980844', "
                        "precision=53), Float('0.10536398467432949', precision=53), "
                        "Float('-0.032567049808429123', precision=53)], [Float('0.049808429118773964', "
                        "precision=53), Float('-0.22605363984674329', precision=53), "
                        "Float('1.0708812260536398', precision=53), Float('0.10536398467432952', "
                        "precision=53)], [Float('-0.27203065134099624', precision=53), "
                        "Float('0.0038314176245210973', precision=53), Float('0.15134099616858235', "
                        "precision=53), Float('1.1168582375478928', precision=53)]])",
         'rows': [['0.5', '0.375', '0.0', '0.125'],
                  ['0.25', '0.375', '0.25', '0.125'],
                  ['0.375', '0.125', '0.25', '0.25'],
                  ['0.125', '0.375', '0.25', '0.25']],
         'stationary': "ImmutableDenseMatrix([[Float('0.33333333333333331', precision=53), "
                       "Float('0.33333333333333331', precision=53), Float('0.16666666666666669', "
                       "precision=53), Float('0.16666666666666666', precision=53)]])"},
 'fl5': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (5, 5) * (0, 5).',
         'classes': '[([Integer(0), Integer(1), Integer(2), Integer(3), Integer(4)], true, Integer(1))]',
         'fundamental': "ImmutableDenseMatrix([[Float('0.84000000000000008', precision=53), Integer(0), "
                        "Float('0.0596707818930041', precision=53), Float('0.0020576131687242974', "
                        "precision=53), Float('0.098271604938271598', precision=53)], "
                        "[Float('0.040000000000000036', precision=53), Float('1.3333333333333333', "
                        "precision=53), Float('-0.28600823045267487', precision=53), "
                        "Float('0.051440329218107005', precision=53), Float('-0.13876543209876546', "
                        "precision=53)], [Float('0.040000000000000036', precision=53), Integer(0), "
                        "Float('1.0967078189300412', precision=53), Float('-0.14609053497942384', "
                        "precision=53), Float('0.0093827160493827055', precision=53)], "
                        "[Float('0.040000000000000036', precision=53), Integer(0), "
                        "Float('-0.23662551440329219', precision=53), Float('1.1872427983539096', "
                        "precision=53), Float('0.0093827160493827038', precision=53)], "
                        "[Float('0.040000000000000042', precision=53), Integer(0), "
                        "Float('0.0596707818930041', precision=53), Float('0.0020576131687242987', "
                        "precision=53), Float('0.89827160493827174', precision=53)]])",
         'rows': [['0.0', '0.0', '0.375', '0.25', '0.375'],
                  ['0.25', '0.25', '0.125', '0.25', '0.125'],
                  ['0.25', '0.0', '0.375', '0.125', '0.25'],
                  ['0.25', '0.0', '0.125', '0.375', '0.25'],
                  ['0.25', '0.0', '0.375', '0.25', '0.125']],
         'stationary': "ImmutableDenseMatrix([[Float('0.19999999999999996', precision=53), Integer(0), "
                       "Float('0.31481481481481483', precision=53), Float('0.24074074074074073', "
                       "precision=53), Float('0.24444444444444446', precision=53)]])"},
 'flabs': {'absorbing': "ImmutableDenseMatrix([[Float('0.5', precision=53), Float('0.5', precision=53)]])",
           'classes': '[([Integer(0)], true, Integer(1)), ([Integer(2)], true, Integer(1)), ([Integer(1)], '
                      'false, Integer(1))]',
           'fundamental': "ImmutableDenseMatrix([[Float('2.0', precision=53)]])",
           'rows': [['1.0', '0', '0'], ['0.25', '0.5', '0.25'], ['0', '0', '1.0']],
           'stationary': "ImmutableDenseMatrix([[Add(Float('1.0', precision=53), Mul(Integer(-1), "
                         "Float('1.0', precision=53), Symbol('tau0'))), Integer(0), Symbol('tau0')]])"},
 'flabs2': {'absorbing': "ImmutableDenseMatrix([[Float('0.2857142857142857', precision=53), "
                         "Float('0.71428571428571441', precision=53)], [Float('0.42857142857142855', "
                         "precision=53), Float('0.5714285714285714', precision=53)]])",
            'classes': '[([Integer(0)], true, Integer(1)), ([Integer(3)], true, Integer(1)), ([Integer(1), '
                       'Integer(2)], false, Integer(1))]',
            'fundamental': "ImmutableDenseMatrix([[Float('1.4285714285714286', precision=53), "
                           "Float('0.5714285714285714', precision=53)], [Float('0.47619047619047616', "
                           "precision=53), Float('1.5238095238095237', precision=53)]])",
            'rows': [['1.0', '0', '0', '0'],
                     ['0.1', '0.2', '0.3', '0.4'],
                     ['0.25', '0.25', '0.25', '0.25'],
                     ['0', '0', '0', '1.0']],
            'stationary': "ImmutableDenseMatrix([[Add(Float('1.0', precision=53), Mul(Integer(-1), "
                          "Float('1.0', precision=53), Symbol('tau0'))), Integer(0), Integer(0), "
                          "Symbol('tau0')]])"},
 'flnonuniq': {'absorbing': "ERR:ValueError:The fundamental matrix doesn't exist.",
               'classes': '[([Integer(0), Integer(1)], true, Integer(1)), ([Integer(2)], true, Integer(1))]',
               'fundamental': "ERR:ValueError:The fundamental matrix doesn't exist.",
               'rows': [['0.5', '0.5', '0'], ['0.5', '0.5', '0'], ['0', '0', '1.0']],
               'stationary': "ImmutableDenseMatrix([[Add(Float('0.5', precision=53), Mul(Integer(-1), "
                             "Float('0.5', precision=53), Symbol('tau0'))), Add(Float('0.5', precision=53), "
                             "Mul(Integer(-1), Float('0.5', precision=53), Symbol('tau0'))), "
                             "Symbol('tau0')]])"},
 'g2': {'limiting': "ImmutableDenseMatrix([[Mul(Symbol('b', positive=True), Pow(Add(Symbol('a', "
                    "positive=True), Symbol('b', positive=True)), Integer(-1))), Mul(Symbol('a', "
                    "positive=True), Pow(Add(Symbol('a', positive=True), Symbol('b', positive=True)), "
                    'Integer(-1)))]])',
        'rows': [['-a', 'a'], ['b', '-b']]},
 'g2s': {'limiting': "ImmutableDenseMatrix([[Mul(Symbol('b', positive=True), Pow(Add(Symbol('a', "
                     "positive=True), Symbol('b', positive=True)), Integer(-1))), Mul(Symbol('a', "
                     "positive=True), Pow(Add(Symbol('a', positive=True), Symbol('b', positive=True)), "
                     'Integer(-1)))]])',
         'rows': [['-a', 'a'], ['b', '-b']]},
 'g3': {'limiting': "ImmutableDenseMatrix([[Mul(Symbol('b', positive=True), Symbol('c', positive=True), "
                    "Pow(Add(Mul(Symbol('a', positive=True), Symbol('b', positive=True)), Mul(Symbol('a', "
                    "positive=True), Symbol('c', positive=True)), Mul(Symbol('b', positive=True), "
                    "Symbol('c', positive=True))), Integer(-1))), Mul(Symbol('a', positive=True), "
                    "Symbol('c', positive=True), Pow(Add(Mul(Symbol('a', positive=True), Symbol('b', "
                    "positive=True)), Mul(Symbol('a', positive=True), Symbol('c', positive=True)), "
                    "Mul(Symbol('b', positive=True), Symbol('c', positive=True))), Integer(-1))), "
                    "Mul(Symbol('a', positive=True), Symbol('b', positive=True), Pow(Add(Mul(Symbol('a', "
                    "positive=True), Symbol('b', positive=True)), Mul(Symbol('a', positive=True), "
                    "Symbol('c', positive=True)), Mul(Symbol('b', positive=True), Symbol('c', "
                    'positive=True))), Integer(-1)))]])',
        'rows': [['-a', 'a', '0'], ['0', '-b', 'b'], ['c', '0', '-c']]},
 'g3b': {'limiting': "ImmutableDenseMatrix([[Mul(Integer(2), Pow(Add(Mul(Integer(2), Symbol('a', "
                     "positive=True)), Mul(Integer(3), Symbol('b', positive=True)), Integer(2)), "
                     "Integer(-1))), Mul(Add(Mul(Integer(2), Symbol('a', positive=True)), Mul(Integer(2), "
                     "Symbol('b', positive=True))), Pow(Add(Mul(Integer(2), Symbol('a', positive=True)), "
                     "Mul(Integer(3), Symbol('b', positive=True)), Integer(2)), Integer(-1))), "
                     "Mul(Symbol('b', positive=True), Pow(Add(Mul(Integer(2), Symbol('a', positive=True)), "
                     "Mul(Integer(3), Symbol('b', positive=True)), Integer(2)), Integer(-1)))]])",
         'rows': [['-a-b', 'a', 'b'], ['1', '-1', '0'], ['0', '2', '-2']]},
 'g3c': {'limiting': "ImmutableDenseMatrix([[Mul(Symbol('b', positive=True), Pow(Add(Mul(Symbol('a', "
                     "positive=True), Symbol('b', positive=True)), Symbol('a', positive=True), Symbol('b', "
                     "positive=True)), Integer(-1))), Mul(Symbol('a', positive=True), Symbol('b', "
                     "positive=True), Pow(Add(Mul(Symbol('a', positive=True), Symbol('b', positive=True)), "
                     "Symbol('a', positive=True), Symbol('b', positive=True)), Integer(-1))), "
                     "Mul(Symbol('a', positive=True), Pow(Add(Mul(Symbol('a', positive=True), Symbol('b', "
                     "positive=True)), Symbol('a', positive=True), Symbol('b', positive=True)), "
                     'Integer(-1)))]])',
         'rows': [['-a', 'a', '0'], ['1', '-2', '1'], ['0', 'b', '-b']]},
 'g3d': {'limiting': "ImmutableDenseMatrix([[Mul(Integer(3), Symbol('a', positive=True), "
                     "Pow(Add(Mul(Integer(4), Symbol('a', positive=True)), Mul(Integer(2), Symbol('b', "
                     'positive=True)), Integer(6)), Integer(-1))), Mul(Integer(3), Pow(Add(Mul(Integer(2), '
                     "Symbol('a', positive=True)), Symbol('b', positive=True), Integer(3)), Integer(-1))), "
                     "Mul(Add(Symbol('a', positive=True), Mul(Integer(2), Symbol('b', positive=True))), "
                     "Pow(Add(Mul(Integer(4), Symbol('a', positive=True)), Mul(Integer(2), Symbol('b', "
                     'positive=True)), Integer(6)), Integer(-1)))]])',
         'rows': [['-2', '1', '1'], ['a', '-a-b', 'b'], ['0', '3', '-3']]},
 'gf': {'limiting': "ImmutableDenseMatrix([[Float('0.33333333333333337', precision=53), "
                    "Float('0.66666666666666663', precision=53)]])",
        'rows': [['-0.5', '0.5'], ['0.25', '-0.25']]},
 'gf2': {'limiting': "ImmutableDenseMatrix([[Float('0.1428571428571429', precision=53), "
                     "Float('0.8571428571428571', precision=53)]])",
         'rows': [['-0.75', '0.75'], ['0.125', '-0.125']]},
 'gf3': {'limiting': 'ERR:IndexError:list index out of range',
         'rows': [['-1.5', '1.0', '0.5'], ['0.25', '-0.5', '0.25'], ['0', '2.0', '-2.0']]},
 'gf3b': {'limiting': "ImmutableDenseMatrix([[Float('0.14285714285714285', precision=53), "
                      "Float('0.2857142857142857', precision=53), Float('0.57142857142857151', "
                      'precision=53)]])',
          'rows': [['-1.0', '0.5', '0.5'], ['0.25', '-0.5', '0.25'], ['0.125', '0.125', '-0.25']]},
 'gf4': {'limiting': 'ERR:IndexError:list index out of range',
         'rows': [['-1.5', '0.5', '0.5', '0.5'],
                  ['0.25', '-0.75', '0.25', '0.25'],
                  ['0', '0.5', '-1.0', '0.5'],
                  ['0.375', '0.125', '0', '-0.5']]},
 'gmix': {'limiting': "ImmutableDenseMatrix([[Float('0.2857142857142857', precision=53), "
                      "Float('0.5714285714285714', precision=53), Float('0.14285714285714285', "
                      'precision=53)]])',
          'rows': [['-1', '1', '0'], ['0.5', '-1.0', '0.5'], ['0', '2', '-2']]},
 'gnonuniq': {'rows': [['-1', '1', '0', '0'],
                       ['1', '-1', '0', '0'],
                       ['0', '0', '-2', '2'],
                       ['0', '0', '3', '-3']]},
 'gp': {'limiting': "ImmutableDenseMatrix([[Mul(Symbol('q'), Pow(Add(Symbol('p'), Symbol('q')), "
                    "Integer(-1))), Mul(Symbol('p'), Pow(Add(Symbol('p'), Symbol('q')), Integer(-1)))]])",
        'rows': [['-p', 'p'], ['q', '-q']]},
 'nonuniq': {'absorbing': 'ImmutableDenseMatrix([[Rational(1, 3), Rational(1, 3), Rational(1, 3)]])',
             'classes': '[([Integer(0), Integer(1)], true, Integer(1)), ([Integer(2)], true, Integer(1)), '
                        '([Integer(3)], false, Integer(1))]',
             'fundamental': 'ImmutableDenseMatrix([[Rational(4, 3)]])',
             'rows': [['1/2', '1/2', '0', '0'],
                      ['1/3', '2/3', '0', '0'],
                      ['0', '0', '1', '0'],
                      ['1/4', '1/4', '1/4', '1/4']],
             'stationary': 'ImmutableDenseMatrix([[Add(Rational(2, 5), Mul(Integer(-1), Rational(2, 5), '
                           "Symbol('tau0'))), Add(Rational(3, 5), Mul(Integer(-1), Rational(3, 5), "
                           "Symbol('tau0'))), Symbol('tau0'), Integer(0)]])"},
 'nonuniq2': {'absorbing': 'ImmutableDenseMatrix([[Rational(1, 2), Rational(1, 2)]])',
              'classes': '[([Integer(0)], true, Integer(1)), ([Integer(1)], true, Integer(1)), '
                         '([Integer(2)], false, Integer(1))]',
              'fundamental': 'ImmutableDenseMatrix([[Rational(3, 2)]])',
              'rows': [['1', '0', '0'], ['0', '1', '0'], ['1/3', '1/3', '1/3']],
              'stationary': "ImmutableDenseMatrix([[Add(Integer(1), Mul(Integer(-1), Symbol('tau0'))), "
                            "Symbol('tau0'), Integer(0)]])"},
 'nonuniqsym': {'absorbing': "ERR:ValueError:The fundamental matrix doesn't exist.",
                'classes': '[([Integer(0), Integer(1)], true, Integer(1)), ([Integer(2)], true, Integer(1))]',
                'fundamental': "ERR:ValueError:The fundamental matrix doesn't exist.",
                'rows': [['p', '1-p', '0'], ['1-p', 'p', '0'], ['0', '0', '1']],
                'stationary': 'ImmutableDenseMatrix([[Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), '
                              "Symbol('tau0'))), Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), "
                              "Symbol('tau0'))), Symbol('tau0')]])"},
 'rf0': {'E2': "Float('0.546875', precision=53)",
         'E2sq': "Float('0.6875', precision=53)",
         'P2': "Float('0.546875', precision=53)",
         'P3': "Float('0.390625', precision=53)",
         'Pand': "Float('0.24609375', precision=53)",
         'Pgt': "Float('0.546875', precision=53)",
         'Pne': "Float('0.546875', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.1944444444444444', precision=53), "
                        "Float('-0.19444444444444436', precision=53)], [Float('-0.13888888888888901', "
                        "precision=53), Float('1.1388888888888891', precision=53)]])",
         'rows': [['0.5625', '0.4375'], ['0.3125', '0.6875']],
         'stationary': "ImmutableDenseMatrix([[Float('0.41666666666666674', precision=53), "
                       "Float('0.58333333333333326', precision=53)]])"},
 'rf1': {'E2': "Float('2.75', precision=53)",
         'E2sq': "Float('5.75', precision=53)",
         'P2': "Float('0.25', precision=53)",
         'P3': "Float('0.125', precision=53)",
         'Pand': 'Integer(0)',
         'Pgt': "Float('1.0', precision=53)",
         'Pne': "Float('1.0', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (5, 5) * (0, 5).',
         'fundamental': "ImmutableDenseMatrix([[Float('0.89867124174996926', precision=53), "
                        "Float('-0.046620332277095985', precision=53), Float('-0.10764867561842406', "
                        "precision=53), Float('0.12408746345477144', precision=53), "
                        "Float('0.13151030269077921', precision=53)], [Float('0.14553316643616182', "
                        "precision=53), Float('0.99940477232541436', precision=53), "
                        "Float('-0.16622608147616463', precision=53), Float('-0.29432257838623266', "
                        "precision=53), Float('0.31561072110082106', precision=53)], "
                        "[Float('0.12879676476252166', precision=53), Float('-0.088461336461196416', "
                        "precision=53), Float('0.85469442061588563', precision=53), "
                        "Float('0.14082386512841161', precision=53), Float('-0.035853714045622459', "
                        "precision=53)], [Float('-0.055303653647520171', precision=53), "
                        "Float('-0.054988533113916077', precision=53), Float('0.08481994362843788', "
                        "precision=53), Float('0.92743474378949942', precision=53), "
                        "Float('0.098037499343498882', precision=53)], [Float('-0.021830850300239839', "
                        "precision=53), Float('0.12074368445930568', precision=53), "
                        "Float('0.042978939444337463', precision=53), Float('0.057141856760210767', "
                        "precision=53), Float('0.8009663696363859', precision=53)]])",
         'rows': [['0.0', '0.0', '0.0', '0.5', '1/2'],
                  ['0.25', '0.0', '0.0', '0.0', '3/4'],
                  ['0.25', '0', '0.0', '0.5', '0.25'],
                  ['0.0', '0', '0.25', '0.25', '0.5'],
                  ['0.0', '0.25', '0.25', '1/2', '0']],
         'stationary': "ImmutableDenseMatrix([[Float('0.062761506276150625', precision=53), "
                       "Float('0.079497907949790794', precision=53), Float('0.17154811715481172', "
                       "precision=53), Float('0.3682008368200837', precision=53), "
                       "Float('0.31799163179916318', precision=53)]])"},
 'rf10': {'E2': "Float('1.3369140625', precision=53)",
          'E2sq': "Float('2.96875', precision=53)",
          'P2': "Float('0.2294921875', precision=53)",
          'P3': "Float('0.291015625', precision=53)",
          'Pand': "Float('0.0751953125', precision=53)",
          'Pgt': "Float('0.689453125', precision=53)",
          'Pne': "Float('0.689453125', precision=53)",
          'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
          'fundamental': "ImmutableDenseMatrix([[Float('0.93165283281101696', precision=53), "
                         "Float('0.00030039992908803808', precision=53), Float('-0.024067854739763531', "
                         "precision=53), Float('0.092114621999658489', precision=53)], "
                         "[Float('0.063068454350508632', precision=53), Float('0.98096018206192004', "
                         "precision=53), Float('0.059164614163947415', precision=53), "
                         "Float('-0.10319325057637614', precision=53)], [Float('-0.036645691563838137', "
                         "precision=53), Float('0.061282604868440609', precision=53), "
                         "Float('0.94917517273828234', precision=53), Float('0.026187913957115148', "
                         "precision=53)], [Float('0.039994106535618482', precision=53), "
                         "Float('-0.025246199437397504', precision=53), Float('0.00065466077619021654', "
                         "precision=53), Float('0.98459743212558892', precision=53)]])",
          'rows': [['0.21875', '0.28125', '0.15625', '0.34375'],
                   ['3/8', '0.25', '0.25', '0.125'],
                   ['1/4', '0.34375', '0.125', '0.28125'],
                   ['0.34375', '0.25', '0.1875', '7/32']],
          'stationary': "ImmutableDenseMatrix([[Float('0.29790631196724271', precision=53), "
                        "Float('0.27655738971440341', precision=53), Float('0.18397671963122247', "
                        "precision=53), Float('0.24155957868713143', precision=53)]])"},
 'rf11': {'E2': "Float('0.630859375', precision=53)",
          'E2sq': "Float('0.65625', precision=53)",
          'P2': "Float('0.630859375', precision=53)",
          'P3': "Float('0.365234375', precision=53)",
          'Pand': "Float('0.2412109375', precision=53)",
          'Pgt': "Float('0.630859375', precision=53)",
          'Pne': "Float('0.630859375', precision=53)",
          'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
          'fundamental': "ImmutableDenseMatrix([[Float('1.0422222222222222', precision=53), "
                         "Float('-0.042222222222222189', precision=53)], [Float('-0.024444444444444477', "
                         "precision=53), Float('1.0244444444444445', precision=53)]])",
          'rows': [['0.40625', '0.59375'], ['0.34375', '0.65625']],
          'stationary': "ImmutableDenseMatrix([[Float('0.3666666666666667', precision=53), "
                        "Float('0.6333333333333333', precision=53)]])"},
 'rf2': {'E2': "Float('2.3916015625', precision=53)",
         'E2sq': "Float('6.71875', precision=53)",
         'P2': "Float('0.1494140625', precision=53)",
         'P3': "Float('0.1796875', precision=53)",
         'Pand': "Float('0.0234375', precision=53)",
         'Pgt': "Float('0.8095703125', precision=53)",
         'Pne': "Float('0.8095703125', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (6, 6) * (0, 6).',
         'fundamental': "ImmutableDenseMatrix([[Float('0.95453819738501899', precision=53), "
                        "Float('-0.085020642922166781', precision=53), Float('0.046212790873952733', "
                        "precision=53), Float('0.14687392132837615', precision=53), "
                        "Float('-0.10952813467321897', precision=53), Float('0.046923868008037661', "
                        "precision=53)], [Float('-0.04804349925339698', precision=53), "
                        "Float('1.0286437321343387', precision=53), Float('-0.027579173560277011', "
                        "precision=53), Float('-0.13694259946154008', precision=53), "
                        "Float('0.16048883611577855', precision=53), Float('0.023432704025096781', "
                        "precision=53)], [Float('0.028192784456531617', precision=53), "
                        "Float('-0.023621601107977622', precision=53), Float('0.97693286362919873', "
                        "precision=53), Float('0.085661780868601239', precision=53), "
                        "Float('0.0075363605306805739', precision=53), Float('-0.074702188377034642', "
                        "precision=53)], [Float('-0.022967567536151978', precision=53), "
                        "Float('0.0961860608497727', precision=53), Float('-0.094740033691017952', "
                        "precision=53), Float('0.98597339263847084', precision=53), "
                        "Float('-0.047132285635487152', precision=53), Float('0.082680433374413609', "
                        "precision=53)], [Float('0.0014853596663788856', precision=53), "
                        "Float('0.031615536677866755', precision=53), Float('-0.011163978042113089', "
                        "precision=53), Float('-0.05565511195965555', precision=53), "
                        "Float('1.0920197881847074', precision=53), Float('-0.058301594527184329', "
                        "precision=53)], [Float('0.10114661312956', precision=53), "
                        "Float('-0.053206528000466401', precision=53), Float('0.12564002583268297', "
                        "precision=53), Float('-0.046177655779566361', precision=53), "
                        "Float('-0.099604840887287452', precision=53), Float('0.97220238570507722', "
                        'precision=53)]])',
         'rows': [['0.125', '1/16', '0.21875', '5/16', '0.09375', '0.1875'],
                  ['0.125', '0.1875', '1/8', '0.0625', '5/16', '0.1875'],
                  ['0.21875', '0.125', '0.15625', '0.25', '0.1875', '0.0625'],
                  ['0.15625', '0.25', '0.0625', '3/16', '0.125', '0.21875'],
                  ['0.1875', '0.1875', '0.15625', '0.125', '0.25', '0.09375'],
                  ['0.28125', '0.125', '0.28125', '0.09375', '0.09375', '1/8']],
         'stationary': "ImmutableDenseMatrix([[Float('0.17987040776010543', precision=53), "
                       "Float('0.1565619471519385', precision=53), Float('0.1644487949303981', "
                       "precision=53), Float('0.17590562893096681', precision=53), "
                       "Float('0.17648835717724934', precision=53), Float('0.14672486404934193', "
                       'precision=53)]])'},
 'rf3': {'E2': "Float('1.4765625', precision=53)",
         'E2sq': "Float('2.6875', precision=53)",
         'P2': "Float('0.22265625', precision=53)",
         'P3': "Float('0.2890625', precision=53)",
         'Pand': "Float('0.09375', precision=53)",
         'Pgt': "Float('0.7578125', precision=53)",
         'Pne': "Float('0.7578125', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
         'fundamental': "ImmutableDenseMatrix([[Float('0.94629382716049371', precision=53), "
                        "Float('-0.13688559670781891', precision=53), Float('0.030864197530864178', "
                        "precision=53), Float('0.15972757201646084', precision=53)], "
                        "[Float('0.14629382716049377', precision=53), Float('0.92978106995884768', "
                        "precision=53), Float('0.030864197530864178', precision=53), "
                        "Float('-0.10693909465020578', precision=53)], [Float('-0.058150617283950642', "
                        "precision=53), Float('0.13570699588477364', precision=53), "
                        "Float('0.91975308641975306', precision=53), Float('0.0026905349794238411', "
                        "precision=53)], [Float('-0.013706172839506202', precision=53), "
                        "Float('0.076447736625514393', precision=53), Float('0.030864197530864182', "
                        "precision=53), Float('0.90639423868312752', precision=53)]])",
         'rows': [['0.25', '0.0625', '0.3125', '0.375'],
                  ['0.4375', '0.1875', '0.3125', '0.0625'],
                  ['3/16', '0.375', '0.1875', '0.25'],
                  ['0.25', '0.3125', '0.3125', '0.125']],
         'stationary': "ImmutableDenseMatrix([[Float('0.27611111111111114', precision=53), "
                       "Float('0.23185185185185186', precision=53), Float('0.27777777777777779', "
                       "precision=53), Float('0.21425925925925929', precision=53)]])"},
 'rf4': {'E2': "Float('2.568359375', precision=53)",
         'E2sq': "Float('10.25', precision=53)",
         'P2': "Float('0.1982421875', precision=53)",
         'P3': "Float('0.150390625', precision=53)",
         'Pand': "Float('0.0390625', precision=53)",
         'Pgt': "Float('0.8505859375', precision=53)",
         'Pne': "Float('0.8505859375', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (6, 6) * (0, 6).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.0003449717897777', precision=53), "
                        "Float('0.010550033633742321', precision=53), Float('-0.047642691177425589', "
                        "precision=53), Float('-0.020887380132470983', precision=53), "
                        "Float('-0.0034897851129662115', precision=53), Float('0.061124850999342749', "
                        "precision=53)], [Float('-0.02501209200923948', precision=53), "
                        "Float('1.0145583245602519', precision=53), Float('0.013653995970423429', "
                        "precision=53), Float('-0.027350983157583105', precision=53), "
                        "Float('-0.0013722760450611739', precision=53), Float('0.025523030681208338', "
                        "precision=53)], [Float('0.095097223007943363', precision=53), "
                        "Float('0.045129956395024123', precision=53), Float('0.94924261025695145', "
                        "precision=53), Float('-0.092268742082975003', precision=53), "
                        "Float('0.029893548117811333', precision=53), Float('-0.027094595694755225', "
                        "precision=53)], [Float('-0.057034241537806243', precision=53), "
                        "Float('0.022028205104832711', precision=53), Float('0.047775478523393894', "
                        "precision=53), Float('1.0655371507102396', precision=53), "
                        "Float('-0.034458733898650312', precision=53), Float('-0.043847858902009623', "
                        "precision=53)], [Float('-0.054192720070431556', precision=53), "
                        "Float('-0.041260778366801081', precision=53), Float('0.079056449389184186', "
                        "precision=53), Float('0.066439420869672763', precision=53), "
                        "Float('0.99774060049842361', precision=53), Float('-0.047782972320047752', "
                        "precision=53)], [Float('0.0012755157946574159', precision=53), "
                        "Float('-0.048555430492001811', precision=53), Float('-0.015366077039840198', "
                        "precision=53), Float('0.043185192770109064', precision=53), "
                        "Float('-0.0044834088243722939', precision=53), Float('1.0239442077914478', "
                        'precision=53)]])',
         'rows': [['0.15625', '0.1875', '0.15625', '0.09375', '5/32', '0.25'],
                  ['1/8', '0.1875', '0.21875', '0.09375', '5/32', '0.21875'],
                  ['0.25', '0.21875', '0.15625', '0.03125', '3/16', '0.15625'],
                  ['3/32', '0.1875', '0.25', '0.1875', '0.125', '0.15625'],
                  ['0.09375', '0.125', '0.28125', '0.1875', '5/32', '5/32'],
                  ['0.15625', '1/8', '0.1875', '0.15625', '0.15625', '0.21875']],
         'stationary': "ImmutableDenseMatrix([[Float('0.15263241018912407', precision=53), "
                       "Float('0.1718589406490646', precision=53), Float('0.20406780484672402', "
                       "precision=53), Float('0.11914999227675926', precision=53), "
                       "Float('0.15890368164281138', precision=53), Float('0.19338717039551673', "
                       'precision=53)]])'},
 'rf5': {'E2': "Float('2.5390625', precision=53)",
         'E2sq': "Float('10.0625', precision=53)",
         'P2': "Float('0.16015625', precision=53)",
         'P3': "Float('0.1015625', precision=53)",
         'Pand': "Float('0.03515625', precision=53)",
         'Pgt': "Float('0.89453125', precision=53)",
         'Pne': "Float('0.89453125', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (6, 6) * (0, 6).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.0809236154038993', precision=53), "
                        "Float('0.12333920639383088', precision=53), Float('-0.18432653332870633', "
                        "precision=53), Float('-0.10678828094154243', precision=53), "
                        "Float('-0.019215588586591677', precision=53), Float('0.10606758105911059', "
                        "precision=53)], [Float('-0.037340211974395542', precision=53), "
                        "Float('0.90256721081527336', precision=53), Float('0.10633695146099453', "
                        "precision=53), Float('-0.00098599164938805578', precision=53), "
                        "Float('-0.024486182333617816', precision=53), Float('0.053908223681133363', "
                        "precision=53)], [Float('0.034496028725812602', precision=53), "
                        "Float('-0.12420301325816351', precision=53), Float('1.1286174104906685', "
                        "precision=53), Float('0.12935383964281025', precision=53), "
                        "Float('-0.075259568763302989', precision=53), Float('-0.093004696837824904', "
                        "precision=53)], [Float('0.019621241928649936', precision=53), "
                        "Float('-0.081159830990783347', precision=53), Float('0.086503902398110996', "
                        "precision=53), Float('0.95888673327197993', precision=53), "
                        "Float('0.10481905092675685', precision=53), Float('-0.088671097534714513', "
                        "precision=53)], [Float('-0.059710954322884278', precision=53), "
                        "Float('0.14840380776857742', precision=53), Float('-0.13810147409552892', "
                        "precision=53), Float('0.049728832627552372', precision=53), "
                        "Float('1.0652383559404093', precision=53), Float('-0.065558567918125837', "
                        "precision=53)], [Float('0.015248601190376399', precision=53), "
                        "Float('0.12046966090933885', precision=53), Float('-0.11485229945586911', "
                        "precision=53), Float('-0.16209047385015371', precision=53), "
                        "Float('-0.03122082989926233', precision=53), Float('1.1724453411055697', "
                        'precision=53)]])',
         'rows': [['0.1875', '0.3125', '0.0625', '0.0625', '0.1875', '0.1875'],
                  ['0.0625', '0.125', '0.3125', '0.125', '0.1875', '0.1875'],
                  ['0.125', '0.125', '0.3125', '0.25', '1/8', '1/16'],
                  ['1/8', '1/8', '0.3125', '1/16', '0.3125', '1/16'],
                  ['0.0625', '3/8', '0.0625', '0.1875', '0.25', '0.0625'],
                  ['1/8', '0.3125', '0.125', '0.0', '0.1875', '0.25']],
         'stationary': "ImmutableDenseMatrix([[Float('0.10512638444589234', precision=53), "
                       "Float('0.21941774579633547', precision=53), Float('0.21151185517579627', "
                       "precision=53), Float('0.13340458395250607', precision=53), "
                       "Float('0.20368648751528104', precision=53), Float('0.12685294311418888', "
                       'precision=53)]])'},
 'rf6': {'E2': "Float('0.8125', precision=53)",
         'E2sq': "Float('1.5', precision=53)",
         'P2': "Float('0.25', precision=53)",
         'P3': "Float('0.4375', precision=53)",
         'Pand': "Float('0.125', precision=53)",
         'Pgt': "Float('0.5625', precision=53)",
         'Pne': "Float('0.5625', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.0625', precision=53), Float('-0.0625', "
                        "precision=53), Integer(0)], [Float('0.0625', precision=53), Float('0.9375', "
                        "precision=53), Integer(0)], [Float('-0.1875', precision=53), Float('0.1875', "
                        "precision=53), Float('1.0', precision=53)]])",
         'rows': [['0.5', '0.25', '0.25'], ['0.5', '0.25', '0.25'], ['0.25', '0.5', '0.25']],
         'stationary': "ImmutableDenseMatrix([[Float('0.4375', precision=53), Float('0.3125', precision=53), "
                       "Float('0.25', precision=53)]])"},
 'rf7': {'E2': "Float('0.5625', precision=53)",
         'E2sq': "Float('0.5', precision=53)",
         'P2': "Float('0.5625', precision=53)",
         'P3': "Float('0.375', precision=53)",
         'Pand': "Float('0.1875', precision=53)",
         'Pgt': "Float('0.5625', precision=53)",
         'Pne': "Float('0.5625', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
         'fundamental': "ImmutableDenseMatrix([[Float('0.88000000000000012', precision=53), "
                        "Float('0.11999999999999993', precision=53)], [Float('0.080000000000000057', "
                        "precision=53), Float('0.91999999999999982', precision=53)]])",
         'rows': [['0.25', '0.75'], ['0.5', '0.5']],
         'stationary': "ImmutableDenseMatrix([[Float('0.39999999999999991', precision=53), "
                       "Float('0.60000000000000009', precision=53)]])"},
 'rf8': {'E2': "Float('1.013671875', precision=53)",
         'E2sq': "Float('1.5625', precision=53)",
         'P2': "Float('0.3427734375', precision=53)",
         'P3': "Float('0.32421875', precision=53)",
         'Pand': "Float('0.12890625', precision=53)",
         'Pgt': "Float('0.6708984375', precision=53)",
         'Pne': "Float('0.6708984375', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.0270631539656621', precision=53), "
                        "Float('-0.058790432598663914', precision=53), Float('0.031727278633001836', "
                        "precision=53)], [Float('-0.081016594198031477', precision=53), "
                        "Float('1.0828674897518082', precision=53), Float('-0.0018508955537767555', "
                        "precision=53)], [Float('0.053296102549082902', precision=53), "
                        "Float('-0.025212258411885315', precision=53), Float('0.97191615586280244', "
                        'precision=53)]])',
         'rows': [['0.34375', '0.28125', '3/8'], ['0.25', '0.40625', '11/32'], ['0.375', '0.3125', '5/16']],
         'stationary': "ImmutableDenseMatrix([[Float('0.32318992654774392', precision=53), "
                       "Float('0.33368310598111223', precision=53), Float('0.3431269674711438', "
                       'precision=53)]])'},
 'rf9': {'E2': "Float('0.375', precision=53)",
         'E2sq': 'Rational(3, 4)',
         'P2': "Float('0.375', precision=53)",
         'P3': "Float('0.375', precision=53)",
         'Pand': "Float('0.1875', precision=53)",
         'Pgt': "Float('0.375', precision=53)",
         'Pne': "Float('0.375', precision=53)",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
         'fundamental': "ImmutableDenseMatrix([[Float('1.5', precision=53), Float('-0.5', precision=53)], "
                        "[Float('-0.5', precision=53), Float('1.5', precision=53)]])",
         'rows': [['0.75', '0.25'], ['0.25', '3/4']],
         'stationary': "ImmutableDenseMatrix([[Float('0.5', precision=53), Float('0.5', precision=53)]])"},
 's2b': {'E2': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Symbol('p'), "
               "Symbol('q')), Mul(Integer(2), Symbol('p')))",
         'E2sq': "Add(Integer(1), Mul(Integer(-1), Symbol('q')))",
         'P2': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Symbol('p'), "
               "Symbol('q')), Mul(Integer(2), Symbol('p')))",
         'P3': "Add(Mul(Integer(-1), Symbol('p'), Symbol('q')), Mul(Integer(-1), Pow(Symbol('q'), "
               "Integer(2))), Mul(Integer(2), Symbol('q')))",
         'Pand': "Mul(Symbol('p'), Add(Integer(1), Mul(Integer(-1), Symbol('p'))))",
         'Pgt': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Symbol('p'), "
                "Symbol('q')), Mul(Integer(2), Symbol('p')))",
         'Pne': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Symbol('p'), "
                "Symbol('q')), Mul(Integer(2), Symbol('p')))",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
         'fundamental': "ImmutableDenseMatrix([[Mul(Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(2), "
                        "Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2))), Integer(-1)), "
                        "Add(Mul(Symbol('p'), Symbol('q')), Symbol('p'), Pow(Symbol('q'), Integer(2)))), "
                        "Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Symbol('p'), Symbol('q')), "
                        "Mul(Integer(-1), Symbol('p'))), Pow(Add(Pow(Symbol('p'), Integer(2)), "
                        "Mul(Integer(2), Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2))), "
                        "Integer(-1)))], [Mul(Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(2), "
                        "Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2))), Integer(-1)), "
                        "Add(Mul(Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), "
                        "Symbol('q')))), Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Symbol('p'), "
                        "Symbol('q')), Symbol('q')), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(2), "
                        "Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2))), Integer(-1)))]])",
         'rows': [['1-p', 'p'], ['q', '1-q']],
         'stationary': "ImmutableDenseMatrix([[Mul(Symbol('q'), Pow(Add(Symbol('p'), Symbol('q')), "
                       "Integer(-1))), Mul(Symbol('p'), Pow(Add(Symbol('p'), Symbol('q')), Integer(-1)))]])"},
 's3c': {'E2': "Add(Rational(5, 4), Mul(Integer(-1), Symbol('p')))",
         'E2sq': "Add(Integer(4), Mul(Integer(-1), Integer(3), Symbol('q')))",
         'P2': "Add(Rational(1, 2), Mul(Integer(-1), Rational(1, 2), Symbol('p')))",
         'P3': "Mul(Symbol('p'), Symbol('q'))",
         'Pand': 'Integer(0)',
         'Pgt': "Add(Rational(3, 4), Mul(Integer(-1), Rational(1, 2), Symbol('p')))",
         'Pne': "Add(Rational(3, 4), Mul(Integer(-1), Rational(1, 2), Symbol('p')))",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                        "Integer(2)), Symbol('q')), Mul(Integer(2), Symbol('p'), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(2), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(-1), Integer(2), Symbol('p')), Mul(Integer(2), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(4), Symbol('q')), Integer(2)), Pow(Add(Mul(Integer(4), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                        "Integer(4), Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), "
                        "Integer(2)), Mul(Integer(4), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), "
                        "Symbol('p')), Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), "
                        "Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(2), Symbol('p'), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), "
                        "Mul(Integer(-1), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Symbol('q')), "
                        "Integer(1)), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('q')), Pow(Symbol('p'), Integer(2)), Mul(Integer(4), Symbol('p'), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(-1), Integer(2), Symbol('p')), Pow(Symbol('q'), Integer(2)), "
                        "Mul(Integer(2), Symbol('q')), Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                        "Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), "
                        "Integer(2)), Mul(Integer(3), Symbol('p'), Symbol('q')), Symbol('p'), "
                        "Mul(Integer(-1), Symbol('q')), Integer(-2)), Pow(Add(Mul(Integer(4), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                        "Integer(4), Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), "
                        "Integer(2)), Mul(Integer(4), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), "
                        "Symbol('p')), Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), "
                        "Integer(1)), Integer(-1)))], [Mul(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('q')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2))), "
                        "Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(-1), Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(2), "
                        "Symbol('p'))), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(4), Pow(Symbol('p'), "
                        "Integer(2)), Symbol('q')), Pow(Symbol('p'), Integer(2)), Mul(Integer(4), "
                        "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), "
                        "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Pow(Symbol('q'), "
                        "Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), Integer(-1))), "
                        "Mul(Add(Mul(Integer(6), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(-1), Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), "
                        "Pow(Symbol('q'), Integer(2)), Symbol('q'), Integer(1)), Pow(Add(Mul(Integer(4), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                        "Integer(4), Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), "
                        "Integer(2)), Mul(Integer(4), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), "
                        "Symbol('p')), Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), "
                        "Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(6), Pow(Symbol('p'), "
                        "Integer(2)), Symbol('q')), Mul(Integer(3), Pow(Symbol('p'), Integer(2))), "
                        "Mul(Integer(5), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(3), "
                        "Symbol('p')), Symbol('q')), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(4), Pow(Symbol('p'), "
                        "Integer(2)), Symbol('q')), Pow(Symbol('p'), Integer(2)), Mul(Integer(4), "
                        "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), "
                        "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Pow(Symbol('q'), "
                        "Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), Integer(-1)))], "
                        "[Mul(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(4), Symbol('p'), Symbol('q'))), "
                        "Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(4), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('q')), Pow(Symbol('p'), Integer(2)), Mul(Integer(4), Symbol('p'), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(-1), Integer(2), Symbol('p')), Pow(Symbol('q'), Integer(2)), "
                        "Mul(Integer(2), Symbol('q')), Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(6), "
                        "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(3), "
                        "Symbol('p'), Symbol('q')), Pow(Symbol('q'), Integer(2))), Pow(Add(Mul(Integer(4), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                        "Integer(4), Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), "
                        "Integer(2)), Mul(Integer(4), Symbol('p'), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), "
                        "Symbol('p')), Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), "
                        "Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                        "Integer(2)), Symbol('q')), Pow(Symbol('p'), Integer(2)), Mul(Integer(9), "
                        "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                        "Mul(Integer(2), Symbol('q')), Integer(1)), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), "
                        "Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(4), "
                        "Pow(Symbol('p'), Integer(2)), Symbol('q')), Pow(Symbol('p'), Integer(2)), "
                        "Mul(Integer(4), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), "
                        "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                        "Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), "
                        'Integer(-1)))]])',
         'rows': [['1/2', '1/2', '0'], ['p', '0', '1-p'], ['0', 'q', '1-q']],
         'stationary': "ImmutableDenseMatrix([[Mul(Integer(2), Symbol('p'), Symbol('q'), "
                       "Pow(Add(Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), "
                       "Symbol('q'), Integer(1)), Integer(-1))), Mul(Symbol('q'), Pow(Add(Mul(Integer(2), "
                       "Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), Symbol('q'), Integer(1)), "
                       "Integer(-1))), Mul(Add(Integer(1), Mul(Integer(-1), Symbol('p'))), "
                       "Pow(Add(Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), "
                       "Symbol('q'), Integer(1)), Integer(-1)))]])"},
 's3d': {'E2': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(3), "
               "Symbol('p'), Symbol('q')), Mul(Integer(2), Symbol('p')), Mul(Integer(-1), Integer(2), "
               "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(2), Symbol('q'), Symbol('r')), "
               "Mul(Integer(4), Symbol('q')))",
         'E2sq': "Add(Integer(4), Mul(Integer(-1), Integer(4), Symbol('r')))",
         'P2': "Add(Mul(Integer(-1), Symbol('p'), Symbol('q')), Mul(Rational(1, 3), Symbol('p')), "
               "Mul(Integer(-1), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Symbol('q'), Symbol('r')), "
               "Mul(Integer(2), Symbol('q')))",
         'P3': "Add(Mul(Integer(-1), Symbol('p'), Symbol('r')), Mul(Integer(-1), Symbol('q'), Symbol('r')), "
               "Mul(Integer(-1), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), Symbol('r')))",
         'Pand': "Mul(Symbol('q'), Add(Mul(Integer(-1), Symbol('p')), Mul(Integer(-1), Symbol('q')), "
                 'Integer(1)))',
         'Pgt': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(2), "
                "Symbol('p'), Symbol('q')), Mul(Rational(5, 3), Symbol('p')), Mul(Integer(-1), "
                "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Symbol('q'), Symbol('r')), Mul(Integer(2), "
                "Symbol('q')))",
         'Pne': "Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(2), "
                "Symbol('p'), Symbol('q')), Mul(Rational(5, 3), Symbol('p')), Mul(Integer(-1), "
                "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Symbol('q'), Symbol('r')), Mul(Integer(2), "
                "Symbol('q')))",
         'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
         'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(15), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(5), Symbol('p'), Symbol('r')), Mul(Integer(2), "
                        "Symbol('p')), Mul(Integer(4), Symbol('q'), Symbol('r')), Mul(Integer(4), "
                        "Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), "
                        "Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                        "Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(12), "
                        "Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('p'), "
                        "Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), Mul(Integer(8), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), "
                        "Integer(-1))), Mul(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('r'), Integer(2))), Mul(Integer(3), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Mul(Integer(3), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), "
                        "Integer(3), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(6), Symbol('p'), "
                        "Symbol('q')), Mul(Integer(-1), Integer(3), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2)))), Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('r'), Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), "
                        "Symbol('p'), Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(8), Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), "
                        "Integer(2)))), Integer(-1))), Mul(Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2))), "
                        "Mul(Integer(15), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                        "Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p'), Symbol('r')), "
                        "Mul(Integer(-1), Integer(2), Symbol('p')), Mul(Integer(4), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(4), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                        "Integer(4), Symbol('q'))), Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('r'), Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), "
                        "Symbol('p'), Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(8), Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), "
                        "Integer(2)))), Integer(-1)))], [Mul(Add(Mul(Integer(6), Symbol('p'), "
                        "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Symbol('p'), Symbol('r')), "
                        "Symbol('p'), Mul(Integer(-1), Integer(2), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(2), Symbol('q')), Mul(Integer(-1), Integer(2), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(-1), Integer(2), Symbol('r'))), Pow(Add(Mul(Integer(9), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(6), "
                        "Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), Integer(2)), "
                        "Mul(Integer(12), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(4), "
                        "Symbol('p'), Symbol('q')), Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(4), Symbol('p'), Symbol('r')), Mul(Integer(4), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(8), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), Integer(-1))), "
                        "Mul(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(3), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                        "Mul(Integer(3), Pow(Symbol('p'), Integer(2))), Mul(Integer(6), Symbol('p'), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(9), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(6), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(3), "
                        "Symbol('p'), Symbol('r')), Mul(Integer(6), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(12), Symbol('q'), Symbol('r')), Mul(Integer(6), Pow(Symbol('r'), "
                        "Integer(2)))), Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), "
                        "Pow(Symbol('r'), Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                        "Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), "
                        "Symbol('p'), Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(8), Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), "
                        "Integer(2)))), Integer(-1))), Mul(Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2))), "
                        "Mul(Integer(6), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                        "Integer(5), Symbol('p'), Symbol('q')), Mul(Integer(2), Symbol('p'), Symbol('r')), "
                        "Mul(Integer(-1), Symbol('p')), Mul(Integer(-1), Integer(2), Pow(Symbol('q'), "
                        "Integer(2))), Mul(Integer(-1), Integer(2), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(2), Symbol('r'))), "
                        "Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                        "Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(12), "
                        "Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('p'), "
                        "Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), Mul(Integer(8), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), "
                        "Integer(-1)))], [Mul(Add(Mul(Integer(15), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(-1), Symbol('p'), Symbol('r')), Mul(Integer(4), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                        "Mul(Integer(-1), Integer(4), Symbol('r'))), Pow(Add(Mul(Integer(9), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(6), "
                        "Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), Integer(2)), "
                        "Mul(Integer(12), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(4), "
                        "Symbol('p'), Symbol('q')), Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(4), Symbol('p'), Symbol('r')), Mul(Integer(4), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(8), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), Integer(-1))), "
                        "Mul(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(-1), Integer(6), Pow(Symbol('p'), Integer(2)), "
                        "Symbol('r')), Mul(Integer(-1), Integer(3), Symbol('p'), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(-1), Integer(3), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                        "Mul(Integer(-1), Integer(6), Symbol('p'), Symbol('r'))), Pow(Add(Mul(Integer(9), "
                        "Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(6), "
                        "Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), Integer(2)), "
                        "Mul(Integer(12), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(4), "
                        "Symbol('p'), Symbol('q')), Mul(Integer(12), Symbol('p'), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(4), Symbol('p'), Symbol('r')), Mul(Integer(4), "
                        "Pow(Symbol('q'), Integer(2))), Mul(Integer(8), Symbol('q'), Symbol('r')), "
                        "Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), Integer(-1))), "
                        "Mul(Add(Mul(Integer(12), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                        "Pow(Symbol('p'), Integer(2)), Mul(Integer(15), Symbol('p'), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(11), "
                        "Symbol('p'), Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), "
                        "Mul(Integer(4), Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('r'))), "
                        "Pow(Add(Mul(Integer(9), Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                        "Integer(2))), Mul(Integer(6), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                        "Pow(Symbol('p'), Integer(2)), Mul(Integer(12), Symbol('p'), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(12), "
                        "Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('p'), "
                        "Symbol('r')), Mul(Integer(4), Pow(Symbol('q'), Integer(2))), Mul(Integer(8), "
                        "Symbol('q'), Symbol('r')), Mul(Integer(4), Pow(Symbol('r'), Integer(2)))), "
                        'Integer(-1)))]])',
         'rows': [['1-p-q', 'p', 'q'], ['1/3', '1/3', '1/3'], ['r', '0', '1-r']],
         'stationary': "ImmutableDenseMatrix([[Mul(Integer(2), Symbol('r'), Pow(Add(Mul(Integer(3), "
                       "Symbol('p'), Symbol('r')), Symbol('p'), Mul(Integer(2), Symbol('q')), "
                       "Mul(Integer(2), Symbol('r'))), Integer(-1))), Mul(Integer(3), Symbol('p'), "
                       "Symbol('r'), Pow(Add(Mul(Integer(3), Symbol('p'), Symbol('r')), Symbol('p'), "
                       "Mul(Integer(2), Symbol('q')), Mul(Integer(2), Symbol('r'))), Integer(-1))), "
                       "Mul(Add(Symbol('p'), Mul(Integer(2), Symbol('q'))), Pow(Add(Mul(Integer(3), "
                       "Symbol('p'), Symbol('r')), Symbol('p'), Mul(Integer(2), Symbol('q')), "
                       "Mul(Integer(2), Symbol('r'))), Integer(-1)))]])"},
 's4': {'E2': "Add(Integer(2), Mul(Integer(-1), Integer(2), Symbol('p')))",
        'E2sq': 'Integer(4)',
        'P2': 'Integer(0)',
        'P3': 'Integer(0)',
        'Pand': 'Integer(0)',
        'Pgt': "Add(Integer(1), Mul(Integer(-1), Symbol('p')))",
        'Pne': "Add(Integer(1), Mul(Integer(-1), Symbol('p')))",
        'absorbing': 'ERR:ShapeError:Matrix size mismatch: (4, 4) * (0, 4).',
        'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), "
                       "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(3), Pow(Symbol('p'), "
                       "Integer(2)), Symbol('q')), Mul(Integer(3), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(4), Symbol('p')), Integer(4)), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(3), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(-1), "
                       "Integer(3), Symbol('q')), Integer(4)), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), "
                       "Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(8), "
                       "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), "
                       "Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                       "Integer(8), Symbol('p')), Integer(4)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                       "Integer(3), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(3), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(3), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(2), Symbol('p')), Integer(-1)), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), "
                       "Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), "
                       "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), "
                       "Integer(-1)), Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Symbol('p'), Symbol('q')), Mul(Integer(2), "
                       "Symbol('p')), Mul(Integer(3), Symbol('q')), Integer(-3)))], "
                       "[Mul(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(3), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(-1), "
                       "Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(3), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(4), Symbol('p'))), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Integer(3), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(-1), "
                       "Symbol('q')), Integer(4)), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), "
                       "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(8), Pow(Symbol('p'), "
                       "Integer(2)), Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), Integer(2))), "
                       "Mul(Integer(8), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(8), "
                       "Symbol('p')), Integer(4)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(5), "
                       "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(5), Pow(Symbol('p'), "
                       "Integer(2))), Mul(Integer(5), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                       "Integer(6), Symbol('p')), Integer(1)), Pow(Add(Mul(Integer(4), Pow(Symbol('p'), "
                       "Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(8), "
                       "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), "
                       "Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                       "Integer(8), Symbol('p')), Integer(4)), Integer(-1))), Mul(Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1)), "
                       "Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(6), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(3), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), Integer(3), "
                       "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(5), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Symbol('q'), "
                       "Integer(-1)))], [Mul(Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), "
                       "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(3), Pow(Symbol('p'), "
                       "Integer(2)), Symbol('q')), Mul(Integer(-1), Symbol('p'), Symbol('q'))), "
                       "Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Integer(5), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(5), Symbol('p'), Symbol('q')), Symbol('q')), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Integer(-1), Integer(3), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(3), Pow(Symbol('p'), Integer(2))), Mul(Integer(7), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(6), Symbol('p')), Integer(3)), "
                       "Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), "
                       "Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), "
                       "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), "
                       "Integer(-1)), Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(5), Symbol('p'), "
                       "Pow(Symbol('q'), Integer(2))), Mul(Integer(7), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(2), Symbol('p')), Mul(Integer(-1), Symbol('q')), "
                       "Integer(1)))], [Mul(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(-1), Integer(3), Symbol('p'), Symbol('q'))), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Integer(3), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(3), Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('q'))), "
                       "Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                       "Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                       "Pow(Symbol('p'), Integer(2)), Mul(Integer(5), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(2), Symbol('p')), Integer(1)), Pow(Add(Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), "
                       "Integer(8), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Integer(4), "
                       "Pow(Symbol('p'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q')), "
                       "Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), Integer(-1))), "
                       "Mul(Pow(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Integer(8), Pow(Symbol('p'), Integer(2)), "
                       "Symbol('q')), Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(8), "
                       "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(8), Symbol('p')), Integer(4)), "
                       "Integer(-1)), Add(Mul(Integer(3), Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                       "Integer(2))), Mul(Integer(-1), Integer(6), Pow(Symbol('p'), Integer(2)), "
                       "Symbol('q')), Mul(Integer(3), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), "
                       "Integer(3), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(9), Symbol('p'), "
                       "Symbol('q')), Mul(Integer(-1), Integer(6), Symbol('p')), Symbol('q'), "
                       'Integer(3)))]])',
        'rows': [['0', '1', '0', '0'], ['p', '0', '1-p', '0'], ['0', 'q', '0', '1-q'], ['0', '0', '1', '0']],
        'stationary': "ImmutableDenseMatrix([[Mul(Symbol('p'), Symbol('q'), Pow(Add(Mul(Integer(2), "
                      "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Integer(2)), "
                      "Integer(-1))), Mul(Symbol('q'), Pow(Add(Mul(Integer(2), Symbol('p'), Symbol('q')), "
                      "Mul(Integer(-1), Integer(2), Symbol('p')), Integer(2)), Integer(-1))), "
                      "Mul(Add(Integer(1), Mul(Integer(-1), Symbol('p'))), Pow(Add(Mul(Integer(2), "
                      "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Integer(2)), "
                      "Integer(-1))), Mul(Pow(Add(Mul(Integer(2), Symbol('p'), Symbol('q')), "
                      "Mul(Integer(-1), Integer(2), Symbol('p')), Integer(2)), Integer(-1)), "
                      "Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p')), Mul(Integer(-1), "
                      "Symbol('q')), Integer(1)))]])"},
 'sym2': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (2, 2) * (0, 2).',
          'classes': '[([Integer(0), Integer(1)], true, Integer(1))]',
          'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(-1), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Symbol('p')), Pow(Symbol('q'), Integer(2)), Symbol('q'), "
                         "Integer(1)), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), "
                         "Integer(-1))), Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Symbol('p'), "
                         "Symbol('q')), Mul(Integer(-1), Symbol('p')), Symbol('q')), "
                         "Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Pow(Symbol('q'), "
                         "Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), Integer(-1)))], "
                         "[Mul(Add(Mul(Integer(-1), Symbol('p'), Symbol('q')), Pow(Symbol('q'), "
                         "Integer(2))), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), "
                         "Integer(-1))), Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Symbol('p'), "
                         "Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), Mul(Integer(2), "
                         "Symbol('q')), Integer(1)), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), "
                         "Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(2), Symbol('q')), Integer(1)), "
                         'Integer(-1)))]])',
          'rows': [['p', '1-p'], ['q', '1-q']],
          'stationary': "ImmutableDenseMatrix([[Mul(Integer(-1), Symbol('q'), Pow(Add(Symbol('p'), "
                        "Mul(Integer(-1), Symbol('q')), Integer(-1)), Integer(-1))), Mul(Add(Symbol('p'), "
                        "Integer(-1)), Pow(Add(Symbol('p'), Mul(Integer(-1), Symbol('q')), Integer(-1)), "
                        'Integer(-1)))]])'},
 'sym3': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
          'classes': '[([Integer(0), Integer(1), Integer(2)], true, Integer(1))]',
          'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2))), Mul(Symbol('p'), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(3), Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(2), "
                         "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Symbol('r')), Mul(Integer(-1), Symbol('p')), Mul(Pow(Symbol('q'), Integer(2)), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), Integer(3), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(3), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(3), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(2), Symbol('r')), Integer(1)), "
                         "Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), "
                         "Mul(Integer(-1), Integer(4), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(4), Symbol('r')), Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                         "Pow(Symbol('p'), Integer(2)), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('p'), Integer(2)), Mul(Symbol('p'), Symbol('q'), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Integer(-1), Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(-1), Symbol('q')), Pow(Symbol('r'), "
                         "Integer(2)), Symbol('r'), Integer(1)), Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                         "Integer(2)), Symbol('q'), Symbol('r')), Mul(Integer(-1), Integer(2), "
                         "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('r')), Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), "
                         "Symbol('p'), Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(8), "
                         "Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p')), Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), Integer(4), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(6), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(4), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('r')), Integer(1)), "
                         "Integer(-1))), Mul(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                         "Integer(2))), Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                         "Mul(Pow(Symbol('p'), Integer(2)), Symbol('r')), Mul(Integer(-1), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(3), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), Symbol('r')), Symbol('p'), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), "
                         "Symbol('q'), Symbol('r')), Symbol('q'), Symbol('r'), Integer(-1)), "
                         "Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), "
                         "Mul(Integer(-1), Integer(4), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(4), Symbol('r')), Integer(1)), Integer(-1)))], "
                         "[Mul(Add(Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Mul(Symbol('p'), Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Symbol('p'), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(-1), Symbol('p'), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Symbol('q'), Symbol('r')), Pow(Symbol('r'), Integer(2)), Mul(Integer(-1), "
                         "Symbol('r'))), Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                         "Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                         "Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), "
                         "Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p')), Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), Integer(4), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(6), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(4), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('r')), Integer(1)), "
                         "Integer(-1))), Mul(Add(Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                         "Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), Integer(2)), "
                         "Mul(Symbol('p'), Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(3), "
                         "Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(2), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Integer(3), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(3), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p')), Mul(Integer(-1), Integer(2), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(2), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Symbol('q')), Mul(Integer(3), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Symbol('r')), Integer(1)), "
                         "Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), "
                         "Mul(Integer(-1), Integer(4), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(4), Symbol('r')), Integer(1)), Integer(-1))), "
                         "Mul(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(4), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(3), Symbol('p'), Symbol('r')), Mul(Integer(2), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), "
                         "Integer(5), Symbol('q'), Symbol('r')), Mul(Integer(-1), Symbol('q')), "
                         "Mul(Integer(3), Symbol('r'))), Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                         "Integer(2)), Symbol('q'), Symbol('r')), Mul(Integer(-1), Integer(2), "
                         "Pow(Symbol('p'), Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('r')), Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), "
                         "Symbol('p'), Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(8), "
                         "Symbol('p'), Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p')), Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), Integer(4), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(6), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(4), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('r')), Integer(1)), "
                         "Integer(-1)))], [Mul(Add(Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Symbol('p'), Symbol('q'), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Symbol('p'), Symbol('r')), Mul(Pow(Symbol('q'), Integer(2)), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(3), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(3), Pow(Symbol('r'), Integer(2)))), "
                         "Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), "
                         "Mul(Integer(-1), Integer(4), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(4), Symbol('r')), Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), "
                         "Pow(Symbol('p'), Integer(2)), Symbol('q'), Symbol('r')), Mul(Pow(Symbol('p'), "
                         "Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Symbol('p'), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Symbol('p'), Symbol('r')), Mul(Integer(-1), Symbol('q'), Pow(Symbol('r'), "
                         "Integer(2))), Pow(Symbol('r'), Integer(2)), Mul(Integer(-1), Symbol('r'))), "
                         "Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                         "Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), "
                         "Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), Integer(2))), Mul(Integer(2), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), Integer(2)), "
                         "Mul(Integer(-1), Integer(4), Symbol('q'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('q')), Mul(Integer(4), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(4), Symbol('r')), Integer(1)), Integer(-1))), "
                         "Mul(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), Integer(2))), "
                         "Mul(Integer(-1), Pow(Symbol('p'), Integer(2)), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('q')), "
                         "Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), Pow(Symbol('p'), "
                         "Integer(2)), Mul(Integer(-1), Symbol('p'), Pow(Symbol('q'), Integer(2)), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p'), Pow(Symbol('q'), "
                         "Integer(2))), Mul(Integer(5), Symbol('p'), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(4), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(6), "
                         "Symbol('p'), Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('p')), "
                         "Mul(Integer(2), Pow(Symbol('q'), Integer(2)), Symbol('r')), Pow(Symbol('q'), "
                         "Integer(2)), Mul(Integer(-1), Integer(6), Symbol('q'), Symbol('r')), "
                         "Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(5), Symbol('r')), "
                         "Integer(1)), Pow(Add(Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('q'), "
                         "Integer(2))), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), Integer(2)), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(-1), Integer(2), Pow(Symbol('p'), "
                         "Integer(2)), Symbol('q')), Mul(Pow(Symbol('p'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('p'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(2), Symbol('p'), "
                         "Pow(Symbol('q'), Integer(2)), Symbol('r')), Mul(Integer(-1), Integer(2), "
                         "Symbol('p'), Pow(Symbol('q'), Integer(2))), Mul(Integer(2), Symbol('p'), "
                         "Symbol('q'), Pow(Symbol('r'), Integer(2))), Mul(Integer(8), Symbol('p'), "
                         "Symbol('q'), Symbol('r')), Mul(Integer(4), Symbol('p'), Symbol('q')), "
                         "Mul(Integer(-1), Integer(4), Symbol('p'), Pow(Symbol('r'), Integer(2))), "
                         "Mul(Integer(-1), Integer(6), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                         "Integer(2), Symbol('p')), Mul(Pow(Symbol('q'), Integer(2)), Pow(Symbol('r'), "
                         "Integer(2))), Mul(Integer(2), Pow(Symbol('q'), Integer(2)), Symbol('r')), "
                         "Pow(Symbol('q'), Integer(2)), Mul(Integer(-1), Integer(4), Symbol('q'), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(-1), Integer(6), Symbol('q'), "
                         "Symbol('r')), Mul(Integer(-1), Integer(2), Symbol('q')), Mul(Integer(4), "
                         "Pow(Symbol('r'), Integer(2))), Mul(Integer(4), Symbol('r')), Integer(1)), "
                         'Integer(-1)))]])',
          'rows': [['p', '1-p', '0'], ['0', 'q', '1-q'], ['r', '0', '1-r']],
          'stationary': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(-1), Symbol('q'), Symbol('r')), "
                        "Symbol('r')), Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p'), "
                        "Symbol('r')), Mul(Integer(-1), Symbol('p')), Mul(Integer(-1), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(-1), Symbol('q')), Mul(Integer(2), Symbol('r')), "
                        "Integer(1)), Integer(-1))), Mul(Add(Mul(Integer(-1), Symbol('p'), Symbol('r')), "
                        "Symbol('r')), Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('p'), "
                        "Symbol('r')), Mul(Integer(-1), Symbol('p')), Mul(Integer(-1), Symbol('q'), "
                        "Symbol('r')), Mul(Integer(-1), Symbol('q')), Mul(Integer(2), Symbol('r')), "
                        "Integer(1)), Integer(-1))), Mul(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                        "Symbol('p')), Mul(Integer(-1), Symbol('q')), Integer(1)), Pow(Add(Mul(Symbol('p'), "
                        "Symbol('q')), Mul(Integer(-1), Symbol('p'), Symbol('r')), Mul(Integer(-1), "
                        "Symbol('p')), Mul(Integer(-1), Symbol('q'), Symbol('r')), Mul(Integer(-1), "
                        "Symbol('q')), Mul(Integer(2), Symbol('r')), Integer(1)), Integer(-1)))]])"},
 'sym3b': {'absorbing': 'ERR:ShapeError:Matrix size mismatch: (3, 3) * (0, 3).',
           'classes': '[([Integer(0), Integer(1), Integer(2)], true, Integer(1))]',
           'fundamental': "ImmutableDenseMatrix([[Mul(Add(Mul(Integer(-1), Integer(8), Symbol('p'), "
                          "Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(28), "
                          "Pow(Symbol('q'), Integer(2))), Mul(Integer(22), Symbol('q')), Integer(8)), "
                          "Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(10), Symbol('p'), "
                          "Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(25), "
                          "Pow(Symbol('q'), Integer(2))), Mul(Integer(20), Symbol('q')), Integer(4)), "
                          "Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(3), Pow(Symbol('p'), "
                          "Integer(2))), Mul(Integer(12), Symbol('p'), Symbol('q')), Mul(Integer(6), "
                          "Symbol('p')), Mul(Integer(-1), Integer(3), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(-1), Integer(15), Symbol('q'))), Pow(Add(Pow(Symbol('p'), "
                          "Integer(2)), Mul(Integer(-1), Integer(10), Symbol('p'), Symbol('q')), "
                          "Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), "
                          "Integer(2))), Mul(Integer(20), Symbol('q')), Integer(4)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(4), Pow(Symbol('p'), Integer(2))), Mul(Integer(-1), "
                          "Integer(14), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(6), "
                          "Symbol('p')), Mul(Integer(13), Symbol('q')), Integer(-4)), "
                          "Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(10), Symbol('p'), "
                          "Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(25), "
                          "Pow(Symbol('q'), Integer(2))), Mul(Integer(20), Symbol('q')), Integer(4)), "
                          "Integer(-1)))], [Mul(Add(Mul(Integer(-1), Integer(2), Symbol('p'), Symbol('q')), "
                          "Mul(Integer(-1), Integer(2), Symbol('p')), Mul(Integer(-1), Integer(2), "
                          "Pow(Symbol('q'), Integer(2))), Integer(4)), Pow(Add(Pow(Symbol('p'), Integer(2)), "
                          "Mul(Integer(-1), Integer(10), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(4), Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(20), Symbol('q')), Integer(4)), Integer(-1))), "
                          "Mul(Add(Mul(Integer(-1), Integer(9), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(3), Symbol('p')), Mul(Integer(27), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(12), Symbol('q')), Integer(6)), Pow(Add(Pow(Symbol('p'), Integer(2)), "
                          "Mul(Integer(-1), Integer(10), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(4), Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(20), Symbol('q')), Integer(4)), Integer(-1))), "
                          "Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Symbol('p'), Symbol('q')), Symbol('p'), "
                          "Mul(Integer(8), Symbol('q')), Integer(-6)), Pow(Add(Pow(Symbol('p'), Integer(2)), "
                          "Mul(Integer(-1), Integer(10), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(4), Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(20), Symbol('q')), Integer(4)), Integer(-1)))], "
                          "[Mul(Add(Mul(Integer(-1), Integer(2), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(2), Pow(Symbol('q'), Integer(2))), Mul(Integer(-1), Integer(10), "
                          "Symbol('q'))), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), "
                          "Integer(10), Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(4), "
                          "Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), Integer(2))), Mul(Integer(20), "
                          "Symbol('q')), Integer(4)), Integer(-1))), Mul(Add(Mul(Integer(-1), Integer(9), "
                          "Symbol('p'), Symbol('q')), Mul(Integer(27), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(-1), Integer(3), Symbol('q'))), Pow(Add(Pow(Symbol('p'), Integer(2)), "
                          "Mul(Integer(-1), Integer(10), Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                          "Integer(4), Symbol('p')), Mul(Integer(25), Pow(Symbol('q'), Integer(2))), "
                          "Mul(Integer(20), Symbol('q')), Integer(4)), Integer(-1))), "
                          "Mul(Add(Pow(Symbol('p'), Integer(2)), Mul(Symbol('p'), Symbol('q')), "
                          "Mul(Integer(-1), Integer(4), Symbol('p')), Mul(Integer(33), Symbol('q')), "
                          "Integer(4)), Pow(Add(Pow(Symbol('p'), Integer(2)), Mul(Integer(-1), Integer(10), "
                          "Symbol('p'), Symbol('q')), Mul(Integer(-1), Integer(4), Symbol('p')), "
                          "Mul(Integer(25), Pow(Symbol('q'), Integer(2))), Mul(Integer(20), Symbol('q')), "
                          'Integer(4)), Integer(-1)))]])',
           'rows': [['1/2', 'p/2', '(1-p)/2'], ['1/3', '1/3', '1/3'], ['0', 'q', '1-q']],
           'stationary': "ImmutableDenseMatrix([[Mul(Integer(-1), Integer(2), Symbol('q'), "
                         "Pow(Add(Symbol('p'), Mul(Integer(-1), Integer(5), Symbol('q')), Integer(-2)), "
                         "Integer(-1))), Mul(Integer(-1), Integer(3), Symbol('q'), Pow(Add(Symbol('p'), "
                         "Mul(Integer(-1), Integer(5), Symbol('q')), Integer(-2)), Integer(-1))), "
                         "Mul(Add(Symbol('p'), Integer(-2)), Pow(Add(Symbol('p'), Mul(Integer(-1), "
                         "Integer(5), Symbol('q')), Integer(-2)), Integer(-1)))]])"},
 'symabs': {'absorbing': "ImmutableDenseMatrix([[Symbol('p'), Add(Integer(1), Mul(Integer(-1), "
                         "Symbol('p')))]])",
            'classes': '[([Integer(0)], true, Integer(1)), ([Integer(2)], true, Integer(1)), ([Integer(1)], '
                       'false, Integer(1))]',
            'fundamental': 'ImmutableDenseMatrix([[Integer(1)]])',
            'rows': [['1', '0', '0'], ['p', '0', '1-p'], ['0', '0', '1']],
            'stationary': "ImmutableDenseMatrix([[Add(Integer(1), Mul(Integer(-1), Symbol('tau0'))), "
                          "Integer(0), Symbol('tau0')]])"},
 'symabs2': {'absorbing': "ImmutableDenseMatrix([[Mul(Symbol('p'), Pow(Add(Mul(Symbol('p'), Symbol('q')), "
                          "Mul(Integer(-1), Symbol('q')), Integer(1)), Integer(-1))), Mul(Add(Integer(1), "
                          "Mul(Integer(-1), Symbol('p'))), Add(Integer(1), Mul(Integer(-1), Symbol('q'))), "
                          "Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('q')), "
                          "Integer(1)), Integer(-1)))], [Mul(Symbol('p'), Symbol('q'), "
                          "Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('q')), "
                          "Integer(1)), Integer(-1))), Mul(Add(Integer(1), Mul(Integer(-1), Symbol('q'))), "
                          "Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), Symbol('q')), "
                          'Integer(1)), Integer(-1)))]])',
             'classes': '[([Integer(0)], true, Integer(1)), ([Integer(3)], true, Integer(1)), ([Integer(1), '
                        'Integer(2)], false, Integer(1))]',
             'fundamental': "ImmutableDenseMatrix([[Mul(Integer(2), Pow(Add(Mul(Symbol('p'), Symbol('q')), "
                            "Mul(Integer(-1), Symbol('q')), Integer(1)), Integer(-1))), Mul(Add(Integer(1), "
                            "Mul(Integer(-1), Symbol('p'))), Pow(Add(Mul(Symbol('p'), Symbol('q')), "
                            "Mul(Integer(-1), Symbol('q')), Integer(1)), Integer(-1)))], [Mul(Integer(2), "
                            "Symbol('q'), Pow(Add(Mul(Symbol('p'), Symbol('q')), Mul(Integer(-1), "
                            "Symbol('q')), Integer(1)), Integer(-1))), Pow(Add(Mul(Symbol('p'), "
                            "Symbol('q')), Mul(Integer(-1), Symbol('q')), Integer(1)), Integer(-1))]])",
             'rows': [['1', '0', '0', '0'],
                      ['p/2', '1/2', '(1-p)/2', '0'],
                      ['0', 'q', '0', '1-q'],
                      ['0', '0', '0', '1']],
             'stationary': "ImmutableDenseMatrix([[Add(Integer(1), Mul(Integer(-1), Symbol('tau0'))), "
                           "Integer(0), Integer(0), Symbol('tau0')]])"}}


def _rec(f):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return S.srepr(f())
    except Exception as e:
        return 'ERR:%s' % type(e).__name__


def _same(expected, got):
    if expected.startswith('ERR:'):
        return got == ':'.join(expected.split(':')[:2])
    return got == expected


def _dtmc(rows):
    T = S.Matrix([[S.sympify(x, locals=LOCALS) for x in row] for row in rows])
    return DiscreteMarkovChain('X', list(range(T.rows)), T), T.rows


DTMC = sorted(k for k, v in CASES.items() if 'stationary' in v)
CTMC = sorted(k for k, v in CASES.items() if 'limiting' in v)


@pytest.mark.parametrize('name', DTMC)
def test_dtmc(name):
    case = CASES[name]
    X, n = _dtmc(case['rows'])
    Eq, Ne = S.Eq, S.Ne
    calls = {'stationary': lambda: X.stationary_distribution(),
             'fundamental': lambda: X.fundamental_matrix(),
             'absorbing': lambda: X.absorbing_probabilities(),
             'P2': lambda: P(Eq(X[2], n - 1), Eq(X[0], 0)),
             'P3': lambda: P(Eq(X[3], 0), Eq(X[1], n - 1)),
             'Pne': lambda: P(Ne(X[2], 0), Eq(X[0], 0)),
             'Pgt': lambda: P(X[2] > 0, Eq(X[0], 0)),
             'E2': lambda: E(X[2], Eq(X[0], 0)),
             'E2sq': lambda: E(X[1]**2, Eq(X[0], n - 1)),
             'Pand': lambda: P(Eq(X[1], 0) & Eq(X[2], n - 1), Eq(X[0], 0)),
             'classes': lambda: X.communication_classes()}
    for key, expected in case.items():
        if key in calls:
            assert _same(expected, _rec(calls[key])), (name, key)


@pytest.mark.parametrize('name', CTMC)
def test_ctmc_limiting(name):
    case = CASES[name]
    G = S.Matrix([[S.sympify(x, locals=LOCALS) for x in row] for row in case['rows']])
    C = ContinuousMarkovChain('C', list(range(G.rows)), G)
    assert _same(case['limiting'], _rec(lambda: C.limiting_distribution()))


@pytest.mark.parametrize('name', ['ctmcq', 'ctmcq2'])
def test_ctmc_query_forms(name):
    T2 = S.Matrix([[-1, 1, 0], [1, -1, 0], [0, 1, -1]])
    C2 = ContinuousMarkovChain('C', [0, 1, 2], T2)
    Eq, Ne = S.Eq, S.Ne
    calls = {'q612': lambda: P(Eq(C2(1), 1), Eq(C2(0), 1)),
             'q2': lambda: P(Eq(C2(2), 0), Eq(C2(0), 1)),
             'qE': lambda: E(C2(1), Eq(C2(0), 1)),
             'qne': lambda: P(Ne(C2(1), 1), Eq(C2(0), 1)),
             'qgt': lambda: P(C2(1) > 0, Eq(C2(0), 2))}
    for key, expected in CASES[name].items():
        if key in calls:
            assert _same(expected, _rec(calls[key])), (name, key)


def test_process_variance_is_sympys_central_moment():
    # SymPy: variance(X, c) = cmoment(X, 2, c), the unsimplified
    # E((X - E(X))**2) (statement 619 of test_ContinuousMarkovChain).
    from symox.stats import variance
    T2 = S.Matrix([[-1, 1, 0], [1, -1, 0], [0, 1, -1]])
    C2 = ContinuousMarkovChain('C', [0, 1, 2], T2)
    h, e3 = S.S.Half, S.exp(-3)
    got = variance(C2(S.Rational(3, 2)), S.Eq(C2(0), 1))
    assert got == (h - e3/2)**2*(e3/2 + h) + (-h - e3/2)**2*(h - e3/2)
    assert S.srepr(got) == S.srepr((h - e3/2)**2*(e3/2 + h) + (-h - e3/2)**2*(h - e3/2))
    T = S.Matrix([[S.Rational(1, 2), S.Rational(1, 2)], [S.Rational(1, 4), S.Rational(3, 4)]])
    X = DiscreteMarkovChain('X', [0, 1], T)
    assert variance(X[1], S.Eq(X[0], 0)) == S.Rational(1, 4)


def test_matrix_rows_are_checked_in_queries_only():
    # SymPy validates only the dimensions at construction; P/E queries
    # raise when (sum(row) - 1) != 0 structurally, the other methods accept
    # the matrix (Float rows summing to 0.99999999999999989, negative
    # entries).
    F = S.Float
    T = S.Matrix([[F(0.1), F(0.2), F(0.7)], [F(0.3), F(0.6), F(0.1)], [F(0.6), F(0.1), F(0.3)]])
    assert (sum(T.row(1)) - 1) != 0
    X = DiscreteMarkovChain('X', [0, 1, 2], T)
    assert X.stationary_distribution().shape == (1, 3)
    with pytest.raises(ValueError, match='Values in a row must sum to 1'):
        P(S.Eq(X[1], 0), S.Eq(X[0], 0))
    with pytest.raises(ValueError, match='Values in a row must sum to 1'):
        E(X[1], S.Eq(X[0], 0))
    Y = DiscreteMarkovChain('Y', [0, 1], S.Matrix([[2, -1], [0, 1]]))
    assert Y.communication_classes() is not None
    assert P(S.Eq(Y[1], 1), S.Eq(Y[0], 0)) == -1
    with pytest.raises(ValueError):
        DiscreteMarkovChain('Z', [0, 1, 2], S.Matrix([[1, 0], [0, 1]]))
