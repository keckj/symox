"""Wave-4 parity of the sympy.polys.matrices DomainMatrix API (WP-C5).

The expected values were computed with the reference SymPy 1.14 (pure Python
ground types) and are compared in a canonical form: the storage format, the
domain, the shape and the nonzero entries as SymPy expressions.
"""
import pytest

import symox as sp
from symox.polys.matrices import DomainMatrix, DM
from symox.polys.matrices.domainscalar import DomainScalar
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.sdm import SDM
from symox.polys.matrices.exceptions import DMNonInvertibleMatrixError

ZZ, QQ, GF = sp.ZZ, sp.QQ, sp.GF
x = sp.Symbol('x')


def canon(o, K):
    if isinstance(o, DomainMatrix):
        L = o.domain
        return ('DM', o.rep.fmt, str(L), o.shape,
                sorted((k, str(L.to_sympy(v))) for k, v in o.to_dok().items()))
    if isinstance(o, DomainScalar):
        return ('DS', str(o.domain), str(o.domain.to_sympy(o.element)))
    if isinstance(o, (list, tuple)):
        return [canon(e, K) for e in o]
    if isinstance(o, bool) or o is None or isinstance(o, str):
        return o
    if type(o) is dict:
        return sorted((k, canon(v, K)) for k, v in o.items())
    try:
        return str(K.to_sympy(K.convert(o)))
    except Exception:
        return str(QQ.to_sympy(QQ.convert(o)))


EXPECTED = [('DM([[2, 4], [6, 8]], ZZ)', 'content()', '2'),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'primitive()',
  ['2', ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '1'), ((0, 1), '2'), ((1, 0), '3'), ((1, 1), '4')])]),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'cancel_denom(ZZ(4))',
  [('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '1'), ((0, 1), '2'), ((1, 0), '3'), ((1, 1), '4')]), '2']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'cancel_denom_elementwise(ZZ(4))',
  [('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '1'), ((0, 1), '1'), ((1, 0), '3'), ((1, 1), '2')]),
   ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '2'), ((0, 1), '1'), ((1, 0), '2'), ((1, 1), '1')])]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'charpoly()', ['1', '-10', '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)', 'charpoly_factor_list()', [[['1', '-10', '-8'], '1']]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'charpoly_factor_blocks()', [[['1', '-10', '-8'], '1']]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'charpoly_base()', ['1', '-10', '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)', 'charpoly_berk()', ['1', '-10', '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)', 'adj_poly_det()', [['-1', '10'], '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'fflu()',
  [('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '1'), ((1, 1), '1')]),
   ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '2'), ((1, 0), '6'), ((1, 1), '-8')]),
   ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '2'), ((1, 1), '-16')]),
   ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '2'), ((0, 1), '4'), ((1, 1), '-8')])]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'scc()', [['0', '1']]),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'pow(3)',
  ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '296'), ((0, 1), '432'), ((1, 0), '648'), ((1, 1), '944')])),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'eval_poly([ZZ(1), ZZ(0), ZZ(2)])',
  ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '30'), ((0, 1), '40'), ((1, 0), '60'), ((1, 1), '90')])),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'solve_den(DM([[1], [1]], ZZ))',
  [('DM', 'dense', 'ZZ', (2, 1), [((0, 0), '4'), ((1, 0), '-4')]), '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'solve_den_charpoly(DM([[1], [1]], ZZ))',
  [('DM', 'dense', 'ZZ', (2, 1), [((0, 0), '4'), ((1, 0), '-4')]), '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'solve_den_rref(DM([[1], [1]], ZZ))',
  [('DM', 'dense', 'ZZ', (2, 1), [((0, 0), '4'), ((1, 0), '-4')]), '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)', 'nullspace()', ('DM', 'dense', 'ZZ', (0, 2), [])),
 ('DM([[2, 4], [6, 8]], ZZ)', 'to_dok()', [((0, 0), '2'), ((0, 1), '4'), ((1, 0), '6'), ((1, 1), '8')]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'to_dod()', [(0, [(0, '2'), (1, '4')]), (1, [(0, '6'), (1, '8')])]),
 ('DM([[2, 4], [6, 8]], ZZ)', 'to_flat_nz()[0]', ['2', '4', '6', '8']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'inv_den()',
  [('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '8'), ((0, 1), '-4'), ((1, 0), '-6'), ((1, 1), '2')]), '-8']),
 ('DM([[2, 4], [6, 8]], ZZ)',
  'adj_det()',
  [('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '8'), ((0, 1), '-4'), ((1, 0), '-6'), ((1, 1), '2')]), '-8']),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'clear_denoms()',
  [('DS', 'ZZ', '12'),
   ('DM', 'dense', 'QQ', (2, 2), [((0, 0), '6'), ((0, 1), '4'), ((1, 0), '3'), ((1, 1), '12')])]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'clear_denoms_rowwise()',
  [('DM', 'sparse', 'ZZ', (2, 2), [((0, 0), '6'), ((1, 1), '4')]),
   ('DM', 'dense', 'QQ', (2, 2), [((0, 0), '3'), ((0, 1), '2'), ((1, 0), '1'), ((1, 1), '4')])]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'clear_denoms(convert=True)',
  [('DS', 'ZZ', '12'),
   ('DM', 'dense', 'ZZ', (2, 2), [((0, 0), '6'), ((0, 1), '4'), ((1, 0), '3'), ((1, 1), '12')])]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'lu()',
  [('DM', 'dense', 'QQ', (2, 2), [((0, 0), '1'), ((1, 0), '1/2'), ((1, 1), '1')]),
   ('DM', 'dense', 'QQ', (2, 2), [((0, 0), '1/2'), ((0, 1), '1/3'), ((1, 1), '5/6')]),
   []]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'lu_solve(DM([[1], [2]], QQ))',
  ('DM', 'dense', 'QQ', (2, 1), [((0, 0), '4/5'), ((1, 0), '9/5')])),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'qr()',
  [('DM', 'dense', 'QQ', (2, 2), [((0, 0), '1/2'), ((0, 1), '-1/3'), ((1, 0), '1/4'), ((1, 1), '2/3')]),
   ('DM', 'dense', 'QQ', (2, 2), [((0, 0), '1'), ((0, 1), '4/3'), ((1, 1), '1')])]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)', 'charpoly()', ['1', '-3/2', '5/12']),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)', 'charpoly_factor_list()', [[['12', '-18', '5'], '1']]),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'inv()',
  ('DM', 'dense', 'QQ', (2, 2), [((0, 0), '12/5'), ((0, 1), '-4/5'), ((1, 0), '-3/5'), ((1, 1), '6/5')])),
 ('DM([[QQ(1,2), QQ(1,3)], [QQ(1,4), 1]], QQ)',
  'solve_den(DM([[1], [2]], QQ))',
  [('DM', 'dense', 'QQ', (2, 1), [((0, 0), '4/5'), ((1, 0), '9/5')]), '1']),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'nullspace()',
  ('DM', 'dense', 'QQ', (1, 3), [((0, 0), '-1'), ((0, 1), '-1'), ((0, 2), '1')])),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'nullspace(divide_last=True)',
  ('DM', 'dense', 'QQ', (1, 3), [((0, 0), '-1'), ((0, 1), '-1'), ((0, 2), '1')])),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'rref()',
  [('DM', 'dense', 'QQ', (3, 3), [((0, 0), '1'), ((0, 2), '1'), ((1, 1), '1'), ((1, 2), '1')]), ['0', '1']]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'rref_den()',
  [('DM', 'dense', 'QQ', (3, 3), [((0, 0), '1'), ((0, 2), '1'), ((1, 1), '1'), ((1, 2), '1')]),
   '1',
   ['0', '1']]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'rowspace()',
  # Reference: rows 0 and 1, which are dependent; the pivot rows 0 and 2 are kept.
  ('DM', 'dense', 'QQ', (2, 3), [((0, 0), '1'), ((0, 1), '2'), ((0, 2), '3'), ((1, 0), '1'), ((1, 2), '1')])),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'fflu()',
  [('DM', 'dense', 'QQ', (3, 3), [((0, 0), '1'), ((1, 2), '1'), ((2, 1), '1')]),
   ('DM',
    'dense',
    'QQ',
    (3, 3),
    [((0, 0), '1'), ((1, 0), '1'), ((1, 1), '-2'), ((2, 0), '2'), ((2, 2), '1')]),
   ('DM', 'dense', 'QQ', (3, 3), [((0, 0), '1'), ((1, 1), '-2'), ((2, 2), '-2')]),
   ('DM',
    'dense',
    'QQ',
    (3, 3),
    [((0, 0), '1'), ((0, 1), '2'), ((0, 2), '3'), ((1, 1), '-2'), ((1, 2), '-2')])]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'lu()',
  [('DM', 'dense', 'QQ', (3, 3), [((0, 0), '1'), ((1, 0), '1'), ((1, 1), '1'), ((2, 0), '2'), ((2, 2), '1')]),
   ('DM',
    'dense',
    'QQ',
    (3, 3),
    [((0, 0), '1'), ((0, 1), '2'), ((0, 2), '3'), ((1, 1), '-2'), ((1, 2), '-2')]),
   [['1', '2']]]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)', 'rank()', '2'),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], QQ)',
  'charpoly_factor_list()',
  [[['1', '0'], '1'], [['1', '-6', '2'], '1']]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ)',
  'nullspace()',
  ('DM', 'dense', 'ZZ', (1, 3), [((0, 0), '-2'), ((0, 1), '-2'), ((0, 2), '2')])),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ)',
  'rref_den()',
  [('DM', 'dense', 'ZZ', (3, 3), [((0, 0), '-2'), ((0, 2), '-2'), ((1, 1), '-2'), ((1, 2), '-2')]),
   '-2',
   ['0', '1']]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ)',
  'fflu()',
  [('DM', 'dense', 'ZZ', (3, 3), [((0, 0), '1'), ((1, 2), '1'), ((2, 1), '1')]),
   ('DM',
    'dense',
    'ZZ',
    (3, 3),
    [((0, 0), '1'), ((1, 0), '1'), ((1, 1), '-2'), ((2, 0), '2'), ((2, 2), '1')]),
   ('DM', 'dense', 'ZZ', (3, 3), [((0, 0), '1'), ((1, 1), '-2'), ((2, 2), '-2')]),
   ('DM',
    'dense',
    'ZZ',
    (3, 3),
    [((0, 0), '1'), ((0, 1), '2'), ((0, 2), '3'), ((1, 1), '-2'), ((1, 2), '-2')])]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ)', 'lu()', ('EXC', 'DMNotAField')),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ).to_sparse()',
  'nullspace()',
  ('DM', 'sparse', 'ZZ', (1, 3), [((0, 0), '-2'), ((0, 1), '-2'), ((0, 2), '2')])),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ).to_sparse()',
  'rref_den()',
  [('DM', 'sparse', 'ZZ', (3, 3), [((0, 0), '-2'), ((0, 2), '-2'), ((1, 1), '-2'), ((1, 2), '-2')]),
   '-2',
   ['0', '1']]),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ).to_sparse()', 'lu()', ('EXC', 'DMNotAField')),
 ('DM([[1, 2, 3], [2, 4, 6], [1, 0, 1]], ZZ).to_sparse()', 'charpoly()', ['1', '-6', '2', '0']),
 ('DM([[0, 1, 0], [0, 0, 1], [1, 0, 0]], QQ)',
  'charpoly_factor_list()',
  [[['1', '-1'], '1'], [['1', '1', '1'], '1']]),
 ('DM([[0, 1, 0], [0, 0, 1], [1, 0, 0]], QQ)', 'scc()', [['0', '1', '2']]),
 ('DM([[0, 1, 0], [0, 0, 1], [1, 0, 0]], QQ)',
  'pow(4)',
  ('DM', 'dense', 'QQ', (3, 3), [((0, 1), '1'), ((1, 2), '1'), ((2, 0), '1')])),
 ('DM([[1, 0, 0, 0, -20160], [0, 1, 0, 0, 33768], [0, 0, 1, 0, 39578], [0, 0, 0, 1, 47757]], ZZ)',
  'lll()',
  ('DM',
   'dense',
   'ZZ',
   (4, 5),
   [((0, 0), '10'),
    ((0, 1), '-3'),
    ((0, 2), '-2'),
    ((0, 3), '8'),
    ((0, 4), '-4'),
    ((1, 0), '3'),
    ((1, 1), '-9'),
    ((1, 2), '8'),
    ((1, 3), '1'),
    ((1, 4), '-11'),
    ((2, 0), '-3'),
    ((2, 1), '13'),
    ((2, 2), '-9'),
    ((2, 3), '-3'),
    ((2, 4), '-9'),
    ((3, 0), '-12'),
    ((3, 1), '-7'),
    ((3, 2), '-11'),
    ((3, 3), '9'),
    ((3, 4), '-1')])),
 ('DM([[1, 0, 0, 0, -20160], [0, 1, 0, 0, 33768], [0, 0, 1, 0, 39578], [0, 0, 0, 1, 47757]], ZZ)',
  'lll_transform()',
  [('DM',
    'dense',
    'ZZ',
    (4, 5),
    [((0, 0), '10'),
     ((0, 1), '-3'),
     ((0, 2), '-2'),
     ((0, 3), '8'),
     ((0, 4), '-4'),
     ((1, 0), '3'),
     ((1, 1), '-9'),
     ((1, 2), '8'),
     ((1, 3), '1'),
     ((1, 4), '-11'),
     ((2, 0), '-3'),
     ((2, 1), '13'),
     ((2, 2), '-9'),
     ((2, 3), '-3'),
     ((2, 4), '-9'),
     ((3, 0), '-12'),
     ((3, 1), '-7'),
     ((3, 2), '-11'),
     ((3, 3), '9'),
     ((3, 4), '-1')]),
   ('DM',
    'dense',
    'ZZ',
    (4, 4),
    [((0, 0), '10'),
     ((0, 1), '-3'),
     ((0, 2), '-2'),
     ((0, 3), '8'),
     ((1, 0), '3'),
     ((1, 1), '-9'),
     ((1, 2), '8'),
     ((1, 3), '1'),
     ((2, 0), '-3'),
     ((2, 1), '13'),
     ((2, 2), '-9'),
     ((2, 3), '-3'),
     ((3, 0), '-12'),
     ((3, 1), '-7'),
     ((3, 2), '-11'),
     ((3, 3), '9')])]),
 ('DM([[1, x], [x, 1]], ZZ[x])', 'charpoly()', ['1', '-2', '1 - x**2']),
 ('DM([[1, x], [x, 1]], ZZ[x])', 'det()', '1 - x**2'),
 ('DM([[1, x], [x, 1]], ZZ[x])',
  'adj_det()',
  [('DM', 'dense', 'ZZ[x]', (2, 2), [((0, 0), '1'), ((0, 1), '-x'), ((1, 0), '-x'), ((1, 1), '1')]),
   '1 - x**2']),
 ('DM([[1, x], [x, 1]], ZZ[x])',
  'rref_den()',
  [('DM', 'dense', 'ZZ[x]', (2, 2), [((0, 0), '1 - x**2'), ((1, 1), '1 - x**2')]), '1 - x**2', ['0', '1']]),
 ('DM([[1, x], [x, 1]], ZZ[x])',
  'solve_den(DM([[1], [0]], ZZ[x]))',
  [('DM', 'dense', 'ZZ[x]', (2, 1), [((0, 0), '1'), ((1, 0), '-x')]), '1 - x**2']),
 ('DM([[1, x], [x, 1]], QQ.frac_field(x))', 'charpoly()', ['1', '-2', '1 - x**2']),
 ('DM([[1, x], [x, 1]], QQ.frac_field(x))',
  'inv()',
  ('DM',
   'dense',
   'QQ(x)',
   (2, 2),
   [((0, 0), '-1/(x**2 - 1)'),
    ((0, 1), 'x/(x**2 - 1)'),
    ((1, 0), 'x/(x**2 - 1)'),
    ((1, 1), '-1/(x**2 - 1)')])),
 ('DM([[1, x], [x, 1]], QQ.frac_field(x))',
  'solve_den_charpoly(DM([[1], [0]], QQ.frac_field(x)))',
  [('DM', 'dense', 'QQ(x)', (2, 1), [((0, 0), '1'), ((1, 0), '-x')]), '1 - x**2']),
 ('DM([[3, 1], [1, 3]], GF(7))', 'charpoly()', ['1', '1', '1']),
 ('DM([[3, 1], [1, 3]], GF(7))', 'charpoly_factor_list()', [[['1', '3'], '1'], [['1', '-2'], '1']]),
 ('DM([[3, 1], [1, 3]], GF(7))',
  'inv()',
  ('DM', 'dense', 'GF(7)', (2, 2), [((0, 0), '3'), ((0, 1), '-1'), ((1, 0), '-1'), ((1, 1), '3')])),
 ('DM([[3, 1], [1, 3]], GF(7))',
  'lu()',
  [('DM', 'dense', 'GF(7)', (2, 2), [((0, 0), '1'), ((1, 0), '-2'), ((1, 1), '1')]),
   ('DM', 'dense', 'GF(7)', (2, 2), [((0, 0), '3'), ((0, 1), '1'), ((1, 1), '-2')]),
   []]),
 ('DM([[3, 1], [1, 3]], GF(7))', 'det()', '1'),
 ('DomainMatrix.from_dok({(0, 1): ZZ(2), (1, 0): ZZ(3)}, (2, 2), ZZ)',
  'to_dense()',
  ('DM', 'dense', 'ZZ', (2, 2), [((0, 1), '2'), ((1, 0), '3')])),
 ('DomainMatrix.from_dok({(0, 1): ZZ(2), (1, 0): ZZ(3)}, (2, 2), ZZ)',
  'to_dod()',
  [(0, [(1, '2')]), (1, [(0, '3')])]),
 ('DomainMatrix.from_dok({(0, 1): ZZ(2), (1, 0): ZZ(3)}, (2, 2), ZZ)', 'charpoly()', ['1', '0', '-6']),
 ('DomainMatrix.from_dod({0: {1: QQ(1,2)}, 1: {0: QQ(3)}}, (2, 2), QQ)',
  'to_dok()',
  [((0, 1), '1/2'), ((1, 0), '3')]),
 ('DomainMatrix.from_dod({0: {1: QQ(1,2)}, 1: {0: QQ(3)}}, (2, 2), QQ)',
  'inv()',
  ('DM', 'sparse', 'QQ', (2, 2), [((0, 1), '1/3'), ((1, 0), '2')])),
 ('DomainMatrix.from_list_flat([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], (2, 2), ZZ)',
  'to_list_flat()',
  ['1', '2', '3', '4']),
 ('DomainMatrix.from_list_flat([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], (2, 2), ZZ)', 'det()', '-2'),
 ('DomainMatrix.from_dict_sympy(2, 2, {0: {0: sp.Rational(1, 2)}, 1: {1: sp.Integer(3)}})',
  'to_dok()',
  [((0, 0), '1/2'), ((1, 1), '3')]),
 ('DomainMatrix.from_dict_sympy(2, 2, {0: {0: sp.Rational(1, 2)}, 1: {1: sp.Integer(3)}})',
  'charpoly()',
  ['1', '-7/2', '3/2'])]


@pytest.mark.parametrize('source, operation, expected', EXPECTED,
                         ids=[f'{i}-{op}' for i, (_, op, _) in enumerate(EXPECTED)])
def test_domainmatrix_methods_match_reference(source, operation, expected):
    M = eval(source)
    if isinstance(expected, tuple) and expected[0] == 'EXC':
        with pytest.raises(Exception) as info:
            eval('M.' + operation)
        assert type(info.value).__name__ == expected[1]
        return
    assert canon(eval('M.' + operation), M.domain) == expected


def test_sdm_and_ddm_helpers_match_reference():
    from symox.polys.matrices.sdm import (sdm_irref, sdm_nullspace_from_rref,
                                          sdm_particular_from_rref, sdm_rref_den)
    from symox.polys.matrices.dense import (ddm_irref_den, ddm_irref, ddm_iinv,
                                            ddm_ilu_split, ddm_berk)
    from symox.polys.matrices.ddm import ddm_idet
    q = QQ
    assert sdm_irref({0: {0: q(1), 1: q(2)}, 1: {0: q(3), 1: q(4)}}) == \
        ({0: {0: q(1)}, 1: {1: q(1)}}, [0, 1], {})
    r, p, nz = sdm_irref({0: {0: q(1), 1: q(2), 2: q(3)}, 1: {0: q(2), 1: q(4), 2: q(6)}})
    assert (r, p, nz) == ({0: {0: q(1), 1: q(2), 2: q(3)}}, [0], {1: {0}, 2: {0}})
    assert sdm_nullspace_from_rref(r, q.one, 3, p, nz) == \
        ([{1: q(1), 0: q(-2)}, {2: q(1), 0: q(-3)}], [1, 2])
    r, p, nz = sdm_irref({0: {0: q(1), 1: q(2), 2: q(5)}, 1: {0: q(3), 1: q(4), 2: q(6)}})
    assert sdm_particular_from_rref(r, 3, p) == {0: q(-4), 1: q(9, 2)}
    assert sdm_rref_den({0: {0: ZZ(2), 1: ZZ(4)}, 1: {0: ZZ(1), 1: ZZ(3)}}, ZZ) == \
        ({0: {0: ZZ(2)}, 1: {1: ZZ(2)}}, ZZ(2), [0, 1])
    a = [[ZZ(2), ZZ(4)], [ZZ(1), ZZ(3)]]
    assert ddm_irref_den(a, ZZ) == (ZZ(2), [0, 1]) and a == [[2, 0], [0, 2]]
    a = [[q(2), q(4)], [q(1), q(3)]]
    assert ddm_irref(a) == [0, 1] and a == [[1, 0], [0, 1]]
    assert ddm_idet([[ZZ(2), ZZ(4)], [ZZ(1), ZZ(3)]], ZZ) == 2
    ainv = [[None, None], [None, None]]
    ddm_iinv(ainv, [[q(2), q(4)], [q(1), q(3)]], q)
    assert ainv == [[q(3, 2), q(-2)], [q(-1, 2), q(1)]]
    L = [[q(0), q(0)], [q(0), q(0)]]
    a = [[q(0), q(4)], [q(1), q(3)]]
    assert ddm_ilu_split(L, a, q) == [(0, 1)]
    assert a == [[1, 3], [0, 4]] and L == [[0, 0], [0, 0]]
    assert ddm_berk([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], ZZ) == [[1], [-5], [-2]]


def test_ddm_and_sdm_nullspace_over_ZZ_stay_exact():
    # The reference DDM/SDM.nullspace over ZZ returns floats (-2.0, 1.0); the
    # exact fraction-free basis is kept and equals nullspace_from_rref.
    rows = [[ZZ(1), ZZ(2), ZZ(3)], [ZZ(2), ZZ(4), ZZ(6)]]
    expected = [[-2, 1, 0], [-3, 0, 1]]
    N, nonpivots = DDM(rows, (2, 3), ZZ).nullspace()
    assert nonpivots == [1, 2] and N.to_list() == expected
    assert all(ZZ.of_type(e) for row in N.to_list() for e in row)
    N, nonpivots = SDM.from_list(rows, (2, 3), ZZ).nullspace()
    assert nonpivots == [1, 2] and N.to_list() == expected


def test_rowspace_returns_independent_rows():
    # The reference keeps the first rank rows, which is wrong for [[0, 0], [1, 2]].
    M = DM([[0, 0], [1, 2]], QQ)
    assert M.rowspace().to_Matrix() == sp.Matrix([[1, 2]])
    assert DM([[1, -1], [2, -2]], QQ).rowspace().to_Matrix() == sp.Matrix([[1, -1]])


def test_negative_powers_use_the_exact_inverse():
    M = DM([[1, 2], [3, 4]], QQ)
    assert (M**-2).to_Matrix() == (sp.Matrix([[1, 2], [3, 4]])**-2)
    with pytest.raises(DMNonInvertibleMatrixError):
        DM([[1, 2], [2, 4]], QQ)**-1


def test_expression_domains():
    s = sp.sin(x)**2 + sp.cos(x)**2
    M = DomainMatrix.from_Matrix(sp.Matrix([[s, 0], [0, 2]])).convert_to(sp.EXRAW)
    # The reference raises TypeError while sorting the EXRAW factors.
    cp = [sp.EXRAW.to_sympy(c) for c in M.charpoly()]
    assert sp.expand(cp[1] + s + 2) == 0 and sp.expand(cp[2] - 2*s) == 0
    with pytest.raises(sp.polys.polyerrors.DomainError):
        DM([[1, 2], [3, 4]], sp.EX).charpoly_factor_list()


def test_eval_poly_over_polynomial_rings():
    # The reference raises TypeError (DomainMatrix + PolyElement) here.
    K = ZZ[x]
    M = DM([[1, 2], [3, 4]], K)
    P = M.eval_poly([K.one, K.zero, K(2)])
    assert P.to_Matrix() == sp.Matrix([[9, 10], [15, 24]])


def test_dfm_requires_flint_like_the_reference():
    with pytest.raises(NotImplementedError):
        DM([[1, 2], [3, 4]], ZZ).to_dfm()
