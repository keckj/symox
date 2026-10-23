"""Wave-4 parity regressions for audit batches c3/c4 (c-rest-b).

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

from symox import *
del test
from symox.abc import x, y, z, s
from symox.categories import *
from symox.physics.control import *
from symox.physics.optics import jones_2_stokes, jones_vector
from symox.vector import CoordSys3D

OA, OB, OC = Object('A'), Object('B'), Object('C')
f = NamedMorphism(OA, OB, 'f')
g = NamedMorphism(OB, OC, 'g')
D = Diagram([f, g], {g*f: 'unique'})
N = CoordSys3D('N')
A, B, C, D2, E = [Object(n) for n in 'ABCDE']
h = NamedMorphism(D2, A, 'h')
k = NamedMorphism(D2, B, 'k')
fA, gB = NamedMorphism(A, B, 'f'), NamedMorphism(B, C, 'g')
dg = Diagram([fA, gB])
dgfk = Diagram([fA, gB, h, k])
dseq = Diagram([fA, gB, NamedMorphism(C, D2, 'h2'), NamedMorphism(D2, E, 'i')])
dloop = Diagram([fA, NamedMorphism(B, A, 'g'), NamedMorphism(A, A, 'u'), NamedMorphism(B, A, 'g')*fA])
tf1 = TransferFunction(s, 1 - s, s)
tf2 = TransferFunction(1, s, s)
tf3 = TransferFunction(5, 1, s)
tf4 = TransferFunction(s - 1, s, s)
tf5 = TransferFunction(0, 1, s)
F1 = MIMOFeedback(TransferFunctionMatrix([[tf1, tf2], [tf3, tf4]]),
                  TransferFunctionMatrix([[tf3, tf5], [tf5, tf5]]), 1)

CASES = {
    'list_frac': ('pretty([x/2, Rational(1,3), x**2], use_unicode=False)',
        ' x        2 \n[-, 1/3, x ]\n 2          '),
    'tuple_one': ('pretty((x/y,), use_unicode=False)',
        ' x  \n(-,)\n y  '),
    'tuple_mat': ('pretty((Matrix([[1, x], [y/2, 3]]), x), use_unicode=True)',
        '⎛⎡1  x⎤   ⎞\n⎜⎢    ⎥   ⎟\n⎜⎢y   ⎥, x⎟\n⎜⎢─  3⎥   ⎟\n⎝⎣2   ⎦   ⎠'),
    'dict': ('pretty({x: x**2/3, y: 1}, use_unicode=False)',
        '     2       \n    x        \n{x: --, y: 1}\n    3        '),
    'set': ('pretty({x**2, y}, use_unicode=False)',
        '  2    \n{x , y}'),
    'empty_set': ('pretty(set())',
        'set()'),
    'frozenset': ('pretty(frozenset([x, y]), use_unicode=False)',
        'frozenset({x, y})'),
    'nested': ('pretty([[x, y/2], (1,)], use_unicode=False)',
        '     y        \n[[x, -], (1,)]\n     2        '),
    'jones': ('pretty([jones_2_stokes(e) for e in [jones_vector(0, 0), jones_vector(pi/2, 0), jones_vector(pi/4, pi/8)]], use_unicode=True)',
        '⎡           ⎡             1             ⎤⎤\n⎢           ⎢                           ⎥⎥\n⎢           ⎢             0             ⎥⎥\n⎢⎡1⎤  ⎡1 ⎤  ⎢                           ⎥⎥\n⎢⎢ ⎥  ⎢  ⎥  ⎢            √2             ⎥⎥\n⎢⎢1⎥  ⎢-1⎥  ⎢            ──             ⎥⎥\n⎢⎢ ⎥, ⎢  ⎥, ⎢            2              ⎥⎥\n⎢⎢0⎥  ⎢0 ⎥  ⎢                           ⎥⎥\n⎢⎢ ⎥  ⎢  ⎥  ⎢      ________     ________⎥⎥\n⎢⎣0⎦  ⎣0 ⎦  ⎢     ╱ 1   √2     ╱ √2   1 ⎥⎥\n⎢           ⎢2⋅  ╱  ─ - ── ⋅  ╱  ── + ─ ⎥⎥\n⎣           ⎣  ╲╱   2   4   ╲╱   4    2 ⎦⎦'),
    'tf': ('pretty(TransferFunction(s + 1, s**2 + 2, s), use_unicode=False)',
        's + 1 \n------\n 2    \ns  + 2'),
    'tf_one': ('pretty(TransferFunction(1, s**2 + 2, s), use_unicode=False)',
        '  1   \n------\n 2    \ns  + 2'),
    'series': ('pretty(Series(TransferFunction(s + 1, s**2 + 2, s), TransferFunction(3, s, s)), use_unicode=False)',
        '/s + 1 \\ /3\\\n|------|*|-|\n| 2    | \\s/\n\\s  + 2/    '),
    'parallel': ('pretty(Parallel(TransferFunction(s + 1, s**2 + 2, s), TransferFunction(3, s, s)), use_unicode=False)',
        's + 1     \n------   3\n 2     + -\ns  + 2   s'),
    'feedback': ('pretty(Feedback(TransferFunction(s + 1, s**2 + 2, s), TransferFunction(3, s, s)), use_unicode=False)',
        '    /s + 1 \\    \n    |------|    \n    | 2    |    \n    \\s  + 2/    \n----------------\n1   /s + 1 \\ /3\\\n- + |------|*|-|\n1   | 2    | \\s/\n    \\s  + 2/    '),
    'feedback_pos': ('pretty(Feedback(TransferFunction(s + 1, s**2 + 2, s), TransferFunction(3, s, s), 1), use_unicode=True)',
        '    ⎛s + 1 ⎞    \n    ⎜──────⎟    \n    ⎜ 2    ⎟    \n    ⎝s  + 2⎠    \n────────────────\n1   ⎛s + 1 ⎞ ⎛3⎞\n─ - ⎜──────⎟⋅⎜─⎟\n1   ⎜ 2    ⎟ ⎝s⎠\n    ⎝s  + 2⎠    '),
    'tfm': ('pretty(TransferFunctionMatrix([[TransferFunction(s, s + 1, s)], [TransferFunction(1, s**2, s)]]), use_unicode=False)',
        '[  s  ]   \n[-----]   \n[s + 1]   \n[     ]   \n[ 1   ]   \n[ --  ]   \n[  2  ]   \n[ s   ]{t}'),
    'tfm_u': ('pretty(TransferFunctionMatrix([[TransferFunction(s, s + 1, s), TransferFunction(2, 1, s)]]), use_unicode=True)',
        '⎡  s    2⎤ \n⎢─────  ─⎥ \n⎣s + 1  1⎦τ'),
    'mimo_series': ('pretty(MIMOSeries(TransferFunctionMatrix([[TransferFunction(s, s + 1, s)]]), TransferFunctionMatrix([[TransferFunction(1, s, s)]])), use_unicode=False)',
        '[1]    [  s  ]   \n[-]   *[-----]   \n[s]{t} [s + 1]{t}'),
    'mimo_parallel': ('pretty(MIMOParallel(TransferFunctionMatrix([[TransferFunction(s, s + 1, s)]]), TransferFunctionMatrix([[TransferFunction(1, s, s)]])), use_unicode=False)',
        '[  s  ]      [1]   \n[-----]    + [-]   \n[s + 1]{t}   [s]{t}'),
    'mimo_feedback': ('pretty(MIMOFeedback(TransferFunctionMatrix([[TransferFunction(s, s + 1, s)]]), TransferFunctionMatrix([[TransferFunction(1, s, s)]])), use_unicode=False)',
        '/    [  s  ]    [1]   \\-1   [  s  ]   \n|I + [-----]   *[-]   |   * [-----]   \n\\    [s + 1]{t} [s]{t}/     [s + 1]{t}'),
    'list_tf': ('pretty([TransferFunction(1, s, s), x], use_unicode=False)',
        ' 1    \n[-, x]\n s    '),
    'diag_keys': ('pretty(sorted(D.premises.keys(), key=default_sort_key), use_unicode=False)',
        '[g*f:A-->C, id:A-->A, id:B-->B, id:C-->C, f:A-->B, g:B-->C]'),
    'diag_prem': ('pretty(D.premises, use_unicode=False, num_columns=200)',
        '{g*f:A-->C: EmptySet, id:A-->A: EmptySet, id:B-->B: EmptySet, id:C-->C: EmptySet, f:A-->B: EmptySet, g:B-->C: EmptySet}'),
    'diag_concl': ('pretty(D.conclusions, use_unicode=True)',
        '{g∘f:A——▶C: {unique}}'),
    'diag_hom': ('pretty(D.hom(OA, OC), use_unicode=False)',
        '({g*f:A-->C}, {g*f:A-->C})'),
    'diag_hom_empty': ('pretty(D.hom(OC, OA), use_unicode=False)',
        '(EmptySet, EmptySet)'),
    'diag_hom_str': ('str(D.hom(OA, OC))',
        '({CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")))}, {CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")))})'),
    'diag_concl_str': ('str(D.conclusions)',
        '{CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))): {unique}}'),
    'vec_add_kind': ('str(Add(N.i*2 + N.j, N.j*N.x).kind)',
        'VectorKind(NumberKind)'),
    'vec_add_val': ('str(Add(N.i*2 + N.j, N.j*N.x))',
        '2*N.i + (N.x + 1)*N.j'),
    'vec_add_one': ('str(Add(N.i))',
        'N.i'),
    'mat_list_wrap': ('pretty([Matrix([[x**10+y**10+x*y*z*1000, 1]]), Matrix([[x**10+y**10+x*y*z*1000, 1]])], use_unicode=False, num_columns=40)',
        ' [ 10                 10   ]  [ 10     >\n[[x   + 1000*x*y*z + y    1], [x   + 1 >\n\n>              10   ] \n> 000*x*y*z + y    1]]'),
    'x_grid_fg': ('pretty(DiagramGrid(dg))',
        'A  B\n    \n   C'),
    'x_grid_fghk': ('pretty(DiagramGrid(dgfk))',
        'A  B  D\n       \n   C   '),
    'x_grid_groups': ('pretty(DiagramGrid(dgfk, groups=[[A, D2], B, C]))',
        'B     C\n       \nA  D   '),
    'x_grid_seq': ("pretty(DiagramGrid(dseq, layout='sequential'))",
        'A  B  C  D  E'),
    'x_grid_seq_t': ("pretty(DiagramGrid(dseq, layout='sequential', transpose=True))",
        'A\n \nB\n \nC\n \nD\n \nE'),
    'x_grid_dims': ('str([(G.width, G.height) for G in (DiagramGrid(dg), DiagramGrid(dgfk), DiagramGrid(dseq))])',
        '[(2, 2), (3, 2), (3, 3)]'),
    'x_grid_morphisms': ('str(sorted(map(str, DiagramGrid(dgfk).morphisms)))',
        '[\'NamedMorphism(Object("A"), Object("B"), "f")\', \'NamedMorphism(Object("B"), Object("C"), "g")\', \'NamedMorphism(Object("D"), Object("A"), "h")\', \'NamedMorphism(Object("D"), Object("B"), "k")\']'),
    'x_xypic_fghk': ('xypic_draw_diagram(dgfk)',
        '\\xymatrix{\nA \\ar[r]_{f} & B \\ar[d]^{g} & D \\ar[l]^{k} \\ar@/_3mm/[ll]_{h} \\\\\n& C & \n}\n'),
    'x_xypic_loop': ('xypic_draw_diagram(dloop)',
        '\\xymatrix{\nA \\ar[r]^{f} \\ar@(u,l)[]^{u} & B \\ar@/^3mm/[l]^{g} \n}\n'),
    'x_xypic_seq': ("xypic_draw_diagram(dseq, layout='sequential')",
        '\\xymatrix{\nA \\ar[r]^{f} & B \\ar[r]^{g} & C \\ar[r]^{h_{2}} & D \\ar[r]^{i} & E \n}\n'),
    'x_mimo_doit': ('pretty(F1.doit(), use_unicode=False)',
        '[  -s           s - 1       ]   \n[-------     -----------    ]   \n[6*s - 1     s*(6*s - 1)    ]   \n[                           ]   \n[5*s - 5  (s - 1)*(6*s + 24)]   \n[-------  ------------------]   \n[6*s - 1     s*(6*s - 1)    ]{t}'),
    'x_mimo_doit_nocancel': ('pretty(F1.doit(cancel=False), use_unicode=False)',
        '[             s*(s - 1)                              s - 1               ]   \n[         -----------------                       -----------            ]   \n[         (1 - s)*(6*s - 1)                       s*(6*s - 1)            ]   \n[                                                                        ]   \n[s*(25*s - 25) + 5*(1 - s)*(6*s - 1)  s*(s - 1)*(6*s - 1) + s*(25*s - 25)]   \n[-----------------------------------  -----------------------------------]   \n[         (1 - s)*(6*s - 1)                       2                      ]   \n[                                                s *(6*s - 1)            ]{t}'),
    'x_mimo_doit_expand': ('pretty(F1.doit(expand=True), use_unicode=False)',
        '[  -s          s - 1      ]   \n[-------      --------    ]   \n[6*s - 1         2        ]   \n[             6*s  - s    ]   \n[                         ]   \n[            2            ]   \n[5*s - 5  6*s  + 18*s - 24]   \n[-------  ----------------]   \n[6*s - 1         2        ]   \n[             6*s  - s    ]{t}'),
    'x_cancel_1': ('str(cancel((s - 1)/(s*(6*s - 1)), expand=False))',
        '(s - 1)/(s*(6*s - 1))'),
    'x_cancel_2': ('str(cancel(Mul(s*(s - 1)*(6*s - 1) + s*(25*s - 25), 1/(s**2*(6*s - 1)), evaluate=False), expand=False))',
        '(s - 1)*(6*s + 24)/(s*(6*s - 1))'),
    'x_cancel_3': ('str(cancel((x + 1)*(x + 2)/(x*(x + 3)), expand=False))',
        '(x + 1)*(x + 2)/(x*(x + 3))'),
    'x_cancel_4': ('str(cancel((x**2 - 1)/(x + 1)/(x - 1), expand=False))',
        '(x**2 - 1)/((x - 1)*(x + 1))'),
    'x_cancel_5': ('str(cancel(Rational(2, 3)*(x - 1)/(Rational(4, 5)*(1 - x)), expand=False))',
        '-5/6'),
    'x_cancel_6': ('str(cancel((x + 1)*(y - 1)/(2*(1 - y)*(x + 2)), expand=False))',
        '(-x - 1)/(2*x + 4)'),
    'x_cancel_7': ('str(cancel((x + I)*(x - I)/(x**2 + 1), expand=False))',
        '(x - I)*(x + I)/(x**2 + 1)'),
    'x_jones_sym': ('str(jones_2_stokes(Matrix([x, I*y])))',
        'Matrix([[Abs(x)**2 + Abs(y)**2], [Abs(x)**2 - Abs(y)**2], [2*im(x*conjugate(y))], [2*re(x*conjugate(y))]])'),
    'x_sort_keys': ('str(sorted([gB*fA, fA, IdentityMorphism(A), B, A, gB], key=default_sort_key))',
        '[Object("A"), Object("B"), CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g"))), IdentityMorphism(Object("A")), NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")]'),
}


@pytest.mark.parametrize("key", sorted(CASES))
def test_matches_sympy(key):
    code, want = CASES[key]
    assert eval(code) == want


def test_morphism_base_class_is_abstract():
    with pytest.raises(NotImplementedError):
        Morphism(A, B)


def test_diagram_grid_error_behavior():
    # SymPy 1.14 fails on an empty diagram and on groups missing an object.
    with pytest.raises(IndexError):
        xypic_draw_diagram(Diagram())
    with pytest.raises(KeyError):
        DiagramGrid(D, groups=[[OB, OA]])


def test_vector_add_kind():
    v = Add(N.i*2 + N.j, N.j*N.x)
    assert str(v.kind) == 'VectorKind(NumberKind)'
    with pytest.raises(TypeError):
        Add(N.i, 1)
