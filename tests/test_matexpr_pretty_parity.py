"""LaTeX and pretty printing of matrix expressions against SymPy 1.14.

The expected strings were produced by the reference SymPy printers
(latex, pretty with and without unicode).
"""

import pytest

import symox as S


def _namespace():
    ns = {k: getattr(S, k) for k in dir(S) if not k.startswith("_")}
    n, x = S.symbols("n x")
    ns.update(n=n, x=x)
    for name in "ABCD":
        ns[name] = S.MatrixSymbol(name, n, n)
    ns["X"] = S.MatrixSymbol("X", 3, 3)
    ns["Y"] = S.MatrixSymbol("Y", 2, 2)
    ns["P"] = S.MatrixSymbol("P", 2, 2)
    ns["Q"] = S.MatrixSymbol("Q", 2, 2)
    ns["v"] = S.MatrixSymbol("v", 3, 1)
    ns["w"] = S.MatrixSymbol("w", 3, 1)
    from symox.matrices.expressions.hadamard import hadamard_power
    ns["hadamard_power"] = hadamard_power
    return ns


NS = _namespace()


# (source, latex, ascii pretty, unicode pretty)
CASES = [
    ('Inverse(A*B)', '\\left(A B\\right)^{-1}', '     -1\n(A*B)  ', '     -1\n(A⋅B)  '),
    ('Inverse(A + B)', '\\left(A + B\\right)^{-1}', '       -1\n(A + B)  ', '       -1\n(A + B)  '),
    ('Inverse(A**2)', '\\left(A^{2}\\right)^{-1}', '    -1\n/ 2\\  \n\\A /  ', '    -1\n⎛ 2⎞  \n⎝A ⎠  '),
    ('Inverse(A.T)', '\\left(A^{T}\\right)^{-1}', '    -1\n/ T\\  \n\\A /  ', '    -1\n⎛ T⎞  \n⎝A ⎠  '),
    ('Inverse(Inverse(A))', '\\left(A^{-1}\\right)^{-1}', '     -1\n/ -1\\  \n\\A  /  ', '     -1\n⎛ -1⎞  \n⎝A  ⎠  '),
    ('Inverse(A**Rational(1, 2))', '\\left(A^{\\frac{1}{2}}\\right)^{-1}', '      -1\n/ 1/2\\  \n\\A   /  ', '      -1\n⎛ 1/2⎞  \n⎝A   ⎠  '),
    ('Inverse(Identity(n) + A)', '\\left(\\mathbb{I} + A\\right)^{-1}', '       -1\n(I + A)  ', '       -1\n(𝕀 + A)  '),
    ('Inverse(-A + B)', '\\left(- A + B\\right)^{-1}', '        -1\n(-A + B)  ', '        -1\n(-A + B)  '),
    ('Inverse(KroneckerProduct(A, B))', '\\left(A \\otimes B\\right)^{-1}', '       -1\n(A x B)  ', '       -1\n(A ⨂ B)  '),
    ('Inverse(HadamardPower(A, 2))', '\\left(A^{\\circ {2}}\\right)^{-1}', '     -1\n/ .2\\  \n\\A  /  ', '     -1\n⎛ ∘2⎞  \n⎝A  ⎠  '),
    ('(A.T*B).I', 'B^{-1} \\left(A^{T}\\right)^{-1}', '        -1\n -1 / T\\  \nB  *\\A /  ', '        -1\n -1 ⎛ T⎞  \nB  ⋅⎝A ⎠  '),
    ('C.T*C.I', 'C^{T} C^{-1}', ' T  -1\nC *C  ', ' T  -1\nC ⋅C  '),
    ('MatPow(A, Rational(1, 2))', 'A^{\\frac{1}{2}}', ' 1/2\nA   ', ' 1/2\nA   '),
    ('MatPow(A, Rational(-1, 2))', 'A^{- \\frac{1}{2}}', ' -1/2\nA    ', ' -1/2\nA    '),
    ('A**Rational(3, 2)', 'A^{\\frac{3}{2}}', ' 3/2\nA   ', ' 3/2\nA   '),
    ('(A + B)**Rational(1, 2)', '\\left(A + B\\right)^{\\frac{1}{2}}', '       1/2\n(A + B)   ', '       1/2\n(A + B)   '),
    ('MatPow(A, x)', 'A^{x}', ' x\nA ', ' x\nA '),
    ('MatPow(A, 2*x)', 'A^{2 x}', ' 2*x\nA   ', ' 2⋅x\nA   '),
    ('MatPow(A, x**2)', 'A^{x^{2}}', ' / 2\\\n \\x /\nA    ', ' ⎛ 2⎞\n ⎝x ⎠\nA    '),
    ('MatPow(A, sin(x))', 'A^{\\sin{\\left(x \\right)}}', ' sin(x)\nA      ', ' sin(x)\nA      '),
    ('MatPow(A, -2)', 'A^{-2}', ' -2\nA  ', ' -2\nA  '),
    ('MatPow(A, 0.5)', 'A^{0.5}', ' 0.5\nA   ', ' 0.5\nA   '),
    ('MatPow(A, -0.5)', 'A^{-0.5}', ' -0.5\nA    ', ' -0.5\nA    '),
    ('MatPow(A, I)', 'A^{i}', ' I\nA ', ' ⅈ\nA '),
    ('(A*B)**2', '\\left(A B\\right)^{2}', '     2\n(A*B) ', '     2\n(A⋅B) '),
    ('(A**2).T', '\\left(A^{T}\\right)^{2}', '    2\n/ T\\ \n\\A / ', '    2\n⎛ T⎞ \n⎝A ⎠ '),
    ('Transpose(Inverse(A))', '\\left(A^{-1}\\right)^{T}', '     T\n/ -1\\ \n\\A  / ', '     T\n⎛ -1⎞ \n⎝A  ⎠ '),
    ('Transpose(HadamardProduct(A, B))', '\\left(A \\circ B\\right)^{T}', '      T\n(A.*B) ', '     T\n(A∘B) '),
    ('Adjoint(A*B)', '\\left(A B\\right)^{\\dagger}', '     +\n(A*B) ', '     †\n(A⋅B) '),
    ('Adjoint(Inverse(A))', '\\left(A^{-1}\\right)^{\\dagger}', '     +\n/ -1\\ \n\\A  / ', '     †\n⎛ -1⎞ \n⎝A  ⎠ '),
    ('A**Rational(-1, 2)*B', 'A^{- \\frac{1}{2}} B', ' -1/2  \nA    *B', ' -1/2  \nA    ⋅B'),
    ('(A*B.T + C.I)**Rational(1, 3)', '\\left(C^{-1} + A B^{T}\\right)^{\\frac{1}{3}}', '            1/3\n/   T    -1\\   \n\\A*B  + C  /   ', '            1/3\n⎛   T    -1⎞   \n⎝A⋅B  + C  ⎠   '),
    ('Rational(1, 2)*A', '\\left(\\frac{1}{2}\\right) A', '1/2*A', '1/2⋅A'),
    ('x*A', 'x A', 'x*A', 'x⋅A'),
    ('(x + 1)*A', '\\left(x + 1\\right) A', '(x + 1)*A', '(x + 1)⋅A'),
    ('-x*A', '- x A', '-x*A', '-x⋅A'),
    ('-(x + 1)*A', '\\left(- x - 1\\right) A', '(-x - 1)*A', '(-x - 1)⋅A'),
    ('-I*A', '- i A', '-I*A', '-ⅈ⋅A'),
    ('-2*I*A', '- 2 i A', '-2*I*A', '-2⋅ⅈ⋅A'),
    ('(1 + I)*A', '\\left(1 + i\\right) A', '(1 + I)*A', '(1 + ⅈ)⋅A'),
    ('-A*B', '- A B', '-A*B', '-A⋅B'),
    ('A*HadamardProduct(B, C)', 'A \\left(B \\circ C\\right)', 'A*(B.*C)', 'A⋅(B∘C)'),
    ('KroneckerProduct(A, B)*C', '\\left(A \\otimes B\\right) C', '(A x B)*C', '(A ⨂ B)⋅C'),
    ('A + B - C', 'A + B - C', '-C + A + B', '-C + A + B'),
    ('A - B*D.I*C', 'A - B D^{-1} C', '    -1      \n-B*D  *C + A', '    -1      \n-B⋅D  ⋅C + A'),
    ('X[0:2, 1:3]', 'X\\left[:2, 1:\\right]', 'X[:2, 1:]', 'X[:2, 1:]'),
    ('X[::2, 1:]', 'X\\left[::2, 1:\\right]', 'X[::2, 1:]', 'X[::2, 1:]'),
    ('ImmutableMatrix([[1, 2], [3, 4]])*Y', '\\left[\\begin{matrix}1 & 2\\\\3 & 4\\end{matrix}\\right] Y', '[1  2]  \n[    ]*Y\n[3  4]  ', '⎡1  2⎤  \n⎢    ⎥⋅Y\n⎣3  4⎦  '),
    ('BlockMatrix([[A, B], [C, D]])', '\\left[\\begin{matrix}A & B\\\\C & D\\end{matrix}\\right]', '[A  B]\n[    ]\n[C  D]', '⎡A  B⎤\n⎢    ⎥\n⎣C  D⎦'),
    ('BlockMatrix([[A]])', '\\left[\\begin{matrix}A\\end{matrix}\\right]', 'A', 'A'),
    ('BlockMatrix([[A, B, C]])', '\\left[\\begin{matrix}A & B & C\\end{matrix}\\right]', '[A  B  C]', '[A  B  C]'),
    ('BlockMatrix([[A], [B]])', '\\left[\\begin{matrix}A\\\\B\\end{matrix}\\right]', '[A]\n[ ]\n[B]', '⎡A⎤\n⎢ ⎥\n⎣B⎦'),
    ('BlockMatrix([[A.I, B*C], [C + D, D**2]])', '\\left[\\begin{matrix}A^{-1} & B C\\\\C + D & D^{2}\\end{matrix}\\right]', '[  -1      ]\n[ A     B*C]\n[          ]\n[        2 ]\n[C + D  D  ]', '⎡  -1      ⎤\n⎢ A     B⋅C⎥\n⎢          ⎥\n⎢        2 ⎥\n⎣C + D  D  ⎦'),
    ('BlockMatrix([[A, ZeroMatrix(n, n)], [Identity(n), D]])', '\\left[\\begin{matrix}A & 0\\\\\\mathbb{I} & D\\end{matrix}\\right]', '[A  0]\n[    ]\n[I  D]', '⎡A  𝟘⎤\n⎢    ⎥\n⎣𝕀  D⎦'),
    ('BlockMatrix([[A, B], [C, D]]).T', '\\left[\\begin{matrix}A^{T} & C^{T}\\\\B^{T} & D^{T}\\end{matrix}\\right]', '[ T   T]\n[A   C ]\n[      ]\n[ T   T]\n[B   D ]', '⎡ T   T⎤\n⎢A   C ⎥\n⎢      ⎥\n⎢ T   T⎥\n⎣B   D ⎦'),
    ('Inverse(BlockMatrix([[A, B], [C, D]]))', '\\left(\\left[\\begin{matrix}A & B\\\\C & D\\end{matrix}\\right]\\right)^{-1}', '        -1\n/[A  B]\\  \n|[    ]|  \n\\[C  D]/  ', '        -1\n⎛⎡A  B⎤⎞  \n⎜⎢    ⎥⎟  \n⎝⎣C  D⎦⎠  '),
    ('BlockDiagMatrix(A, B)', '\\left[\\begin{matrix}A & 0\\\\0 & B\\end{matrix}\\right]', '[A  0]\n[    ]\n[0  B]', '⎡A  𝟘⎤\n⎢    ⎥\n⎣𝟘  B⎦'),
    ('BlockDiagMatrix(A, B, C)', '\\left[\\begin{matrix}A & 0 & 0\\\\0 & B & 0\\\\0 & 0 & C\\end{matrix}\\right]', '[A  0  0]\n[       ]\n[0  B  0]\n[       ]\n[0  0  C]', '⎡A  𝟘  𝟘⎤\n⎢       ⎥\n⎢𝟘  B  𝟘⎥\n⎢       ⎥\n⎣𝟘  𝟘  C⎦'),
    ('BlockMatrix([[A, B]])*BlockMatrix([[C], [D]])', '\\left[\\begin{matrix}A & B\\end{matrix}\\right] \\left[\\begin{matrix}C\\\\D\\end{matrix}\\right]', '       [C]\n[A  B]*[ ]\n       [D]', '       ⎡C⎤\n[A  B]⋅⎢ ⎥\n       ⎣D⎦'),
    ('A', 'A', 'A', 'A'),
    ('A*B', 'A B', 'A*B', 'A⋅B'),
    ('A*B + C', 'A B + C', 'A*B + C', 'A⋅B + C'),
    ('A.T', 'A^{T}', ' T\nA ', ' T\nA '),
    ('Adjoint(A)', 'A^{\\dagger}', ' +\nA ', ' †\nA '),
    ('A.I', 'A^{-1}', ' -1\nA  ', ' -1\nA  '),
    ('A**2', 'A^{2}', ' 2\nA ', ' 2\nA '),
    ('A**-1', 'A^{-1}', ' -1\nA  ', ' -1\nA  '),
    ('2*A', '2 A', '2*A', '2⋅A'),
    ('Identity(n)', '\\mathbb{I}', 'I', '𝕀'),
    ('ZeroMatrix(n, n)', '0', '0', '𝟘'),
    ('OneMatrix(n, n)', '1', '1', '𝟙'),
    ('HadamardProduct(A, B)', 'A \\circ B', 'A.*B', 'A∘B'),
    ('HadamardProduct(A, B, C)', 'A \\circ B \\circ C', 'A.*B.*C', 'A∘B∘C'),
    ('HadamardPower(A, 2)', 'A^{\\circ {2}}', ' .2\nA  ', ' ∘2\nA  '),
    ('HadamardPower(A, x)', 'A^{\\circ {x}}', ' .x\nA  ', ' ∘x\nA  '),
    ('KroneckerProduct(A, B)', 'A \\otimes B', 'A x B', 'A ⨂ B'),
    ('KroneckerProduct(A, B, C)', 'A \\otimes B \\otimes C', 'A x B x C', 'A ⨂ B ⨂ C'),
    ('Trace(A)', '\\operatorname{tr}\\left(A \\right)', 'tr(A)', 'tr(A)'),
    ('Trace(A*B)', '\\operatorname{tr}\\left(A B \\right)', 'tr(A*B)', 'tr(A⋅B)'),
    ('Determinant(A)', '\\left|{A}\\right|', '|A|', '│A│'),
    ('Determinant(A + B)', '\\left|{A + B}\\right|', '|A + B|', '│A + B│'),
    ('A[0, 1]', '{A}_{0,1}', 'A_01', 'A₀₁'),
    ('X[0, 1]', '{X}_{0,1}', 'X_01', 'X₀₁'),
    ('(A + B)[0, 1]', '{A}_{0,1} + {B}_{0,1}', 'A_01 + B_01', 'A₀₁ + B₀₁'),
    ('A[x, n - 1]', '{A}_{x,n - 1}', '(A)[x, n - 1]', '(A)[x, n - 1]'),
    ('x*Trace(A)', 'x \\operatorname{tr}\\left(A \\right)', 'x*tr(A)', 'x⋅tr(A)'),
    ('Trace(A)**2', '\\operatorname{tr}\\left(A \\right)^{2}', '     2\ntr(A) ', '     2\ntr(A) '),
    ('Trace(A) + Determinant(B)', '\\left|{B}\\right| + \\operatorname{tr}\\left(A \\right)', '|B| + tr(A)', '│B│ + tr(A)'),
    ("MatrixSymbol('alpha', n, n)", '\\alpha', 'alpha', 'α'),
    ("MatrixSymbol('A_1', n, n)", 'A_{1}', 'A_1', 'A₁'),
    ("MatrixSymbol('Ab^2', n, n)", 'Ab^{2}', 'Ab^2', 'Ab²'),
    ('A + B', 'A + B', 'A + B', 'A + B'),
    ('A*B*C', 'A B C', 'A*B*C', 'A⋅B⋅C'),
    ('A.T*B', 'A^{T} B', ' T  \nA *B', ' T  \nA ⋅B'),
    ('Adjoint(A)*B', 'A^{\\dagger} B', ' +  \nA *B', ' †  \nA ⋅B'),
    ('(A + B).T', 'A^{T} + B^{T}', ' T    T\nA  + B ', ' T    T\nA  + B '),
    ('(A + B)*C', '\\left(A + B\\right) C', '(A + B)*C', '(A + B)⋅C'),
    ('A*(B + C)', 'A \\left(B + C\\right)', 'A*(B + C)', 'A⋅(B + C)'),
    ('2*A*B', '2 A B', '2*A*B', '2⋅A⋅B'),
    ('x*A*B + C', 'x A B + C', 'x*A*B + C', 'x⋅A⋅B + C'),
    ('A - B', 'A - B', '-B + A', '-B + A'),
    ('-A', '- A', '-A', '-A'),
    ('A + 2*B', 'A + 2 B', '2*B + A', '2⋅B + A'),
    ('(A + B).I*C', '\\left(A + B\\right)^{-1} C', '       -1  \n(A + B)  *C', '       -1  \n(A + B)  ⋅C'),
    ('Transpose(A*B)', '\\left(A B\\right)^{T}', '     T\n(A*B) ', '     T\n(A⋅B) '),
    ('Adjoint(A + B)', '\\left(A + B\\right)^{\\dagger}', '       +\n(A + B) ', '       †\n(A + B) '),
    ('Adjoint(A).T', '\\left(A^{T}\\right)^{\\dagger}', '    +\n/ T\\ \n\\A / ', '    †\n⎛ T⎞ \n⎝A ⎠ '),
    ('HadamardProduct(A + B, C)', '\\left(A + B\\right) \\circ C', '(A + B).*C', '(A + B)∘C'),
    ('HadamardProduct(A*B, C)', '\\left(A B\\right) \\circ C', '(A*B).*C', '(A⋅B)∘C'),
    ('KroneckerProduct(A + B, C)', '\\left(A + B\\right) \\otimes C', '(A + B) x C', '(A + B) ⨂ C'),
    ('HadamardPower(A + B, 2)', '\\left(A + B\\right)^{\\circ {2}}', '     .2\nA + B  ', '     ∘2\nA + B  '),
    ('HadamardPower(A.T, 2)', '\\left(A^{T}\\right)^{\\circ {2}}', '    .2\n/ T\\  \n\\A /  ', '    ∘2\n⎛ T⎞  \n⎝A ⎠  '),
    ('X[1:2, :]', 'X\\left[1:2, :\\right]', 'X[1:2, :]', 'X[1:2, :]'),
    ('X[:, 1]', 'X\\left[:, 1:2\\right]', 'X[:, 1:2]', 'X[:, 1:2]'),
    ('Identity(3)', '\\mathbb{I}', 'I', '𝕀'),
    ('ZeroMatrix(2, 3)', '0', '0', '𝟘'),
    ('BlockMatrix([[Identity(2), ZeroMatrix(2, 2)], [ZeroMatrix(2, 2), Identity(2)]])', '\\left[\\begin{matrix}\\mathbb{I} & 0\\\\0 & \\mathbb{I}\\end{matrix}\\right]', '[I  0]\n[    ]\n[0  I]', '⎡𝕀  𝟘⎤\n⎢    ⎥\n⎣𝟘  𝕀⎦'),
    ('MatPow(A + B, 2)', '\\left(A + B\\right)^{2}', '       2\n(A + B) ', '       2\n(A + B) '),
    ('MatPow(A.T, 2)', '\\left(A^{T}\\right)^{2}', '    2\n/ T\\ \n\\A / ', '    2\n⎛ T⎞ \n⎝A ⎠ '),
    ('MatPow(A*B, x)', '\\left(A B\\right)^{x}', '     x\n(A*B) ', '     x\n(A⋅B) '),
    ('x**2*A', 'x^{2} A', ' 2  \nx *A', ' 2  \nx ⋅A'),
    ('sqrt(2)*A', '\\sqrt{2} A', '  ___  \n\\/ 2 *A', '√2⋅A'),
    ('A*B - C*D', 'A B - C D', 'A*B -C*D', 'A⋅B -C⋅D'),
    ('Trace(A).T if False else Trace(A.T)', '\\operatorname{tr}\\left(A^{T} \\right)', '  / T\\\ntr\\A /', '  ⎛ T⎞\ntr⎝A ⎠'),
    ('-sqrt(2)*A', '- \\sqrt{2} A', '   ___  \n-\\/ 2 *A', '-√2⋅A'),
    ('-x*HadamardProduct(A,B)', '- x \\left(A \\circ B\\right)', '-x*(A.*B)', '-x⋅(A∘B)'),
    ('-2*HadamardProduct(A,B)', '- 2 A \\circ B', '-2*(A.*B)', '-2⋅(A∘B)'),
    ('A - x*B', '- x B + A', '-x*B + A', '-x⋅B + A'),
    ('A - I*B', 'A + - i B', '-I*B + A', '-ⅈ⋅B + A'),
    ('A-(x+1)*B', '\\left(- x - 1\\right) B + A', '(-x - 1)*B + A', '(-x - 1)⋅B + A'),
    ('Rational(-1,2)*A', '- \\frac{1}{2} A', '-1/2*A', '-1/2⋅A'),
    ('-2*A + B', '- 2 A + B', '-2*A + B', '-2⋅A + B'),
    ('A.applyfunc(sin)', '{\\left( d \\mapsto \\sin{\\left(d \\right)} \\right)}_{\\circ}\\left({A}\\right)', '(d -> sin(d)).(A)', '(d ↦ sin(d))˳(A)'),
    ('A.applyfunc(Lambda(x, x**2))', '{\\left( x \\mapsto x^{2} \\right)}_{\\circ}\\left({A}\\right)', '/      2\\    \n\\x -> x /.(A)', '⎛     2⎞    \n⎝x ↦ x ⎠˳(A)'),
    ('A.conjugate()', '\\left(A^{T}\\right)^{\\dagger}', '    +\n/ T\\ \n\\A / ', '    †\n⎛ T⎞ \n⎝A ⎠ '),
    ('DotProduct(v, w)', 'v \\cdot w', 'v*w', 'v⋅w'),
    ('Determinant(BlockMatrix([[P,Q],[Q,P]]))', '\\left|{\\begin{matrix}P & Q\\\\Q & P\\end{matrix}}\\right|', '|P  Q|\n|    |\n|Q  P|', '│P  Q│\n│    │\n│Q  P│'),
    ('Trace(BlockMatrix([[P,Q],[Q,P]]))', '\\operatorname{tr}\\left(\\left[\\begin{matrix}P & Q\\\\Q & P\\end{matrix}\\right] \\right)', '  /[P  Q]\\\ntr|[    ]|\n  \\[Q  P]/', '  ⎛⎡P  Q⎤⎞\ntr⎜⎢    ⎥⎟\n  ⎝⎣Q  P⎦⎠'),
    ('Trace(A)*Determinant(B)', '\\left|{B}\\right| \\operatorname{tr}\\left(A \\right)', '|B|*tr(A)', '│B│⋅tr(A)'),
    ('x*A[0,1]', 'x {A}_{0,1}', 'x*A_01', 'x⋅A₀₁'),
    ('A[0,1]**2', '{A}_{0,1}^{2}', '    2\nA_01 ', '   2\nA₀₁ '),
    ('Trace(A**2)**2', '\\operatorname{tr}\\left(A^{2} \\right)^{2}', '      2\n  / 2\\ \ntr\\A / ', '      2\n  ⎛ 2⎞ \ntr⎝A ⎠ '),
    ('sin(Trace(A))', '\\sin{\\left(\\operatorname{tr}\\left(A \\right) \\right)}', 'sin(tr(A))', 'sin(tr(A))'),
    ('BlockDiagMatrix(P, Q).T', '\\left[\\begin{matrix}P^{T} & 0\\\\0 & Q^{T}\\end{matrix}\\right]', '[ T    ]\n[P   0 ]\n[      ]\n[     T]\n[0   Q ]', '⎡ T    ⎤\n⎢P   𝟘 ⎥\n⎢      ⎥\n⎢     T⎥\n⎣𝟘   Q ⎦'),
    ('Inverse(BlockDiagMatrix(P, Q))', '\\left(\\left[\\begin{matrix}P & 0\\\\0 & Q\\end{matrix}\\right]\\right)^{-1}', '        -1\n/[P  0]\\  \n|[    ]|  \n\\[0  Q]/  ', '        -1\n⎛⎡P  𝟘⎤⎞  \n⎜⎢    ⎥⎟  \n⎝⎣𝟘  Q⎦⎠  '),
    ('BlockMatrix([[P,Q],[Q,P]]).T', '\\left[\\begin{matrix}P^{T} & Q^{T}\\\\Q^{T} & P^{T}\\end{matrix}\\right]', '[ T   T]\n[P   Q ]\n[      ]\n[ T   T]\n[Q   P ]', '⎡ T   T⎤\n⎢P   Q ⎥\n⎢      ⎥\n⎢ T   T⎥\n⎣Q   P ⎦'),
    ('Adjoint(BlockMatrix([[P,Q],[Q,P]]))', '\\left[\\begin{matrix}P & Q\\\\Q & P\\end{matrix}\\right]^{\\dagger}', '      +\n[P  Q] \n[    ] \n[Q  P] ', '      †\n⎡P  Q⎤ \n⎢    ⎥ \n⎣Q  P⎦ '),
    ('HadamardPower(A, x+1)', 'A^{\\circ \\left({x + 1}\\right)}', ' .(x + 1)\nA        ', ' ∘(x + 1)\nA        '),
    ('HadamardPower(A, -x)', 'A^{\\circ {- x}}', ' .(-x)\nA     ', ' ∘(-x)\nA     '),
    ('HadamardProduct(A, B.T)', 'A \\circ B^{T}', '    T\nA.*B ', '   T\nA∘B '),
    ('HadamardProduct(A.T, A**2)', 'A^{T} \\circ A^{2}', ' T   2\nA .*A ', ' T  2\nA ∘A '),
    ('KroneckerProduct(A*B, C)', '\\left(A B\\right) \\otimes C', '(A*B) x C', '(A⋅B) ⨂ C'),
    ('KroneckerProduct(A.T, C)', 'A^{T} \\otimes C', ' T    \nA  x C', ' T    \nA  ⨂ C'),
    ('A*B*KroneckerProduct(C,D)', 'A B \\left(C \\otimes D\\right)', 'A*B*(C x D)', 'A⋅B⋅(C ⨂ D)'),
    ('(x+1)*(A+B)', '\\left(x + 1\\right) \\left(A + B\\right)', '(x + 1)*(A + B)', '(x + 1)⋅(A + B)'),
    ('x*(A*B)', 'x A B', 'x*A*B', 'x⋅A⋅B'),
    ('X[0,:]', 'X\\left[:1, :\\right]', 'X[:1, :]', 'X[:1, :]'),
    ('X[1:,1:]**2', '\\left(X\\left[1:, 1:\\right]\\right)^{2}', '           2\n(X[1:, 1:]) ', '           2\n(X[1:, 1:]) '),
    ('X.T[0:2, 0:1]', '\\left(X^{T}\\right)\\left[:2, :1\\right]', '/ T\\        \n\\X /[:2, :1]', '⎛ T⎞        \n⎝X ⎠[:2, :1]'),
    ('Y[0,0]*Y[1,1] - Y[0,1]*Y[1,0]', '{Y}_{0,0} {Y}_{1,1} - {Y}_{0,1} {Y}_{1,0}', 'Y_00*Y_11 - Y_01*Y_10', 'Y₀₀⋅Y₁₁ - Y₀₁⋅Y₁₀'),
    ('A[x, 0]', '{A}_{x,0}', '(A)[x, 0]', '(A)[x, 0]'),
    ('(A.T)[0, 1]', '{A}_{1,0}', 'A_10', 'A₁₀'),
    ('(A+B).T[0,1]', '{A}_{1,0} + {B}_{1,0}', 'A_10 + B_10', 'A₁₀ + B₁₀'),
    ("MatrixSymbol('alpha_1', n, n)*A", '\\alpha_{1} A', 'alpha_1*A', 'α₁⋅A'),
    ("Transpose(MatrixSymbol('A^2', n, n))", '\\left(A^{2}\\right)^{T}', '   T\nA^2 ', '  T\nA² '),
    ("MatPow(MatrixSymbol('A^2', n, n), 2)", '\\left(A^{2}\\right)^{2}', '   2\nA^2 ', '  2\nA² '),
    ('Eq(A*B, C)', 'A B = C', 'A*B = C', 'A⋅B = C'),
    ('A*B + 2*Identity(n)', '2 \\mathbb{I} + A B', '2*I + A*B', '2⋅𝕀 + A⋅B'),
    ('A + ZeroMatrix(n,n)*B + OneMatrix(n,n)', 'A + 1', 'A + 1', 'A + 𝟙'),
    ('3*A*B - 2*C*D + A', 'A + 3 A B - 2 C D', '-2*C*D + 3*A*B + A', '-2⋅C⋅D + 3⋅A⋅B + A'),
    ('2*A.T*B.I', '2 A^{T} B^{-1}', '   T  -1\n2*A *B  ', '   T  -1\n2⋅A ⋅B  '),
    ('A**-2', 'A^{-2}', ' -2\nA  ', ' -2\nA  '),
    ('(2*A)**2', '\\left(2 A\\right)^{2}', '     2\n(2*A) ', '     2\n(2⋅A) '),
    ('(-A)**2', '\\left(- A\\right)^{2}', '    2\n(-A) ', '    2\n(-A) '),
    ('Inverse(-A)', '\\left(- A\\right)^{-1}', '    -1\n(-A)  ', '    -1\n(-A)  '),
    ('Transpose(-A)', '\\left(- A\\right)^{T}', '    T\n(-A) ', '    T\n(-A) '),
    ('Transpose(2*A)', '\\left(2 A\\right)^{T}', '     T\n(2*A) ', '     T\n(2⋅A) '),
]

# (source, latex settings, latex)
SETTINGS_CASES = [
    ('A*B', {'mat_symbol_style': 'bold'}, '\\mathbf{A} \\mathbf{B}'),
    ('Identity(n) + A', {'mat_symbol_style': 'bold'}, '\\mathbf{I} + \\mathbf{A}'),
    ('ZeroMatrix(n,n)', {'mat_symbol_style': 'bold'}, '\\mathbf{0}'),
    ('OneMatrix(n,n)', {'mat_symbol_style': 'bold'}, '\\mathbf{1}'),
    ("MatrixSymbol('A_1', n, n)", {'mat_symbol_style': 'bold'}, '\\mathbf{A}_{1}'),
    ('Adjoint(A)', {'adjoint_style': 'star'}, 'A^{\\ast}'),
    ('Adjoint(A*B)', {'adjoint_style': 'hermitian'}, '\\left(A B\\right)^{\\mathsf{H}}'),
    ('BlockMatrix([[P,Q],[Q,P]])', {'mat_delim': '('}, '\\left(\\begin{matrix}P & Q\\\\Q & P\\end{matrix}\\right)'),
    ('BlockMatrix([[P,Q],[Q,P]])', {'mat_delim': ''}, '\\begin{matrix}P & Q\\\\Q & P\\end{matrix}'),
    ('BlockMatrix([[P,Q],[Q,P]])', {'mat_str': 'array'}, '\\left[\\begin{array}{cc}P & Q\\\\Q & P\\end{array}\\right]'),
    ('BlockMatrix([[P,Q],[Q,P]])', {'mode': 'inline'}, '$\\left[\\begin{smallmatrix}P & Q\\\\Q & P\\end{smallmatrix}\\right]$'),
    ('A*B + C', {'mode': 'equation'}, '\\begin{equation}A B + C\\end{equation}'),
    ('Matrix([[1,2],[3,4]])', {'mat_delim': '('}, '\\left(\\begin{matrix}1 & 2\\\\3 & 4\\end{matrix}\\right)'),
]


@pytest.mark.parametrize("source, tex, ascii_form, unicode_form", CASES)
def test_matrix_expression_printing(source, tex, ascii_form, unicode_form):
    expr = eval(source, NS)
    assert S.latex(expr) == tex
    assert S.pretty(expr, use_unicode=False) == ascii_form
    assert S.pretty(expr, use_unicode=True) == unicode_form


@pytest.mark.parametrize("source, settings, tex", SETTINGS_CASES)
def test_matrix_expression_latex_settings(source, settings, tex):
    assert S.latex(eval(source, NS), **settings) == tex
