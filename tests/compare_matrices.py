"""Differential testing of symox matrices against SymPy.

Every case is a Python expression evaluated in a SymPy namespace and in a
symox namespace. We check

* ``str`` equality (Python ``str`` of the result, so tuples / lists / dicts
  of matrices use the multi-line ``repr`` of matrices like SymPy), and
* semantic equality when the strings differ (the symox output is parsed
  back with SymPy and compared entry-wise with ``simplify``).

Run directly for a report (``-v`` lists every case), or through pytest
(``test_matrices_against_sympy``).
"""

from __future__ import annotations

import sys

import sympy
import symox

CASES = r"""
# --- construction and printing ---
Matrix([[1, 2], [3, 4]])
repr(Matrix([[1, 2], [3, 4]]))
Matrix([1, 2, 3])
repr(Matrix([1, -2, 30]))
Matrix(2, 3, [1, 2, 3, 4, 5, 6])
Matrix(3, 3, lambda i, j: i + j)
Matrix(2, 2, lambda i, j: x**i*y**j)
Matrix([[x, y], [1, x + y]])
Matrix([[Rational(1, 2), x], [x**2, 1/(x + 1)]])
Matrix([])
Matrix(0, 3, [])
Matrix([[]])
eye(3)
zeros(2, 3)
ones(2)
ones(2, 3)
diag(1, 2, 3)
diag(Matrix([[1, 2], [3, 4]]), 5)
diag(x, y, z)
repr(Matrix([[x**2 + 1, -1], [Rational(-3, 4), x]]))
repr(Matrix([[1, 2, 3]]))
Matrix([[1.5, 2], [3, 4.25]])
Matrix([[I, 1 + I], [2, -I]])
Matrix(Matrix([[1, 2], [3, 4]]))
Matrix([[sin(x), cos(x)], [exp(x), log(x)]])
# --- shape and indexing ---
Matrix([[1, 2, 3], [4, 5, 6]]).shape
Matrix([[1, 2, 3], [4, 5, 6]]).rows
Matrix([[1, 2, 3], [4, 5, 6]]).cols
Matrix([[1, 2, 3], [4, 5, 6]])[0, 1]
Matrix([[1, 2, 3], [4, 5, 6]])[4]
Matrix([[1, 2, 3], [4, 5, 6]])[-1]
Matrix([[1, 2, 3], [4, 5, 6]])[-1, -1]
Matrix([[1, 2, 3], [4, 5, 6]])[:, 0]
Matrix([[1, 2, 3], [4, 5, 6]])[0, :]
Matrix([[1, 2, 3], [4, 5, 6]])[:, 1:]
Matrix([[1, 2, 3], [4, 5, 6]])[1:, ::2]
Matrix([[1, 2, 3], [4, 5, 6]])[:3]
Matrix([[1, 2, 3], [4, 5, 6]]).row(1)
Matrix([[1, 2, 3], [4, 5, 6]]).col(-1)
Matrix([[1, 2, 3], [4, 5, 6]]).tolist()
list(Matrix([[1, 2], [3, 4]]))
len(Matrix([[1, 2, 3], [4, 5, 6]]))
Matrix([[1, 2, 3], [4, 5, 6]]).extract([0, 1], [0, 2])
Matrix([[1, 2, 3], [4, 5, 6]]).reshape(3, 2)
Matrix([[1, 2], [3, 4]]).row_insert(1, Matrix([[0, 0]]))
Matrix([[1, 2], [3, 4]]).col_insert(0, Matrix([7, 8]))
Matrix.hstack(Matrix([[1], [2]]), Matrix([[3, 4], [5, 6]]))
Matrix.vstack(Matrix([[1, 2]]), Matrix([[3, 4], [5, 6]]))
Matrix([[1, 2]]).row_join(Matrix([[3]]))
Matrix([[1, 2]]).col_join(Matrix([[3, 4]]))
Matrix([[1, 2], [3, 4]]).T
Matrix([[1, I], [2 + I, x]]).H
Matrix([[1, I], [2 + I, 3]]).C
Matrix([[1, 2], [3, 4]]).trace()
Matrix([[x, 1], [y, x**2]]).trace()
Matrix([[1, 2], [2, 1]]).is_symmetric()
Matrix([[1, 2], [3, 1]]).is_symmetric()
Matrix([[0, 0], [0, 0]]).is_zero_matrix
Matrix([[1, 2], [3, 4]]) == Matrix([[1, 2], [3, 4]])
Matrix([[1, 2], [3, 4]]) == Matrix([[1, 2], [3, 5]])
Matrix([[1, 2], [3, 4]]) != Matrix([[1, 2], [3, 5]])
# --- arithmetic ---
Matrix([[1, 2], [3, 4]]) + Matrix([[x, 1], [0, y]])
Matrix([[1, 2], [3, 4]]) - Matrix([[x, 1], [0, y]])
2*Matrix([[1, 2], [3, 4]])
Matrix([[x, y]])*3
x*Matrix([[1, 2], [3, x]])
Matrix([[x, y]])/2
Matrix([[x + 1, 2]])/x
-Matrix([[1, -x], [y, 0]])
Matrix([[1, 2], [3, 4]])*Matrix([[5, 6], [7, 8]])
Matrix([[1, 2], [3, 4]]) @ Matrix([5, 6])
Matrix([[a, b], [c, d]])*Matrix([[x + 1, 1], [1, x]])
Matrix([[x + 1, 1], [1, x]])*Matrix([[x + 1, 1], [1, x]])
Matrix([[a, b], [c, d]])**2
Matrix([[x + 1, 1], [1, x]])**2
Matrix([[x + 1, 1], [1, x]])**3
Matrix([[1, 1], [1, 0]])**10
Matrix([[1, 2], [3, 4]])**0
Matrix([[1, 2], [3, 4]])**-1
Matrix([[1, 2], [3, 4]])**-2
Matrix([[2, 0], [0, x]])**5
Matrix([[x, 1], [0, x]])**4
Matrix([[Rational(1, 2), 1], [0, 3]])**3
Matrix([[1, 2], [3, 4]]).multiply_elementwise(Matrix([[x, y], [z, 1]]))
# --- element-wise ---
Matrix([[x, x**2], [sin(x), 1]]).diff(x)
Matrix([[x*y, x**2*y**3]]).diff(x, y)
Matrix([[x, y], [x*y, 1]]).subs(x, 2)
Matrix([[x, y], [x*y, 1]]).subs({x: y, y: 3})
Matrix([[(x + 1)**2, x*(x - 1)]]).expand()
Matrix([[x, 2*x]]).applyfunc(lambda e: e**2)
Matrix([[pi, sqrt(2)], [Rational(1, 3), E]]).evalf()
Matrix([[x, y], [1, z]]).free_symbols == {x, y, z}
# --- determinants ---
Matrix([[5]]).det()
Matrix([[1, 2], [3, 4]]).det()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).det()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).det()
Matrix(4, 4, lambda i, j: (i + 1)**j).det()
Matrix(5, 5, lambda i, j: (3*i + 7*j) % 11 - 5).det()
Matrix(6, 6, lambda i, j: Rational(1, i + j + 1)).det()
Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]).det()
Matrix([[Rational(1, 2), Rational(1, 3)], [Rational(1, 4), Rational(1, 5)]]).det()
Matrix([[a, b], [c, d]]).det()
Matrix([[a, b, c], [d, e, f], [g, h, k]]).det()
Matrix([[x, y, 1], [1, x, y], [y, 1, x]]).det()
Matrix(4, 4, lambda i, j: x**(i + j) + i).det()
Matrix([[x, 1, 0, 0], [1, x, 0, 0], [0, 0, x, 2], [0, 0, 3, x]]).det()
Matrix(4, 4, lambda i, j: x + i*j if i != j else y).det()
Matrix(4, 4, lambda i, j: symbols('a%d%d' % (i, j))).det()
Matrix(5, 5, lambda i, j: (x + i)**(j % 2) + j).det()
Matrix([[x, 1, 0, 2], [0, x, 1, 0], [1, 0, x, 1], [2, 1, 0, x]]).det()
Matrix([[1/x, 1], [1, x]]).det()
Matrix([[1/x, y], [x, 1/y]]).det()
Matrix([[sin(x), cos(x)], [-cos(x), sin(x)]]).det()
Matrix([[x, 1, 2, 3], [0, y, 1, 0], [0, 0, z, 1], [0, 0, 0, x*y]]).det()
Matrix([[1, 2], [3, 4]]).det(method='berkowitz')
Matrix(4, 4, lambda i, j: x**(i + j) + i).det(method='berkowitz')
Matrix(4, 4, lambda i, j: (i + 2)**j).det(method='lu')
Matrix([[1.5, 2], [3, 4]]).det()
# --- inverses ---
Matrix([[1, 2], [3, 4]]).inv()
Matrix([[2, 0, 1], [1, 3, 2], [1, 1, 1]]).inv()
Matrix(4, 4, lambda i, j: Rational(1, i + j + 1)).inv()
Matrix([[Rational(1, 2), 3], [Rational(2, 7), 1]]).inv()
Matrix([[a, b], [c, d]]).inv()
Matrix([[x, 1], [2, x]]).inv()
Matrix([[2*x, 1], [2, x]]).inv()
Matrix([[x/2, 1], [3, x]]).inv()
Matrix([[1, x], [1, y]]).inv()
Matrix([[1/x, 1], [1, x + 1]]).inv()
Matrix([[x, y, 1], [1, x, y], [y, 1, x]]).inv()
Matrix([[1, x, x**2], [1, y, y**2], [1, z, z**2]]).inv()
Matrix([[x, 1, 0], [0, x, 1], [0, 0, x]]).inv()
Matrix([[1, 2], [3, 4]]).inv(method='LU')
Matrix([[1, 2], [3, 4]]).inv(method='ADJ')
Matrix([[a, b], [c, d]]).inv(method='LU')
Matrix([[sin(x), 1], [1, x]]).inv()
Matrix([[1, 2], [3, 4]]).adjugate()
Matrix([[a, b], [c, d]]).adjugate()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).adjugate()
Matrix([[x, 1, 0], [y, 2, 1], [0, 1, x]]).adjugate()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).cofactor(0, 1)
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).minor(1, 1)
Matrix([[x, 1, 0], [y, 2, 1], [0, 1, x]]).cofactor_matrix()
# --- row reduction ---
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).rref()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rref()
Matrix([[0, 2, 4], [1, 1, 1], [2, 4, 6]]).rref()
Matrix([[Rational(1, 2), 1, 3], [2, Rational(1, 3), 4]]).rref()
Matrix([[1, 2, 1, 0], [2, 4, 0, 1]]).rref()
Matrix(4, 6, lambda i, j: (i*j + i + 2*j) % 5).rref()
Matrix([[x, 1, 2], [1, y, 3]]).rref()
Matrix([[x, x + 1]]).rref()
Matrix([[a, b], [c, d]]).rref()
Matrix([[1, x], [x, 1]]).rref()
Matrix([[x, 1], [x**2, x]]).rref()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rank()
Matrix([[1, 2], [2, 4]]).rank()
Matrix([[x, 1], [x, 1]]).rank()
Matrix([[x, y], [y, x]]).rank()
Matrix(5, 5, lambda i, j: i*j).rank()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).nullspace()
Matrix([[1, 2, 1, 0], [2, 4, 0, 1]]).nullspace()
Matrix([[1, 1, 1], [1, 1, 1]]).nullspace()
Matrix([[x, 1, 2], [1, y, 3]]).nullspace()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]]).nullspace()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).columnspace()
Matrix([[1, 2, 3], [2, 4, 6]]).columnspace()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rowspace()
Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).echelon_form()
# --- LU and linear systems ---
Matrix([[1, 2], [3, 4]]).LUdecomposition()
Matrix([[0, 1], [1, 0]]).LUdecomposition()
Matrix([[a, b], [c, d]]).LUdecomposition()
Matrix([[2, 1, 1], [4, -6, 0], [-2, 7, 2]]).LUdecomposition()
Matrix([[1, 2], [3, 4]]).LUsolve(Matrix([5, 6]))
Matrix([[a, b], [c, d]]).LUsolve(Matrix([1, 2]))
Matrix([[2, 1, 1], [4, -6, 0], [-2, 7, 2]]).LUsolve(Matrix([5, -2, 9]))
Matrix([[1, 2], [3, 4]]).solve(Matrix([5, 6]))
Matrix([[a, b], [c, d]]).solve(Matrix([1, 2]))
Matrix([[2, 1, 1], [4, -6, 0], [-2, 7, 2]]).solve(Matrix([5, -2, 9]))
Matrix([[1, 2], [3, 4]]).solve(Matrix([[1, 0], [0, 1]]))
Matrix([[x, 1], [1, x]]).solve(Matrix([1, 0]))
Matrix([[1, 2], [2, 4]]).gauss_jordan_solve(Matrix([3, 6]))
# --- characteristic polynomials and eigenvalues ---
Matrix([[1, 2], [3, 4]]).charpoly(x)
Matrix([[1, 2], [3, 4]]).charpoly()
Matrix([[1, 2], [3, 4]]).charpoly(x).as_expr()
Matrix([[x, 1], [2, x]]).charpoly(y)
Matrix([[x, 1], [2, x]]).charpoly(y).as_expr()
Matrix([[x, 1], [2, x]]).charpoly(x).as_expr()
Matrix([[Rational(1, 2), 1], [3, 4]]).charpoly(x)
Matrix([[1, 2, 0], [0, 3, 1], [1, 0, 2]]).charpoly(x).as_expr()
Matrix(4, 4, lambda i, j: (i + 1)*(j + 2) % 5).charpoly(x).as_expr()
Matrix([[a, b], [c, d]]).charpoly(x).as_expr()
Matrix([[1, 2], [3, 4]]).eigenvals()
Matrix([[2, 1], [1, 2]]).eigenvals()
Matrix([[1, 1], [0, 1]]).eigenvals()
Matrix([[0, -1], [1, 0]]).eigenvals()
Matrix([[2, 0, 0], [0, 3, 4], [0, 4, 9]]).eigenvals()
Matrix([[1, 2, 3], [0, 4, 5], [0, 0, 6]]).eigenvals()
Matrix([[6, -2, 2], [-2, 3, -1], [2, -1, 3]]).eigenvals()
Matrix([[1, 1, 0], [1, 1, 0], [0, 0, 2]]).eigenvals()
Matrix([[4, 1], [2, 3]]).eigenvals()
Matrix([[0, 1], [-2, -3]]).eigenvals()
Matrix([[1, 2], [2, 1]]).eigenvals()
Matrix([[3, -1, 0], [-1, 3, 0], [0, 0, 1]]).eigenvals()
Matrix([[Rational(1, 2), 1], [1, Rational(1, 3)]]).eigenvals()
Matrix([[x, 1], [1, x]]).eigenvals()
Matrix([[1, 2], [3, 4]]).eigenvects()
Matrix([[2, 1], [1, 2]]).eigenvects()
Matrix([[1, 1], [0, 1]]).eigenvects()
Matrix([[0, -1], [1, 0]]).eigenvects()
Matrix([[2, 0, 0], [0, 3, 4], [0, 4, 9]]).eigenvects()
Matrix([[6, -2, 2], [-2, 3, -1], [2, -1, 3]]).eigenvects()
Matrix([[1, 1, 0], [1, 1, 0], [0, 0, 2]]).eigenvects()
Matrix([[4, 1], [2, 3]]).eigenvects()
Matrix([[2, 0], [0, 2]]).eigenvects()
# --- calculus and vectors ---
Matrix([x*y, x + y]).jacobian([x, y])
Matrix([x**2*y, sin(x), exp(x*y)]).jacobian(Matrix([x, y]))
jacobian([x*y*z, x**2 + y], [x, y, z])
hessian(x**2*y + y**3*z, (x, y, z))
hessian(sin(x*y), [x, y])
wronskian([1, x, x**2], x)
wronskian([sin(x), cos(x)], x)
wronskian([exp(x), exp(2*x)], x)
Matrix([3, 4]).norm()
Matrix([[1, 2], [3, 4]]).norm()
Matrix([x, y]).norm()
Matrix([1, 2, 3]).dot(Matrix([4, 5, 6]))
Matrix([x, y, z]).dot(Matrix([1, 2, 3]))
Matrix([1, 0, 0]).cross(Matrix([0, 1, 0]))
Matrix([[x, y, z]]).cross(Matrix([[1, 2, 3]]))
# --- harder cases ---
Matrix(5, 5, lambda i, j: x**((i*j) % 3) + i - j).det()
Matrix(5, 5, lambda i, j: (x + y)**((i + j) % 2) - i*j).det()
Matrix(4, 4, lambda i, j: 1/(x + i + j)).det()
Matrix([[1/x, 1, 0, 0], [1, 1/y, 1, 0], [0, 1, 1/x, 1], [0, 0, 1, 1/y]]).det()
Matrix([[x, 1, 0], [1, x, 1], [0, 1, x]]).inv()
Matrix([[1, 1/x], [1/y, 1]]).inv()
Matrix([[x + 1, x], [x, x - 1]]).inv()
Matrix([[x**2, x], [1, 1]]).inv()
Matrix([[2, x], [x, 2]]).inv()
Matrix([[Rational(1, 2)*x, 1], [1, x]]).inv()
Matrix([[a, 0, 0], [0, b, 0], [0, 0, c]]).inv()
Matrix([[a, b, 0], [0, c, d], [0, 0, e]]).inv()
Matrix(3, 3, lambda i, j: x**(i*j)).inv()
Matrix([[x, y], [y, x]]).inv()
Matrix([[x - 1, 2], [3, x + 1]]).inv()
Matrix([[sqrt(2), 1], [1, sqrt(2)]]).inv()
Matrix([[1, 2], [3, 4]]).evalf().inv()
Matrix([[I, 1], [1, I]]).inv()
Matrix([[1, I], [-I, 2]]).det()
Matrix([[x, 0, 1], [0, x, 0], [1, 0, x]]).rref()
Matrix([[x, y, 1], [1, x, y]]).rref()
Matrix([[1, x, 0], [0, 0, 1], [0, 0, x]]).rref()
Matrix([[x, 0], [0, 0], [0, y]]).rref()
Matrix([[x, x**2], [1, x]]).rref()
Matrix([[x - 1, 1], [1, x - 1]]).rref()
Matrix([[1, 2.5], [3, 4]]).rref()
Matrix([[x, 1, 1], [1, x, 1], [1, 1, x]]).nullspace()
Matrix([[x, x], [y, y]]).nullspace()
Matrix([[1, x, x**2], [1, x, x**2]]).nullspace()
Matrix([[1, x], [x, x**2]]).nullspace()
Matrix([[x, y], [1, 1]]).columnspace()
Matrix([[x, 1], [x**2, x]]).rank()
Matrix([[x, 1, 2], [1, y, 3], [2, 3, 4]]).rank()
Matrix([[x, y, 1], [x, y, 1], [1, 2, 3]]).rank()
Matrix([[x, 1], [1, y]]).LUdecomposition()
Matrix([[x, 1, 0], [y, 2, 1], [0, 1, x]]).LUdecomposition()
Matrix([[0, x], [y, 1]]).LUdecomposition()
Matrix([[x, 1, 0], [y, 2, 1], [0, 1, x]]).LUsolve(Matrix([1, 0, 0]))
Matrix([[x, 1, 0], [y, 2, 1], [0, 1, x]]).solve(Matrix([1, 0, 0]))
Matrix([[x, 1], [1, y]]).solve(Matrix([a, b]))
Matrix([[1, 1, 1], [1, 2, 4], [1, 3, 9]]).solve(Matrix([a, b, c]))
Matrix([[3, 1], [1, 3]]).eigenvects()
Matrix([[1, 4], [1, 1]]).eigenvects()
Matrix([[5, 4, 2], [4, 5, 2], [2, 2, 2]]).eigenvects()
Matrix([[2, 1, 0], [0, 2, 0], [0, 0, 3]]).eigenvects()
Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]]).eigenvals()
Matrix([[1, 2, 0], [2, 1, 0], [0, 0, 5]]).eigenvects()
Matrix([[1, 1], [-1, 1]]).eigenvects()
Matrix([[0, 2], [1, 0]]).eigenvects()
Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]]).eigenvals()
Matrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]]).eigenvects()
Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]]).eigenvals()
Matrix([[4, -2], [1, 1]]).eigenvals()
Matrix([[Rational(1, 2), Rational(1, 3)], [Rational(1, 4), Rational(1, 5)]]).eigenvals()
Matrix([[Rational(1, 2), Rational(1, 3)], [Rational(1, 4), Rational(1, 5)]]).eigenvects()
Matrix([[x, 0], [0, y]]).eigenvals()
Matrix([[a, b], [b, a]]).eigenvals()
Matrix([[a, b], [c, d]]).eigenvals()
Matrix([[x, 1], [0, y]]).eigenvects()
Matrix([[1, x], [x, 1]]).eigenvals()
Matrix([[0, 1], [x, 0]]).eigenvals()
Matrix([[x, 1], [1, x]]).eigenvects()
Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).eigenvals()
Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).eigenvects()
Matrix([[0, 0], [0, 0]]).eigenvects()
Matrix([[x, y], [z, x]]).charpoly(a).as_expr()
Matrix([[1/x, 1], [1, x]]).charpoly(y)
Matrix([[1/x, 1], [1, x]]).charpoly(y).as_expr()
Matrix([[x, 1], [1, y]]).charpoly(z)
Matrix([[sin(x), 1], [1, x]]).charpoly(y).as_expr()
Matrix(3, 3, lambda i, j: x**i + j).charpoly(y).as_expr()
Matrix([[1, 2], [3, 4]]).charpoly(x).all_coeffs()
Matrix([[x, 1], [1, x]])**-1
Matrix([[1, 2], [3, 4]])**-3
Matrix([[x, 1], [1, 0]])**5
Matrix([[sin(x), cos(x)], [-cos(x), sin(x)]])**2
Matrix([[x/2, 1], [1, 1/x]])**2
Matrix([[1, x], [0, 1]])**20
Matrix([[1, 2], [3, 4]]).adjugate()*Matrix([[1, 2], [3, 4]])
Matrix(4, 4, lambda i, j: x**(i + j) + i).adjugate()
Matrix([[x, 1, 0, 0], [0, y, 1, 0], [0, 0, z, 1], [1, 0, 0, x]]).det()
Matrix([[x, 1, 0, 0], [0, y, 1, 0], [0, 0, z, 1], [1, 0, 0, x]]).det(method='berkowitz')
Matrix([[x, 1, 0, 0], [0, y, 1, 0], [0, 0, z, 1], [1, 0, 0, x]]).det(method='lu')
Matrix([[a, b, c, d], [b, a, d, c], [c, d, a, b], [d, c, b, a]]).det()
Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]]).det()
Matrix([[x, 2], [3, 4]]).det()
Matrix([[x**2, y], [y, 1/x]]).det()
Matrix([[sqrt(2), 1, 0, 0], [1, sqrt(2), 1, 0], [0, 1, sqrt(2), 1], [0, 0, 1, sqrt(2)]]).det()
Matrix([[1.0, 2.0, 3.0, 1.0], [4.0, 5.0, 6.0, 2.0], [7.0, 8.0, 10.0, 1.0], [1.0, 1.0, 1.0, 5.0]]).det()
Matrix([[1, 2], [3, 4]]).norm()
Matrix([[x, 1], [y, 2]]).norm()
Matrix([[x, y, z]]).dot(Matrix([[a, b, c]]))
Matrix([1, I]).dot(Matrix([1, I]))
Matrix([[1, 2], [3, 4]]).dot(Matrix([1, 2]))
Matrix([[x, 2*x], [y, x*y]]).jacobian([x, y])
Matrix([[1, 2], [3, 4]]).row(-1)
Matrix([[1, 2, 3], [4, 5, 6]])[[0, 1], [0, 2]]
Matrix([[1, 2, 3], [4, 5, 6]])[1, [0, 2]]
Matrix([[1, 2], [3, 4]]).diagonal()
Matrix([[1, 2], [3, 4]]).minor_submatrix(0, 1)
Matrix([[1, 2], [3, 4]]).is_square
Matrix([[1, 2, 3], [4, 5, 6]]).is_square
Matrix([[1, 2], [3, 4]]).strongly_connected_components()
Matrix([[1, 0, 1], [0, 2, 0], [1, 0, 3]]).strongly_connected_components()
Matrix([[1, 2], [3, 4]]).T.T == Matrix([[1, 2], [3, 4]])
Matrix([[x, y], [y, x]]).subs(x, 1).det()
Matrix([[exp(x), exp(2*x)], [exp(x), 2*exp(2*x)]]).det()
Matrix([[exp(x), exp(-x)], [exp(x), -exp(-x)]]).det()
Matrix([[cos(x), -sin(x)], [sin(x), cos(x)]]).inv()
Matrix([[cos(x), -sin(x)], [sin(x), cos(x)]]).det()
Matrix([[x, sin(x)], [cos(x), 1]]).det()
Matrix([[x, sin(x)], [cos(x), 1]]).rref()
Matrix([[1, x], [x, sin(x)]]).inv()
Matrix(3, 3, lambda i, j: sin(x)**(i + j)).det()
Matrix(4, 4, lambda i, j: Rational(i + 1, j + 2)).rref()
Matrix(3, 4, lambda i, j: (i + 1)**j).rref()
Matrix(3, 3, lambda i, j: x**(i + j)).rref()
Matrix(2, 3, lambda i, j: x*i + y*j + 1).rref()
Matrix(2, 3, lambda i, j: x*i + y*j + 1).nullspace()
""".strip().splitlines()

CASES = [c for c in CASES if c.strip() and not c.lstrip().startswith("#")]

# Cases where symox deliberately or knowingly differs (documented).
KNOWN_DEVIATIONS: set[str] = set()


def namespace(mod):
    ns = {name: getattr(mod, name) for name in dir(mod) if not name.startswith("_")}
    names = "x y z a b c d e f g h k"
    ns.update(zip(names.split(), mod.symbols(names)))
    ns["Rational"] = mod.Rational
    if mod is sympy:
        # SymPy has no top-level jacobian
        ns["jacobian"] = lambda f, v: sympy.Matrix(f).jacobian(v)
    return ns


def _to_sympy(obj):
    """Best effort conversion of a symox result to SymPy."""
    if isinstance(obj, symox.Matrix):
        return obj.to_sympy()
    if isinstance(obj, symox.Expr):
        return obj.to_sympy()
    if isinstance(obj, tuple):
        return tuple(_to_sympy(o) for o in obj)
    if isinstance(obj, list):
        return [_to_sympy(o) for o in obj]
    if isinstance(obj, dict):
        return {_to_sympy(k): _to_sympy(v) for k, v in obj.items()}
    if isinstance(obj, symox.PurePoly):
        return obj.as_expr().to_sympy()
    return obj


def _equiv(a, b):
    if isinstance(a, sympy.PurePoly):
        a = a.as_expr()
    if isinstance(a, (tuple, list)):
        return type(a) is type(b) and len(a) == len(b) and all(_equiv(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        if not isinstance(b, dict) or len(a) != len(b):
            return False
        # keys may differ structurally: match by simplification
        rest = list(b.items())
        for k, v in a.items():
            for i, (k2, v2) in enumerate(rest):
                if _equiv(k, k2) and _equiv(v, v2):
                    rest.pop(i)
                    break
            else:
                return False
        return True
    if isinstance(a, sympy.MatrixBase):
        if not isinstance(b, sympy.MatrixBase) or a.shape != b.shape:
            return False
        return all(sympy.simplify(p - q) == 0 for p, q in zip(a, b))
    if isinstance(a, sympy.Basic) and isinstance(b, sympy.Basic):
        return a == b or sympy.simplify(a - b) == 0
    return a == b


def compare(case):
    ns_s = namespace(sympy)
    ns_x = namespace(symox)
    try:
        exp_s = eval(case, ns_s)
    except Exception as e:  # noqa: BLE001
        exp_s = e
    try:
        exp_x = eval(case, ns_x)
    except Exception as e:  # noqa: BLE001
        if isinstance(exp_s, Exception):
            return "same-error", type(e).__name__, type(exp_s).__name__
        return "error", f"{type(e).__name__}: {e}", str(exp_s)
    if isinstance(exp_s, Exception):
        return "DIFFERENT", str(exp_x), f"{type(exp_s).__name__}: {exp_s}"
    s_s, s_x = str(exp_s), str(exp_x)
    if s_s == s_x:
        return "same", s_x, s_s
    try:
        if _equiv(exp_s, _to_sympy(exp_x)):
            return "equivalent", s_x, s_s
    except Exception:  # noqa: BLE001
        pass
    return "DIFFERENT", s_x, s_s


def run(verbose=True, show_all=False):
    stats = {}
    bad = []
    for case in CASES:
        status, got, want = compare(case)
        stats[status] = stats.get(status, 0) + 1
        if status not in ("same", "same-error") or show_all:
            if verbose:
                g = got.replace("\n", "")
                w = want.replace("\n", "")
                print(f"[{status:10s}] {case}\n      symox: {g}\n      sympy: {w}")
            if status in ("DIFFERENT", "error") and case not in KNOWN_DEVIATIONS:
                bad.append((case, got, want))
    if verbose:
        print("\nsummary:", stats, f"(total {len(CASES)})")
    return bad


def test_matrices_against_sympy():
    bad = run(verbose=False)
    assert not bad, "\n".join(f"{c}: symox={g} sympy={w}" for c, g, w in bad)


if __name__ == "__main__":
    sys.exit(1 if run(show_all="-v" in sys.argv) else 0)
