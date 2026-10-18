"""Finite crystallographic Cartan types, with exact root arithmetic."""
from collections import deque
from fractions import Fraction
from functools import lru_cache
import operator
import re

from .. import Matrix, Rational


class Standard_Cartan:
    def __init__(self, series, rank):
        self.series, self.n = series, rank

    def __repr__(self):
        return '%s%s' % (self.series,self.n)

    def __hash__(self):
        return hash((self.series,self.n))

    def __eq__(self,other):
        return isinstance(other,Standard_Cartan) and (self.series,self.n) == (other.series,other.n)

    def rank(self):
        return self.n

    def dimension(self):
        return self.n+1 if self.series == 'A' else 8 if self.series == 'E' else 3 if self.series == 'G' else self.n

    def _simple(self,index):
        index = operator.index(index)
        if not 1 <= index <= self.n:
            raise ValueError('simple-root index must lie between 1 and the rank')
        i,letter,n = index-1,self.series,self.n
        root = [Fraction(0)]*self.dimension()
        if letter == 'A' or letter in 'BCD' and index < n:
            root[i],root[i+1] = Fraction(1),Fraction(-1)
        elif letter in 'BC':
            root[-1] = Fraction(1 if letter == 'B' else 2)
        elif letter == 'D':
            root[-2:] = [Fraction(1),Fraction(1)]
        elif letter == 'E':
            if index == 1:
                root = [Fraction(1,2)]+[Fraction(-1,2)]*6+[Fraction(1,2)]
            elif index == 2:
                root[:2] = [Fraction(1),Fraction(1)]
            else:
                root[index-3],root[index-2] = Fraction(-1),Fraction(1)
        elif letter == 'F':
            if index < 3:
                root[index],root[index+1] = Fraction(1),Fraction(-1)
            elif index == 3:
                root[3] = Fraction(1)
            else:
                root = [Fraction(1,2)]+[Fraction(-1,2)]*3
        elif letter == 'G':
            root = list(map(Fraction,[0,1,-1] if index == 1 else [1,-2,1]))
        return tuple(root)

    def simple_root(self,i):
        return [Rational(v.numerator,v.denominator) for v in self._simple(i)]

    def basic_root(self,i,j):
        if i == j or not (0 <= i < self.dimension() and 0 <= j < self.dimension()):
            raise ValueError('invalid distinct ambient root indices')
        root = [0]*self.dimension()
        root[i],root[j] = (-1,1) if self.series == 'E' else (1,-1)
        return root

    def cartan_matrix(self):
        simple = [self._simple(i) for i in range(1,self.n+1)]
        norms = [sum(x*x for x in root) for root in simple]
        return Matrix([[int(2*sum(a*b for a,b in zip(row,col))/norms[j])
                        for j,col in enumerate(simple)] for row in simple])

    @lru_cache(maxsize=32)
    def _positive_coefficients(self):
        n = self.n
        matrix = [[int(v) for v in self.cartan_matrix().row(i)] for i in range(n)]
        roots = {tuple(int(i == j) for j in range(n)) for i in range(n)}
        queue = deque(roots)
        while queue:
            root = queue.popleft()
            for i in range(n):
                reflected = list(root)
                reflected[i] -= sum(root[j]*matrix[j][i] for j in range(n))
                reflected = tuple(reflected)
                if reflected not in roots:
                    roots.add(reflected)
                    queue.append(reflected)
        if len(roots) != self.roots():
            raise ArithmeticError('reflection closure does not have the expected finite root count')
        return tuple(sorted((r for r in roots if all(v >= 0 for v in r)),key=lambda r:(sum(r),r)))

    def positive_roots(self):
        simple = [self._simple(i) for i in range(1,self.n+1)]
        values = []
        for coefficients in self._positive_coefficients():
            ambient = [sum(c*root[j] for c,root in zip(coefficients,simple)) for j in range(self.dimension())]
            values.append([Rational(v.numerator,v.denominator) for v in ambient])
        if self.series in 'ABCD':
            ordered = []
            for i in range(self.dimension()-1):
                for j in range(i+1,self.dimension()):
                    root = self.basic_root(i,j)
                    ordered.append(root)
                    if self.series != 'A':
                        root = root[:]
                        root[j] = 1
                        ordered.append(root)
            if self.series in 'BC':
                for i in range(self.n):
                    root = [0]*self.n
                    root[i] = 1 if self.series == 'B' else 2
                    ordered.append(root)
            if {tuple(r) for r in ordered} != {tuple(r) for r in values}:
                raise ArithmeticError('classical root ordering does not match reflection closure')
            values = ordered
        return dict(enumerate(values,1))

    def highest_root(self):
        coefficients = self._positive_coefficients()[-1]
        return [sum(c*self.simple_root(i+1)[j] for i,c in enumerate(coefficients)) for j in range(self.dimension())]

    def roots(self):
        n,letter = self.n,self.series
        if letter == 'A':
            return n*(n+1)
        if letter in 'BC':
            return 2*n*n
        if letter == 'D':
            return 2*n*(n-1)
        return {'E':{6:72,7:126,8:240},'F':{4:48},'G':{2:12}}[letter][n]

    def basis(self):
        return self.roots()+self.n

    def lie_algebra(self):
        n = self.n
        return {'A':'su(%s)' % (n+1),'B':'so(%s)' % (2*n+1),
                'C':'sp(%s)' % (2*n),'D':'so(%s)' % (2*n)}.get(self.series,str(self))

    def dynkin_diagram(self):
        matrix = self.cartan_matrix()
        if self.series == 'E':
            return '        2\n        0\n        |\n        |\n'+'---'.join('0' for _ in range(self.n-1))+'\n1   '+'   '.join(str(i) for i in range(3,self.n+1))
        if self.series == 'D':
            indent = ' '*(4*(self.n-3))
            return indent+str(self.n-1)+'\n'+indent+'0\n'+indent+'|\n'+indent+'|\n'+'---'.join('0' for _ in range(self.n-1))+'\n'+'   '.join(map(str,range(1,self.n-1)))+'   '+str(self.n)
        nodes = '0'
        for i in range(self.n-1):
            a,b = -int(matrix[i,i+1]),-int(matrix[i+1,i])
            nodes += ('---' if a == b == 1 else '≡>≡' if a == 3 else '≡<≡' if b == 3 else '=>=' if a > b else '=<=')+'0'
        return nodes+'\n'+'   '.join(map(str,range(1,self.n+1)))


class CartanType_generator:
    def __call__(self, specification):
        if isinstance(specification,Standard_Cartan):
            return specification
        if isinstance(specification,str):
            match = re.fullmatch(r'([A-G])(\d+)',specification)
            if not match:
                raise ValueError('expected a finite Cartan type such as A3 or E8')
            letter,rank = match.group(1),int(match.group(2))
        elif isinstance(specification,(tuple,list)) and len(specification) == 2:
            letter,rank = specification[0],operator.index(specification[1])
        else:
            raise TypeError('Cartan type must be a string or (series, rank)')
        minimum = {'A':1,'B':2,'C':2,'D':3}
        valid = rank >= minimum[letter] if letter in minimum else letter == 'E' and rank in (6,7,8) or letter == 'F' and rank == 4 or letter == 'G' and rank == 2
        if not valid:
            raise ValueError('invalid rank for finite Cartan type %s%s' % (letter,rank))
        return Standard_Cartan(letter,rank)


CartanType = CartanType_generator()
