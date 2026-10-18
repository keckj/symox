"""Deferred vectors, indexed contractions and matrix scalar entry points."""
from collections import Counter
import operator

from . import _core
from ._compat import _ExprType

_indexed = {}
_indices = {}
_deferred = {}


class DeferredVector:
    def __init__(self,name,**assumptions):
        from . import Symbol
        self._symbol = Symbol(name,**assumptions)
        self.name = str(name)

    def __repr__(self):
        return "DeferredVector(%r)" % self.name

    def __str__(self):
        return self.name

    def _symox_(self):
        return self._symbol

    def __hash__(self):
        return hash((DeferredVector,self._symbol))

    def __eq__(self,other):
        return isinstance(other,DeferredVector) and self._symbol == other._symbol

    def __getitem__(self,index):
        from . import Symbol
        index = operator.index(index)
        if index < 0:
            raise IndexError('deferred vector index out of range')
        value = Symbol('%s[%s]' % (self.name,index))
        _deferred[value] = self,index
        return value


class Idx:
    def __init__(self,label,range=None,**kwargs):
        from . import Symbol,sympify,Integer
        if kwargs:
            raise TypeError('unexpected index options')
        self.label = Symbol(label,integer=True) if isinstance(label,str) else sympify(label)
        if self.label.is_integer is False:
            raise TypeError('index labels must be integer-valued')
        if range is None:
            self.lower = self.upper = None
        elif isinstance(range,(tuple,list)) or getattr(range,'is_Tuple',False):
            values = range.args if getattr(range,'is_Tuple',False) else range
            if len(values) != 2:
                raise ValueError('index bounds must be (lower, upper)')
            self.lower,self.upper = map(sympify,values)
        else:
            self.lower,self.upper = Integer(0),sympify(range)-1
        if self.lower is not None:
            if self.lower.is_integer is False or self.upper.is_integer is False:
                raise TypeError('index bounds must be integers')
            if self.lower.is_number and self.upper.is_number and self.upper < self.lower:
                raise ValueError('index upper bound is below its lower bound')
        _indices[self.label] = self

    @property
    def args(self):
        from . import Tuple
        return (self.label,) if self.lower is None else (self.label,Tuple(self.lower,self.upper))

    @property
    def free_symbols(self):
        return {self}

    def _symox_(self):
        return self.label

    def __repr__(self):
        return str(self.label)

    def __hash__(self):
        return hash((Idx,self.args))

    def __eq__(self,other):
        return isinstance(other,Idx) and self.args == other.args

    def __getattr__(self,name):
        return getattr(self.label,name)


class IndexConformanceException(ValueError):
    """Outer Einstein indices disagree or occur more than twice."""


def _index_result(expr):
    from .codegen.ast import Indexed
    from . import sympify
    indexed = expr if isinstance(expr,Indexed) else _indexed.get(expr)
    if indexed is not None:
        counts = Counter(indexed.indices)
    else:
        expr = sympify(expr)
        if expr.is_Add:
            results = [_index_result(a)[0] for a in expr.args]
            nonscalar = [r for r in results if r]
            if nonscalar and any(r != nonscalar[0] for r in nonscalar):
                raise IndexConformanceException('indices are not consistent across sum terms')
            return nonscalar[0] if nonscalar else set(),()
        if expr.is_Mul:
            counts = Counter(i for a in expr.args for i in _index_result(a)[0])
        elif expr.is_Pow:
            base,exponent = expr.args
            first,second = _index_result(base)[0],_index_result(exponent)[0]
            return first | second,()
        elif expr.is_Function:
            results = [_index_result(a)[0] for a in expr.args if isinstance(a,_core.Expr)]
            return set().union(*results),()
        else:
            return set(),()
    if any(count > 2 for count in counts.values()):
        raise IndexConformanceException('an Einstein index occurs more than twice')
    return {i for i,n in counts.items() if n == 1},tuple(sorted((i for i,n in counts.items() if n == 2),key=str))


def get_indices(expr):
    return {_indices.get(i,i) for i in _index_result(expr)[0]},{}


def get_contraction_structure(expr):
    from . import sympify
    from .codegen.ast import Indexed
    node = expr if isinstance(expr,Indexed) else sympify(expr)
    terms = node.args if not isinstance(node,Indexed) and node.is_Add else (node,)
    result = {}
    for term in terms:
        dummy = _index_result(term)[1]
        key = tuple(_indices.get(i,i) for i in dummy) or None
        result.setdefault(key,set()).add(term)
        nested = []
        if not isinstance(term,Indexed):
            for child in term.args:
                if isinstance(child,_core.Expr) and child.args:
                    structure = get_contraction_structure(child)
                    if any(k is not None for k in structure):
                        nested.append(structure)
        if nested:
            result[term] = nested
    return result


class Permanent(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr.class_name == 'Permanent')

    def __new__(cls,mat):
        from . import MatrixExpr
        return _core._apply('Permanent',MatrixExpr(mat))


def per(matexpr):
    from . import Matrix
    if isinstance(matexpr,Matrix) or getattr(matexpr,'_symox_sparse',False):
        return matexpr.per()
    return Permanent(matexpr)


def diagonalize_vector(vector):
    from . import Matrix,DiagMatrix,diag
    if isinstance(vector,Matrix):
        if 1 not in vector.shape:
            raise ValueError('expected a row or column vector')
        return diag(*tuple(vector))
    return DiagMatrix(vector)


def matrix_symbols(expr):
    from . import MatrixSymbol
    result,seen = set(),set()
    def visit(value):
        if isinstance(value,MatrixSymbol):
            result.add(value)
        elif isinstance(value,(tuple,list)) or getattr(value,'_symox_blockgrid',False):
            for part in value:
                visit(part)
        else:
            key = id(value)
            if key in seen:
                return
            seen.add(key)
            for part in getattr(value,'args',()):
                visit(part)
    visit(expr)
    return result


def blockcut(expr,rowsizes,colsizes):
    from . import MatrixExpr,BlockMatrix,sympify,simplify
    matrix = MatrixExpr(expr)
    rows,cols = tuple(map(sympify,rowsizes)),tuple(map(sympify,colsizes))
    if not rows or not cols or any(n.is_integer is False or n.is_nonpositive is True for n in rows+cols):
        raise ValueError('block sizes must be positive integers')
    if simplify(sum(rows)-matrix.rows) != 0 or simplify(sum(cols)-matrix.cols) != 0:
        raise ValueError('block sizes do not cover the matrix')
    blocks,start = [],sympify(0)
    for size in rows:
        row,left = [],sympify(0)
        for width in cols:
            row.append(matrix[start:start+size,left:left+width])
            left += width
        blocks.append(row)
        start += size
    return BlockMatrix(blocks)


def install(namespace):
    from functools import wraps
    from .codegen.ast import Indexed
    if getattr(Indexed,'_index_registry_installed',False):
        return
    original = Indexed._symox_
    def expression(self):
        try:
            value = original(self)
        except ValueError:
            if self.base.shape is not None or len(self.indices) == 1:raise
            from . import Symbol
            value = Symbol(str(self))
        _indexed[value] = self
        return value
    Indexed._symox_ = expression
    Indexed._index_registry_installed = True
    old_lambdify = namespace['lambdify']
    @wraps(old_lambdify)
    def lambdify(args,expr,*rest,**kwargs):
        arguments = [args] if isinstance(args,DeferredVector) else list(args) if isinstance(args,(tuple,list)) else [args]
        if not any(isinstance(a,DeferredVector) for a in arguments):
            return old_lambdify(args,expr,*rest,**kwargs)
        def symbols(value):
            if isinstance(value,(list,tuple)):
                return set().union(*(symbols(v) for v in value))
            return value.free_symbols
        free = symbols(expr)
        flat,selectors = [],[]
        for position,arg in enumerate(arguments):
            if isinstance(arg,DeferredVector):
                components = sorted(((value,owner[1]) for value,owner in _deferred.items() if owner[0] == arg and value in free),key=lambda item:item[1])
                for value,index in components:
                    flat.append(value)
                    selectors.append((position,index))
            else:
                flat.append(arg)
                selectors.append((position,None))
        compiled = old_lambdify(flat,expr,*rest,**kwargs)
        def call(*values):
            if len(values) != len(arguments):
                raise TypeError('wrong number of numerical arguments')
            return compiled(*(values[p] if i is None else values[p][i] for p,i in selectors))
        return call
    namespace['lambdify'] = lambdify
