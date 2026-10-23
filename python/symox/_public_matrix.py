"""Deferred vectors, indexed contractions and matrix scalar entry points."""
import operator

from . import _core
from ._compat import _ExprType

_indexed = {}
_indices = {}
_deferred = {}


class DeferredVector:
    # SymPy's ``DeferredVector(Symbol, NotIterable)``: list() would never end.
    _iterable = False

    def __iter__(self):
        raise TypeError("'DeferredVector' object is not iterable")

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
    def __new__(cls,label,range=None,**kwargs):
        # Idx.__new__ returns an integer label itself.
        from . import sympify
        if not isinstance(label,str):
            value = sympify(label)
            if value.is_Number:
                if not value.is_integer:
                    raise TypeError('Index is not an integer number.')
                return value
        return super().__new__(cls)

    def __init__(self,label,range=None,**kwargs):
        from . import Symbol,sympify,Integer
        if kwargs:
            raise TypeError('unexpected index options')
        self.label = Symbol(label,integer=True) if isinstance(label,str) else sympify(label)
        if not self.label.is_integer:
            raise TypeError('Idx object requires an integer label.')
        # Idx.__new__: integer (or infinite) bounds, a range ``n`` meaning
        # ``(0, n - 1)``; the reference does not order the bounds.
        from . import S
        if range is None:
            self.lower = self.upper = None
        elif isinstance(range,(tuple,list)) or getattr(range,'is_Tuple',False):
            values = range.args if getattr(range,'is_Tuple',False) else range
            if len(values) != 2:
                raise ValueError('\nIdx range tuple must have length 2, but got %s' % len(values))
            values = tuple(map(sympify,values))
            for bound in values:
                if bound.is_integer is False and bound is not S.Infinity and bound is not S.NegativeInfinity:
                    raise TypeError('Idx object requires integer bounds.')
            self.lower,self.upper = values
        else:
            range = sympify(range)
            if range is not S.Infinity and range.is_integer is False:
                raise TypeError('Idx object requires an integer dimension.')
            self.lower,self.upper = Integer(0),range-1
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

    # Idx is an Expr in the reference: arithmetic acts on its label, so
    # ``x[i + 1]`` (apply_finite_diff stencils) is the shifted index.
    def __add__(self,other): return self.label+_idx_operand(other)
    def __radd__(self,other): return _idx_operand(other)+self.label
    def __sub__(self,other): return self.label-_idx_operand(other)
    def __rsub__(self,other): return _idx_operand(other)-self.label
    def __mul__(self,other): return self.label*_idx_operand(other)
    def __rmul__(self,other): return _idx_operand(other)*self.label
    def __truediv__(self,other): return self.label/_idx_operand(other)
    def __rtruediv__(self,other): return _idx_operand(other)/self.label
    def __pow__(self,other): return self.label**_idx_operand(other)
    def __rpow__(self,other): return _idx_operand(other)**self.label
    def __neg__(self): return -self.label
    def __pos__(self): return self.label


def _idx_operand(other):
    from . import sympify
    return other.label if isinstance(other,Idx) else sympify(other)


class IndexConformanceException(ValueError):
    """Outer Einstein indices disagree or occur more than twice."""


def _as_indexed(expr):
    from .codegen.ast import Indexed
    if isinstance(expr,Indexed):return expr
    from .codegen.ast import _INDEXED_SYMBOLS
    try:return _indexed.get(expr) or _INDEXED_SYMBOLS.get(expr)
    except TypeError:return None


def _unique_and_repeated(inds):
    # SymPy's index_methods._unique_and_repeated: first-occurrence order.
    uniq = {}
    for i in inds:
        uniq[i] = 0 if i in uniq else 1
    return [i for i,n in uniq.items() if n],[i for i,n in uniq.items() if not n]


def _remove_repeated(inds):
    unique,repeated = _unique_and_repeated(inds)
    return set(unique),tuple(repeated)


def _index_object(i):
    try:return _indices.get(i,i)
    except TypeError:return i


def _indexed_indices(indexed):
    return [_index_object(i) for i in indexed.indices]


def _get_indices_Mul(expr,return_dummies=False):
    inds = [i for arg in expr.args for i in get_indices(arg)[0]]
    inds,dummies = _remove_repeated(inds)
    if return_dummies:return inds,{},dummies
    return inds,{}


def _index_atom(expr):
    return isinstance(expr,_core.Expr) and not expr.args


def get_indices(expr):
    """Outer indices of ``expr``, following sympy.tensor.index_methods."""
    from . import sympify
    indexed = _as_indexed(expr)
    if indexed is not None:
        unique,_ = _unique_and_repeated(_indexed_indices(indexed))
        return set(unique),{}
    if expr is None:return set(),{}
    if isinstance(expr,Idx):return {expr},{}
    expr = sympify(expr)
    if _index_atom(expr):return set(),{}
    if expr.is_Mul:
        return _get_indices_Mul(expr)
    if expr.is_Add:
        results = [get_indices(a)[0] for a in expr.args]
        nonscalar = [r for r in results if r != set()]
        if not nonscalar:return set(),{}
        if any(r != nonscalar[0] for r in nonscalar[1:]):
            raise IndexConformanceException('Indices are not consistent: %s' % expr)
        return nonscalar[0],{}
    if expr.is_Pow or getattr(expr,'class_name','') == 'exp':
        base,exponent = expr.as_base_exp()
        return get_indices(base)[0] | get_indices(exponent)[0],{}
    if getattr(expr,'class_name','') == 'Piecewise':
        return set(),{}
    if expr.is_Function:
        found = set()
        for arg in expr.args:
            found |= get_indices(arg)[0]
        return found,{}
    return set(),{}


def _trivial(structure):
    return None in structure and len(structure) == 1


def get_contraction_structure(expr):
    """Dummy-index structure of ``expr``, following sympy.tensor.index_methods."""
    from . import sympify
    indexed = _as_indexed(expr)
    if indexed is not None:
        _,key = _remove_repeated(_indexed_indices(indexed))
        return {key or None:{expr}}
    expr = sympify(expr)
    if _index_atom(expr):return {None:{expr}}
    if expr.is_Mul:
        _,_,key = _get_indices_Mul(expr,return_dummies=True)
        result = {key or None:{expr}}
        nested = []
        for factor in expr.args:
            structure = get_contraction_structure(factor)
            if not _trivial(structure):nested.append(structure)
        if nested:result[expr] = nested
        return result
    if expr.is_Pow or getattr(expr,'class_name','') == 'exp':
        base,exponent = expr.as_base_exp()
        dicts = [d for d in (get_contraction_structure(base),get_contraction_structure(exponent)) if not _trivial(d)]
        result = {None:{expr}}
        if dicts:result[expr] = dicts
        return result
    if expr.is_Add:
        result = {}
        for term in expr.args:
            for key,value in get_contraction_structure(term).items():
                if key in result:result[key] |= value
                else:result[key] = value
        return result
    if getattr(expr,'class_name','') == 'Piecewise':
        return {None:expr}
    if expr.is_Function:
        deep = [d for d in (get_contraction_structure(a) for a in expr.args) if not _trivial(d)]
        result = {None:{expr}}
        if deep:result[expr] = deep
        return result
    return {None:{expr}}


class Permanent(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr.class_name == 'Permanent')

    def __new__(cls,mat):
        from . import MatrixExpr
        return _core._apply('Permanent',MatrixExpr(mat))


def per(matexpr):
    from . import DenseMatrix as Matrix
    if isinstance(matexpr,Matrix) or getattr(matexpr,'_symox_sparse',False):
        return matexpr.per()
    return Permanent(matexpr)


def diagonalize_vector(vector):
    from . import DenseMatrix as Matrix,DiagMatrix,diag
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
            from .codegen.ast import _INDEXED_SYMBOLS
            _INDEXED_SYMBOLS.setdefault(value,self)
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
