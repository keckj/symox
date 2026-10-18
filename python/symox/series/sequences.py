"""Formula and periodic sequences using native scalar mathematics."""
import itertools
import math
import operator

from .. import _core

MAX_TERMS=100_000
MAX_PERIOD=4096


def _scalar(value):
    from .. import sympify
    value=sympify(value)
    if not isinstance(value,_core.Expr):
        raise TypeError('Sequence entries must be scalar expressions')
    return value


def _items(value):
    if isinstance(value,_core.Expr) and getattr(value.func,'__name__','')=='Tuple':
        return value.args
    if isinstance(value,(tuple,list)):
        return tuple(value)
    raise ValueError('Expected a tuple or list')


def _limits(values, limits, periodic=False):
    from .. import Dummy, oo
    free=set().union(*(value.free_symbols for value in values))
    def variable():
        if len(free)==1:
            return next(iter(free))
        if free and not periodic:
            raise ValueError('Supply an index variable for a formula with multiple free symbols')
        return Dummy('k')
    if limits is None:
        limits=(variable(),0,oo)
    else:
        limits=_items(limits)
        if len(limits)==2:
            limits=(variable(),*limits)
        if len(limits)!=3:
            raise ValueError('Sequence limits require (index, start, stop) or (start, stop)')
    index,start,stop=map(_scalar,limits)
    if not index.is_Symbol:
        raise ValueError('Sequence index must be a symbol')
    if start==-oo and stop==oo:
        raise ValueError('Both sequence endpoints cannot be unbounded')
    if start.is_integer is False and start not in (-oo,oo) or stop.is_integer is False and stop not in (-oo,oo):
        raise ValueError('Sequence endpoints must be integers or symbolic bounds')
    return index,start,stop


def _empty_limits(limits):
    from .. import oo
    _,start,stop=limits
    return start==oo or stop==-oo or (stop-start).is_negative is True


def _compare_bound(condition):
    try:
        return bool(condition)
    except TypeError:
        raise NotImplementedError('Sequence bounds cannot be decided for this coefficient') from None


class SeqBase:
    __slots__=('_args','_limits')
    is_Sequence=True
    is_Atom=is_Symbol=is_Number=is_Add=is_Mul=is_Pow=is_Function=False
    is_commutative=None

    @property
    def args(self):return self._args
    @property
    def func(self):return type(self)
    @property
    def variables(self):return (self._limits[0],)
    @property
    def gen(self):return self.args[0]
    @property
    def start(self):return self._limits[1]
    @property
    def stop(self):return self._limits[2]
    @property
    def length(self):return self.stop-self.start+1
    @property
    def interval(self):
        from .. import Interval
        return Interval(self.start,self.stop)
    @property
    def free_symbols(self):
        return (set().union(*(value.free_symbols for value in self._expressions()))
                -set(self.variables))|self.start.free_symbols|self.stop.free_symbols

    def _ith_point(self, index):
        from .. import oo
        return self.stop-index if self.start==-oo else self.start+index

    def coeff(self, pt):
        pt=_scalar(pt)
        if _compare_bound(pt<self.start) or _compare_bound(pt>self.stop):
            raise IndexError('Coefficient index is outside the sequence interval')
        return self._eval_coeff(pt)

    def __getitem__(self, index):
        from .. import oo
        if not isinstance(index,slice):
            index=operator.index(index)
            if index<0:
                raise IndexError('Sequence positions must be nonnegative')
            return self.coeff(self._ith_point(index))
        start=0 if index.start is None else operator.index(index.start)
        if index.stop is None:
            if self.length==oo or self.length.is_integer is not True:
                raise ValueError('An unbounded or symbolic sequence slice needs a finite stop')
            stop=int(self.length)
        else:
            stop=operator.index(index.stop)
        step=1 if index.step is None else operator.index(index.step)
        if not step:
            raise ValueError('Slice step cannot be zero')
        distance=stop-start if step>0 else start-stop
        count=max(0,(distance+abs(step)-1)//abs(step))
        if count>MAX_TERMS:
            raise ValueError('Sequence slices are limited to 100000 terms')
        if any(value<0 for value in (start,stop)):
            raise IndexError('Sequence positions must be nonnegative')
        return [self[position] for position in range(start,stop,step)]

    def __iter__(self):
        from .. import oo
        if self.length==oo:
            for index in itertools.count():
                yield self[index]
        elif self.length.is_integer is True:
            for index in range(int(self.length)):
                yield self[index]
        else:
            raise NotImplementedError('Iteration requires a concrete sequence length')

    def __add__(self, other):
        return SeqAdd(self,other)
    def __sub__(self, other):
        if not isinstance(other,SeqBase):
            return NotImplemented
        return SeqAdd(self,other.coeff_mul(-1))
    def __mul__(self, other):
        return SeqMul(self,other)
    def __rmul__(self, other):
        return SeqMul(other,self)
    def __neg__(self):return self.coeff_mul(-1)

    def __eq__(self, other):
        return type(self) is type(other) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))
    def __repr__(self):return '%s(%s)' % (type(self).__name__,', '.join(map(str,self.args)))
    __str__=__repr__
    def __reduce__(self):return type(self),self.args
    def copy(self):return self
    def doit(self, **options):return self

    def sort_key(self, order=None):
        from .. import S, Tuple, default_sort_key
        def key(value):
            return Tuple(*value).sort_key(order=order) if isinstance(value,tuple) else default_sort_key(value,order)
        name='EmptySequence' if self is EmptySequence else type(self).__name__
        return (5,0,name),(len(self.args),tuple(key(arg) for arg in self.args)),S.One.sort_key(),S.One

    def xreplace(self, mapping):
        if self in mapping:
            return mapping[self]
        return self._transform(lambda value:value.xreplace(mapping),replace_index=True,mapping=mapping)

    def subs(self, *args, **options):
        # The index is bound in the generating expression. Its appearances
        # in interval endpoints remain free and can still be substituted.
        from collections.abc import Mapping
        from .. import Dummy
        index=self.variables[0]
        pairs=(args[0].items() if len(args)==1 and isinstance(args[0],Mapping) else
               args[0] if len(args)==1 and isinstance(args[0],(list,tuple)) else
               (args,) if len(args)==2 else ())
        replacement_symbols=set()
        for old,new in pairs:
            if old==index:
                continue
            if not hasattr(new,'free_symbols'):
                try:new=_scalar(new)
                except TypeError:continue
            replacement_symbols.update(new.free_symbols)
        result_index=Dummy(index.name,integer=True) if index in replacement_symbols else index
        def substitute(value):
            token=Dummy('sequence_index')
            return value.xreplace({index:token}).subs(*args,**options).xreplace({token:result_index})
        return self._transform(substitute,boundary=lambda value:value.subs(*args,**options),
                               replace_index=True,mapping={index:result_index})

    def simplify(self, **options):
        from .. import simplify
        return self._transform(lambda value:simplify(value,**options))
    def expand(self, **options):
        from .. import expand
        return self._transform(lambda value:expand(value,**options))

    def find_linear_recurrence(self, n, d=None, gfvar=None):
        from .. import Matrix, S, cancel, simplify
        n=operator.index(n)
        if not 0<=n<=256:
            raise ValueError('Recurrence searches require between 0 and 256 terms')
        degree=min(n//2,64 if d is None else operator.index(d))
        if degree<0 or degree>64:
            raise ValueError('Recurrence search order must lie between 0 and 64')
        values=[cancel(simplify(value)) for value in self[:n]]
        coefficients=[]
        if any(value!=0 for value in values):
            for order in range(1,degree+1):
                matrix=Matrix([[values[k-j-1] for j in range(order)] for k in range(order,n)])
                rhs=Matrix(values[order:])
                try:
                    solution,parameters=matrix.gauss_jordan_solve(rhs)
                except ValueError:
                    continue
                replacements={parameter:S.Zero for parameter in parameters}
                candidate=[cancel(value.subs(replacements)) for value in solution]
                if all(cancel(values[k]-sum(candidate[j]*values[k-j-1] for j in range(order)))==0
                       for k in range(order,n)):
                    coefficients=candidate
                    break
        if gfvar is None:
            return coefficients
        if not coefficients:
            return [],None
        variable=_scalar(gfvar)
        numerator=sum((values[k]-sum(coefficients[j-1]*values[k-j] for j in range(1,k+1)))
                      *variable**k for k in range(len(coefficients)))
        denominator=1-sum(value*variable**(j+1) for j,value in enumerate(coefficients))
        return coefficients,cancel(numerator/denominator)


class SeqFormula(SeqBase):
    __slots__=()
    def __new__(cls, formula, limits=None):
        formula=_scalar(formula)
        limits=_limits((formula,),limits)
        if _empty_limits(limits):return EmptySequence
        self=object.__new__(cls)
        self._limits=limits
        self._args=(formula,limits)
        return self
    @property
    def formula(self):return self.args[0]
    def _expressions(self):return (self.formula,)
    def _eval_coeff(self, pt):return self.formula.subs(self.variables[0],pt)
    def coeff_mul(self, coeff):return SeqFormula(self.formula*_scalar(coeff),self._limits)
    def _transform(self, function, boundary=None, replace_index=False, mapping=None):
        boundary=function if boundary is None else boundary
        index=self.variables[0]
        if replace_index:index=mapping.get(index,index)
        return SeqFormula(function(self.formula),(index,boundary(self.start),boundary(self.stop)))
    def to_sympy(self):
        import sympy
        return sympy.SeqFormula(self.formula.to_sympy(),tuple(v.to_sympy() for v in self._limits))


class SeqPer(SeqBase):
    __slots__=()
    def __new__(cls, periodical, limits=None):
        values=[]
        active=set()
        def collect(items,depth=0):
            if depth>32 or id(items) in active:
                raise ValueError('Periodic sequence input is cyclic or nested deeper than 32')
            active.add(id(items))
            try:
                for value in items:
                    if isinstance(value,(list,tuple)) or (isinstance(value,_core.Expr)
                            and getattr(value.func,'__name__','')=='Tuple'):
                        collect(_items(value),depth+1)
                    else:
                        if len(values)>=MAX_PERIOD:
                            raise ValueError('Sequence periods are limited to 4096 entries')
                        values.append(_scalar(value))
            finally:
                active.remove(id(items))
        collect(_items(periodical))
        values=tuple(values)
        limits=_limits(values,limits,periodic=True)
        if not values or _empty_limits(limits):return EmptySequence
        self=object.__new__(cls)
        self._limits=limits
        self._args=(values,limits)
        return self
    @property
    def periodical(self):return self.args[0]
    @property
    def period(self):return len(self.periodical)
    def _expressions(self):return self.periodical
    def _phase(self, pt):
        from .. import oo
        difference=(self.stop-pt if self.start==-oo else pt-self.start)
        if difference.is_integer is not True or not difference.is_number:
            raise NotImplementedError('Periodic coefficients need a concrete integer offset')
        return int(difference)%self.period
    def _eval_coeff(self, pt):
        return self.periodical[self._phase(pt)].subs(self.variables[0],pt)
    def coeff_mul(self, coeff):
        coeff=_scalar(coeff)
        return SeqPer(tuple(value*coeff for value in self.periodical),self._limits)
    def _transform(self, function, boundary=None, replace_index=False, mapping=None):
        boundary=function if boundary is None else boundary
        index=self.variables[0]
        if replace_index:index=mapping.get(index,index)
        return SeqPer(tuple(function(value) for value in self.periodical),
                      (index,boundary(self.start),boundary(self.stop)))
    def to_sympy(self):
        import sympy
        return sympy.SeqPer(tuple(v.to_sympy() for v in self.periodical),
                            tuple(v.to_sympy() for v in self._limits))


class _SeqOperation(SeqBase):
    __slots__=()
    def __new__(cls, *arguments, evaluate=True):
        from .. import Max, Min, Dummy, S
        if not all(isinstance(arg,SeqBase) for arg in arguments):
            raise TypeError('Sequence arithmetic requires sequence operands; use coeff_mul for scalars')
        additive=cls is SeqAdd
        if not additive and any(arg is EmptySequence for arg in arguments):return EmptySequence
        arguments=tuple(arg for arg in arguments if arg is not EmptySequence)
        if not arguments:return EmptySequence
        if len(arguments)==1:return arguments[0]
        if evaluate:
            arguments=tuple(child for arg in arguments for child in
                            (arg.args if isinstance(arg,cls) else (arg,)))
        free=set().union(*(arg.free_symbols for arg in arguments))
        index=arguments[0].variables[0]
        if index in free:index=Dummy('n',integer=True)
        limits=(index,Max(*(arg.start for arg in arguments)),Min(*(arg.stop for arg in arguments)))
        if _empty_limits(limits):return EmptySequence
        def combine(values):
            if additive:return sum(values,S.Zero)
            result=S.One
            for value in values:result*=value
            return result
        if evaluate and all(isinstance(arg,SeqFormula) for arg in arguments):
            return SeqFormula(combine(arg.formula.subs(arg.variables[0],index) for arg in arguments),limits)
        if evaluate and all(isinstance(arg,SeqPer) for arg in arguments):
            period=math.lcm(*(arg.period for arg in arguments))
            if period<=MAX_PERIOD:
                from .. import oo
                anchor=limits[2] if limits[1]==-oo else limits[1]
                direction=-1 if limits[1]==-oo else 1
                try:
                    patterns=[combine(arg.periodical[arg._phase(anchor+direction*i)]
                                      .subs(arg.variables[0],index) for arg in arguments)
                              for i in range(period)]
                except NotImplementedError:
                    pass
                else:
                    return SeqPer(patterns,limits)
        self=object.__new__(cls)
        self._limits=limits
        self._args=tuple(sorted(arguments,key=lambda arg:str(arg)))
        return self
    @property
    def free_symbols(self):return set().union(*(arg.free_symbols for arg in self.args))
    @property
    def gen(self):return tuple(arg.gen for arg in self.args)
    def _eval_coeff(self, pt):
        from .. import S
        values=[arg.coeff(pt) for arg in self.args]
        if isinstance(self,SeqAdd):return sum(values,S.Zero)
        result=S.One
        for value in values:result*=value
        return result
    def coeff_mul(self, coeff):
        if isinstance(self,SeqAdd):return SeqAdd(*(arg.coeff_mul(coeff) for arg in self.args))
        return SeqMul(self.args[0].coeff_mul(coeff),*self.args[1:])
    def subs(self, *args, **options):return type(self)(*(arg.subs(*args,**options) for arg in self.args))
    def xreplace(self, mapping):
        return mapping[self] if self in mapping else type(self)(*(arg.xreplace(mapping) for arg in self.args))
    def _transform(self, function, **options):
        return type(self)(*(arg._transform(function,**options) for arg in self.args))
    def to_sympy(self):
        import sympy
        return getattr(sympy,type(self).__name__)(*(arg.to_sympy() for arg in self.args),evaluate=False)


class SeqAdd(_SeqOperation):__slots__=()
class SeqMul(_SeqOperation):__slots__=()


class _EmptySequence(SeqBase):
    __slots__=()
    @property
    def start(self):
        from .. import oo
        return oo
    @property
    def stop(self):
        from .. import oo
        return -oo
    @property
    def variables(self):return ()
    @property
    def free_symbols(self):return set()
    @property
    def length(self):
        from .. import S
        return S.Zero
    @property
    def interval(self):
        from .. import S
        return S.EmptySet
    def __iter__(self):return iter(())
    def __getitem__(self,index):
        if isinstance(index,slice):return []
        raise IndexError('EmptySequence has no coefficients')
    def coeff(self,pt):raise IndexError('EmptySequence has no coefficients')
    def coeff_mul(self,coeff):return self
    def subs(self,*args,**options):return self
    def xreplace(self,mapping):return mapping.get(self,self)
    def _transform(self,*args,**options):return self
    def __repr__(self):return 'EmptySequence'
    __str__=__repr__
    def __reduce__(self):return _empty_sequence,()
    def to_sympy(self):
        from sympy import EmptySequence
        return EmptySequence


def _empty_sequence():return EmptySequence

EmptySequence=object.__new__(_EmptySequence)
EmptySequence._args=()


def sequence(seq, limits=None):
    if isinstance(seq,SeqBase):
        if limits is not None:
            raise ValueError('An existing sequence already has limits')
        return seq
    if isinstance(seq,(list,tuple)) or (isinstance(seq,_core.Expr) and
            getattr(seq.func,'__name__','')=='Tuple'):
        return SeqPer(seq,limits)
    return SeqFormula(seq,limits)


__all__=['SeqFormula','SeqPer','SeqAdd','SeqMul','EmptySequence','sequence']
