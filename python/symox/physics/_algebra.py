"""Ordered operator polynomials with native commuting scalar coefficients."""
import symox as S
from ..vector._basis import _BasisDependent


def _terms(value):return value._terms if isinstance(value,OperatorExpr) else {():S.sympify(value)}


def _word(word):
    coefficient = S.Integer(1);stack = []
    for atom in word:
        if getattr(atom,'is_identity',False):continue
        if stack and hasattr(stack[-1],'_multiply_atom'):
            reduction = stack[-1]._multiply_atom(atom)
            if reduction is not None:
                stack.pop();factor,replacement = reduction;coefficient *= factor;stack.extend(replacement);continue
        stack.append(atom)
    return tuple(stack),coefficient


def _word_repr(word):
    """Noncommutative product text; equal neighbours print as powers like SymPy's Mul."""
    parts = [];index = 0
    while index<len(word):
        run = index
        while run+1<len(word) and word[run+1]==word[index]:run += 1
        text = word[index]._atom_repr();count = run-index+1
        if count>1:
            simple = text.replace('_','').isalnum() or (text.endswith(')') and text.count('(')==1 and text[0].isalpha())
            text = '%s**%d'%(text if simple else '(%s)'%text,count)
        parts.append(text);index = run+1
    return '*'.join(parts)


class OperatorExpr(_BasisDependent):
    is_commutative = False
    def __init__(self,terms=None):
        self._terms = {word:S.sympify(value) for word,value in (terms or {}).items() if value!=0}
    @staticmethod
    def _make(terms):
        terms = {word:S.simplify(value) for word,value in terms.items() if value!=0}
        terms = {word:value for word,value in terms.items() if value!=0}
        if not terms:return S.Integer(0)
        if set(terms)=={()}:return terms[()]
        if len(terms)==1:
            word,value = next(iter(terms.items()))
            if len(word)==1 and value==1:return word[0]
        return OperatorExpr(terms)
    @property
    def free_symbols(self):return set().union(*(c.free_symbols for c in self._terms.values()),*(a._parameter_symbols() for w in self._terms for a in w))
    @property
    def args(self):return tuple(self._terms.items())
    def _symox_(self):
        terms = [coefficient*S._core._apply('OrderedProduct',*(a._native_atom() for a in word),evaluate=False) if word else coefficient for word,coefficient in self._terms.items()]
        return S.Add(*terms)
    def __repr__(self):
        # Print like the reference Add/Mul printer: each operator word stands
        # in as a placeholder symbol sorting after commuting factors, so
        # coefficients, signs, denominators and term order follow it.
        words = {}
        for word in self._terms:
            if word:words[word] = _word_repr(word)
        if not words:return str(self._terms.get((),0))
        placeholders = {word:'~~%s~%d'%(text,i) for i,(word,text) in enumerate(words.items())}
        try:
            native = S.Add(*(c*S.Symbol(placeholders[w],commutative=False) if w else c for w,c in self._terms.items()))
            text = str(native)
        except Exception:
            return ' + '.join((words[w] if c==1 else str(c)+'*'+words[w]) if w else str(c) for w,c in self._terms.items())
        for word,name in sorted(placeholders.items(),key=lambda item:-len(item[1])):text = text.replace(name,words[word])
        return text
    __str__ = __repr__
    def __eq__(self,other):
        try:return self._terms==_terms(other)
        except (TypeError,ValueError):return False
    # Terms never change after construction, so operator polynomials can be
    # arguments of commutators and tensor products, like SymPy's Add and Mul.
    def __hash__(self):return hash((OperatorExpr,frozenset(self._terms.items())))
    def __bool__(self):return bool(self._terms)
    def __add__(self,other):
        terms = dict(self._terms)
        for word,value in _terms(other).items():terms[word] = terms.get(word,S.Integer(0))+value
        return self._make(terms)
    __radd__ = __add__
    def __neg__(self):return self*-1
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        terms = {}
        for left,a in self._terms.items():
            for right,b in _terms(other).items():
                word,factor = _word(left+right);terms[word] = terms.get(word,S.Integer(0))+a*b*factor
        return self._make(terms)
    def __rmul__(self,other):
        if isinstance(other,OperatorExpr):return other*self
        scalar = S.sympify(other);return self._make({word:scalar*c for word,c in self._terms.items()})
    def __truediv__(self,other):
        if isinstance(other,OperatorExpr):return self*other**-1
        return self*(1/S.sympify(other))
    def __pow__(self,power):
        power = S.sympify(power)
        if not power.is_Integer:return PowerAtom(self,power)
        if power<0:return PowerAtom(self,power)
        result = S.Integer(1)
        for _ in range(int(power)):result = self.__rmul__(result)
        return result
    def expand(self,**kwargs):
        # Only the quantum hints rewrite operator atoms; the algebra itself is
        # always kept expanded over sums.
        if not any(kwargs.get(hint) for hint in ('commutator','tensorproduct')):return self
        result = S.Integer(0)
        for word,coefficient in self._terms.items():
            term = coefficient
            for atom in word:term = term*atom.expand(**kwargs)
            result = result+term
        return result
    def simplify(self,**kwargs):return self._make({w:S.simplify(v,**kwargs) for w,v in self._terms.items()})
    def trigsimp(self,**kwargs):return self._make({w:S.trigsimp(v,**kwargs) for w,v in self._terms.items()})
    def evalf(self,*args,**kwargs):return self._make({w:v.evalf(*args,**kwargs) for w,v in self._terms.items()})
    def doit(self,**hints):
        # SymPy's Mul/Add.doit: evaluate every factor (commutators, inner products).
        if not any(type(atom).doit is not Atom.doit for word in self._terms for atom in word) and \
                not any(hasattr(c,'doit') and c.has(S.Function) for c in self._terms.values()):
            return self
        result = S.Integer(0)
        for word,coefficient in self._terms.items():
            term = coefficient.doit(**hints) if hasattr(coefficient,'doit') else coefficient
            for atom in word:term = term*atom.doit(**hints)
            result = result+term
        return result
    def _diff_once(self,var):
        # Product rule: coefficients and atoms with their own derivative
        # (such as wavefunctions); other operators are constant.
        result=S.Integer(0)
        for word,coefficient in self._terms.items():
            term=S.diff(coefficient,var)
            for atom in word:term=term*atom
            result+=term
            for i,atom in enumerate(word):
                derivative=getattr(atom,'_atom_diff',None)
                if derivative is None:continue
                term=coefficient
                for j,other in enumerate(word):term=term*(derivative(var) if i==j else other)
                result+=term
        return result
    def diff(self,*symbols,**kwargs):
        if not symbols:
            free=self.free_symbols
            if len(free)!=1:raise ValueError('specify differentiation variables to differentiate %s'%self)
            symbols=tuple(free)
        pending=list(symbols);steps=[]
        while pending:
            var=pending.pop(0)
            if isinstance(var,(tuple,list)):var,count=var
            elif pending and isinstance(pending[0],int) or (pending and isinstance(pending[0],S.Expr) and pending[0].is_Integer):count=pending.pop(0)
            else:count=1
            steps.extend([var]*int(count))
        result=self
        for var in steps:result=result._diff_once(var) if isinstance(result,OperatorExpr) else S.diff(result,var)
        return result
    def _eval_symox_diff(self,*symbols):return self.diff(*symbols)
    def integrate(self,*args,**kwargs):return self._make({w:S.integrate(v,*args,**kwargs) for w,v in self._terms.items()})
    def subs(self,*args,**kwargs):
        mapping = args[0] if len(args)==1 else {args[0]:args[1]}
        scalar_mapping = {k:v for k,v in mapping.items() if not isinstance(k,OperatorExpr)}
        result = S.Integer(0)
        for word,value in self._terms.items():
            term = value.subs(scalar_mapping,**kwargs)
            for atom in word:term = term*mapping.get(atom,atom._subs_parameters(scalar_mapping,**kwargs))
            result += term
        return result
    def xreplace(self,rule):return self.subs(rule,simultaneous=True)
    def has(self,*patterns):return any(pattern==a for w in self._terms for a in w for pattern in patterns)


class Atom(OperatorExpr):
    def __init__(self,*labels):
        self.labels = tuple(labels);self._terms = {(self,):S.Integer(1)}
    def _key(self):return type(self),self.labels
    def __hash__(self):return hash(self._key())
    def __eq__(self,other):return isinstance(other,Atom) and self._key()==other._key()
    def _parameter_symbols(self):return set().union(*(v.free_symbols for v in self.labels if isinstance(v,S.Expr)))
    def expand(self,**kwargs):return self
    def doit(self,**hints):return self
    @property
    def args(self):return self.labels
    def _atom_repr(self):return '%s(%s)'%(type(self).__name__,', '.join(map(str,self.labels)))
    def __repr__(self):return self._atom_repr()
    __str__ = __repr__
    def _native_atom(self):
        return S._core._apply(type(self).__name__,*(v._symox_() if isinstance(v,OperatorExpr) else S.Symbol(v) if isinstance(v,str) else S.sympify(v) for v in self.labels),evaluate=False)
    def _subs_parameters(self,mapping,**kwargs):return type(self)(*(v.subs(mapping,**kwargs) if isinstance(v,S.Expr) else v for v in self.labels))


class PowerAtom(Atom):
    def __init__(self,base,power):self.base,self.power = base,power;super().__init__(base,power)
    def _atom_repr(self):return '(%s)**%s'%(self.base,self.power)
    def _subs_parameters(self,mapping,**kwargs):return PowerAtom(self.base.subs(mapping,**kwargs),self.power.subs(mapping,**kwargs))
