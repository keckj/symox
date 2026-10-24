"""Ordered operator polynomials with native commuting scalar coefficients."""
import symox as S
from ..vector._basis import _BasisDependent


def _terms(value):return value._terms if isinstance(value,OperatorExpr) else {():S.sympify(value)}


def _word(word):
    coefficient = S.Integer(1);stack = []
    for atom in word:
        if getattr(atom,'is_identity',False):continue
        if stack:
            merged = _merge_powers(stack[-1],atom)
            if merged is not None:
                stack.pop();stack.extend(merged);continue
        if stack and hasattr(stack[-1],'_multiply_atom'):
            reduction = stack[-1]._multiply_atom(atom)
            if reduction is not None:
                stack.pop();factor,replacement = reduction;coefficient *= factor;stack.extend(replacement);continue
        stack.append(atom)
    return tuple(stack),coefficient


def _power_parts(atom):
    if isinstance(atom,PowerAtom) and isinstance(atom.base,Atom):return atom.base,atom.power
    return atom,S.Integer(1)


def _merge_powers(left,right):
    # SymPy's Mul combines adjacent noncommutative powers of one base, so
    # A*A**-1 is 1 and A**-1*A**-1 is A**-2.
    if not (isinstance(left,PowerAtom) or isinstance(right,PowerAtom)):return None
    (b1,e1),(b2,e2) = _power_parts(left),_power_parts(right)
    if b1!=b2:return None
    power = e1+e2
    if power==0:return ()
    if power.is_Integer and power>0:return (b1,)*int(power)
    return (PowerAtom(b1,power),)


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
        # Native form (matrix entries, scalar functions): a noncommutative
        # product of atom nodes, so native products keep the operator order
        # and print, conjugate, transpose and adjoint like the reference
        # noncommutative Mul of operators.
        terms = [coefficient*S.Mul(*(_atom_node(a) for a in word)) if word else coefficient for word,coefficient in self._terms.items()]
        return S.Add(*terms)
    def __repr__(self):
        # Print like the reference Add/Mul printer: each operator atom stands
        # in as a noncommutative placeholder symbol, so coefficients, signs,
        # denominators and term order follow it.  The placeholder name starts
        # with the atom class, then its text, so the noncommutative factor
        # keys order like the reference ones (class key, then arguments).
        if not any(self._terms):return str(self._terms.get((),0))
        names,texts = {},{}
        def placeholder(atom):
            if atom not in names:
                text = atom._atom_repr()
                name = '~~%s~%s~%d'%(type(atom).__name__,text,len(names))
                names[atom] = name;texts[name] = text
            return S.Symbol(names[atom],commutative=False)
        def factor(atom):
            if isinstance(atom,PowerAtom) and isinstance(atom.base,Atom):return placeholder(atom.base)**atom.power
            return placeholder(atom)
        try:
            native = S.Add(*(c*S.Mul(*(factor(a) for a in w)) if w else c for w,c in self._terms.items()))
            text = str(native)
        except Exception:
            words = {w:_word_repr(w) for w in self._terms if w}
            return ' + '.join((words[w] if c==1 else str(c)+'*'+words[w]) if w else str(c) for w,c in self._terms.items())
        for name in sorted(texts,key=len,reverse=True):
            atom = texts[name]
            simple = atom.replace('_','').isalnum() or (atom.endswith(')') and atom.count('(')==1 and atom[0].isalpha()) or (atom.startswith('[') and atom.endswith(']')) or (atom.startswith('{') and atom.endswith('}'))
            text = text.replace(name+'**',(atom if simple else '(%s)'%atom)+'**').replace(name,atom)
        return text
    __str__ = __repr__
    def __eq__(self,other):
        if isinstance(other,S.Expr) and not isinstance(other,OperatorExpr) and other.has(_AtomNode):
            # A native form taken out of a matrix or a scalar function.
            return self._symox_()==other
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
    def _atom_repr(self):
        # Reference StrPrinter._print_Pow for a noncommutative base:
        # parenthesize by precedence, so A**(-1), A**(1/2), (A + B)**n.
        base = self.base._atom_repr() if isinstance(self.base,Atom) else str(self.base)
        simple = base.replace('_','').isalnum() or (base.endswith(')') and base.count('(')==1 and base[0].isalpha())
        if not simple:base = '(%s)'%base
        power = self.power
        if power==S.Rational(1,2):return 'sqrt(%s)'%(self.base._atom_repr() if isinstance(self.base,Atom) else self.base)
        if power==S.Rational(-1,2):return '1/sqrt(%s)'%(self.base._atom_repr() if isinstance(self.base,Atom) else self.base)
        text = str(power)
        if not (power.is_Atom and not power.is_negative and not power.is_Rational or power.is_Integer and power>=0):text = '(%s)'%text
        return '%s**%s'%(base,text)
    def _subs_parameters(self,mapping,**kwargs):return PowerAtom(self.base.subs(mapping,**kwargs),self.power.subs(mapping,**kwargs))


# Python atoms of the native atom nodes.
_NODE_ATOMS = {}


def _atom_node(atom):
    if isinstance(atom,PowerAtom):
        base = atom.base
        base = _atom_node(base) if isinstance(base,Atom) else base._symox_() if isinstance(base,OperatorExpr) else S.sympify(base)
        return base**atom.power
    native = atom._native_atom()
    try:_NODE_ATOMS.setdefault(native,atom)
    except TypeError:pass
    return _AtomNode(native)


class _AtomNode(S.Function):
    """Native stand-in of an operator atom (a noncommutative symbol-like node)."""
    is_commutative = False
    _symox_user_hooks = True
    @property
    def atom(self):
        try:return _NODE_ATOMS.get(self.args[0])
        except TypeError:return None
    def _sympystr(self,printer,*args):
        atom = self.atom
        return atom._atom_repr() if atom is not None else printer._print(self.args[0])
    _sympyrepr = _sympystr
    def _symox_adjoint(self):
        # Reference Operator/State adjoint: Dagger(a), Dagger(|k>) = <k|.
        atom = self.atom
        adjoint = getattr(atom,'adjoint',None) if atom is not None else None
        if adjoint is None:
            from .quantum._core import Dagger
            if atom is None:return None
            value = Dagger(atom)
        else:value = adjoint()
        return value._symox_() if isinstance(value,OperatorExpr) else S.sympify(value)
    def _symox_transpose(self):
        # Expr._eval_transpose of a noncommutative operator.
        atom = self.atom
        if getattr(atom,'is_hermitian',None):return S.conjugate(self)
        return None
