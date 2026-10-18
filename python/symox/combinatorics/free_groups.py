"""Reduced words in finitely generated free groups."""
import operator


def free_group(symbols):
    group=FreeGroup(symbols)
    return (group,)+group.generators


class FreeGroup:
    def __init__(self,symbols):
        from ..polys.rings import _symbols
        self.symbols=_symbols(symbols)
        if len(set(self.symbols))!=len(self.symbols):raise ValueError('duplicate free generators')
        self.rank=len(self.symbols)
        self.identity=FreeGroupElement(self,())
        self.generators=tuple(FreeGroupElement(self,((s,1),)) for s in self.symbols)

    def __eq__(self,other):return isinstance(other,FreeGroup) and self.symbols==other.symbols
    def __hash__(self):return hash(self.symbols)
    def __contains__(self,word):return isinstance(word,FreeGroupElement) and word.group==self
    def __repr__(self):return '<free group on %s>'%', '.join(map(str,self.symbols))
    def order(self):
        from .. import oo
        return oo if self.rank else 1

    def dtype(self,syllables):
        return FreeGroupElement(self,syllables)


class FreeGroupElement:
    def __init__(self,group,syllables):
        self.group=group;stack=[]
        for s,e in syllables:
            if s not in group.symbols:raise ValueError('unknown free generator')
            e=operator.index(e)
            if not e:continue
            if stack and stack[-1][0]==s:
                e+=stack.pop()[1]
            if e:stack.append((s,e))
        self.array_form=tuple(stack)

    def __eq__(self,other):return isinstance(other,FreeGroupElement) and self.group==other.group and self.array_form==other.array_form
    def __hash__(self):return hash((self.group,self.array_form))
    def __bool__(self):return bool(self.array_form)
    @property
    def is_identity(self):return not self.array_form
    @property
    def letter_form(self):
        return tuple(s if e>0 else -s for s,e in self.array_form for _ in range(abs(e)))
    def __len__(self):return sum(abs(e) for s,e in self.array_form)
    def number_syllables(self):return len(self.array_form)
    def syllable_exponent(self,i):return self.array_form[i][1]
    def syllable_generator(self,i):return self.array_form[i][0]
    def exponent_sum(self,gen):
        if len(gen.array_form)!=1 or gen.array_form[0][1]!=1:raise ValueError('expected a generator')
        return sum(e for s,e in self.array_form if s==gen.array_form[0][0])
    def __mul__(self,other):
        if not isinstance(other,FreeGroupElement) or self.group!=other.group:return NotImplemented
        return FreeGroupElement(self.group,self.array_form+other.array_form)
    def inverse(self):return FreeGroupElement(self.group,((s,-e) for s,e in reversed(self.array_form)))
    __invert__=inverse
    def __truediv__(self,other):return self*other.inverse()
    def __pow__(self,n):
        if isinstance(n,FreeGroupElement):return n.inverse()*self*n
        n=operator.index(n);base=self if n>=0 else self.inverse();n=abs(n);result=self.group.identity
        if len(self)*n>100000:raise NotImplementedError('free word expansion exceeds 100000 letters')
        while n:
            if n&1:result=result*base
            base=base*base;n//=2
        return result
    def commutator(self,other):return self.inverse()*other.inverse()*self*other
    def subword(self,from_i,to_j,strict=True):
        if strict and not 0<=from_i<=to_j<=len(self):raise ValueError('invalid subword bounds')
        letters=[]
        for s,e in self.array_form:letters.extend([(s,1 if e>0 else -1)]*abs(e))
        return FreeGroupElement(self.group,letters[from_i:to_j])
    def __repr__(self):
        return '*'.join(str(s) if e==1 else '%s**%d'%(s,e) for s,e in self.array_form) or '<identity>'
    def to_sympy(self):
        from sympy.combinatorics import free_group
        group,*gens=free_group(','.join(map(str,self.group.symbols)))
        value=group.identity
        for s,e in self.array_form:value*=gens[self.group.symbols.index(s)]**e
        return value
