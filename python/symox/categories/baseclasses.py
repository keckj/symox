"""Structural categories and finite diagrams, without a symbolic dependency."""
from types import MappingProxyType


class _Value:
    def _set_args(self,*args):
        object.__setattr__(self,'args',tuple(args))

    def __setattr__(self,name,value):
        raise AttributeError('category values are immutable')

    def __hash__(self):
        return hash((type(self),self.args))

    def __eq__(self,other):
        return type(self) is type(other) and self.args == other.args

    def _plain_repr(self):
        return '%s(%s)' % (type(self).__name__,', '.join(map(repr,self.args)))

    def __repr__(self):
        # The reference prints category values with its StrPrinter, which
        # quotes names with double quotes: Object("A").
        from ..printing.str import sstr,StrPrinter
        if hasattr(type(self),'_sympystr') or any(hasattr(StrPrinter,'_print_'+c.__name__) for c in type(self).__mro__):return sstr(self)
        return '%s(%s)' % (type(self).__name__,', '.join(sstr(a) for a in self.args))

    __str__ = __repr__


class Object(_Value):
    def __init__(self,name):
        if not isinstance(name,str) or not name:
            raise ValueError('an object requires a nonempty string name')
        self._set_args(name)

    @property
    def name(self):
        return self.args[0]

    def _latex(self,printer):
        from .. import Symbol
        return printer.doprint(Symbol(self.name))

    def _pretty(self,printer):
        from ..printing.pretty.pretty_symbology import pretty_symbol
        from ..printing.pretty.stringpict import prettyForm
        return prettyForm(pretty_symbol(self.name))


def _pretty_arrow(name,morphism):
    from ..printing.pretty.pretty_symbology import xsym
    from ..printing.pretty.stringpict import prettyForm
    return prettyForm(name + morphism.domain._pretty(None).s + xsym('-->') + morphism.codomain._pretty(None).s)


class Morphism(_Value):
    def __init__(self,domain,codomain):
        if not isinstance(domain,Object) or not isinstance(codomain,Object):
            raise TypeError('morphism endpoints must be category objects')
        self._set_args(domain,codomain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    def compose(self,other):
        return CompositeMorphism(other,self)

    __mul__ = compose

    def _latex(self,printer):
        return r'%s\rightarrow %s' % (printer._print(self.domain),printer._print(self.codomain))


class IdentityMorphism(Morphism):
    def __init__(self,domain):
        if not isinstance(domain,Object):
            raise TypeError('identity domain must be an Object')
        self._set_args(domain)

    @property
    def codomain(self):
        return self.domain

    def _latex(self,printer):
        from .. import Symbol
        return '%s:%s' % (printer._print(Symbol('id')),Morphism._latex(self,printer))

    def _pretty(self,printer):
        return _pretty_arrow('id:',self)


class NamedMorphism(Morphism):
    def __init__(self,domain,codomain,name):
        super().__init__(domain,codomain)
        if not isinstance(name,str) or not name:
            raise ValueError('a named morphism requires a nonempty name')
        self._set_args(domain,codomain,name)

    @property
    def name(self):
        return self.args[2]

    def _latex(self,printer):
        from .. import Symbol
        return '%s:%s' % (printer._print(Symbol(self.name)),Morphism._latex(self,printer))

    def _pretty(self,printer):
        from ..printing.pretty.pretty_symbology import pretty_symbol
        return _pretty_arrow(pretty_symbol(self.name)+':',self)


class CompositeMorphism(Morphism):
    def __new__(cls,*components):
        if len(components) == 1 and not isinstance(components[0],Morphism):
            components = tuple(components[0])
        if not components or any(not isinstance(c,Morphism) for c in components):
            raise TypeError('composition requires at least one morphism')
        if any(a.codomain != b.domain for a,b in zip(components,components[1:])):
            raise ValueError('uncomposable morphisms')
        normalized = []
        for morphism in components:
            if isinstance(morphism,CompositeMorphism):
                normalized.extend(morphism.components)
            elif not isinstance(morphism,IdentityMorphism):
                normalized.append(morphism)
        if not normalized:
            return components[0]
        if len(normalized) == 1:
            return normalized[0]
        obj = super().__new__(cls)
        obj._set_args(tuple(normalized))
        return obj

    def __init__(self,*components):
        pass

    @property
    def components(self):
        return self.args[0]

    @property
    def domain(self):
        return self.components[0].domain

    @property
    def codomain(self):
        return self.components[-1].codomain

    def flatten(self,new_name):
        return NamedMorphism(self.domain,self.codomain,new_name)

    def _latex(self,printer):
        from .. import Symbol
        names = [printer._print(Symbol(c.name)) for c in reversed(self.components)]
        return r'\circ '.join(names) + ':' + Morphism._latex(self,printer)

    def _pretty(self,printer):
        from ..printing.pretty.pretty_symbology import pretty_symbol,xsym
        names = [pretty_symbol(c.name) for c in reversed(self.components)]
        return _pretty_arrow(xsym('.').join(names)+':',self)

    def _sympystr(self,printer):
        return 'CompositeMorphism((%s))' % ', '.join(printer._print(c) for c in self.components)


def _properties(value):
    if value is None:
        return frozenset()
    if isinstance(value,str):
        return frozenset((value,))
    try:
        return frozenset(value)
    except TypeError:
        return frozenset((value,))


def _entries(values):
    return list(values.items()) if hasattr(values,'items') else [(v,()) for v in values]


def _closure(values,identities=True,components=True):
    result = {}
    def add(morphism,properties):
        if not isinstance(morphism,Morphism):
            raise TypeError('diagram entries must be morphisms')
        properties = _properties(properties)
        if isinstance(morphism,IdentityMorphism) and properties:
            raise ValueError('identity morphisms cannot carry additional properties')
        result[morphism] = result.get(morphism,frozenset()) | properties
    for morphism,properties in _entries(values):
        add(morphism,properties)
        if components and isinstance(morphism,CompositeMorphism):
            for part in morphism.components:
                add(part,())
    if identities:
        for morphism in list(result):
            add(IdentityMorphism(morphism.domain),())
            add(IdentityMorphism(morphism.codomain),())
    # All simple composable paths are finite, including in diagrams with loops.
    # Explicitly supplied cyclic words are retained, but their unbounded powers
    # are not enumerated as part of a finite diagram.
    edges = [(m,p) for m,p in result.items() if not isinstance(m,(IdentityMorphism,CompositeMorphism))]
    outgoing = {}
    for morphism,properties in edges:
        outgoing.setdefault(morphism.domain,[]).append((morphism,properties))
    def visit(path,properties,seen):
        last = path[-1]
        for morphism,props in outgoing.get(last.codomain,()):
            if morphism.codomain in seen:
                continue
            word = CompositeMorphism(*path,morphism)
            add(word,properties & props)
            if len(result) > 10000:
                raise ValueError('diagram exceeds 10000 finite path entries')
            visit(path+[morphism],properties & props,seen | {morphism.codomain})
    for morphism,properties in edges:
        visit([morphism],properties,{morphism.domain,morphism.codomain})
    return result


class Diagram(_Value):
    __repr__ = __str__ = _Value._plain_repr

    def __init__(self,premises=(),conclusions=()):
        premises = _closure(premises)
        objects = frozenset(o for m in premises for o in (m.domain,m.codomain))
        conclusions = {m:p for m,p in _closure(conclusions,False,False).items()
                       if m.domain in objects and m.codomain in objects}
        self._set_args(frozenset(premises.items()),frozenset(conclusions.items()),objects)

    @property
    def premises(self):
        return MappingProxyType(dict(self.args[0]))

    @property
    def conclusions(self):
        return MappingProxyType(dict(self.args[1]))

    @property
    def objects(self):
        return self.args[2]

    def hom(self,A,B):
        return tuple(frozenset(m for m in mapping if m.domain == A and m.codomain == B)
                     for mapping in (self.premises,self.conclusions))

    def is_subdiagram(self,diagram):
        if not isinstance(diagram,Diagram):
            raise TypeError('expected a Diagram')
        return all(m in own and own[m] == props for own,sub in
                   ((self.premises,diagram.premises),(self.conclusions,diagram.conclusions)) for m,props in sub.items())

    def subdiagram_from_objects(self,objects):
        objects = frozenset(objects)
        if not objects <= self.objects:
            raise ValueError('objects must belong to the original diagram')
        mappings = [frozenset((m,p) for m,p in mapping.items() if m.domain in objects and m.codomain in objects)
                    for mapping in (self.premises,self.conclusions)]
        result = object.__new__(Diagram)
        result._set_args(*mappings,objects)
        return result


class Category(_Value):
    def _pretty(self,printer):
        return Object._pretty(self,printer)

    def _latex(self,printer):
        from .. import Symbol
        return r'\mathbf{%s}' % printer._print(Symbol(self.name))

    def __init__(self,name,objects=(),commutative_diagrams=()):
        if not isinstance(name,str) or not name:
            raise ValueError('a category requires a nonempty string name')
        objects,diagrams = frozenset(objects),frozenset(commutative_diagrams)
        if any(not isinstance(o,Object) for o in objects) or any(not isinstance(d,Diagram) for d in diagrams):
            raise TypeError('expected category objects and diagrams')
        self._set_args(name,objects,diagrams)

    name = property(lambda self:self.args[0])
    objects = property(lambda self:self.args[1])
    commutative_diagrams = property(lambda self:self.args[2])

    def hom(self,A,B):
        raise NotImplementedError('an abstract category does not specify its hom sets')

    def all_morphisms(self):
        raise NotImplementedError('an abstract category does not enumerate its morphisms')
