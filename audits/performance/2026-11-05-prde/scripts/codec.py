"""Encoding of prde call arguments, DifferentialExtension states and results
that replays identically in SymPy 1.14 and symox.

Symbols are encoded with srepr; Dummies are renamed by name only (the
dummy_index differs between processes), and every Dummy of one record with
the same name is the same object on replay. Domains are encoded
structurally (ZZ, QQ, ZZ[gens], QQ[gens], ZZ(gens), QQ(gens), EX, other).
"""
import re

_DUMMY = re.compile(r"Dummy\('([^']*)'(?:, dummy_index=\d+)?\)")


def mod():
    import importlib
    import os
    return importlib.import_module(os.environ.get('PRDE_BACKEND_MODULE', 'symox'))


def srepr(S, e):
    return _DUMMY.sub(lambda m: "Dummy('%s')" % m.group(1), S.srepr(e))


def dom(S, d):
    s = str(d)
    if s in ('ZZ', 'QQ', 'EX', 'EXRAW'):
        return [s]
    if getattr(d, 'is_PolynomialRing', False) or getattr(d, 'is_FractionField', False):
        kind = 'poly' if d.is_PolynomialRing else 'frac'
        return [kind, dom(S, d.dom), [srepr(S, g) for g in d.symbols]]
    return ['other', s]


def poly(S, p):
    return ['Poly', srepr(S, p.as_expr()), [srepr(S, g) for g in p.gens], dom(S, p.get_domain())]


def de(S, DE):
    return {'D': [poly(S, d) for d in DE.D], 'T': [srepr(S, t) for t in DE.T],
            'exts': list(DE.exts) if DE.exts is not None else None,
            'extargs': [srepr(S, a) if a is not None else None for a in DE.extargs] if DE.extargs is not None else None,
            'level': DE.level, 'x': srepr(S, DE.x)}


def value(S, v):
    """Result/argument encoding (recursive)."""
    from importlib import import_module
    name = S.__name__
    PM = import_module(name + '.polys.polymatrix').PolyMatrix
    if v is None or isinstance(v, (bool, str)):
        return v
    if isinstance(v, int):
        return ['int', str(v)]
    if isinstance(v, S.Poly):
        return poly(S, v)
    if isinstance(v, PM):
        return ['PM', v.rows, v.cols, [srepr(S, e) for e in v._items] if hasattr(v, '_items') else
                [srepr(S, e) for e in v.to_Matrix()], dom(S, v.domain), [srepr(S, g) for g in v.gens]]
    if isinstance(v, tuple):
        return ['tuple', [value(S, i) for i in v]]
    if isinstance(v, list):
        return ['list', [value(S, i) for i in v]]
    if isinstance(v, S.Basic):
        return ['expr', srepr(S, v)]
    return ['repr', type(v).__name__, repr(v)]


class Decoder:
    def __init__(self, S):
        self.S = S
        self.dummies = {}
        from importlib import import_module
        self.PM = import_module(S.__name__ + '.polys.polymatrix').PolyMatrix
        self.ns = {k: getattr(S, k) for k in dir(S) if not k.startswith('_')}
        self.ns['Dummy'] = self.dummy
        self.ns['Integer'] = S.Integer
        self.ns['Rational'] = S.Rational

    def dummy(self, name, **kw):
        if name not in self.dummies:
            self.dummies[name] = self.S.Dummy(name)
        return self.dummies[name]

    def expr(self, s):
        return eval(s, self.ns)

    def dom(self, d):
        S = self.S
        if d[0] in ('ZZ', 'QQ', 'EX', 'EXRAW'):
            return getattr(S, d[0])
        if d[0] in ('poly', 'frac'):
            base = self.dom(d[1])
            gens = [self.expr(g) for g in d[2]]
            return base.poly_ring(*gens) if d[0] == 'poly' else base.frac_field(*gens)
        raise ValueError('unsupported domain %r' % (d,))

    def poly(self, p):
        _, e, gens, d = p
        return self.S.Poly(self.expr(e), *[self.expr(g) for g in gens], domain=self.dom(d))

    def value(self, v):
        if v is None or isinstance(v, (bool, str)):
            return v
        tag = v[0]
        if tag == 'int':
            return int(v[1])
        if tag == 'Poly':
            return self.poly(v)
        if tag == 'tuple':
            return tuple(self.value(i) for i in v[1])
        if tag == 'list':
            return [self.value(i) for i in v[1]]
        if tag == 'expr':
            return self.expr(v[1])
        raise ValueError(tag)

    def de(self, d):
        from importlib import import_module
        risch = import_module(self.S.__name__ + '.integrals.risch')
        ext = {'D': [self.poly(p) for p in d['D']], 'T': [self.expr(t) for t in d['T']]}
        if d['exts'] is not None:
            ext['exts'] = list(d['exts'])
        if d['extargs'] is not None:
            ext['extargs'] = [self.expr(a) if a is not None else None for a in d['extargs']]
        DE = risch.DifferentialExtension(extension=ext)
        DE.x = self.expr(d['x'])
        while DE.level > d['level']:
            DE.decrement_level()
        return DE
