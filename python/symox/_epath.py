"""Expression paths over native scalar trees and Python sequences."""
import re


class EPath:
    def __new__(cls, path):
        if isinstance(path, cls):
            return path
        return super().__new__(cls)

    def __init__(self, path):
        if isinstance(path, EPath):
            return
        if not path:
            raise ValueError('empty EPath')
        if not path.startswith('/'):
            raise NotImplementedError('only root paths are supported')
        self._path = path
        self._selectors = []
        for part in path[1:].split('/'):
            m = re.fullmatch(r'([A-Za-z_][A-Za-z_0-9?]*(?:\|[A-Za-z_][A-Za-z_0-9?]*)*)?(\*|\[-?\d*(?::-?\d*){0,2}\])?', part.strip())
            if m is None or not part.strip():
                raise ValueError('invalid path selector: %r' % part)
            query, span = m.groups()
            if span and span != '*':
                values = span[1:-1].split(':')
                if len(values) == 1:
                    span = int(values[0])
                else:
                    span = slice(*(int(v) if v else None for v in values))
                    if span.step == 0:
                        raise ValueError('slice step cannot be zero')
            else:
                span = slice(None)
            self._selectors.append(((query or '').split('|') if query else [], span))

    def __repr__(self):
        return 'EPath(%r)' % self._path

    @staticmethod
    def _children(expr):
        from . import Expr
        if isinstance(expr, Expr):
            if expr.is_Add:
                return list(expr.as_ordered_terms())
            if expr.is_Mul:
                return list(expr.as_ordered_factors())
            return list(expr.args)
        if isinstance(expr, (list, tuple)):
            return list(expr)
        return []

    @staticmethod
    def _matches(expr, query):
        if not query:
            return True
        attrs = [q[:-1] for q in query if q.endswith('?')]
        types = [q for q in query if not q.endswith('?')]
        names = {c.__name__ for c in type(expr).__mro__}
        name = getattr(expr, 'class_name', None)
        names.add(name)
        if name is not None:
            # Native Expr is a single Python class for all scalar AST heads.
            names.update(('Expr', 'Basic'))
            if getattr(expr, 'is_Number', False):
                names.add('Number')
            if getattr(expr,'is_Atom',False):
                names.update(('Atom','AtomicExpr'))
            if getattr(expr,'is_Integer',False):
                names.add('Rational')
            if getattr(expr,'is_Function',False):
                names.add('Function')
            if name=='Dummy':
                names.add('Symbol')
            names.add(getattr(getattr(expr, 'func', None), '__name__', None))
        return bool(set(types) & names) or bool(attrs and all(hasattr(expr, a) for a in attrs))

    def _walk(self, expr, depth, func, selected):
        if depth == len(self._selectors):
            if selected is not None:
                selected.append(expr)
                return expr
            return func(expr)
        children = self._children(expr)
        query, span = self._selectors[depth]
        indices = range(len(children))[span] if isinstance(span, slice) else [span]
        changed = False
        for i in indices:
            if not -len(children) <= i < len(children):
                continue
            value = children[i]
            if self._matches(value, query):
                children[i] = self._walk(value, depth+1, func, selected)
                changed |= children[i] != value
        if selected is not None or not changed:
            return expr
        if isinstance(expr, (tuple, list)):
            return type(expr)(children)
        return expr.func(*children)

    def select(self, expr):
        result = []
        self._walk(expr, 0, None, result)
        return result

    def apply(self, expr, func, args=None, kwargs=None):
        return self._walk(expr, 0, lambda e: func(e, *(args or ()), **(kwargs or {})), None)


def epath(path, expr=None, func=None, args=None, kwargs=None):
    path = EPath(path)
    if expr is None:
        return path
    return path.select(expr) if func is None else path.apply(expr, func, args, kwargs)
