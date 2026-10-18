"""General functional composition using native expressions and polynomials."""
from ..calculus._analysis import scalar, variable, walk


def decompogen(f, symbol):
    """Decompose a scalar expression into outermost-first compositions.

    Unary functions, powers with one varying argument, a common polynomial
    generator and native univariate polynomial decomposition are supported.
    An irreducible or unsupported decomposition returns the original function.
    Every returned chain reconstructs the input by composition.
    """
    from .. import Poly
    from ..polys.polyerrors import PolynomialError, CoercionFailed
    f=scalar(f)
    symbol=variable(symbol)
    list(walk(f))
    def decompose(expression,depth):
        if depth>256:raise NotImplementedError('Functional decomposition exceeds 256 levels')
        if symbol not in expression.free_symbols:return [expression]
        argument=None
        if expression.is_Pow:
            base,exponent=expression.args
            if symbol not in exponent.free_symbols:argument=base
            elif symbol not in base.free_symbols:argument=exponent
        elif expression.is_Function or expression.is_UndefinedFunction:
            varying={arg for arg in expression.args if symbol in arg.free_symbols}
            if len(varying)==1:argument=next(iter(varying))
        if argument is not None:
            if argument==symbol:return [expression]
            return [expression.xreplace({argument:symbol})]+decompose(argument,depth+1)
        head=getattr(expression.func,'__name__','')
        if head in ('Min','Max'):
            arguments=[]
            shared=None
            for arg in expression.args:
                if symbol not in arg.free_symbols:
                    arguments.append(arg)
                    continue
                chain=decompose(arg,depth+1)
                if len(chain)==1:chain=[symbol]+chain
                if shared is None:shared=chain[1:]
                elif shared!=chain[1:]:return [expression]
                arguments.append(chain[0])
            if shared is not None:
                outer=expression.func(*arguments)
                return shared if outer==symbol else [outer]+shared
        try:
            polynomial=Poly(expression)
            generators=[generator for generator in polynomial.gens if symbol in generator.free_symbols]
            if len(generators)==1 and generators[0]!=symbol:
                generator=generators[0]
                if generator==expression:return [expression]
                return [expression.xreplace({generator:symbol})]+decompose(generator,depth+1)
            polynomial=Poly(expression,symbol)
            return [part.as_expr() for part in polynomial.decompose()]
        except (ValueError,NotImplementedError,PolynomialError,CoercionFailed):
            return [expression]
    return decompose(f,0)


def compogen(g_s, symbol):
    """Compose a nonempty finite sequence, with the outer function first."""
    symbol=variable(symbol)
    expressions=[]
    for expression in g_s:
        if len(expressions)>=256:raise NotImplementedError('Functional composition exceeds 256 functions')
        expression=scalar(expression)
        list(walk(expression))
        expressions.append(expression)
    if not expressions:raise ValueError('Supply at least one function to compose')
    result=expressions[-1]
    for expression in reversed(expressions[:-1]):
        result=expression.subs(symbol,result)
        list(walk(result))
    return result


__all__=['decompogen','compogen']
