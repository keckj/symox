"""Leading powers expressed in the original variable along a chosen ray."""
from . import _core


def _data(expression, variable, logx, cdir, origin=0):
    from . import S, Dummy, log, sympify
    variable=sympify(variable)
    origin=sympify(origin)
    direction=sympify(cdir)
    if not variable.is_Symbol:raise ValueError('Leading term variable must be a symbol')
    if origin.has(variable):raise ValueError('Expansion center must be independent of the variable')
    if origin.is_finite is False:
        raise NotImplementedError('Use asymptotic series for an infinite expansion center')
    if direction==0:direction=S.One
    if direction.has(variable) or direction.is_zero is not False or direction.is_finite is False:
        raise ValueError('Approach direction must be a finite, nonzero constant')
    if direction.is_positive is True:direction=S.One
    elif direction.is_negative is True:direction=S.NegativeOne
    if expression==0:return S.Zero,S.Zero
    # SymPy tests free occurrences (has_free): a bound x, as in
    # Derivative(1, x), leaves the expression as its own leading term
    if variable not in expression.free_symbols:return expression,S.Zero
    small=Dummy('leading_t',positive=True)
    marker=Dummy('leading_log')
    translated=expression.subs(variable,origin+direction*small)
    coefficient,power=_core._leadterm(translated,small)
    # The native kernel's exponent is in a positive local parameter t.
    # Divide out d**power to return the coefficient of (x-origin)**power.
    coefficient=coefficient.xreplace({log(small):marker-log(direction)})/direction**power
    if small in coefficient.free_symbols:
        raise ValueError('Leading coefficient still depends on the expansion variable')
    formal_log=log(variable-origin) if logx is None else sympify(logx)
    return coefficient.xreplace({marker:formal_log}),power


def leadterm(self, x, logx=None, cdir=0, **options):
    origin=options.pop('x0',0)
    direction=options.pop('dir',None)
    if options:raise TypeError('Unknown leadterm options: '+', '.join(sorted(options)))
    if direction is not None:
        if direction not in ('+','-'):raise ValueError('Direction must be + or -')
        if cdir!=0:raise TypeError('Supply the approach direction only once')
        cdir=1 if direction=='+' else -1
    return _data(self,x,logx,cdir,origin)


def as_leading_term(self, *symbols, logx=None, cdir=0):
    from . import sympify
    result=self
    for symbol in symbols:
        symbol=sympify(symbol)
        coefficient,power=_data(result,symbol,logx,cdir)
        result=coefficient*symbol**power
    return result
