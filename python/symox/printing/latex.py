"""LaTeX output using the native printer, including custom _latex hooks."""
from .printer import Printer

class LatexPrinter(Printer):
    printmethod='_latex'
    _default_settings = dict(Printer._default_settings,
                             symbol_names={}, disable_split_super_sub=False)

    def _print_Symbol(self, expr, style='plain'):
        name = self._settings['symbol_names'].get(expr)
        if name is not None:
            return name
        return self._deal_with_super_sub(expr.name, style=style)

    _print_RandomSymbol = _print_Symbol

    def _split_super_sub(self, name):
        from .conventions import split_super_sub
        if name is None or '{' in name:
            return name, [], []
        if self._settings['disable_split_super_sub']:
            name, supers, subs = name.replace('_', r'\_').replace('^', r'\^'), [], []
        else:
            name, supers, subs = split_super_sub(name)
        return translate(name), [translate(s) for s in supers], [translate(s) for s in subs]

    def _deal_with_super_sub(self, string, style='plain'):
        name, supers, subs = self._split_super_sub(string)
        if style == 'bold':
            name = r'\mathbf{' + name + '}'
        if supers:
            name += '^{' + ' '.join(supers) + '}'
        if subs:
            name += '_{' + ' '.join(subs) + '}'
        return name

    def emptyPrinter(self,expr):
        from symox import _core
        settings = dict(self._settings, mode='plain')
        if settings.get('fold_short_frac') is None:
            settings['fold_short_frac'] = self._settings.get('mode') == 'inline'
        return _core.latex(expr, **settings)

    def doprint(self, expr):
        mode = self._settings.get('mode', 'plain')
        if mode not in ('plain', 'inline', 'equation', 'equation*'):
            raise ValueError("'mode' must be one of 'inline', 'plain', 'equation' or 'equation*'")
        result = self._print(expr)
        if mode == 'plain':
            return result
        if mode == 'inline':
            return '$' + result + '$'
        if self._settings.get('itex', False):
            return '$$' + result + '$$'
        return r'\begin{' + mode + '}' + result + r'\end{' + mode + '}'


def latex(expr,**settings):
    return LatexPrinter(settings).doprint(expr)

def print_latex(expr,**settings):
    print(latex(expr,**settings))


def multiline_latex(lhs,rhs,terms_per_line=1,environment='align*',use_dots=False,**settings):
    from symox import sympify
    import operator
    terms_per_line=operator.index(terms_per_line)
    if terms_per_line<1:raise ValueError('terms_per_line must be positive')
    lhs,rhs=sympify(lhs),sympify(rhs)
    if environment not in ('align*','eqnarray','IEEEeqnarray'):raise ValueError('Unknown environment: '+environment)
    array=environment!='align*'
    lines=[];terms=rhs.as_ordered_terms()
    for start in range(0,len(terms),terms_per_line):
        chunk=[]
        for index,term in enumerate(terms[start:start+terms_per_line],start):
            minus=term.could_extract_minus_sign()
            text=latex(-term if minus else term,**settings)
            chunk.append(('- ' if minus else '+ ' if index else '')+text)
        prefix=(latex(lhs,**settings)+(' & = & ' if array else ' = & ')) if start==0 else ('& & ' if array else '& ')
        suffix=''
        if start+terms_per_line<len(terms):suffix=(' \\dots' if use_dots else '')+('\\nonumber' if array else '')+r'\\'
        lines.append(prefix+' '.join(chunk)+suffix)
    return '\\begin{'+environment+'}'+('{rCl}' if environment=='IEEEeqnarray' else '')+'\n'+'\n'.join(lines)+'\n\\end{'+environment+'}'

accepted_latex_functions='arcsin arccos arctan sin cos tan sinh cosh tanh sqrt ln log sec csc cot coth re im frac root arg'.split()


def latex_escape(s):
    replacements={'\\':r'\textbackslash','~':r'\textasciitilde','^':r'\textasciicircum'}
    replacements.update({char:'\\'+char for char in '&%$#_{}'})
    return ''.join(replacements.get(char,char) for char in str(s))


def translate(s):
    from symox import _core
    # Translate a single token without applying symbol subscript grammar twice.
    return _core._latex_translate(str(s))
