"""LaTeX output using the native printer, including custom _latex hooks."""
from .printer import Printer

class LatexPrinter(Printer):
    printmethod='_latex'
    def emptyPrinter(self,expr):
        from symox import _core
        return _core.latex(expr,**self._settings)

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
