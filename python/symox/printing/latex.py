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

    def _print_NDimArray(self, expr):
        # SymPy's LaTeX layout of N-dim arrays: nested matrix blocks.
        import itertools
        if expr.rank() == 0:
            return self._print(expr[()])
        settings = self._settings
        mat_str = settings.get('mat_str')
        if mat_str is None:
            if settings.get('mode', 'plain') == 'inline':
                mat_str = 'smallmatrix'
            elif expr.shape[-1] <= 10:
                mat_str = 'matrix'
            else:
                mat_str = 'array'
        block_str = r'\begin{%MATSTR%}%s\end{%MATSTR%}'.replace('%MATSTR%', mat_str)
        if mat_str == 'array':
            block_str = block_str.replace('%s', '{' + 'c'*expr.shape[0] + '}%s')
        left_delim = settings.get('mat_delim', '[')
        if left_delim:
            right_delim = {'(': ')', '[': ']'}[left_delim]
            block_str = r'\left' + left_delim + block_str + r'\right' + right_delim
        level_str = [[] for i in range(expr.rank() + 1)]
        for outer_i in itertools.product(*(range(i) for i in expr.shape)):
            level_str[-1].append(self._print(expr[outer_i]))
            even = True
            for back_outer_i in range(expr.rank()-1, -1, -1):
                if len(level_str[back_outer_i+1]) < expr.shape[back_outer_i]:
                    break
                if even:
                    level_str[back_outer_i].append(r' & '.join(level_str[back_outer_i+1]))
                else:
                    level_str[back_outer_i].append(block_str % (r'\\'.join(level_str[back_outer_i+1])))
                    if len(level_str[back_outer_i+1]) == 1:
                        level_str[back_outer_i][-1] = r'\left[' + level_str[back_outer_i][-1] + r'\right]'
                even = not even
                level_str[back_outer_i+1] = []
        out_str = level_str[0][0]
        if expr.rank() % 2 == 1:
            out_str = block_str % out_str
        return out_str

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
