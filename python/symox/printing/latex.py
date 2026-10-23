"""LaTeX output using the native printer, including custom _latex hooks."""
from ._legacy import LegacyPrinter as Printer

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

    # SymPy's printers of Python containers, explicit matrices and objects
    # that are not expressions (the native printer takes expressions only).
    def _print(self, expr, **kwargs):
        if (getattr(expr, 'is_Matrix', False) and hasattr(expr, 'tolist')
                and hasattr(expr, 'rows') and not isinstance(expr, type)):
            if not (self.printmethod and hasattr(type(expr), self.printmethod)):
                return self._print_MatrixBase(expr)
        return super()._print(expr, **kwargs)

    def _separator(self):
        separator = self._settings.get('decimal_separator', 'period')
        if separator == 'comma':
            return ';'
        if separator == 'period':
            return ','
        raise ValueError('Unknown Decimal Separator')

    def _add_parens_lspace(self, s):
        return r"\left( {}\right)".format(s)

    def _print_tuple(self, expr):
        sep = self._separator()
        if len(expr) == 1:
            # 1-tuple needs a trailing separator
            return self._add_parens_lspace(self._print(expr[0]) + sep)
        return self._add_parens_lspace(
            (sep + r" \  ").join([self._print(i) for i in expr]))

    _print_Tuple = _print_tuple

    def _print_Float(self, expr):
        text = self.emptyPrinter(expr)
        if self._settings.get('decimal_separator', 'period') == 'comma':
            text = text.replace('.', '{,}')
        return text

    def _print_list(self, expr):
        sep = self._separator()
        return r"\left[ %s\right]" % (sep + r" \  ").join([self._print(i) for i in expr])

    def _print_dict(self, d):
        from symox import default_sort_key
        keys = sorted(d.keys(), key=default_sort_key)
        items = ["%s : %s" % (self._print(key), self._print(d[key])) for key in keys]
        return r"\left\{ %s\right\}" % r", \  ".join(items)

    def _print_set(self, s):
        from symox import default_sort_key
        items = sorted(s, key=default_sort_key)
        sep = self._separator()
        return r"\left\{%s\right\}" % (sep + " ").join(map(self._print, items))

    _print_frozenset = _print_set

    def _print_Dict(self, expr):
        return self._print_dict(dict(expr))

    def _print_bool(self, e):
        return r"\text{%s}" % e

    _print_NoneType = _print_bool

    def _print_matrix_contents(self, expr):
        lines = []
        for line in range(expr.rows):
            lines.append(" & ".join([self._print(i) for i in expr[line, :]]))
        mat_str = self._settings.get('mat_str')
        if mat_str is None:
            if self._settings.get('mode', 'plain') == 'inline':
                mat_str = 'smallmatrix'
            elif (expr.cols <= 10) is True:
                mat_str = 'matrix'
            else:
                mat_str = 'array'
        out_str = r'\begin{%MATSTR%}%s\end{%MATSTR%}'
        out_str = out_str.replace('%MATSTR%', mat_str)
        if mat_str == 'array':
            out_str = out_str.replace('%s', '{' + 'c'*expr.cols + '}%s')
        return out_str % r"\\".join(lines)

    def _print_MatrixBase(self, expr):
        out_str = self._print_matrix_contents(expr)
        left_delim = self._settings.get('mat_delim', '[')
        if left_delim:
            right_delim = {'(': ')', '[': ']'}[left_delim]
            out_str = r'\left' + left_delim + out_str + r'\right' + right_delim
        return out_str

    def emptyPrinter(self,expr):
        from symox import _core
        from symox.core.basic import Basic
        if isinstance(expr, (str, type)) or not isinstance(expr, (Basic, int, float, complex)):
            # SymPy prints objects that are not expressions monospaced; a
            # string is not parsed.
            converted = None
            if not isinstance(expr, (str, type)):
                try:
                    converted = _core.sympify(expr)
                except Exception:
                    converted = None
            if converted is None:
                return r"\mathtt{\text{%s}}" % latex_escape(str(expr))
            expr = converted
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
    # Symbol subclasses with their own _latex, such as diffgeom's
    # BaseScalarField (\mathbf{x}), print through symbol_names.
    try:symbols = expr.free_symbols if hasattr(expr,'free_symbols') else ()
    except Exception:symbols = ()
    printer = LatexPrinter(settings)
    own = {s:s._latex(printer) for s in symbols if '_latex' in type(s).__dict__}
    if own:printer = LatexPrinter({**settings,'symbol_names':{**own,**settings.get('symbol_names',{})}})
    return printer.doprint(expr)

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
