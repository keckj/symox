"""Typed SMT-LIB expressions and declarations for native scalar constraints."""
import re
from ._legacy import LegacyPrinter as Printer


class SMTLibPrinter(Printer):
    printmethod='_smtlib'
    def __init__(self,settings=None,symbol_table=None):
        super().__init__(settings)
        self.symbol_table=dict(symbol_table or {})
        self.types={bool:'Bool',int:'Int',float:'Real',**self._settings.get('known_types',{})}
        self.constants=self._settings.get('known_constants',{})
        self.functions=self._settings.get('known_functions',{})
    def name(self,value):
        name=str(value)
        if re.fullmatch(r'[a-zA-Z_+*/=<>!?$%&~^.-][a-zA-Z0-9_+*/=<>!?$%&~^.-]*',name) and name not in ('true','false','let','forall','exists'):
            return name
        if '|' in name or '\\' in name:raise ValueError('SMT-LIB quoted symbols cannot contain vertical bars or backslashes')
        return '|'+name+'|'
    def _print(self,expr,**kwargs):
        from symox import sympify
        if not hasattr(expr,'class_name'):expr=sympify(expr)
        if expr in self.constants:return self.constants[expr]
        name=expr.class_name
        if expr.is_Symbol:return self.name(expr.name)
        if expr.is_Integer:
            return str(int(expr))
        if expr.is_Rational:return '(/ '+str(expr.p)+' '+str(expr.q)+')'
        if expr.is_Float:return self._print_float(str(expr))
        if name in ('BooleanTrue','BooleanFalse'):return name.removeprefix('Boolean').lower()
        if name=='Piecewise':
            from .mathml import _pairs
            pairs=_pairs(expr)
            if str(pairs[-1][1]) not in ('True','true'):raise ValueError('Piecewise requires a default branch')
            result=self._print(pairs[-1][0])
            for value,condition in reversed(pairs[:-1]):result='(ite '+self._print(condition)+' '+self._print(value)+' '+result+')'
            return result
        mapping={'Add':'+','Mul':'*','Pow':'pow','Equality':'=','Unequality':'distinct',
                 'StrictLessThan':'<','LessThan':'<=','StrictGreaterThan':'>','GreaterThan':'>=',
                 'And':'and','Or':'or','Not':'not','Xor':'xor','Implies':'=>','ITE':'ite',
                 'Abs':'abs','Min':'min','Max':'max','exp':'exp','log':'log','sin':'sin','cos':'cos',
                 'tan':'tan','asin':'arcsin','acos':'arccos','atan':'arctan','atan2':'arctan2',
                 'sinh':'sinh','cosh':'cosh','tanh':'tanh'}
        operator=self.functions.get(expr.func,self.functions.get(name,mapping.get(name)))
        if operator is None and expr.is_UndefinedFunction:operator=self.name(expr.name)
        if operator is None:raise NotImplementedError('No SMT-LIB mapping for '+name)
        return '('+operator+' '+' '.join(self._print(a) for a in expr.args)+')'
    @staticmethod
    def _print_float(value):
        from decimal import Decimal
        text=format(Decimal(str(value)),'f')
        if not Decimal(text).is_finite():raise ValueError('SMT-LIB does not represent NaN or infinity')
        if '.' not in text:text+='.0'
        return '(- '+text[1:]+')' if text.startswith('-') else text


def _signature(value):
    """Argument and result types of a ``Callable[[...], T]`` symbol-table entry."""
    args=getattr(value,'__args__',None)
    if args is None and isinstance(value,(tuple,list)):args=tuple(value)
    return tuple(args) if args is not None else None


def _function_key(table,node):
    for key in (node.func,node.name):
        try:
            if key in table:return key
        except TypeError:pass
    return None


def _undefined_calls(expr):
    from symox import preorder_traversal
    return [node for node in preorder_traversal(expr) if getattr(node,'is_UndefinedFunction',False)]


def _is_boolean_node(e):
    return (getattr(e,'is_Relational',False) or getattr(e,'is_Boolean',False)
            or getattr(e,'class_name','') in ('And','Or','Not','Xor','Implies','ITE','BooleanTrue','BooleanFalse','Equivalent','Nand','Nor'))


def _auto_infer_smtlib_types(*exprs,symbol_table=None):
    # sympy.printing.smtlib._auto_infer_smtlib_types
    from symox import preorder_traversal
    symbols=dict(symbol_table) if symbol_table else {}
    def safe_update(syms,inferred):
        for sym in syms:
            old=symbols.setdefault(sym,inferred)
            if old!=inferred:
                raise TypeError(f"Could not infer type of `{sym}`. Apparently both `{old}` and `{inferred}`?")
    nodes=[list(preorder_traversal(e)) for e in exprs]
    safe_update({e for e in exprs if e.is_Symbol},bool)
    safe_update({a for ns in nodes for n in ns if _is_boolean_node(n) and not getattr(n,'is_Relational',False) and not getattr(n,'is_Symbol',False)
                 for a in n.args if getattr(a,'is_Symbol',False)},bool)
    for wanted in (bool,int):
        found=set()
        for e in exprs:
            for call in _undefined_calls(e):
                key=_function_key(symbols,call)
                if key is None:continue
                params=_signature(symbols[key]) or ()
                found|={a for a,param in zip(call.args,params) if a.is_Symbol and param==wanted}
        safe_update(found,wanted)
    atoms={n for ns in nodes for n in ns if getattr(n,'is_Symbol',False)}
    safe_update({a for a in atoms if a.is_integer},int)
    safe_update({a for a in atoms if a.is_real and not a.is_integer},float)
    equalities=[n for ns in nodes for n in ns if getattr(n,'class_name','')=='Equality']
    pairs=[(r.lhs,r.rhs) for r in equalities if r.lhs.is_Symbol]+[(r.rhs,r.lhs) for r in equalities if r.rhs.is_Symbol]
    for infer,related in pairs:
        key=_function_key(symbols,related) if getattr(related,'is_UndefinedFunction',False) else None
        inference=(symbols[infer] if infer in symbols else
                   symbols[related] if related in symbols else
                   (_signature(symbols[key]) or (None,))[-1] if key is not None else
                   bool if getattr(related,'is_Boolean',False) or _is_boolean_node(related) else
                   int if related.is_integer or related.is_Integer else
                   float if related.is_real else None)
        if inference:safe_update({infer},inference)
    return symbols


def smtlib_code(expr,auto_assert=True,auto_declare=True,precision=None,symbol_table=None,
                known_types=None,known_constants=None,known_functions=None,
                prefix_expressions=None,suffix_expressions=None,log_warn=None):
    from symox import sympify
    log_warn=log_warn or (lambda _:None)
    expressions=[sympify(e) for e in expr] if isinstance(expr,(list,tuple)) else [sympify(expr)]
    table=_auto_infer_smtlib_types(*expressions,symbol_table=symbol_table or {})
    printer=SMTLibPrinter({'precision':precision,'known_types':known_types or {},
                          'known_constants':known_constants or {},'known_functions':known_functions or {}},table)
    # [FALLBACK RULES] of the reference: untyped symbols become Real with a
    # warning, an undefined function needs a known name or a signature.
    for e in expressions:
        for symbol in sorted(e.free_symbols,key=lambda s:s.name):
            if symbol not in printer.constants and symbol not in printer.symbol_table:
                log_warn(f"Could not infer type of `{symbol}`. Defaulting to float.")
                printer.symbol_table[symbol]=float
        for call in _undefined_calls(e):
            if _function_key(printer.functions,call) is None and _function_key(printer.symbol_table,call) is None:
                raise TypeError(f"Unknown type of undefined function `{call}`. "
                                "Must be mapped to ``str`` in known_functions or mapped to ``Callable[..]`` in symbol_table.")
    declarations=[]
    if auto_declare:
        symbols=set().union(*(e.free_symbols for e in expressions))
        emitted={}
        for symbol in sorted(symbols,key=lambda s:s.name):
            name=printer.name(symbol.name)
            if name in emitted and emitted[name]!=symbol:raise ValueError('Ambiguous SMT-LIB symbol name: '+name)
            emitted[name]=symbol
            if symbol in printer.constants:continue
            kind=printer.symbol_table[symbol]
            if kind not in printer.types:raise TypeError('Unknown SMT-LIB symbol type: '+repr(kind))
            declarations.append('(declare-const '+name+' '+printer.types[kind]+')')
        functions={node.name:node for e in expressions for node in _undefined_calls(e)}
        for name,node in sorted(functions.items()):
            if _function_key(printer.functions,node) is not None:continue
            signature=_signature(printer.symbol_table[_function_key(printer.symbol_table,node)])
            if not signature or len(signature)!=len(node.args)+1:
                raise TypeError('Undefined SMT-LIB functions need an argument/result type signature in symbol_table')
            types=[printer.types[t] for t in signature]
            declarations.append('(declare-fun '+printer.name(name)+' ('+' '.join(types[:-1])+') '+types[-1]+')')
    body=[]
    for e in expressions:
        text=printer.doprint(e)
        key=_function_key(printer.symbol_table,e) if getattr(e,'is_UndefinedFunction',False) else None
        boolean=(_is_boolean_node(e) or printer.symbol_table.get(e) is bool
                 or key is not None and (_signature(printer.symbol_table[key]) or (None,))[-1] is bool)
        if auto_assert and boolean:text='(assert '+text+')'
        elif auto_assert:log_warn(f"Non-Boolean expression `{e}` will not be asserted. Converting to SMTLib verbatim.")
        body.append(text)
    def render(items):return [e if isinstance(e,str) else printer.doprint(e) for e in items or ()]
    return '\n'.join([*render(prefix_expressions),*sorted(declarations),*body,*render(suffix_expressions)])
