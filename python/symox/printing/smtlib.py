"""Typed SMT-LIB expressions and declarations for native scalar constraints."""
import re
from .printer import Printer


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
            value=int(expr)
            return str(value) if value>=0 else '(- '+str(-value)+')'
        if expr.is_Rational:return '(/ '+self._print_float(expr.p)+' '+self._print_float(expr.q)+')'
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


def smtlib_code(expr,auto_assert=True,auto_declare=True,precision=None,symbol_table=None,
                known_types=None,known_constants=None,known_functions=None,
                prefix_expressions=None,suffix_expressions=None,log_warn=None):
    from symox import sympify,preorder_traversal
    expressions=[sympify(e) for e in expr] if isinstance(expr,(list,tuple)) else [sympify(expr)]
    printer=SMTLibPrinter({'precision':precision,'known_types':known_types or {},
                          'known_constants':known_constants or {},'known_functions':known_functions or {}},symbol_table)
    lines=list(prefix_expressions or ())
    if auto_declare:
        symbols=set().union(*(e.free_symbols for e in expressions))
        emitted={}
        for symbol in sorted(symbols,key=lambda s:s.name):
            name=printer.name(symbol.name)
            if name in emitted and emitted[name]!=symbol:raise ValueError('Ambiguous SMT-LIB symbol name: '+name)
            emitted[name]=symbol
            if symbol in printer.constants:continue
            kind=printer.symbol_table.get(symbol)
            if kind is None:
                kind=int if symbol.is_integer is True else bool if getattr(symbol,'is_Boolean',False) else float
                if log_warn:log_warn('Inferred '+printer.types[kind]+' for '+str(symbol))
            if kind not in printer.types:raise TypeError('Unknown SMT-LIB symbol type: '+repr(kind))
            lines.append('(declare-const '+printer.name(symbol.name)+' '+printer.types[kind]+')')
        functions={node.name:node for e in expressions for node in preorder_traversal(e) if getattr(node,'is_UndefinedFunction',False)}
        for name,node in sorted(functions.items()):
            if name in printer.functions or node.func in printer.functions:continue
            signature=printer.symbol_table.get(node.func,printer.symbol_table.get(name))
            if not isinstance(signature,(tuple,list)) or len(signature)!=len(node.args)+1:
                raise TypeError('Undefined SMT-LIB functions need an argument/result type tuple in symbol_table')
            types=[printer.types[t] for t in signature]
            lines.append('(declare-fun '+printer.name(name)+' ('+' '.join(types[:-1])+') '+types[-1]+')')
    for e in expressions:
        text=printer.doprint(e)
        boolean=getattr(e,'is_Relational',False) or e.class_name in ('And','Or','Not','Xor','Implies','ITE','BooleanTrue','BooleanFalse') or printer.symbol_table.get(e) is bool
        lines.append('(assert '+text+')' if auto_assert and boolean else text)
    lines.extend(suffix_expressions or ())
    return '\n'.join(lines)
