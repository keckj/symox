"""Typed structural selectors and occurrence queries over native graphs."""
from collections import Counter, defaultdict


_RELATIONS={'Eq':'Rel(==)','Ne':'Rel(!=)','Lt':'Rel(<)','Le':'Rel(<=)',
            'Gt':'Rel(>)','Ge':'Rel(>=)'}
_SPECIAL_HEADS={'log','Derivative','Piecewise','Order','Contains','And','Or',
    'Not','Xor','Implies','Equivalent','Nand','Nor','ITE','bell','hyper'}
_NON_FUNCTIONS={'Max','Min','Lambda','Order','And','Or','Not','Xor','Implies',
    'Equivalent','Nand','Nor','ITE','Contains','AppliedPredicate','Predicate'}
_NON_EXPR_HEADS=set(_RELATIONS.values())|(_NON_FUNCTIONS-{'Max','Min','Lambda','Order'})|{
    'BooleanTrue','BooleanFalse','true','false','Equality','Unequality',
    'LessThan','GreaterThan','StrictLessThan','StrictGreaterThan'}


def _head(value):
    return getattr(getattr(value,'func',None),'__name__','')


def _is_expr(value):
    from . import Expr, MatrixExpr, Tuple
    if isinstance(value,MatrixExpr):return True
    if not isinstance(value,Expr) or isinstance(value,Tuple):return False
    if getattr(value,'is_UndefinedFunction',False):return True
    return _head(value) not in _NON_EXPR_HEADS


def _is_function(value):
    if getattr(value,'is_UndefinedFunction',False):return True
    return (_head(value)=='Piecewise' or bool(getattr(value,'is_Function',False))) and _head(value) not in _NON_FUNCTIONS


def matcher(query, matching=True):
    from . import (Expr, Function, Atom, AtomicExpr, Tuple, AlgebraicNumber,
                  sympify, _is_function_head, __dict__ as namespace)
    from . import _core
    from ._node_types import native_module
    from .combinatorics import Permutation
    from ._matrixexpr_methods import _NameAtom, _TupleView
    if type(query) is str:query=sympify(query)
    # The defining module: native classes report SymPy's under the facade.
    foreign=isinstance(query,type) and str(native_module(query) or '').startswith('sympy.')
    undefined=isinstance(query,type) and any(base.__name__=='UndefinedFunction'
        and base.__module__.startswith('sympy.') for base in type(query).__mro__)
    if undefined:
        nargs=query.nargs
        nargs=tuple(map(int,nargs)) if getattr(nargs,'is_FiniteSet',False) else None
        query=Function(query.__name__,nargs=nargs,**dict(query._kwargs))
    elif foreign and query.__name__ in namespace:
        query=namespace[query.__name__]
    if _is_function_head(query):
        return lambda value:bool(getattr(value,'is_UndefinedFunction',False)) and value.func==query
    if isinstance(query,Permutation):return lambda value:isinstance(value,Permutation) and value==query
    if isinstance(query,AlgebraicNumber):return lambda value:isinstance(value,AlgebraicNumber) and value==query
    if isinstance(query,_core.Head):return lambda value:getattr(value,'func',None)==query
    if query is Expr:return _is_expr
    if query is Function:return _is_function
    if query is Atom:return lambda value:bool(getattr(value,'is_Atom',False)) and _head(value) not in {
        'BooleanTrue','BooleanFalse','true','false'}
    if query is AtomicExpr:return lambda value:_is_expr(value) and bool(getattr(value,'is_Atom',False))
    if query is Tuple:return lambda value:isinstance(value,(Tuple,_TupleView))
    if isinstance(query,type):return lambda value:isinstance(value,query)
    constructor=getattr(query,'_symox_head_name',None)
    name=getattr(query,'__name__',None)
    if not constructor and (name in _SPECIAL_HEADS or name in _RELATIONS) and namespace.get(name) is query:
        constructor=_RELATIONS.get(name,name)
    if constructor:
        return lambda value:not getattr(value,'is_UndefinedFunction',False) and _head(value)==constructor
    if callable(query) and not isinstance(query,Expr):return lambda value:bool(query(value))
    if not isinstance(query,(Expr,_NameAtom)):
        try:query=sympify(query)
        except (TypeError,ValueError):pass
    def test(value):
        if isinstance(value,_NameAtom):return isinstance(query,_NameAtom) and value==query
        if value==query:return True
        method=getattr(value,'match',None)
        return matching and method is not None and method(query) is not None
    return test


def _graph(expression):
    from . import Expr, Set, _wrap_set, sympify
    from ._matrixexpr_methods import _children, _TupleView
    owners={}
    def normalize(value):
        if isinstance(value,Expr):return _wrap_set(value)
        if isinstance(value,int):return sympify(value)
        if isinstance(value,_TupleView):return value
        if isinstance(value,(tuple,list)):return _TupleView(value)
        return value
    def key(value):
        underlying=value._e if isinstance(value,Set) else value
        identity=getattr(underlying,'_dag_id',None)
        if isinstance(value,_TupleView):identifier=(type(value),id(value))
        elif identity is not None:identifier=(type(value),identity)
        else:
            try:hash(value);identifier=(type(value),value)
            except TypeError:identifier=(type(value),id(value))
        owners[identifier]=value
        return identifier
    root=normalize(expression)
    root_key=key(root)
    records={}
    active=set()
    order=[]
    pending=[(root,False)]
    while pending:
        value,visited=pending.pop()
        identifier=key(value)
        if visited:
            active.remove(identifier)
            order.append(identifier)
            continue
        if identifier in active:raise ValueError('Cyclic structural-query input')
        if identifier in records:continue
        children=[normalize(child) for child in _children(value)]
        records[identifier]=(value,[key(child) for child in children])
        if len(records)>100000:
            raise NotImplementedError('Structural query exceeds 100000 distinct nodes')
        active.add(identifier)
        pending.append((value,True))
        pending.extend((child,False) for child in reversed(children))
    return root_key,records,order


def find(self, query, group=False):
    test=matcher(query)
    root,records,order=_graph(self)
    if not group:return {value for value,_ in records.values() if test(value)}
    weights=defaultdict(int,{root:1})
    found=Counter()
    for identifier in reversed(order):
        value,children=records[identifier]
        weight=weights[identifier]
        if test(value):found[value]+=weight
        for child in children:weights[child]+=weight
    return dict(found)


def count(self, query):
    return sum(find(self,query,group=True).values())
