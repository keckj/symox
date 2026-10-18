"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'cse')

from collections import defaultdict
from dataclasses import dataclass
import symox as S

@dataclass(frozen=True)
class Unevaluated:
    func:object
    args:tuple
    def __init__(self,func,args):object.__setattr__(self,'func',func);object.__setattr__(self,'args',tuple(args))
    @property
    def free_symbols(self):return set().union(*(getattr(arg,'free_symbols',set()) for arg in self.args))
    def _symox_(self):
        args=[S.sympify(arg) for arg in self.args]
        try:return self.func(*args,evaluate=False)
        except TypeError:return self.func(*args)
    def as_unevaluated_basic(self):return self._symox_()


class FuncArgTracker:
    def __init__(self,funcs):
        self.value_numbers={};self.value_number_to_value=[];self.arg_to_funcset=[];self.func_to_argset=[]
        for func_i,function in enumerate(funcs):
            numbers={self.get_or_add_value_number(arg) for arg in function.args};self.func_to_argset.append(numbers)
            for number in numbers:self.arg_to_funcset[number].add(func_i)
    def get_or_add_value_number(self,value):
        if value not in self.value_numbers:
            self.value_numbers[value]=len(self.value_number_to_value);self.value_number_to_value.append(value);self.arg_to_funcset.append(set())
        return self.value_numbers[value]
    def get_args_in_value_order(self,argset):return [self.value_number_to_value[i] for i in sorted(argset)]
    def get_common_arg_candidates(self,argset,min_func_i=0):
        counts=defaultdict(int)
        for arg in argset:
            for function in self.arg_to_funcset[arg]:
                if function>=min_func_i:counts[function]+=1
        return {function:count for function,count in counts.items() if count>=2}
    def get_subset_candidates(self,argset,restrict_to_funcset=None):
        candidates=set(range(len(self.func_to_argset))) if restrict_to_funcset is None else set(restrict_to_funcset)
        for argument in argset:candidates&=self.arg_to_funcset[argument]
        return candidates
    def update_func_argset(self,func_i,new_argset):
        previous=self.func_to_argset[func_i];new=set(new_argset)
        for argument in previous-new:self.arg_to_funcset[argument].discard(func_i)
        for argument in new-previous:self.arg_to_funcset[argument].add(func_i)
        self.func_to_argset[func_i]=new
    def stop_arg_tracking(self,func_i):self.update_func_argset(func_i,set())


def reps_toposort(r):
    pending=[(symbol,S.sympify(expression)) for symbol,expression in r];keys={symbol for symbol,expression in pending};done=set();ordered=[]
    while pending:
        ready=next((i for i,(symbol,expression) in enumerate(pending) if not (expression.free_symbols&keys)-done),None)
        if ready is None:raise ValueError('cyclic common-subexpression replacements')
        symbol,expression=pending.pop(ready);ordered.append((symbol,expression));done.add(symbol)
    return ordered


def cse_separate(r,e):
    replacements=list(r);expressions=[]
    for expression in e:
        if getattr(expression,'class_name','')=='Equality' and expression.lhs.is_Symbol:replacements.append((expression.lhs,expression.rhs))
        else:expressions.append(expression)
    return [reps_toposort(replacements),expressions]


def cse_release_variables(r,e):
    if not r:return r,e
    occupied={str(s) for s,expr in r}|{str(s) for expr in e for s in expr.free_symbols}
    outputs=[]
    for i in range(len(e)):
        name='_'+str(i)
        while name in occupied:name='_'+name
        occupied.add(name);outputs.append(S.Symbol(name))
    assignments=list(r)+list(zip(outputs,e));last={}
    for i,(symbol,expression) in enumerate(assignments):
        for value in expression.free_symbols:last[value]=i
    temporaries={symbol for symbol,expression in r};result=[]
    for i,(symbol,expression) in enumerate(assignments):
        result.append((symbol,expression))
        result.extend((value,None) for value in sorted(expression.free_symbols&temporaries,key=str) if last[value]==i)
        if symbol in temporaries and symbol not in last:result.append((symbol,None))
    return result,tuple(outputs)


def preprocess_for_cse(expr,optimizations):
    for before,after in optimizations:
        if before is not None:expr=before(expr)
    return expr

def postprocess_for_cse(expr,optimizations):
    for before,after in reversed(optimizations):
        if after is not None:expr=after(expr)
    return expr


def match_common_args(func_class,funcs,opt_subs):
    tracker=FuncArgTracker(funcs)
    for i,function in enumerate(funcs):
        for j,count in sorted(tracker.get_common_arg_candidates(tracker.func_to_argset[i],i+1).items()):
            common=tracker.func_to_argset[i]&tracker.func_to_argset[j]
            if len(common)<2:continue
            grouped=Unevaluated(func_class,tracker.get_args_in_value_order(common));number=tracker.get_or_add_value_number(grouped)
            for k in (i,j):
                remaining=(tracker.func_to_argset[k]-common)|{number};tracker.update_func_argset(k,remaining)
                opt_subs[funcs[k]]=Unevaluated(func_class,tracker.get_args_in_value_order(remaining))


def opt_cse(exprs,order='canonical'):
    families={S.Add:[],S.Mul:[]};seen=set();result={}
    for expression in exprs:
        for node in S.preorder_traversal(expression):
            if node in seen:continue
            seen.add(node)
            if node.is_Add:families[S.Add].append(node)
            elif node.is_Mul and node.is_commutative:families[S.Mul].append(node)
    for cls,functions in families.items():match_common_args(cls,functions,result)
    return result


def tree_cse(exprs,symbols,opt_subs=None,order='canonical',ignore=()):
    substitutions={old:S.sympify(new) for old,new in (opt_subs or {}).items()}
    expressions=[S.sympify(expr).xreplace(substitutions) for expr in exprs]
    return S.cse(expressions,symbols=symbols,order=order,ignore=ignore)

__all__ += ['FuncArgTracker','Unevaluated','cse_release_variables','cse_separate','match_common_args','opt_cse',
 'postprocess_for_cse','preprocess_for_cse','reps_toposort','tree_cse']
