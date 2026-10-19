"""Operation labels over native scalar and matrix expression graphs."""
from collections import Counter, defaultdict


def count_ops(expression, visual=False):
    from . import S, Symbol, Add, fraction, Expr, MatrixExpr, Matrix, Set, sympify
    from ._matrixexpr_methods import _children

    relations={'Rel(==)':'EQ','Rel(!=)':'NE','Rel(<)':'LT','Rel(<=)':'LE',
               'Rel(>)':'GT','Rel(>=)':'GE','Equality':'EQ','Unequality':'NE',
               'StrictLessThan':'LT','LessThan':'LE','StrictGreaterThan':'GT','GreaterThan':'GE'}
    booleans={'And','Or','Not','Xor','Implies','Equivalent','Nand','Nor','ITE','Contains'}
    product_sizes={}

    def multiplication_size(value):
        pending=[(value,False)]
        while pending:
            node,visited=pending.pop()
            if node in product_sizes:continue
            if (getattr(node,'is_Pow',False) and node.base.is_number
                    and node.exp.is_Integer and abs(int(node.exp))>100000):
                raise NotImplementedError('Fraction normalization exceeds the numeric power resource limit')
            if not getattr(node,'is_Mul',False):product_sizes[node]=1
            elif visited:product_sizes[node]=min(10001,sum(product_sizes[arg] for arg in node.args))
            else:
                pending.append((node,True))
                pending.extend((arg,False) for arg in node.args if arg not in product_sizes)
        return product_sizes[value]

    def denominator_possible(value):
        if value.is_Rational:return value.q!=1
        if value.is_Pow:exponent=value.exp
        elif getattr(value.func,'__name__','')=='exp':exponent=value.args[0]
        else:return False
        return exponent.is_negative is True or exponent.is_positive is not True and exponent.is_Mul

    def negative(value):
        if isinstance(value,MatrixExpr):
            return value.is_MatMul and value.as_coeff_matrices()[0].could_extract_minus_sign()
        return isinstance(value,Expr) and value.as_coeff_Mul()[0].is_negative is True

    def instruction(value, mode):
        counts=Counter()
        if isinstance(value,dict):
            return counts,[(item,'top') for pair in value.items() for item in pair]
        if isinstance(value,(tuple,list,set,frozenset)) or getattr(value,'_symox_blockgrid',False):
            return counts,[(item,mode) for item in value]
        if getattr(value,'is_NDimArray',False):
            return counts,[(item,'top') for _,item in value._items()]
        if getattr(value,'_symox_sparse',False):
            return counts,[(item,'top') for item in value.todok().values()]
        if isinstance(value,Matrix):
            return counts,[(item,'top') for item in value]
        if isinstance(value,Set):value=value._e
        if not isinstance(value,(Expr,MatrixExpr)):
            try:value=sympify(value)
            except (TypeError,ValueError):return counts,[]
            if not isinstance(value,(Expr,MatrixExpr)):return counts,[]
        name=getattr(value.func,'__name__',str(value.func))
        if mode=='top':
            if name in relations or name in booleans:
                counts[relations.get(name,name.upper())]+=1
                return counts,[(arg,'top') for arg in value.args]
            if name in ('Tuple','FiniteSet'):
                return counts,[(arg,'top') for arg in value.args]
            mode='expr'
        if getattr(value,'is_Rational',False):
            if value.p<0:counts['NEG']+=1
            if value.q!=1:counts['DIV']+=1
            return counts,[]
        mul=getattr(value,'is_Mul',False) or getattr(value,'is_MatMul',False)
        add=getattr(value,'is_Add',False) or getattr(value,'is_MatAdd',False)
        if mul:
            if negative(value):
                counts['NEG']+=1
                value=value.as_two_terms()[1] if value.args[0]==-1 else -value
            potential=not isinstance(value,MatrixExpr) and any(denominator_possible(arg) for arg in value.args)
            if potential:
                if multiplication_size(value)>10000:
                    raise NotImplementedError('Fraction normalization exceeds 10000 multiplicative factors')
                numerator,denominator=fraction(value)
            else:numerator,denominator=value,S.One
            if getattr(numerator,'is_Integer',False):
                counts['DIV']+=1
                if numerator.is_negative:counts['NEG']+=1
                return counts,[(denominator,'expr')]
            if denominator!=1:
                counts['DIV']+=1
                children=[(numerator,'expr')]
                if not denominator.is_Integer:children.append((denominator,'expr'))
                return counts,children
        elif add:
            signs=[negative(term) for term in value.args]
            children=[]
            for index,(term,minus) in enumerate(zip(value.args,signs)):
                children.append((-term if minus else term,'expr'))
                if index:counts['SUB' if minus else 'ADD']+=1
            if all(signs):counts['NEG']+=1
            elif signs[0]:counts['SUB']+=1;counts['ADD']-=1
            return counts,children
        if getattr(value,'is_Pow',False) and value.exp==-1:
            counts['DIV']+=1
            return counts,[(value.base,'expr')]
        if value==S.Exp1:
            counts['EXP']+=1
            return counts,[]
        if getattr(value,'is_Pow',False) and value.base==S.Exp1:
            counts['EXP']+=1
            return counts,[(value.exp,'expr')]
        mul=getattr(value,'is_Mul',False) or getattr(value,'is_MatMul',False)
        name=getattr(value.func,'__name__',str(value.func))
        if mul or name in ('And','Or','Max','Min'):
            counts[name.upper()]+=len(value.args)-1
        elif value.args and name!='Order' and (
                getattr(value,'is_Pow',False) or getattr(value,'is_Function',False)
                or name in ('Derivative','Integral','Sum','Piecewise')):
            counts[('FUNC_' if getattr(value,'is_UndefinedFunction',False) else '')+name.upper()]+=1
        children=[] if value.is_Symbol else [(child,'expr') for child in _children(value)]
        return counts,children

    # Count each graph instruction once, then propagate tree occurrence
    # multiplicities. Operation kinds are unchanged by alpha renaming.
    owners={}
    def key(value,mode):
        try:hash(value)
        except TypeError:
            identifier=(mode,type(value),id(value))
            owners[identifier]=value
            return identifier
        return mode,type(value),value
    root=key(expression,'top')
    instructions={}
    order=[]
    active=set()
    pending=[(expression,'top',False)]
    while pending:
        value,mode,visited=pending.pop()
        node=key(value,mode)
        if visited:
            active.remove(node)
            order.append(node)
            continue
        if node in active:raise ValueError('Cyclic operation-count input')
        if node in instructions:continue
        local,children=instruction(value,mode)
        instructions[node]=(local,[key(child,child_mode) for child,child_mode in children])
        if len(instructions)>100_000:
            raise NotImplementedError('Operation counting exceeds 100000 distinct nodes')
        active.add(node)
        pending.append((value,mode,True))
        pending.extend((child,child_mode,False) for child,child_mode in reversed(children))
    weights=defaultdict(int,{root:1})
    counts=Counter()
    for node in reversed(order):
        weight=weights[node]
        local,children=instructions[node]
        for name,number in local.items():counts[name]+=weight*number
        for child in children:weights[child]+=weight
    return Add(*(number*Symbol(name) for name,number in counts.items() if number)) if visual else sum(counts.values())
