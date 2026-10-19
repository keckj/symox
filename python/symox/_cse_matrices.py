"""Flatten concrete matrix entries for CSE without densifying sparse inputs."""


def is_concrete_matrix(value):
    from . import _core
    return isinstance(value,_core.Matrix) or getattr(value,'_symox_sparse',False)


def flatten(expressions):
    from . import Matrix
    values=[];plans=[]
    for expression in expressions:
        if getattr(expression,'_symox_sparse',False):
            entries=expression.todok()
            keys=tuple(entries)
            values.extend(entries.values())
            def rebuild(items,template=expression,keys=keys):
                return template._new(template.rows,template.cols,dict(zip(keys,items)))
            plans.append((len(keys),rebuild))
        elif is_concrete_matrix(expression):
            values.extend(expression)
            def rebuild(items,template=expression):
                matrix=Matrix(template.rows,template.cols,items)
                return matrix.as_immutable() if template._immutable else matrix
            plans.append((expression.rows*expression.cols,rebuild))
        else:
            values.append(expression)
            plans.append((1,lambda items:items[0]))
    def restore(reduced):
        result=[];start=0
        for count,rebuild in plans:
            result.append(rebuild(reduced[start:start+count]))
            start+=count
        return result
    return values,restore
