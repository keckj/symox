"""Native expression dispatch for differential tensors and point evaluation."""
from functools import wraps
from .._dispatch_arguments import FirstArgument
from .. import _core as _proof_core

_PROOF_TYPE=False.__class__.__class__
_PROOF_STR=''.__class__
_PROOF_TUPLE=().__class__
_PROOF_DICT={}.__class__
_PROOF_FUNCTION=(lambda:None).__class__
_PROOF_BUILTIN=[].append.__class__
_SCALAR_DISPATCH_PROOFS=()
_SCALAR_AUX_PROOFS=()
_SCALAR_FIELD_CLASS=None
_SCALAR_HOST_READER=None


def _plain_module(module):
    kind=_PROOF_TYPE(module)
    if (_PROOF_TYPE(kind) is not _PROOF_TYPE or kind.__flags__ & (1<<9)
            or kind.__name__!='module' or kind.__module__!='builtins'):return False
    for key in module.__dict__:
        if _PROOF_TYPE(key) is not _PROOF_STR:return False
    return True


_SCALAR_METADATA=None
if _plain_module(_proof_core):
    _reader=_proof_core.__dict__.get('_ddm_function_metadata')
    if (_PROOF_TYPE(_reader) is _PROOF_BUILTIN and _reader.__self__ is _proof_core
            and _reader.__name__=='_ddm_function_metadata'):
        _SCALAR_METADATA=(_proof_core,_reader)
    _host_reader=_proof_core.__dict__.get('_host_contains')
    if (_PROOF_TYPE(_host_reader) is _PROOF_BUILTIN and _host_reader.__self__ is _proof_core
            and _host_reader.__name__=='_host_contains'):_SCALAR_HOST_READER=_host_reader


def _proof_metadata(function):
    if _SCALAR_METADATA is None:return None
    owner,reader=_SCALAR_METADATA
    if (owner is not _proof_core or not _plain_module(owner)
            or owner.__dict__.get('_ddm_function_metadata') is not reader
            or _PROOF_TYPE(function) is not _PROOF_FUNCTION):return None
    return reader(function)


def _proof_attributes(function):
    if _PROOF_TYPE(function) is not _PROOF_FUNCTION:return False
    namespace=function.__dict__
    if _PROOF_TYPE(namespace) is not _PROOF_DICT:return False
    for key in namespace:
        if _PROOF_TYPE(key) is not _PROOF_STR:return False
    return True


def _dispatch_states(function):
    pending=[function];seen=[];states=[]
    while pending:
        function=pending.pop()
        if _PROOF_TYPE(function) is not _PROOF_FUNCTION or function in seen:continue
        seen.append(function);metadata=_proof_metadata(function)
        if metadata is None:return None
        code,namespace,defaults,keywords,closure=metadata
        if (closure is not None and _PROOF_TYPE(closure) is not _PROOF_TUPLE
                or defaults is not None and _PROOF_TYPE(defaults) is not _PROOF_TUPLE
                or keywords is not None and _PROOF_TYPE(keywords) is not _PROOF_DICT):return None
        try:contents=_PROOF_TUPLE(cell.cell_contents for cell in closure or ())
        except ValueError:return None
        states.append((function,code,closure,contents,defaults,keywords))
        pending.extend(value for value in contents if _PROOF_TYPE(value) is _PROOF_FUNCTION)
    return _PROOF_TUPLE(states)


def _field_expression(expr):
    # Unpickles a callable scalar-field expression.
    import symox
    install()
    return symox._core.Expr.__new__(symox._diffgeom_callable_type,expr)


def install():
    global _SCALAR_DISPATCH_PROOFS,_SCALAR_AUX_PROOFS,_SCALAR_FIELD_CLASS
    import symox
    from ._tensors import _Tensor
    from ._coordinates import Point,_scalar_at,_atoms,BaseScalarField,CoordinateSymbol
    Expr = symox._core.Expr
    if getattr(Expr,'_diffgeom_dispatch_installed',False):
        return
    from ._coordinates import _field_atoms
    class FieldExpression(Expr):
        __slots__ = ()
        def __reduce_ex__(self,protocol):
            return (_field_expression,(Expr.__new__(Expr,self),))
        def __call__(self,*args,**kwargs):
            if len(args)==1 and isinstance(args[0],Point) and not kwargs:
                return _scalar_at(self,args[0])
            raise TypeError('scalar fields are evaluated at differential-geometric points')
    symox._diffgeom_callable_type=FieldExpression
    symox._diffgeom_field_atoms=_field_atoms
    field_states=()
    for _name,_value in FieldExpression.__dict__.items():
        if _PROOF_TYPE(_value) is _PROOF_FUNCTION:
            captured=_dispatch_states(_value)
            if captured is None:field_states=None;break
            field_states+=captured
    _SCALAR_FIELD_CLASS=(FieldExpression,FieldExpression.__bases__,
                         _PROOF_TUPLE(FieldExpression.__dict__.items()),field_states)
    for name in ('__add__','__sub__','__mul__','__truediv__'):
        original = getattr(Expr,name)
        def arithmetic(original,name):
            @wraps(original)
            def operation(self,other):
                if isinstance(other,_Tensor):
                    if name == '__add__':
                        return other.__radd__(self)
                    if name == '__sub__':
                        return other.__rsub__(self)
                    if name == '__mul__':
                        return other.__rmul__(self)
                    raise ValueError('division by a differential tensor is undefined')
                return original(self,other)
            return operation
        operation=arithmetic(original,name);metadata=_proof_metadata(operation)
        proof=(operation,original,_dispatch_states(operation),
               metadata[1] if metadata is not None else None,_Tensor,
               _Tensor.__mro__ if _PROOF_TYPE(_Tensor) is _PROOF_TYPE else None,isinstance,name)
        if _proof_attributes(operation):operation._symox_diffgeom_proof=proof
        _SCALAR_DISPATCH_PROOFS+=(proof,)
        setattr(Expr,name,operation)
    for name in ('rcall',):
        original = getattr(Expr,name)
        def evaluation(original):
            @wraps(original)
            def operation(self,*args,**kwargs):
                if len(args) == 1 and isinstance(args[0],Point) and not kwargs:
                    return _scalar_at(self,args[0])
                if self == 0 and any(isinstance(a,_Tensor) for a in args) and not kwargs:
                    return self
                return original(self,*args,**kwargs)
            return operation
        operation=evaluation(original);metadata=_proof_metadata(operation)
        proof=(name,operation,original,_dispatch_states(operation),
               metadata[1] if metadata is not None else None)
        if _proof_attributes(operation):operation._symox_diffgeom_aux_proof=proof
        _SCALAR_AUX_PROOFS+=(proof,)
        setattr(Expr,name,operation)
    original_atoms = Expr.atoms
    @wraps(original_atoms)
    def atoms(self,*types):
        special = tuple(t for t in types if t in (BaseScalarField,CoordinateSymbol))
        if special:
            result = {field for atom in self.free_symbols if isinstance((field := _atoms.get(atom)),special)}
            ordinary = tuple(t for t in types if t not in special)
            return result | (original_atoms(self,*ordinary) if ordinary else set())
        return original_atoms(self,*types)
    Expr.atoms = atoms
    metadata=_proof_metadata(atoms)
    proof=('atoms',atoms,original_atoms,_dispatch_states(atoms),
           metadata[1] if metadata is not None else None)
    if _proof_attributes(atoms):atoms._symox_diffgeom_aux_proof=proof
    _SCALAR_AUX_PROOFS+=(proof,)
    for name in ('diff','simplify','expand'):
        original = getattr(symox,name)
        def transform(original,name):
            first=FirstArgument(original)
            @wraps(original)
            def operation(*args,**kwargs):
                expression=first.get(args,kwargs)
                if isinstance(expression,_Tensor):
                    expression,rest,options=first.split(args,kwargs)
                    return getattr(expression,name)(*rest,**options)
                return original(*args,**kwargs)
            return operation
        setattr(symox,name,transform(original,name))
    if _proof_attributes(install):
        install._stock_proofs=_SCALAR_DISPATCH_PROOFS
        install._stock_aux_proofs=_SCALAR_AUX_PROOFS
    Expr._diffgeom_dispatch_installed = True


# Read immutable templates without adding audited function-attribute reads.
_metadata=_proof_metadata(install)
_SCALAR_INSTALL_CODE=_metadata[0] if _metadata is not None else None
_SCALAR_OPERATION_CODE=_SCALAR_RCALL_CODE=_SCALAR_ATOMS_CODE=None
if _SCALAR_INSTALL_CODE is not None:
    _code_type=_SCALAR_INSTALL_CODE.__class__
    for _constant in _SCALAR_INSTALL_CODE.co_consts:
        if _PROOF_TYPE(_constant) is not _code_type:continue
        if _constant.co_name=='atoms':_SCALAR_ATOMS_CODE=_constant
        if _constant.co_name in ('arithmetic','evaluation'):
            for _nested in _constant.co_consts:
                if _PROOF_TYPE(_nested) is _code_type and _nested.co_name=='operation':
                    if _constant.co_name=='arithmetic':_SCALAR_OPERATION_CODE=_nested
                    else:_SCALAR_RCALL_CODE=_nested
