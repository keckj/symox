"""Dimensional exponent arithmetic and exact linear unit conversion."""
import symox as S

_QUANTITIES = {}
_DIMENSIONS = {}
_PREFIXES = {}
_DIMENSION_DATA = {}


def _native(value):return S.sympify(value)

def _clean(values):return {str(k):S.sympify(v) for k,v in values.items() if v != 0}

def _combine(first,second,multiplier=1):
    result = dict(first)
    for name,power in second.items():result[name] = S.simplify(result.get(name,S.Integer(0))+power*multiplier)
    return _clean(result)


def _dimensional_name(dependencies):return S.Mul(*(S.Symbol(n)**p for n,p in sorted(dependencies.items())))


# Named dimensions with a symbol, such as Dimension(length, L), by name.
_SYMBOLIC_DIMENSIONS = {}


class Dimension:
    is_commutative = True
    def __init__(self,name,symbol=None):
        if isinstance(name,(list,tuple,dict)):raise TypeError('dimension name must be a scalar expression')
        if symbol is not None and not isinstance(symbol,(str,S.Expr)):raise TypeError('dimension symbol must be a string or Symbol')
        self.name = S.Symbol(name) if isinstance(name,str) else _native(name)
        self.symbol = S.Symbol(symbol) if isinstance(symbol,str) else symbol
        self._atom = S._core._apply('Dimension',self.name,evaluate=False)
        _DIMENSIONS[self._atom] = self
        if self.symbol is not None:_SYMBOLIC_DIMENSIONS.setdefault(str(self.name),self)
    def _symox_(self):return self._atom
    def __repr__(self):return 'Dimension(%s%s)'%(self.name,'' if self.symbol is None else ', '+str(self.symbol))
    def __hash__(self):return hash(self.name)
    def __eq__(self,other):return isinstance(other,Dimension) and self.name == other.name
    @property
    def free_symbols(self):return self.name.free_symbols
    @property
    def args(self):return (self.name,) if self.symbol is None else (self.name,self.symbol)
    def __mul__(self,other):
        if isinstance(other,Quantity):raise TypeError('dimensions and quantities cannot be combined')
        if isinstance(other,Dimension):return Dimension(self.name*other.name)
        value = _native(other)
        return self if value.is_number else self._atom*value
    __rmul__ = __mul__
    def __truediv__(self,other):return self*(other**-1 if isinstance(other,Dimension) else 1/_native(other))
    def __rtruediv__(self,other):return self**-1 if _native(other).is_number else _native(other)*self._atom**-1
    def __pow__(self,power):
        if isinstance(power,Dimension):
            if DimensionSystem.get_default().get_dimensional_dependencies(power):raise TypeError('dimension exponent must be dimensionless')
            power = S.Integer(1)
        power = S.Symbol(power) if isinstance(power,str) else _native(power)
        return Dimension(self.name**power)
    def __neg__(self):return self
    def __add__(self,other):
        if isinstance(other,Quantity):raise TypeError('dimensions and quantities cannot be added')
        if isinstance(other,Dimension) and self == other:return self
        return self._atom+_native(other)
    __radd__ = __add__
    def __sub__(self,other):return self+other if isinstance(other,Dimension) else self._atom-_native(other)
    def __rsub__(self,other):return _native(other)+self._atom
    def has_integer_powers(self,dim_sys):return all(p.is_integer is True for p in dim_sys.get_dimensional_dependencies(self).values())


class DimensionSystem:
    _default = None
    def __init__(self,base_dims,derived_dims=(),dimensional_dependencies=None):
        self.base_dims = tuple(v if isinstance(v,Dimension) else Dimension(v) for v in base_dims)
        self.derived_dims = tuple(v if isinstance(v,Dimension) else Dimension(v) for v in derived_dims)
        self._dependencies = {str(k.name if isinstance(k,Dimension) else k):_clean({str(n.name if isinstance(n,Dimension) else n):p for n,p in v.items()}) for k,v in (dimensional_dependencies or {}).items()}
        for d in self.base_dims:self._dependencies.setdefault(str(d.name),{str(d.name):S.Integer(1)})
    @classmethod
    def get_default(cls):
        if cls._default is None:cls._default = cls(('length','mass','time','current','temperature','amount_of_substance','luminous_intensity'),dimensional_dependencies=_DIMENSION_DATA)
        return cls._default
    def _expand(self,name,trail=()):
        key = str(name)
        if key in trail:raise ValueError('circular dimension dependency')
        data = self._dependencies.get(key,_DIMENSION_DATA.get(key))
        if data is None:return {key:S.Integer(1)}
        result = {}
        for n,p in data.items():
            if n == key:result = _combine(result,{n:p})
            else:result = _combine(result,self._expand(n,trail+(key,)),p)
        return result
    def _name_dependencies(self,value):
        value = _native(value)
        if value.is_number:return {}
        if value.is_Symbol:return self._expand(value)
        if value.is_Mul:
            result = {}
            for arg in value.args:result = _combine(result,self._name_dependencies(arg))
            return result
        if value.is_Pow:
            base,power = value.args
            return _clean({n:p*power for n,p in self._name_dependencies(base).items()})
        if value.is_Add:
            dimensions = [self._name_dependencies(a) for a in value.args]
            if any(d != dimensions[0] for d in dimensions):raise TypeError('cannot add incompatible dimensions')
            return dimensions[0]
        return _analyze(value,UnitSystem.get_default_unit_system())[1]
    def get_dimensional_dependencies(self,dimension,mark_dimensionless=False):
        if isinstance(dimension,Dimension):result = self._name_dependencies(dimension.name)
        else:
            value = S.Symbol(dimension) if isinstance(dimension,str) else _native(dimension)
            if any(symbol in _QUANTITIES for symbol in value.free_symbols):
                result = _analyze(value,UnitSystem.get_default_unit_system(),dimension_system=self)[1]
            else:
                # SymPy reads any other argument as a dimension name
                # expression, such as get_dimensional_expr's result.
                result = self._name_dependencies(value)
        # Reference keys are the system's base dimensions, with their symbols,
        # except for a tabled derived dimension, whose stored keys carry none.
        name = dimension.name if isinstance(dimension,Dimension) else None
        derived = name is not None and name.is_Symbol and str(name) not in {str(d.name) for d in self.base_dims} and (str(name) in self._dependencies or str(name) in _DIMENSION_DATA)
        bases = {} if derived else {str(d.name):d for d in self.base_dims}
        def key(n):
            base = bases.get(str(n))
            if base is not None and base.symbol is None:base = _SYMBOLIC_DIMENSIONS.get(str(n),base)
            return base or Dimension(n)
        return {key(n):p for n,p in result.items()} if result else ({Dimension(1):S.Integer(1)} if mark_dimensionless else {})
    def equivalent_dims(self,first,second):return self.get_dimensional_dependencies(first)==self.get_dimensional_dependencies(second)
    def is_dimensionless(self,dimension):return not self.get_dimensional_dependencies(dimension)
    def extend(self,new_base_dims, new_derived_dims=(),new_dim_deps=None):
        dependencies = dict(self._dependencies);dependencies.update(new_dim_deps or {})
        return DimensionSystem(self.base_dims+tuple(new_base_dims),self.derived_dims+tuple(new_derived_dims),dependencies)
    def can_transf_matrix(self):
        names = sorted(str(d.name) for d in self.base_dims)
        return S.Matrix([[self._name_dependencies(d.name).get(n,0) for d in self.base_dims] for n in names])
    @property
    def is_consistent(self):return self.can_transf_matrix().rank() == len(self.base_dims)
    def __repr__(self):return 'DimensionSystem(%s)'%(self.base_dims,)


class _QuantityType(type):
    # Inside expressions a quantity is represented by its registered atom;
    # SymPy's args are Quantity instances, so isinstance, has(Quantity) and
    # atoms(Quantity) recognise the atom as the quantity it stands for.
    def __instancecheck__(cls, value):
        if type.__instancecheck__(cls, value):return True
        if not getattr(value, 'is_Symbol', False):return False
        try:quantity = _QUANTITIES.get(value)
        except TypeError:return False
        return quantity is not None and type.__instancecheck__(cls, quantity)


class Quantity(metaclass=_QuantityType):
    is_commutative = True
    is_positive = True
    is_real = True
    def __init__(self,name,abbrev=None,latex_repr=None,pretty_unicode_repr=None,pretty_ascii_repr=None,mathml_presentation_repr=None,is_prefixed=False,**assumptions):
        self.name = S.Symbol(str(name));self.abbrev = S.Symbol(str(abbrev if abbrev is not None else name))
        self._atom = S.Symbol(str(name),positive=True)
        self._dimension = Dimension(str(name));self._scale_factor = S.Integer(1)
        self.is_prefixed = is_prefixed;self.is_physical_constant = False
        self._latex_repr = latex_repr
        _QUANTITIES[self._atom] = self
    def _symox_(self):return self._atom
    def __repr__(self):
        # SymPy prints a quantity through StrPrinter, whose ``abbrev``
        # setting selects the abbreviation.
        from ...printing.str import abbreviate_units
        return str(self.abbrev if abbreviate_units() else self.name)
    def __hash__(self):return hash(self._atom)
    def __eq__(self,other):return self._atom == getattr(other,'_atom',other)
    @property
    def args(self):return self.name,self.abbrev
    @property
    def free_symbols(self):return {self._atom}
    def has(self, *patterns):
        # Quantity is a named atom with public name/abbreviation arguments.
        # Its dimensional scale is metadata, not a structural child.
        for pattern in patterns:
            if isinstance(pattern,type) and isinstance(self,pattern):return True
            # A plain Symbol equal to the atom is not this quantity.
            if type.__instancecheck__(Quantity,pattern) and self == pattern:return True
            if any(argument.has(pattern) for argument in self.args):return True
        return False
    def as_ordered_factors(self, order=None):return [self]
    def replace(self, query, value, map=False, simultaneous=True, exact=None):
        result = self._atom.replace(query,value,map=map,simultaneous=simultaneous,exact=exact)
        if map:
            expression,changes = result
            return (self if expression == self._atom else expression,changes)
        return self if result == self._atom else result
    @property
    def dimension(self):return self._dimension
    @property
    def scale_factor(self):return self._scale_factor
    def set_global_dimension(self,dimension):self._dimension = dimension if isinstance(dimension,Dimension) else Dimension(dimension)
    def set_global_relative_scale_factor(self,scale_factor,reference_quantity):
        factor,dimension = _analyze(_native(reference_quantity),UnitSystem.get_default_unit_system())
        self._scale_factor = _native(scale_factor)*factor
        self._dimension = Dimension(_dimensional_name(dimension))
    def convert_to(self,other,unit_system='SI'):return convert_to(self,other,unit_system)
    def __add__(self,other):
        if isinstance(other,Dimension):raise TypeError('quantities and dimensions cannot be added')
        return self._atom+other
    __radd__ = __add__
    def __neg__(self):return -self._atom
    def __sub__(self,other):return self._atom-other
    def __rsub__(self,other):return other-self._atom
    def __mul__(self,other):
        if isinstance(other,Dimension):raise TypeError('quantities and dimensions cannot be multiplied')
        return self._atom*other
    __rmul__ = __mul__
    def __truediv__(self,other):return self._atom/other
    def __rtruediv__(self,other):return other/self._atom
    def __pow__(self,power):return self._atom**power
    def __abs__(self):return self._atom


class PhysicalConstant(Quantity):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.is_physical_constant = True


class Prefix:
    def __init__(self,name,abbrev,exponent,base=10,latex_repr=None):
        self.name,self.abbrev = S.Symbol(str(name)),S.Symbol(str(abbrev));self.exponent,self.base = _native(exponent),_native(base)
        self.scale_factor = self.base**self.exponent;self._atom = S.Symbol(str(name),positive=True);self._latex_repr = latex_repr
        _PREFIXES[self._atom] = self
    @property
    def args(self):return self.name,self.abbrev,self.exponent,self.base
    def _symox_(self):return self._atom
    def __repr__(self):return "Prefix(%r, %r, %s, %s)"%(str(self.name),str(self.abbrev),self.exponent,self.base)
    def __hash__(self):return hash(self._atom)
    def __eq__(self,other):return isinstance(other,Prefix) and self.args==other.args
    def _match(self,value):
        if value == 1:return S.Integer(1)
        matches = [p for p in _PREFIXES.values() if p.base == self.base and p.scale_factor == value]
        return matches[0] if matches else value
    def __mul__(self,other):
        if isinstance(other,Prefix):return self._match(self.scale_factor*other.scale_factor)
        if isinstance(other,Quantity):return self.scale_factor*_native(other)
        return self._atom*other
    __rmul__ = __mul__
    def __truediv__(self,other):return self._match(self.scale_factor/other.scale_factor) if isinstance(other,Prefix) else self._atom/other
    def __rtruediv__(self,other):return self._match(1/self.scale_factor) if other == 1 else other/self._atom
    def __pow__(self,power):return self._atom**power


class UnitSystem:
    _registry = {}
    _default = None
    def __init__(self,base_units,units=(),name='',descr='',dimension_system=None,derived_units=None):
        self._base_units = tuple(base_units);self._units = tuple(dict.fromkeys(self._base_units+tuple(units)))
        self.name,self.descr = name,descr;self._dimension_system = dimension_system or DimensionSystem.get_default()
        self.derived_units = dict(derived_units or {});self._quantity_dimensions = {};self._quantity_scale_factors = {}
        if name:self._registry[name] = self
    @classmethod
    def get_unit_system(cls,unit_system):
        if isinstance(unit_system,cls):return unit_system
        if unit_system not in cls._registry:raise ValueError('unknown unit system: '+str(unit_system))
        return cls._registry[unit_system]
    @classmethod
    def get_default_unit_system(cls):return cls._default
    def get_dimension_system(self):return self._dimension_system
    @property
    def dim(self):return len(self._base_units)
    @property
    def is_consistent(self):return self._dimension_system.is_consistent
    def get_quantity_dimension(self,quantity):return self._quantity_dimensions.get(_native(quantity),quantity.dimension)
    def get_quantity_scale_factor(self,quantity):return self._quantity_scale_factors.get(_native(quantity),quantity.scale_factor)
    def set_quantity_dimension(self,quantity,dimension):self._quantity_dimensions[_native(quantity)] = dimension if isinstance(dimension,Dimension) else Dimension(dimension)
    def set_quantity_scale_factor(self,quantity,scale_factor):self._quantity_scale_factors[_native(quantity)] = _analyze(_native(scale_factor),self)[0]
    def get_dimensional_expr(self,expr):return _dimensional_name(_analyze(_native(expr),self)[1])
    def _collect_factor_and_dimension(self,expr):
        factor,dimensions = _analyze(_native(expr),self)
        return factor,Dimension(_dimensional_name(dimensions))
    def extend(self,base,units=(),name='',descr='',dimension_system=None,derived_units=None):
        return UnitSystem(self._base_units+tuple(base),self._units+tuple(units),name,descr,dimension_system or self._dimension_system,derived_units or self.derived_units)
    def __repr__(self):return self.name or 'UnitSystem(%s)'%(self._base_units,)


def _analyze(expr,system,dimension_system=None):
    dims = dimension_system or system.get_dimension_system()
    if expr in _QUANTITIES:
        q = _QUANTITIES[expr]
        return system.get_quantity_scale_factor(q),dims._name_dependencies(system.get_quantity_dimension(q).name)
    if expr in _PREFIXES:return _PREFIXES[expr].scale_factor,{}
    if expr in _DIMENSIONS:return S.Integer(1),dims._name_dependencies(_DIMENSIONS[expr].name)
    if expr.is_Mul:
        factor,dimension = S.Integer(1),{}
        for arg in expr.args:
            f,d = _analyze(arg,system,dims);factor *= f;dimension = _combine(dimension,d)
        return factor,dimension
    if expr.is_Add:
        values = [_analyze(arg,system,dims) for arg in expr.args]
        if any(d != values[0][1] for _,d in values):raise ValueError('addition combines incompatible dimensions')
        return S.Add(*(f for f,_ in values)),values[0][1]
    if expr.is_Pow:
        base,power = expr.args;f,d = _analyze(base,system,dims);p,pdim = _analyze(power,system,dims)
        if pdim:raise ValueError('a power must be dimensionless')
        return f**p,_clean({n:e*p for n,e in d.items()})
    if expr.is_Atom:return expr,{}
    if expr.class_name == 'Abs':
        factor,dimension = _analyze(expr.args[0],system,dims);return S.Abs(factor),dimension
    if expr.class_name == 'atan2':
        a,b = [_analyze(v,system,dims) for v in expr.args]
        if a[1] != b[1]:raise ValueError('atan2 arguments must have matching dimensions')
        return S.atan2(a[0],b[0]),{}
    arguments = [_analyze(v,system,dims) for v in expr.args if isinstance(v,S.Expr)]
    if any(d for _,d in arguments):raise ValueError('function arguments must be dimensionless')
    return expr.xreplace({a:f for a,(f,d) in zip(expr.args,arguments)}),{}


def _target_powers(source,target_dimensions):
    if not source:
        for i,d in enumerate(target_dimensions):
            if not d:return [S.Integer(int(i==j)) for j in range(len(target_dimensions))]
    names = sorted(set(source).union(*(set(d) for d in target_dimensions)))
    if not names:return [S.Integer(0)]*len(target_dimensions)
    A = S.Matrix([[d.get(n,0) for d in target_dimensions] for n in names]);b = S.Matrix([source.get(n,0) for n in names])
    reduced,pivots = A.row_join(b).rref()
    if len(target_dimensions) in pivots:return None
    result = [S.Integer(0)]*len(target_dimensions)
    for row,column in enumerate(pivots):result[column] = reduced[row,-1]
    return result


def convert_to(expr,target_units,unit_system='SI'):
    expr = _native(expr);system = UnitSystem.get_unit_system(unit_system)
    if getattr(target_units,'is_Tuple',False):target_units = target_units.args
    targets = [_native(v) for v in target_units] if isinstance(target_units,(list,tuple)) else [_native(target_units)]
    if not targets:raise ValueError('at least one target unit is required')
    if expr.is_Add:return S.Add(*(convert_to(a,targets,system) for a in expr.args))
    try:
        factor,dimension = _analyze(expr,system)
        target_values = [_analyze(t,system) for t in targets]
        powers = _target_powers(dimension,[d for _,d in target_values])
    except (ValueError,TypeError):return expr
    if powers is None:return expr
    if not any(p != 0 for p in powers) and any(a in _QUANTITIES for a in expr.free_symbols):return expr
    result = factor
    for target,(scale,d),power in zip(targets,target_values,powers):result *= (target/scale)**power
    return S.simplify(result)


def quantity_simplify(expr,across_dimensions=False,unit_system=None):
    expr = _native(expr)
    if across_dimensions:
        if unit_system is None:raise ValueError('across-dimension simplification requires a unit system')
        system = UnitSystem.get_unit_system(unit_system)
        try:dimension = _analyze(expr,system)[1]
        except ValueError:return expr
        for d,q in system.derived_units.items():
            if system.get_dimension_system()._name_dependencies(d.name) == dimension:return convert_to(expr,q,system)
        return convert_to(expr,system._base_units,system)
    system = UnitSystem.get_unit_system(unit_system or 'SI');choices = {};replacements = {}
    for atom in sorted(expr.free_symbols,key=str):
        if atom in _PREFIXES:replacements[atom] = _PREFIXES[atom].scale_factor
        if atom not in _QUANTITIES:continue
        quantity = _QUANTITIES[atom];dimension = tuple(sorted(_analyze(atom,system)[1].items()))
        if dimension in choices:
            target = choices[dimension];replacements[atom] = system.get_quantity_scale_factor(quantity)/system.get_quantity_scale_factor(target)*_native(target)
        else:choices[dimension] = quantity
    return S.simplify(expr.xreplace(replacements))


def check_dimensions(expr,unit_system='SI'):
    expr = _native(expr);system = UnitSystem.get_unit_system(unit_system);stack = [expr]
    while stack:
        node = stack.pop();stack.extend(node.args)
        if node.is_Add:
            dimensions = []
            for part in node.args:
                if part.free_symbols and not any(s in _QUANTITIES or s in _DIMENSIONS for s in part.free_symbols):continue
                dimensions.append(_analyze(part,system)[1])
            if dimensions and any(d != dimensions[0] for d in dimensions):raise ValueError('addition combines incompatible dimensions')
    return expr


def prefix_unit(unit,prefixes):
    result = []
    for prefix in prefixes.values():
        quantity = Quantity(str(prefix.name)+str(unit.name),str(prefix.abbrev)+str(unit.abbrev),is_prefixed=True)
        quantity.set_global_relative_scale_factor(prefix.scale_factor,unit);result.append(quantity)
    return result


def _abbreviated(expr):
    """``expr`` with quantity atoms renamed to their abbreviations.

    An abbreviation shared by two quantities of ``expr`` or equal to the name
    of another free symbol keeps the full name, since a renamed atom would
    then combine with a different one."""
    try:
        symbols = expr.free_symbols
    except Exception:
        return expr
    quantities = [symbol for symbol in symbols if symbol in _QUANTITIES]
    if not quantities:
        return expr
    abbreviations = [str(_QUANTITIES[symbol].abbrev) for symbol in quantities]
    names = {getattr(symbol, "name", None) for symbol in symbols}
    rule = {}
    for symbol, abbreviation in zip(quantities, abbreviations):
        if abbreviation == symbol.name or abbreviations.count(abbreviation) > 1 or abbreviation in names:
            continue
        rule[symbol] = S.Symbol(abbreviation, positive=True)
    return expr.xreplace(rule) if rule else expr


def _install_abbreviated_printing():
    # Quantities are native positive symbols; Basic.__str__ prints them
    # through StrPrinter, so the ``abbrev`` setting reaches every expression.
    from functools import wraps
    from ... import _core
    from ...printing.str import abbreviate_units
    if getattr(_core.Expr.__str__, '_symox_quantities', False):
        return
    native_str, native_repr = _core.Expr.__str__, _core.Expr.__repr__

    def wrap(native):
        @wraps(native)
        def method(self):
            if abbreviate_units():
                renamed = _abbreviated(self)
                if renamed is not self:
                    return native(renamed)
            return native(self)
        method._symox_quantities = True
        return method

    _core.Expr.__str__, _core.Expr.__repr__ = wrap(native_str), wrap(native_repr)


_install_abbreviated_printing()
