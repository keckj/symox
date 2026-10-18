"""Physical dimensions, quantities, prefixes and exact unit conversions."""
import symox as _S
from ._core import (Dimension,DimensionSystem,Quantity,PhysicalConstant,Prefix,UnitSystem,
                    convert_to,quantity_simplify,check_dimensions,prefix_unit,_DIMENSION_DATA)
from ._catalog import QUANTITIES,DIMENSIONS,PREFIXES_DATA,ALIASES,PUBLIC_NAMES

_DIMENSION_DATA.update({name:{n:_S.sympify(e) for n,e in data[1].items()} for name,data in DIMENSIONS.items()})
for _,_,deps,name,_ in QUANTITIES.values():_DIMENSION_DATA.setdefault(name,{n:_S.sympify(e) for n,e in deps.items()})
_DIMENSION_DATA['angle'] = {}
_dimensions = {name:Dimension(name,symbol) for name,(symbol,deps) in DIMENSIONS.items()}
_quantities = {}
for _name,(_abbrev,_factor,_deps,_dimension,_constant) in QUANTITIES.items():
    _q = (PhysicalConstant if _constant else Quantity)(_name,_abbrev)
    _q._scale_factor = _S.sympify(_factor)
    _q._dimension = _dimensions.get(_dimension) or Dimension(_S.sympify(_dimension))
    _quantities[_name] = _q
_prefixes = {name:Prefix(name,abbrev,_S.sympify(exponent),base) for name,(abbrev,exponent,base) in PREFIXES_DATA.items()}
for _alias,(_kind,_key) in ALIASES.items():globals()[_alias] = {'quantity':_quantities,'dimension':_dimensions,'prefix':_prefixes}[_kind][_key]
Unit = Quantity
dHg0 = 13.5951

_SI_BASES = (meter,kilogram,second,ampere,kelvin,mole,candela)
_DERIVED = {length:meter,mass:kilogram,time:second,current:ampere,temperature:kelvin,amount_of_substance:mole,luminous_intensity:candela,
 force:newton,energy:joule,power:watt,pressure:pascal,frequency:hertz,charge:coulomb,voltage:volt,impedance:ohm,
 conductance:siemens,capacitance:farad,inductance:henry,magnetic_density:tesla,magnetic_flux:weber,velocity:meter/second,acceleration:meter/second**2}
SI = UnitSystem(_SI_BASES,tuple(_quantities.values()),'SI',dimension_system=DimensionSystem.get_default(),derived_units=_DERIVED)
UnitSystem._default = SI
MKS = UnitSystem((meter,kilogram,second),tuple(_quantities.values()),'MKS',dimension_system=DimensionSystem(('length','mass','time'),dimensional_dependencies=_DIMENSION_DATA),derived_units=_DERIVED)
MKSA = UnitSystem((meter,kilogram,second,ampere),tuple(_quantities.values()),'MKSA',dimension_system=DimensionSystem(('length','mass','time','current'),dimensional_dependencies=_DIMENSION_DATA),derived_units=_DERIVED)

from . import prefixes
prefixes.PREFIXES = {str(v.abbrev):v for v in _prefixes.values() if v.base == 10}
prefixes.BIN_PREFIXES = {str(v.abbrev):v for v in _prefixes.values() if v.base == 2}
for _name,_prefix in _prefixes.items():setattr(prefixes,_name,_prefix)
from .systems import mks,mksa,si


def find_unit(quantity,unit_system='SI'):
    system = UnitSystem.get_unit_system(unit_system)
    if isinstance(quantity,str):return sorted(n for n in ALIASES if ALIASES[n][0]=='quantity' and (quantity in n or quantity in str(globals()[n])))
    target = quantity if isinstance(quantity,Dimension) else system._collect_factor_and_dimension(quantity)[1]
    return sorted(n for n in ALIASES if ALIASES[n][0]=='quantity' and system.get_dimension_system().equivalent_dims(globals()[n].dimension,target))


__all__ = list(PUBLIC_NAMES)
