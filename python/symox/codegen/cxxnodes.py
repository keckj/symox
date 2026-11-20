"""
AST nodes specific to C++.
"""

from .ast import Attribute, String, Token, Type, none


class using(Token):
    """ Represents a 'using' statement in C++ """
    _fields = ('type', 'alias')
    defaults = {'alias': none}

    def _coerce(self, field, value):
        # SymPy's ``_construct_type = Type`` and ``_construct_alias = String``.
        if field == 'type' and not isinstance(value, Type):
            return Type(value)
        if field == 'alias' and value != none and not isinstance(value, String):
            return String(value)
        return super()._coerce(field, value)


constexpr = Attribute('constexpr')
