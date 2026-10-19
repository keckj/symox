"""Characteristic-polynomial coefficient callbacks over native domains."""


def install(matrix):
    native = matrix.charpoly

    def charpoly(self, x=None, simplify=None):
        if simplify is None:
            return native(self, x)
        from . import Add, Expr, PurePoly, construct_domain, sympify
        if not callable(simplify) or isinstance(simplify, Expr) and not simplify.is_Lambda:
            raise TypeError('simplify must be callable or None')
        # The native constructor handles generator identity and coefficient
        # collisions. DomainMatrix retains unsimplified EXRAW coefficients
        # until the caller has seen them, as required by this callback API.
        generator = native(self, x).gen
        domain_matrix = self.to_DM()
        coefficients = [sympify(simplify(domain_matrix.domain.to_sympy(value)))
                        for value in domain_matrix.charpoly()]
        if not all(isinstance(value,Expr) for value in coefficients):
            raise TypeError('simplify must return scalar expressions')
        while len(coefficients) > 1 and coefficients[0] == 0:
            coefficients.pop(0)
        domain, _ = construct_domain(coefficients)
        degree = len(coefficients)-1
        expression = Add(*(coefficient*generator**(degree-index)
                           for index, coefficient in enumerate(coefficients)))
        return PurePoly._from_coeffs(coefficients, generator, str(domain), expression)

    matrix.charpoly = charpoly
