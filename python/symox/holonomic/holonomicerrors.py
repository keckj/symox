""" Common Exceptions for `holonomic` module. """


class BaseHolonomicError(Exception):
    def new(self, *args):
        raise NotImplementedError("abstract base class")


class NotPowerSeriesError(BaseHolonomicError):
    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        return 'A Power Series does not exists for %s about %s.' % (self.holonomic, self.x0)


class NotHolonomicError(BaseHolonomicError):
    def __init__(self, m):
        self.m = m

    def __str__(self):
        return self.m


class SingularityError(BaseHolonomicError):
    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        return '%s has a singularity at %s.' % (self.holonomic, self.x0)


class NotHyperSeriesError(BaseHolonomicError):
    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        return 'Power series expansion of %s about %s is not hypergeometric' % (self.holonomic, self.x0)
