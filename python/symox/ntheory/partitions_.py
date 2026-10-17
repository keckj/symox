"""``npartitions`` (``sympy.ntheory.partitions_``)."""


def npartitions(n, verbose=False):
    """Number of partitions of ``n`` (the ``partition`` head)."""
    from .. import partition
    return partition(n)
