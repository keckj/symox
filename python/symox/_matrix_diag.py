"""SymPy's ``MatrixBase.diag`` and ``diag``: block placement with the
``strict``, ``unpack``, ``rows``, ``cols`` and ``cls`` options."""


def _is_sequence(value):
    # sympy.utilities.iterables.is_sequence
    return hasattr(value, '__getitem__') and hasattr(value, '__iter__') and not isinstance(value, str)


def _loose_block(value):
    """``SparseMatrix._handle_creation_inputs`` of a list: a flat list is a
    column and a ragged list of rows is padded with zeros."""
    from . import Matrix
    if value and any(isinstance(row, list) for row in value):
        # A scalar among rows is a row of one entry, as the reference.
        value = [row if isinstance(row, list) else [row] for row in value]
        cols = max(len(row) for row in value)
        entries = {(i, j): entry for i, row in enumerate(value) for j, entry in enumerate(row)}
        return len(value), cols, entries
    block = Matrix(value)
    r, c = block.shape
    return r, c, {(i, j): block[i, j] for i in range(r) for j in range(c)}


def diag_entries(args, strict=False, unpack=True, rows=None, cols=None):
    from . import Matrix, sympify
    if unpack and len(args) == 1 and _is_sequence(args[0]) and not getattr(args[0], 'is_Matrix', False):
        args = args[0]
    entries = {}
    rmax = cmax = 0
    for m in args:
        if isinstance(m, list):
            if strict:
                block = Matrix(m)
                r, c = block.shape
                m = block.tolist()
            else:
                r, c, block = _loose_block(m)
                for (i, j), value in block.items():
                    entries[(i + rmax, j + cmax)] = value
                rmax += r
                cmax += c
                continue
        elif hasattr(m, 'shape'):
            r, c = m.shape
            # A matrix expression has no tolist in the reference, which
            # raises AttributeError; its explicit entries are used instead.
            m = m.tolist() if hasattr(m, 'tolist') else m.as_explicit().tolist()
        else:
            entries[(rmax, cmax)] = m
            rmax += 1
            cmax += 1
            continue
        for i, mi in enumerate(m):
            for j, value in enumerate(mi):
                entries[(i + rmax, j + cmax)] = value
        rmax += r
        cmax += c
    if rows is None:
        rows, cols = cols, rows
    if rows is None:
        rows, cols = rmax, cmax
    else:
        cols = rows if cols is None else cols
    for size in (rows, cols):
        value = sympify(size)
        if not getattr(value, 'is_Integer', False):
            raise ValueError('%s is not an integer' % (size,))
    rows, cols = int(rows), int(cols)
    if rows < rmax or cols < cmax:
        raise ValueError('The constructed matrix is {} x {} but a size of {} x {} '
                         'was specified.'.format(rmax, cmax, rows, cols))
    return rows, cols, entries


def _build(klass, rows, cols, entries):
    from . import Matrix, ImmutableMatrix, S, sympify
    flat = [sympify(entries[(i, j)]) if (i, j) in entries else S.Zero
            for i in range(rows) for j in range(cols)]
    result = Matrix(rows, cols, flat)
    if klass is Matrix or klass is None:
        return result
    if klass is ImmutableMatrix:
        return result.as_immutable()
    return klass(result)


def matrix_diag(kls, *args, strict=False, unpack=True, rows=None, cols=None, **kwargs):
    klass = kwargs.get('cls', kls)
    rows, cols, entries = diag_entries(args, strict, unpack, rows, cols)
    return _build(klass, rows, cols, entries)


def install(namespace):
    from . import _core
    Matrix = namespace['Matrix']

    def diag_method(*args, strict=False, unpack=True, rows=None, cols=None, **kwargs):
        """Returns a matrix with the provided values placed on the diagonal
        (SymPy's ``MatrixBase.diag``)."""
        return matrix_diag(Matrix, *args, strict=strict, unpack=unpack, rows=rows, cols=cols, **kwargs)

    # A static method, as the native one: the subclass support binds the
    # class and applies ``cls=``.
    _core.Matrix.diag = staticmethod(diag_method)

    def diag(*values, strict=True, unpack=False, **kwargs):
        """Returns a matrix with the provided values placed on the diagonal
        (SymPy's ``diag``: ``Matrix.diag`` with ``strict=True`` and
        ``unpack=False``)."""
        return Matrix.diag(*values, strict=strict, unpack=unpack, **kwargs)

    return diag
