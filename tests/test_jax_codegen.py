"""Execute standalone JAX code, including JIT and differentiated branches."""

import inspect
import os
import subprocess
import sys

import numpy as np
import pytest
import symox as sm

jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")


@pytest.mark.parametrize("jit",[False,True])
def test_function_matrix_native_lambda_batching_jit_and_gradients(jit):
    i,j=sm.symbols('i j',integer=True)
    a,b=sm.symbols('a b')
    matrix=sm.FunctionMatrix(2,3,sm.Lambda((i,j),a*(i+1)+b*j*j+a*b))
    outputs=(matrix,matrix[::-1,::2],matrix.diff(a),sm.MatrixElement(matrix,0,2))
    function=sm.jax_lambdify((a,b),outputs,jit=jit)
    av=jnp.asarray([[1.],[2.]])
    bv=jnp.asarray([[2.,3.,4.]])

    def reference(x,y):
        rows=jnp.arange(2).reshape(2,1)
        cols=jnp.arange(3).reshape(1,3)
        dense=x[...,None,None]*(rows+1)+y[...,None,None]*cols**2+(x*y)[...,None,None]
        derivative=jnp.broadcast_to((rows+1)+y[...,None,None],dense.shape)
        return dense,dense[...,::-1,::2],derivative,dense[...,0,2]

    for value,target in zip(function(av,bv),reference(av,bv)):
        np.testing.assert_allclose(value,target,rtol=2e-6)
        assert value.dtype==av.dtype
    gradients=jax.grad(lambda x,y:sum(jnp.sum(v**2) for v in function(x,y)),argnums=(0,1))(av,bv)
    references=jax.grad(lambda x,y:sum(jnp.sum(v**2) for v in reference(x,y)),argnums=(0,1))(av,bv)
    for value,target in zip(gradients,references):
        np.testing.assert_allclose(value,target,rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_matrix_applyfunc_free_parameter_total_derivative_and_gradients(jit):
    matrix = sm.MatrixSymbol('A', 2, 3)
    z, t = sm.symbols('z t', real=True)
    applied = (t*matrix).applyfunc(sm.Lambda(z, z*z+t*z+t*t))
    outputs = (applied, applied.diff(t), applied.T,
               applied[::-1, ::2], sm.MatrixElement(applied, 0, 2))
    function = sm.jax_lambdify((matrix, t), outputs, jit=jit)
    av = jnp.linspace(-.8, 1.2, 36).reshape(2, 3, 2, 3)
    tv = jnp.asarray([[.3], [.8]])

    def reference(a, t):
        scalar = t[..., None, None]
        parent = scalar*a
        dense = parent**2+scalar*parent+scalar**2
        derivative = (2*parent+scalar)*a+parent+2*scalar
        return dense, derivative, jnp.swapaxes(dense, -1, -2), dense[..., ::-1, ::2], dense[..., 0, 2]

    for value, target in zip(function(av, tv), reference(av, tv)):
        np.testing.assert_allclose(value, target, rtol=2e-6, atol=1e-6)
        assert value.dtype == av.dtype and value.device == av.device
    gradients = jax.grad(lambda a, t: sum(jnp.sum(v**2) for v in function(a, t)), argnums=(0, 1))(av, tv)
    references = jax.grad(lambda a, t: sum(jnp.sum(v**2) for v in reference(a, t)), argnums=(0, 1))(av, tv)
    for value, target in zip(gradients, references):
        np.testing.assert_allclose(value, target, rtol=2e-6, atol=1e-6)


def test_distinct_dummy_inputs_jit_and_gradients():
    first, second = sm.Dummy('x'), sm.Dummy('x')
    function = sm.jax_lambdify((first, second), first**2 + 2*second, jit=True)
    a, b = jnp.asarray([1., 2.]), jnp.asarray([3., 4.])
    np.testing.assert_allclose(function(a, b), a**2+2*b)
    da, db = jax.grad(lambda x, y: jnp.sum(function(x, y)), argnums=(0, 1))(a, b)
    np.testing.assert_allclose(da, 2*a)
    np.testing.assert_allclose(db, jnp.full_like(b, 2))


@pytest.mark.parametrize("jit", [False, True])
def test_strided_reverse_composed_and_empty_matrix_slices(jit):
    matrix = sm.MatrixSymbol('A', 5, 6)
    outputs = (matrix[::-1, ::-2], matrix[::2, 1::2], matrix[4:1:-2, 5:0:-2],
               matrix[::2, ::-1][::-1, 1::2], matrix[1:1, ::-1], matrix[::-1, 10:20])
    function = sm.jax_lambdify(matrix, outputs, jit=jit)
    a = jnp.arange(180., dtype=jnp.float32).reshape(2, 3, 5, 6)

    def reference(value):
        return (value[..., ::-1, ::-2], value[..., ::2, 1::2], value[..., 4:1:-2, 5:0:-2],
                value[..., ::2, ::-1][..., ::-1, 1::2],
                value[..., 1:1, ::-1], value[..., ::-1, 10:20])

    for value, target in zip(function(a), reference(a)):
        assert value.dtype == a.dtype
        np.testing.assert_allclose(value, target, rtol=2e-6)
    actual = jax.grad(lambda value: sum(jnp.sum(v**2) for v in function(value)))(a)
    expected = jax.grad(lambda value: sum(jnp.sum(v**2) for v in reference(value)))(a)
    np.testing.assert_allclose(actual, expected, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_diagonal_maps_broadcast_and_gradients(jit):
    vector, matrix = sm.MatrixSymbol('v', 3, 1), sm.MatrixSymbol('A', 2, 3)
    function = sm.jax_lambdify((vector, matrix),
                              (sm.DiagMatrix(vector), sm.DiagonalOf(matrix),
                               sm.DiagonalMatrix(matrix)), jit=jit)
    v, a = jnp.arange(6.).reshape(2, 1, 3, 1), jnp.arange(18.).reshape(1, 3, 2, 3)

    def reference(x, y):
        x, y = jnp.broadcast_to(x, (2, 3, 3, 1)), jnp.broadcast_to(y, (2, 3, 2, 3))
        return (x*jnp.eye(3), jnp.diagonal(y, axis1=-2, axis2=-1)[..., None],
                y*jnp.eye(2, 3))

    for value, target in zip(function(v, a), reference(v, a)):
        np.testing.assert_allclose(value, target, rtol=2e-6)
    gradients = jax.grad(lambda x, y: sum(jnp.sum(value**2) for value in function(x, y)),
                         argnums=(0, 1))(v, a)
    references = jax.grad(lambda x, y: sum(jnp.sum(value**2) for value in reference(x, y)),
                         argnums=(0, 1))(v, a)
    for value, target in zip(gradients, references):
        np.testing.assert_allclose(value, target, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_diagonal_maps_complex_adjoint(jit):
    vector, matrix = sm.MatrixSymbol('v', 3, 1), sm.MatrixSymbol('A', 2, 3)
    function = sm.jax_lambdify((vector, matrix),
                              (sm.DiagMatrix(vector).H, sm.DiagonalOf(matrix).H,
                               sm.DiagonalMatrix(matrix).H), jit=jit)
    v = jnp.asarray([[1+2j], [3-1j], [2j]])
    a = jnp.asarray([[1+1j, 2j, 3], [4, 5-2j, 6j]])
    adjoint = lambda x: jnp.swapaxes(jnp.conj(x), -2, -1)
    expected = (adjoint(jnp.diag(v[:, 0])), adjoint(jnp.diag(a)[:, None]),
                adjoint(a*jnp.eye(2, 3)))
    for value, target in zip(function(v, a), expected):
        np.testing.assert_allclose(value, target, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_empty_diagonal_outputs_keep_batch_shape_and_dtype(jit):
    x = sm.Symbol('x')
    outputs = (sm.DiagMatrix(sm.ZeroMatrix(0, 1)), sm.DiagonalOf(sm.ZeroMatrix(0, 4)),
               sm.DiagonalMatrix(sm.ZeroMatrix(0, 4)))
    function = sm.jax_lambdify(x, outputs, jit=jit)
    value = jnp.ones((2, 3), dtype=jnp.float32)
    for result, shape in zip(function(value), [(0, 0), (0, 1), (0, 4)]):
        assert result.shape == (2, 3) + shape
        assert result.dtype == value.dtype


@pytest.mark.parametrize("jit", [False, True])
def test_block_matrix_broadcast_products_and_gradients(jit):
    a, b, d = (sm.MatrixSymbol('A', 2, 2), sm.MatrixSymbol('B', 2, 1),
               sm.MatrixSymbol('D', 1, 1))
    matrix = sm.BlockMatrix([[a, b], [b.T, d]])
    function = sm.jax_lambdify((a, b, d),
                              (matrix, matrix.T, sm.block_collapse(matrix*matrix)), jit=jit)
    av = jnp.arange(8.).reshape(2, 1, 2, 2)
    bv = jnp.arange(6.).reshape(1, 3, 2, 1)
    dv = jnp.ones((1, 1))

    def reference(x, y, z):
        x, y, z = (jnp.broadcast_to(x, (2, 3, 2, 2)),
                   jnp.broadcast_to(y, (2, 3, 2, 1)), jnp.broadcast_to(z, (2, 3, 1, 1)))
        result = jnp.concatenate((jnp.concatenate((x, y), axis=-1),
                                  jnp.concatenate((jnp.swapaxes(y, -2, -1), z), axis=-1)), axis=-2)
        return result, jnp.swapaxes(result, -2, -1), result @ result

    for value, target in zip(function(av, bv, dv), reference(av, bv, dv)):
        np.testing.assert_allclose(value, target, rtol=2e-6)
    gradients = jax.grad(lambda x, y, z: sum(jnp.sum(value**2) for value in function(x, y, z)),
                         argnums=(0, 1, 2))(av, bv, dv)
    references = jax.grad(lambda x, y, z: sum(jnp.sum(value**2) for value in reference(x, y, z)),
                         argnums=(0, 1, 2))(av, bv, dv)
    for value, target in zip(gradients, references):
        np.testing.assert_allclose(value, target, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_rectangular_block_diagonal_and_slice_batched_gradients(jit):
    a, b = sm.MatrixSymbol('A', 2, 3), sm.MatrixSymbol('B', 3, 2)
    function = sm.jax_lambdify((a, b), (sm.BlockDiagMatrix(a, b), a[:, 1:3]), jit=jit)
    av = jnp.arange(12.).reshape(2, 1, 2, 3)
    bv = jnp.arange(18.).reshape(1, 3, 3, 2)

    def reference(x, y):
        x, y = jnp.broadcast_to(x, (2, 3, 2, 3)), jnp.broadcast_to(y, (2, 3, 3, 2))
        top = jnp.concatenate((x, jnp.zeros((2, 3, 2, 2))), axis=-1)
        bottom = jnp.concatenate((jnp.zeros((2, 3, 3, 3)), y), axis=-1)
        return (jnp.concatenate((top, bottom), axis=-2), x[..., :, 1:3])

    for result, expected in zip(function(av, bv), reference(av, bv)):
        np.testing.assert_allclose(result, expected, rtol=2e-6)
    actual = jax.grad(lambda x, y: sum(jnp.sum(value) for value in function(x, y)), argnums=(0, 1))(av, bv)
    expected = jax.grad(lambda x, y: sum(jnp.sum(value) for value in reference(x, y)), argnums=(0, 1))(av, bv)
    for result, expected in zip(actual, expected):
        np.testing.assert_allclose(result, expected, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_hadamard_kronecker_batched_codegen_and_grad(jit):
    a, b = [sm.MatrixSymbol(name, 2, 2) for name in ('A', 'B')]
    outputs = (sm.HadamardProduct(a, b), sm.HadamardPower(a, 2),
               sm.HadamardPower(2, a), sm.HadamardPower(a, b), sm.KroneckerProduct(a, b))
    function = sm.jax_lambdify((a, b), outputs, jit=jit)
    av = jnp.arange(8.).reshape(2, 1, 2, 2)/10+1
    bv = jnp.arange(12.).reshape(1, 3, 2, 2)/10+1

    def reference(x, y):
        x, y = jnp.broadcast_arrays(x, y)
        kron = (x[..., :, None, :, None]*y[..., None, :, None, :]).reshape(2, 3, 4, 4)
        return (x*y, x**2, 2**x, x**y, kron)

    for result, expected in zip(function(av, bv), reference(av, bv)):
        np.testing.assert_allclose(result, expected, rtol=2e-6)
    actual = jax.grad(lambda x, y: sum(jnp.sum(value) for value in function(x, y)), argnums=(0, 1))(av, bv)
    expected = jax.grad(lambda x, y: sum(jnp.sum(value) for value in reference(x, y)), argnums=(0, 1))(av, bv)
    for result, expected in zip(actual, expected):
        np.testing.assert_allclose(result, expected, rtol=2e-6)


@pytest.mark.parametrize("jit",[False,True])
def test_matrix_trace_determinant_codegen_and_gradient(jit):
    a=sm.MatrixSymbol("A",2,2)
    function=sm.jax_lambdify(a,sm.Trace(a)+sm.Determinant(a),jit=jit)
    value=jnp.asarray([[[1.,2.],[3.,4.]],[[2.,1.],[0.,3.]]])
    expected=lambda x:jnp.trace(x,axis1=-2,axis2=-1)+jnp.linalg.det(x)
    np.testing.assert_allclose(function(value),expected(value),rtol=1e-6,atol=1e-6)
    actual=jax.grad(lambda x:jnp.sum(function(x)))(value)
    reference=jax.grad(lambda x:jnp.sum(expected(x)))(value)
    np.testing.assert_allclose(actual,reference,rtol=1e-6,atol=1e-6)


@pytest.mark.parametrize("jit", [False, True])
def test_matrix_symbols_batch_codegen_and_grad(jit):
    a, b = sm.MatrixSymbol("A", 2, 3), sm.MatrixSymbol("B", 3, 2)
    f = sm.jax_lambdify((a, b), a*b, jit=jit)
    ta = jnp.arange(12.).reshape(2, 2, 3)
    tb = jnp.arange(6.).reshape(3, 2)
    np.testing.assert_allclose(f(ta, tb), ta @ tb)
    actual = jax.grad(lambda x: jnp.sum(f(x, tb)))(ta)
    expected = jax.grad(lambda x: jnp.sum(x @ tb))(ta)
    np.testing.assert_allclose(actual, expected)
    with pytest.raises(ValueError, match="Matrix input"):
        f(ta[..., :2], tb)


@pytest.fixture(scope="module", autouse=True)
def _require_requested_gpu():
    if os.environ.get("SYMOX_JAX_REQUIRE_GPU") == "1":
        assert jax.default_backend() in ("gpu", "cuda"), "JAX tests require a CUDA GPU backend"


def _enable_x64(enabled=True):
    context = getattr(jax, "enable_x64", None)
    if context is None:
        context = jax.experimental.enable_x64
    return context(enabled)


@pytest.mark.parametrize("jit", [False, True])
@pytest.mark.parametrize("cse", [False, True])
def test_broadcast_matrix_and_structured_outputs(jit, cse):
    x, y = sm.symbols("x y")
    expression = sm.Matrix([[sm.sin(x+y), 1], [x*y, sm.Rational(1, 2)]])
    f = sm.jax_lambdify((x, y), (expression, [x > y, sm.I*x]), jit=jit, cse=cse)
    tx, ty = jnp.arange(4.0).reshape(4, 1), jnp.arange(3.0).reshape(1, 3)
    bx, by = jnp.broadcast_arrays(tx, ty)
    actual = f(tx, ty)
    expected = jnp.stack([jnp.sin(bx+by), jnp.ones_like(bx), bx*by,
                          jnp.full_like(bx, .5)], axis=-1).reshape(4, 3, 2, 2)
    assert isinstance(actual, tuple) and isinstance(actual[1], list)
    np.testing.assert_allclose(actual[0], expected, rtol=1e-6, atol=1e-7)
    np.testing.assert_array_equal(actual[1][0], bx > by)
    np.testing.assert_allclose(actual[1][1], 1j*bx)


@pytest.mark.parametrize("jit", [False, True])
def test_grad_hessian_vmap_and_symbolic_gradient(jit):
    x, y = sm.symbols("x y")
    expr = sm.sin(x*y) + sm.exp(x)/(1+y*y) + sm.sqrt(x+2)
    f = sm.jax_lambdify((x, y), expr, jit=jit)
    reference = lambda a, b: jnp.sin(a*b) + jnp.exp(a)/(1+b*b) + jnp.sqrt(a+2)
    point = (jnp.asarray(.2), jnp.asarray(.3))
    symbolic = sm.jax_lambdify((x, y), [sm.diff(expr, x), sm.diff(expr, y)])(*point)
    gradient = jax.grad(f, argnums=(0, 1))(*point)
    for actual, expected in zip(gradient, symbolic):
        np.testing.assert_allclose(actual, expected, rtol=2e-6)
    np.testing.assert_allclose(jax.hessian(f)(*point), jax.hessian(reference)(*point), rtol=2e-6)
    values = jnp.linspace(.1, 1, 7)
    np.testing.assert_allclose(jax.vmap(f, in_axes=(0, None))(values, point[1]),
                               reference(values, point[1]), rtol=2e-6)


@pytest.mark.parametrize("operation", ["reciprocal", "sqrt", "log"])
@pytest.mark.parametrize("jit", [False, True])
def test_inactive_piecewise_domains_and_second_derivatives(operation, jit):
    x = sm.Symbol("x")
    branch = {"reciprocal": 1/x, "sqrt": sm.sqrt(x), "log": sm.log(x)}[operation]
    expr = sm.Piecewise((branch, x > 0), (x*x, True))
    f = sm.jax_lambdify(x, expr, jit=jit)

    def reference(value):
        safe = jnp.where(value > 0, value, 1)
        active = {"reciprocal": lambda: 1/safe, "sqrt": lambda: jnp.sqrt(safe),
                  "log": lambda: jnp.log(safe)}[operation]()
        return jnp.where(value > 0, active, value*value)

    values = jnp.array([-2., 0., .25, 1.])
    np.testing.assert_allclose(f(values), reference(values), rtol=2e-6)
    for degree in (1, 2):
        df, expected = f, reference
        for _ in range(degree):
            df, expected = jax.grad(df), jax.grad(expected)
        actual = jax.jit(jax.vmap(df))(values)
        assert np.isfinite(actual).all()
        np.testing.assert_allclose(actual, jax.vmap(expected)(values), rtol=2e-6)


def test_empty_batches_constants_and_numpy_output():
    x = sm.Symbol("x")
    f = sm.jax_lambdify(x, np.array([[x, 1], [sm.pi, sm.Rational(1, 2)]], dtype=object), jit=True)
    assert f(jnp.empty((2, 0, 3))).shape == (2, 0, 3, 2, 2)
    np.testing.assert_allclose(sm.jax_lambdify((), sm.pi, jit=True)(), np.pi, rtol=1e-6)
    np.testing.assert_allclose(sm.jax_lambdify(x, sm.sqrt(x))(jnp.array([1, 4, 9])), [1, 2, 3])


@pytest.mark.parametrize("module", ["jax", jax])
def test_lambdify_route_and_source_inspection(module):
    x = sm.Symbol("x")
    f = sm.lambdify(x, sm.sin(x)+x*x, module, jit=True)
    np.testing.assert_allclose(f(.3), jnp.sin(.3)+.09, rtol=1e-6)
    assert "def " in inspect.getsource(f)
    assert "import symox" not in f.__source__


def test_nested_arguments_identifiers_and_distinct_assumptions():
    x, y = sm.Symbol("lambda"), sm.Symbol("α")
    a, b = sm.Symbol("same", positive=True), sm.Symbol("same", negative=True)
    f = sm.jax_lambdify(((x, y), a, b), x+2*y+3*a+4*b, jit=True)
    np.testing.assert_allclose(f((1., 2.), 3., -4.), -2.)
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        sm.jax_lambdify((a, a), a)
    with pytest.raises(ValueError, match="[Ss]ymbol|[Uu]nbound"):
        sm.jax_lambdify(x, x+y)


def test_custom_functions_and_standalone_source(tmp_path):
    x = sm.Symbol("x")
    special = sm.Function("custom")(x)
    f = sm.jax_lambdify(x, special, user_functions={"custom": jnp.sin}, jit=True)
    np.testing.assert_allclose(f(jnp.arange(3.)), jnp.sin(jnp.arange(3.)))
    source = sm.jaxcode_function("generated", x, (special, x*x),
                                 user_functions={"custom": "jax.numpy.sin"})
    script = tmp_path / "jax_export.py"
    script.write_text("import jax\nimport builtins\n"
                      "original_import = builtins.__import__\n"
                      "def checked_import(name, *args, **kwargs):\n"
                      "    if name.split('.')[0] in {'symox', 'sympy', 'torch'}:\n"
                      "        raise AssertionError(name)\n"
                      "    return original_import(name, *args, **kwargs)\n"
                      "builtins.__import__ = checked_import\n" + source +
                      "\nx = jnp.asarray(.3)\ny = jax.jit(generated)(x)\n"
                      "assert abs(float(y[0]) - float(jnp.sin(x))) < 1e-6\n"
                      "assert abs(float(jax.grad(lambda t: generated(t)[1])(x)) - .6) < 1e-6\n")
    subprocess.run([sys.executable, str(script)], check=True, capture_output=True, text=True)


def test_optional_import_and_generation_errors():
    subprocess.run([sys.executable, "-c", "import sys; import symox; assert 'jax' not in sys.modules; assert 'torch' not in sys.modules"],
                   check=True, capture_output=True, text=True)
    x = sm.Symbol("x")
    with pytest.raises((TypeError, ValueError)):
        sm.jaxcode_function("invalid name", x, x)
    with pytest.raises((ValueError, NotImplementedError)):
        sm.jax_lambdify(x, sm.besselj(3, x))
    with pytest.raises((TypeError, ValueError)):
        sm.jax_lambdify(x, x, jit_options={"invalid": True})


def test_low_level_printer():
    x = sm.Symbol("x")
    source = sm.jaxcode(sm.sin(x)+x*x)
    assert "jnp.sin" in source and "torch" not in source
    with pytest.raises((ValueError, NotImplementedError)):
        sm.jaxcode(sm.Function("unknown")(x))


@pytest.mark.parametrize("name,module,function,values", [
    ("sin", "numpy", "sin", [-.5, .2, 1.]),
    ("cos", "numpy", "cos", [-.5, .2, 1.]),
    ("tan", "numpy", "tan", [-.5, .2, 1.]),
    ("asin", "numpy", "arcsin", [-.5, .2, .9]),
    ("acos", "numpy", "arccos", [-.5, .2, .9]),
    ("atan", "numpy", "arctan", [-2., 0., 2.]),
    ("sinh", "numpy", "sinh", [-.5, .2, 1.]),
    ("cosh", "numpy", "cosh", [-.5, .2, 1.]),
    ("tanh", "numpy", "tanh", [-.5, .2, 1.]),
    ("asinh", "numpy", "arcsinh", [-2., 0., 2.]),
    ("acosh", "numpy", "arccosh", [1.1, 2., 4.]),
    ("atanh", "numpy", "arctanh", [-.5, .2, .9]),
    ("exp", "numpy", "exp", [-.5, .2, 1.]),
    ("log", "numpy", "log", [.5, 1., 2.]),
    ("Abs", "numpy", "abs", [-1., 0., 1.]),
    ("sign", "numpy", "sign", [-1., 0., 1.]),
    ("floor", "numpy", "floor", [-1.5, .5, 2.5]),
    ("ceiling", "numpy", "ceil", [-1.5, .5, 2.5]),
    ("erf", "special", "erf", [-1., 0., 1.]),
    ("erfc", "special", "erfc", [-1., 0., 1.]),
    ("loggamma", "special", "gammaln", [.5, 1.5, 2.5]),
    ("digamma", "special", "digamma", [.5, 1.5, 2.5]),
    ("erfinv", "special", "erfinv", [-.5, 0., .5]),
])
def test_elementary_function_values(name, module, function, values):
    import jax.scipy.special as jsp
    x = sm.Symbol("x", positive=True) if name == "loggamma" else sm.Symbol("x")
    inputs = jnp.asarray(values)
    expected = getattr(jnp if module == "numpy" else jsp, function)(inputs)
    actual = sm.jax_lambdify(x, getattr(sm, name)(x), jit=True)(inputs)
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=1e-7)


def test_min_max_heaviside_and_all_logical_operators():
    x, y = sm.symbols("x y")
    inputs = jnp.array([-2., 0., .5, 2.])
    other = jnp.asarray(.25)
    p, q = x > 0, x < 1
    expressions = [sm.Max(x, y, 1), sm.Min(x, y, -1), sm.Heaviside(x),
                   sm.And(p, q), sm.Or(p, q), sm.Not(p), sm.Xor(p, q),
                   sm.Implies(p, q), sm.Equivalent(p, q), sm.ITE(p, q, sm.false),
                   sm.Eq(x, 0), sm.Ne(x, 0), sm.atan2(x, y)]
    a, b = inputs > 0, inputs < 1
    expected = [jnp.maximum(jnp.maximum(inputs, other), 1),
                jnp.minimum(jnp.minimum(inputs, other), -1), jnp.heaviside(inputs, .5),
                a & b, a | b, ~a, a ^ b, (~a) | b, a == b, a & b,
                inputs == 0, inputs != 0, jnp.arctan2(inputs, other)]
    for actual, reference in zip(sm.jax_lambdify((x, y), expressions, jit=True)(inputs, other), expected):
        np.testing.assert_allclose(actual, reference, rtol=2e-6)


def test_complex_values_and_real_output_dtypes():
    x = sm.Symbol("x")
    inputs = jnp.asarray([1+2j, -2+1j], dtype=jnp.complex64)
    expressions = [sm.sin(x), sm.exp(x), sm.re(x), sm.im(x), sm.conjugate(x), sm.arg(x), sm.Abs(x)]
    expected = [jnp.sin(inputs), jnp.exp(inputs), inputs.real, inputs.imag,
                inputs.conj(), jnp.angle(inputs), jnp.abs(inputs)]
    for actual, reference in zip(sm.jax_lambdify(x, expressions, jit=True)(inputs), expected):
        assert actual.dtype == reference.dtype
        np.testing.assert_allclose(actual, reference, rtol=2e-6)


def test_supported_special_functions():
    import jax.scipy.special as jsp
    x = sm.Symbol("x")
    inputs = jnp.asarray([.5, 1.5, 2.5])
    expressions = [sm.gamma(x), sm.polygamma(2, x), sm.trigamma(x), sm.zeta(x+2),
                   sm.erfcinv(x/3), sm.besseli(0, x), sm.besseli(1, x)]
    expected = [jsp.gamma(inputs), jsp.polygamma(2, inputs), jsp.polygamma(1, inputs),
                jsp.zeta(inputs+2, jnp.ones_like(inputs)), jsp.erfinv(1-inputs/3),
                jsp.i0(inputs), jsp.i1(inputs)]
    for actual, reference in zip(sm.jax_lambdify(x, expressions, jit=True)(inputs), expected):
        np.testing.assert_allclose(actual, reference, rtol=2e-6, atol=1e-7)
    for unsupported in [sm.besselj(0, x), sm.bessely(1, x), sm.besselk(0, x)]:
        with pytest.raises((ValueError, NotImplementedError)):
            sm.jax_lambdify(x, unsupported)


@pytest.mark.parametrize("enabled", [False, True])
def test_dtype_inference_respects_local_x64_configuration(enabled):
    x, y = sm.symbols("x y")
    f = sm.jax_lambdify((x, y), [x+y, sm.pi, sm.I*x, x > y], jit=True)
    original = jax.config.x64_enabled
    with _enable_x64(enabled):
        float_dtype = jnp.float64 if enabled else jnp.float32
        complex_dtype = jnp.complex128 if enabled else jnp.complex64
        inputs = jnp.asarray([.25, .5], dtype=float_dtype)
        result = f(inputs, .125)
        assert [value.dtype for value in result] == [float_dtype, float_dtype, complex_dtype, jnp.bool_]
        np.testing.assert_allclose(result[0], inputs+.125)
        np.testing.assert_allclose(result[1], np.pi, rtol=1e-6)
        np.testing.assert_allclose(result[2], 1j*inputs)
    assert jax.config.x64_enabled == original


@pytest.mark.parametrize("as_string", [False, True])
@pytest.mark.parametrize("jit", [False, True])
def test_explicit_device_and_dtype_for_python_scalars(as_string, jit):
    x = sm.Symbol("x")
    f = sm.jax_lambdify(x, [sm.sin(x), sm.Rational(1, 3)], jit=jit)
    device = jax.devices("cpu")[0]
    with _enable_x64():
        results = f(.5, dtype=jnp.float64, device="cpu" if as_string else device)
        for result in results:
            assert result.dtype == jnp.float64 and result.device == device
        np.testing.assert_allclose(results[0], np.sin(.5), rtol=1e-14)
        np.testing.assert_allclose(results[1], 1/3, rtol=1e-14)


@pytest.mark.parametrize("jit_options", [None, {"out_shardings": (None, None)}])
def test_explicit_device_compiles_constants_on_requested_backend(jit_options):
    device = jax.devices("cpu")[0]
    f = sm.jax_lambdify((), (sm.pi, sm.Rational(1, 3)), jit=True, jit_options=jit_options)
    compiled = f.lower(device="cpu").compile()
    assert compiled.runtime_executable().local_devices() == [device]
    for result, expected in zip(f(device="cpu"), (np.pi, 1/3)):
        assert result.device == device
        np.testing.assert_allclose(result, expected, rtol=1e-6)


@pytest.mark.parametrize("target_backend", ["cpu", "gpu"])
def test_explicit_device_moves_inputs_before_compilation(target_backend):
    try:
        gpu = jax.devices("gpu")[0]
    except RuntimeError:
        pytest.skip("Requires both CPU and GPU devices")
    cpu = jax.devices("cpu")[0]
    target, origin = (cpu, gpu) if target_backend == "cpu" else (gpu, cpu)
    x, y = sm.symbols("x y")
    f = sm.jax_lambdify(((x, y),), (sm.sin(x)+y*y, [sm.pi, x > y]), jit=True)
    inputs = tuple(jax.device_put(np.array([.25, .5], dtype=np.float32)+shift, origin)
                   for shift in (0., .1))
    compiled = f.lower(inputs, device=target).compile()
    assert compiled.runtime_executable().local_devices() == [target]
    result = f(inputs, device=target)
    assert all(leaf.device == target for leaf in jax.tree_util.tree_leaves(result))
    np.testing.assert_allclose(result[0], np.sin(np.asarray(inputs[0]))+np.asarray(inputs[1])**2,
                               rtol=1e-6)
    # Reuse the same generated function in the opposite direction, and allow
    # normal JAX placement inference when no explicit device is supplied.
    back = tuple(jax.device_put(value, target) for value in inputs)
    assert f(back, device=origin)[0].device == origin
    assert f(inputs)[0].device == origin


def test_explicit_device_preserves_nested_autodiff_and_jit_options():
    x, y = sm.symbols("x y")
    cpu = jax.devices("cpu")[0]
    f = sm.jax_lambdify(((x, y),), sm.sin(x)+y*y, jit=True,
                        jit_options={"inline": True})
    wrapped = lambda value: f((value, .3), device="cpu")
    with _enable_x64():
        result = f((.2, .3), device=cpu, dtype=jnp.float64)
        assert result.dtype == jnp.float64 and result.device == cpu
        np.testing.assert_allclose(jax.grad(wrapped)(.2), np.cos(.2), rtol=1e-12)
        np.testing.assert_allclose(jax.hessian(wrapped)(.2), -np.sin(.2), rtol=1e-12)
        values = jax.device_put(np.array([.1, .2, .3]), cpu)
        actual = jax.vmap(wrapped)(values)
        assert actual.device == cpu
        np.testing.assert_allclose(actual, np.sin(np.asarray(values))+.09, rtol=1e-12)
    assert f.eval_shape((.2, .3), device=cpu).shape == ()
    f.clear_cache()


def test_loggamma_preserves_assumptions_and_rejects_ambiguous_domains():
    for x in [sm.Symbol("x"), sm.Symbol("x", negative=True)]:
        with pytest.raises((ValueError, NotImplementedError), match="loggamma|positive|domain"):
            sm.jax_lambdify(x, sm.loggamma(x))
    x = sm.Symbol("x", positive=True)
    f = sm.jax_lambdify(x, sm.loggamma(x), jit=True)
    np.testing.assert_allclose(f(jnp.asarray(2.5)), .28468287, rtol=2e-6)


def test_boolean_input_preserves_numeric_constant_and_boolean_outputs():
    x = sm.Symbol("x")
    inputs = jnp.asarray([False, True])
    constants, truth, values = sm.jax_lambdify(x, (2, sm.true, x), jit=True)(inputs)
    assert jnp.issubdtype(constants.dtype, jnp.number) and constants.dtype != jnp.bool_
    assert truth.dtype == values.dtype == jnp.bool_
    np.testing.assert_array_equal(constants, [2, 2])
    np.testing.assert_array_equal(truth, [True, True])
    np.testing.assert_array_equal(values, inputs)


def test_cse_retains_assumptions_for_domain_sensitive_functions():
    import jax.scipy.special as jsp
    x, y = sm.symbols("x y", positive=True)
    f = sm.jax_lambdify((x, y), [sm.loggamma(x+y), sm.sin(x+y)], cse=True, jit=True)
    assert "_sx_t_0 = " in f.__source__
    inputs = (jnp.asarray([.5, 1.]), jnp.asarray([1., 2.]))
    result = f(*inputs)
    np.testing.assert_allclose(result[0], jsp.gammaln(inputs[0]+inputs[1]), rtol=2e-6)
    np.testing.assert_allclose(result[1], jnp.sin(inputs[0]+inputs[1]), rtol=2e-6)


def test_joint_cse_is_emitted_and_preserves_gradients():
    x, y = sm.symbols("x y")
    shared = sm.sin(x+y) + sm.exp(x*y)
    outputs = (shared**2, shared**3+x)
    cse = sm.jax_lambdify((x, y), outputs, cse=True, jit=True)
    direct = sm.jax_lambdify((x, y), outputs, cse=False, jit=True)
    assert "_sx_t_0 = " in cse.__source__
    assert "_sx_t_0 = " not in direct.__source__
    inputs = (jnp.asarray(.2), jnp.asarray(.4))
    for actual, expected in zip(cse(*inputs), direct(*inputs)):
        np.testing.assert_allclose(actual, expected, rtol=2e-6)
    actual = jax.grad(lambda a, b: sum(cse(a, b)), argnums=(0, 1))(*inputs)
    expected = jax.grad(lambda a, b: sum(direct(a, b)), argnums=(0, 1))(*inputs)
    np.testing.assert_allclose(actual, expected, rtol=2e-6)


def test_piecewise_first_match_default_nan_and_empty_dense():
    x = sm.Symbol("x")
    inputs = jnp.asarray([-1., 0., 1., 2.])
    expression = sm.Piecewise((1, x < 1), (2, x < 2), (3, True))
    np.testing.assert_array_equal(sm.jax_lambdify(x, expression, jit=True)(inputs), [1, 1, 2, 3])
    partial = sm.jax_lambdify(x, sm.Piecewise((x*x, x > 0)), jit=True)
    np.testing.assert_allclose(partial(inputs), [np.nan, np.nan, 1, 4], equal_nan=True)
    empty = sm.jax_lambdify(x, np.empty((2, 0), dtype=object), jit=True)(inputs)
    assert empty.shape == (4, 2, 0)


def test_nested_printer_symbol_and_function_override():
    from symox.printing.jax import JaxPrinter
    x = sm.Symbol("x")

    class Renamed(JaxPrinter):
        def _print_Symbol(self, expr):
            return "v_" + expr.name

    class Cosine(Renamed):
        def _print_sin(self, expr):
            return "jnp.cos(" + self._print(expr.args[0]) + ")"

    inputs = jnp.asarray([.2, .7])
    helper = lambda value, reference: jnp.asarray(value)
    expression = sm.sin(x) + x*x
    renamed = Renamed().doprint(expression)
    actual = eval(renamed, {"jnp": jnp, "_sx_tensor": helper, "v_x": inputs})
    np.testing.assert_allclose(actual, jnp.sin(inputs)+inputs*inputs, rtol=2e-6)
    overridden = Cosine().doprint(expression)
    actual = eval(overridden, {"jnp": jnp, "_sx_tensor": helper, "v_x": inputs})
    np.testing.assert_allclose(actual, jnp.cos(inputs)+inputs*inputs, rtol=2e-6)


@pytest.mark.parametrize("jit", [False, True])
@pytest.mark.parametrize("pivot", ["A", "D"])
def test_block_schur_rectangular_complex_decomposition_reconstruction_and_gradients(jit, pivot):
    # Complex64 reconstruction must use full float32 products on GPUs, rather
    # than the default reduced-precision matmul kernels.
    with jax.default_matmul_precision("highest"):
        _check_block_schur_rectangular_complex(jit, pivot)


def _check_block_schur_rectangular_complex(jit, pivot):
    t = sm.Symbol("t", real=True)
    if pivot == "A":
        a = sm.ImmutableMatrix([[2+sm.I+t]])
        b = sm.ImmutableMatrix([[1+sm.I, 2, -sm.I]])
        c = sm.ImmutableMatrix([[2], [1-sm.I]])
        d = sm.ImmutableMatrix([[3+t, 1+sm.I, 2-t], [-1, 4-sm.I+t, 1+2*t]])
    else:
        a = sm.ImmutableMatrix([[1+t, 2+sm.I, 3], [-sm.I, 4-t, 1+2*t]])
        b = sm.ImmutableMatrix([[2, 1+sm.I], [1, -1]])
        c = sm.ImmutableMatrix([[1, 2, -sm.I], [3, 0, 2+sm.I]])
        d = sm.ImmutableMatrix([[3+sm.I+t, 0], [0, 4-2*sm.I+t]])
    matrix = sm.BlockMatrix([[a, b], [c, d]])
    factors = matrix.LDUdecomposition() if pivot == "A" else matrix.UDLdecomposition()
    function = sm.jax_lambdify(t, (matrix.schur(pivot), *factors), jit=jit)
    tv = jnp.asarray([[-.3, .2], [.4, .8]])

    def reference(value):
        scale = value[..., None, None]
        dtype = jnp.complex128 if value.dtype == jnp.float64 else jnp.complex64
        def constant(values):
            array = jnp.asarray(values, dtype=dtype)
            return jnp.broadcast_to(array, value.shape+array.shape)
        if pivot == "A":
            aa = constant([[2+1j]])+scale
            bb = constant([[1+1j, 2, -1j]])
            cc = constant([[2], [1-1j]])
            dd = constant([[3, 1+1j, 2], [-1, 4-1j, 1]])+scale*constant([[1, 0, -1], [0, 1, 2]])
            complement = dd-cc@jnp.linalg.solve(aa, bb)
        else:
            aa = constant([[1, 2+1j, 3], [-1j, 4, 1]])+scale*constant([[1, 0, 0], [0, -1, 2]])
            bb = constant([[2, 1+1j], [1, -1]])
            cc = constant([[1, 2, -1j], [3, 0, 2+1j]])
            dd = constant([[3+1j, 0], [0, 4-2j]])+scale*constant([[1, 0], [0, 1]])
            complement = aa-bb@jnp.linalg.solve(dd, cc)
        dense = jnp.concatenate((jnp.concatenate((aa, bb), axis=-1), jnp.concatenate((cc, dd), axis=-1)), axis=-2)
        return complement, dense

    def actual_outputs(value):
        complement, left, diagonal, right = function(value)
        return complement, left@diagonal@right

    complement, left, diagonal, right = function(tv)
    expected_complement, expected_dense = reference(tv)
    assert left.shape[-2:] == (expected_dense.shape[-2],)*2
    assert diagonal.shape[-2:] == expected_dense.shape[-2:]
    assert right.shape[-2:] == (expected_dense.shape[-1],)*2
    for actual, expected in zip(actual_outputs(tv), (expected_complement, expected_dense)):
        assert actual.dtype == expected.dtype
        np.testing.assert_allclose(actual, expected, rtol=3e-5, atol=3e-6)
    # Differentiate the real squared norm of complex outputs; backend solve
    # and direct dense assembly provide an independent reconstruction oracle.
    loss = lambda outputs: sum(jnp.sum(jnp.abs(value)**2) for value in outputs)
    actual_gradient = jax.grad(lambda value: loss(actual_outputs(value)))(tv)
    expected_gradient = jax.grad(lambda value: loss(reference(value)))(tv)
    np.testing.assert_allclose(actual_gradient, expected_gradient, rtol=3e-5, atol=3e-5)
