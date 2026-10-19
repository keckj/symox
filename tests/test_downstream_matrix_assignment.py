"""Matrix assignments preserve scalar printer overrides and row-major offsets."""
import ctypes
import shutil
import subprocess

import pytest
import sympy as sp
import symox as sx
from symox.printing.c import C99CodePrinter
from sympy.printing.c import C99CodePrinter as ReferencePrinter


@pytest.mark.parametrize('shape',[(1,1),(1,3),(2,2),(3,1)])
def test_matrix_assignment_matches_reference(shape):
    x=sx.Symbol('x')
    matrix=sx.Matrix(*shape,lambda i,j:(i+1)*x+j)
    output=sx.MatrixSymbol('out',*shape)
    assert C99CodePrinter().doprint(matrix,output)==ReferencePrinter().doprint(matrix.to_sympy(),output.to_sympy())


def test_matrix_assignment_executes_with_nested_input_overrides(tmp_path):
    compiler=shutil.which('cc')
    if compiler is None:pytest.skip('C compiler unavailable')
    x=sx.Symbol('x');f=sx.Function('f');time=sx.Symbol('t')
    class Printer(C99CodePrinter):
        def _print_Symbol(self,value):
            return 'input[0]' if value==x else super()._print_Symbol(value)
        def _print_Function(self,value):
            return 'input[1]' if value==f(time) else super()._print_Function(value)
    matrix=sx.Matrix([[x+f(time),x*f(time)],[x**2,f(time)**2]])
    body=Printer().doprint(matrix,sx.MatrixSymbol('output',2,2))
    source=tmp_path/'matrix.c';library=tmp_path/'matrix.so'
    source.write_text('#include <math.h>\nvoid evaluate(double *input,double *output){\n'+body+'\n}\n')
    subprocess.run([compiler,'-shared','-fPIC',str(source),'-lm','-o',str(library)],check=True,capture_output=True)
    function=ctypes.CDLL(str(library)).evaluate
    function.argtypes=[ctypes.POINTER(ctypes.c_double)]*2
    values=(ctypes.c_double*2)(2,3);result=(ctypes.c_double*4)()
    function(values,result)
    assert list(result)==[5,6,4,9]


def test_matrix_assignment_rejects_mismatched_shapes():
    with pytest.raises(ValueError):C99CodePrinter().doprint(sx.ones(2,2),sx.MatrixSymbol('out',3,1))


def test_octave_matrix_assignment_dispatches_symbolic_target():
    from symox.printing.octave import OctaveCodePrinter
    from sympy.printing.octave import OctaveCodePrinter as Reference
    x=sx.Symbol('x');matrix=sx.Matrix([[x,2],[3,4]]);target=sx.MatrixSymbol('out',2,2)
    assert OctaveCodePrinter().doprint(matrix,target)==Reference().doprint(matrix.to_sympy(),target.to_sympy())
