"""Compile and execute generated routines, rather than checking export names."""
import ctypes
import importlib
import math
from pathlib import Path
import shutil
import subprocess
import sys
import pytest
import numpy as np
import symox as S
from symox.codegen.ast import (Variable,Declaration,FunctionDefinition,Return,CodeBlock,For,Assignment,aug_assign,real,integer)


def compile_c(tmp_path,source,name='f',arguments=1):
    path=tmp_path/(name+'.c');lib=tmp_path/(name+'.so');path.write_text(source)
    result=subprocess.run(['cc','-shared','-fPIC','-std=c99',str(path),'-lm','-o',str(lib)],text=True,capture_output=True)
    assert result.returncode==0,result.stderr+'\n'+source
    library=ctypes.CDLL(str(lib));function=getattr(library,name);function.argtypes=[ctypes.c_double]*arguments;function.restype=ctypes.c_double
    function._library=library;return function


def test_typed_ast_emits_working_c_and_python(tmp_path):
    from symox.codegen.cutils import render_as_source_file
    from symox.codegen.pyutils import render_as_module
    x=S.Symbol('x');node=FunctionDefinition(real,'f',(Variable(x,real),),[Return(x*x+S.sin(x))])
    compiled=compile_c(tmp_path,render_as_source_file(node))
    namespace={};exec(render_as_module(node),namespace)
    for value in (-2.,0.,1.,3.):assert compiled(value)==pytest.approx(namespace['f'](value),rel=1e-14)


def test_newton_iteration_compiles_and_converges(tmp_path):
    from symox.codegen.algorithms import newtons_method_function
    from symox.codegen.cutils import render_as_source_file
    from symox.codegen.pyutils import render_as_module
    x=S.Symbol('x');node=newtons_method_function(x*x-2,x,func_name='f',itermax=30)
    cfunction=compile_c(tmp_path,render_as_source_file(node));namespace={};exec(render_as_module(node),namespace)
    for start in (.5,1.,2.,10.):
        assert cfunction(start)==pytest.approx(math.sqrt(2),rel=1e-13)
        assert namespace['f'](start)==pytest.approx(math.sqrt(2),rel=1e-13)


def test_code_generator_c_output_buffer_and_fortran(tmp_path):
    from symox.utilities.codegen import codegen
    x,y=S.symbols('x y')
    files=codegen(('f',x*x+S.sin(x)),'C',prefix=str(tmp_path/'cfunc'))
    assert len(files)==2 and compile_c(tmp_path,files[0][1])(1.)==pytest.approx(1+math.sin(1))
    matrix=S.Matrix([[x,y],[x+y,x*y]])
    csource=codegen(('matrix_fn',matrix),'C')[0][1]
    cpath=tmp_path/'matrix.c';cpath.write_text(csource);lib=tmp_path/'matrix.so'
    subprocess.run(['cc','-shared','-fPIC',str(cpath),'-lm','-o',str(lib)],check=True)
    library=ctypes.CDLL(str(lib));fn=library.matrix_fn;fn.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double)];fn.restype=None
    out=(ctypes.c_double*4)();fn(2,3,out);assert list(out)==[2,3,5,6]
    if shutil.which('gfortran'):
        source=codegen(('f',x*x+S.sin(x)),'F95')[0][1];path=tmp_path/'fortran.f90';path.write_text(source)
        subprocess.run(['gfortran','-shared','-fPIC',str(path),'-o',str(tmp_path/'fortran.so')],check=True)
        flib=ctypes.CDLL(str(tmp_path/'fortran.so'));ff=flib.f_;ff.argtypes=[ctypes.POINTER(ctypes.c_double)];ff.restype=ctypes.c_double
        value=ctypes.c_double(1);assert ff(ctypes.byref(value))==pytest.approx(1+math.sin(1))


def test_rust_routine_compiles_and_executes(tmp_path):
    from symox.utilities.codegen import codegen
    if not shutil.which('rustc'):pytest.skip('Rust compiler unavailable')
    x=S.Symbol('x');text=codegen(('f',x*x+S.sin(x)+1),'Rust')[0][1]
    path=tmp_path/'routine.rs';path.write_text(text+'\nfn main() { println!("{}", f(1.0)); }\n')
    result=subprocess.run(['rustc',str(path),'-o',str(tmp_path/'routine')],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    assert float(subprocess.check_output([str(tmp_path/'routine')],text=True))==pytest.approx(2+math.sin(1))


def test_stable_codegen_rewriting_calls_real_c_intrinsics(tmp_path):
    from symox.codegen.rewriting import optimize,expm1_opt,log1p_opt,exp2_opt
    from symox.codegen.cutils import render_as_source_file
    x=S.Symbol('x');expr=optimize(S.exp(x)-1,[expm1_opt]);assert expr.class_name=='expm1'
    assert log1p_opt(S.log(1+x)).class_name=='log1p'
    assert exp2_opt(2**x).class_name=='exp2'
    node=FunctionDefinition(real,'f',(Variable(x,real),),[Return(expr)])
    function=compile_c(tmp_path,render_as_source_file(node))
    for value in (1e-16,-1e-16,1.,-1.):assert function(value)==pytest.approx(math.expm1(value),rel=1e-15,abs=1e-30)


def test_fortran_program_nodes_compile(tmp_path):
    from symox.codegen.fnodes import Program,Do
    from symox.printing.fortran import FCodePrinter
    from symox.codegen.ast import Print
    if not shutil.which('gfortran'):pytest.skip('Fortran compiler unavailable')
    i=S.Symbol('i',integer=True)
    program=Program('demo',[Declaration(Variable(i,integer)),Do([Print((i,))],i,1,3)])
    text=FCodePrinter().doprint(program);path=tmp_path/'demo.f90';path.write_text(text)
    result=subprocess.run(['gfortran',str(path),'-o',str(tmp_path/'demo')],text=True,capture_output=True)
    assert result.returncode==0,result.stderr+'\n'+text
    assert list(map(int,subprocess.check_output([str(tmp_path/'demo')],text=True).split()))==[1,2,3]


def test_custom_numerical_functions_and_lambda_sources():
    from symox.utilities.lambdify import implemented_function,lambdastr
    x,y=S.symbols('x y');f=implemented_function('user_numeric',lambda value:value+2)
    function=S.lambdify(x,f(f(x))+S.sin(x),'numpy')
    values=np.array([0.,1.,2.]);assert function(values)==pytest.approx(values+4+np.sin(values))
    nested=eval(lambdastr((x,(y,)),x+y),{'math':math})
    assert nested(2,(3,))==5
    assert eval(lambdastr(x,S.sin(x)+x*x),{'math':math})(1)==pytest.approx(1+math.sin(1))


def test_recurrence_decorators_and_warning_metadata():
    from symox.utilities.memoization import recurrence_memo,assoc_recurrence_memo
    from symox.utilities.decorator import deprecated,conserve_mpmath_dps
    from symox.utilities.exceptions import SymPyDeprecationWarning
    import mpmath
    @recurrence_memo([1])
    def fac(n,prev):return n*prev[-1]
    assert fac(6)==720 and fac.cache_length()==7 and fac.fetch_item(slice(2,4))==[2,6]
    @assoc_recurrence_memo(lambda n:1)
    def triangle(n,m,prev):return prev[n][m-1]+1
    assert triangle(4,3)==4
    @deprecated('old',deprecated_since_version='0.1',active_deprecations_target='test')
    def old():return 3
    with pytest.warns(SymPyDeprecationWarning):assert old()==3
    previous=mpmath.mp.dps
    @conserve_mpmath_dps
    def change():mpmath.mp.dps=40
    change();assert mpmath.mp.dps==previous


def test_string_helpers_preserve_content():
    from symox.utilities.misc import rawlines,strlines,replace,translate,ordinal
    for text in ('hello','a\nb\n','space  \nend','a\\b\'\"'):
        assert eval(rawlines(text))==text
        assert eval(strlines(text,c=20))==text
    assert replace('foo',{'oo':'ar','f':'b'})=='bar'
    mapping={'abc':'x','a':'y'};assert translate('abc',mapping)=='x' and mapping=={'abc':'x','a':'y'}
    assert ordinal(11)=='11th' and ordinal(22)=='22nd'


def test_typed_rust_integer_arithmetic_and_mixed_loop_compile(tmp_path):
    from symox.printing.rust import RustCodePrinter
    if not shutil.which('rustc'):pytest.skip('Rust compiler unavailable')
    n,i,total=S.symbols('n i total')
    integer_function=FunctionDefinition(integer,'square',(Variable(n,integer),),[Return(n*n+2)])
    float_function=FunctionDefinition(real,'accumulate',(Variable(n,integer),),[
        Declaration(Variable(total,real,S.Integer(0))),
        For(i,S.Range(0,n),[aug_assign(total,'+',i)]),Return(total+1)])
    printer=RustCodePrinter()
    source=printer.doprint(integer_function)+'\n'+printer.doprint(float_function)
    source+='\nfn main() { println!("{} {}", square(4), accumulate(5)); }\n'
    path=tmp_path/'typed.rs';path.write_text(source)
    result=subprocess.run(['rustc',str(path),'-o',str(tmp_path/'typed')],text=True,capture_output=True)
    assert result.returncode==0,result.stderr+'\n'+source
    assert subprocess.check_output([str(tmp_path/'typed')],text=True).split()==['18','11']


def test_fortran_concurrent_do_compiles(tmp_path):
    from symox.codegen.fnodes import Program,Do
    from symox.printing.fortran import FCodePrinter
    if not shutil.which('gfortran'):pytest.skip('Fortran compiler unavailable')
    i=S.Symbol('i',integer=True)
    program=Program('concurrent_demo',[Declaration(Variable(i,integer)),Do([],i,1,3,concurrent=True)])
    source=FCodePrinter().doprint(program)
    path=tmp_path/'concurrent.f90';path.write_text(source)
    result=subprocess.run(['gfortran',str(path),'-o',str(tmp_path/'concurrent')],text=True,capture_output=True)
    assert result.returncode==0,result.stderr+'\n'+source
