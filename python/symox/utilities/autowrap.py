"""Compile numeric routines and genuine NumPy double-precision ufuncs.

Compilation uses the requested optional backend. No symbolic reference package
is imported, and failed compiler commands report their captured diagnostics.
"""
import importlib.util
import itertools
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import tempfile
from .codegen import CCodeGen, FCodeGen, OutputArgument, InOutArgument
from .lambdify import implemented_function


class CodeWrapError(RuntimeError):pass
_counter=itertools.count()

def _run(command,cwd,verbose=False,env=None):
    result=subprocess.run(command,cwd=cwd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
    if verbose:print(result.stdout)
    if result.returncode:raise CodeWrapError(' '.join(map(str,command))+'\n'+result.stdout)

def _load(name,directory):
    paths=list(Path(directory).glob(name+'*'+sysconfig.get_config_var('EXT_SUFFIX')))
    if not paths:raise CodeWrapError('compiler produced no Python extension '+name)
    spec=importlib.util.spec_from_file_location(name,paths[0]);module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CodeWrapper:
    def __init__(self,generator,filepath=None,flags=None,verbose=False,**kwargs):
        self.generator,self.filepath,self.flags,self.verbose=generator,filepath,list(flags or []),verbose
        self.options=kwargs;self.module_name='symox_wrapped_'+str(next(_counter));self.filename=self.module_name
    def wrap_code(self,routine,helpers=None):
        directory=Path(self.filepath or tempfile.mkdtemp(prefix='symox-autowrap-')).resolve()
        directory.mkdir(parents=True,exist_ok=True)
        self.generator.write([*(helpers or []),routine],str(directory/self.filename),to_files=True)
        return self._compile(routine,directory)
    def _compile(self,routine,directory):raise NotImplementedError('choose a compilation backend')


class DummyWrapper(CodeWrapper):
    def _compile(self,routine,directory):
        def wrapped():return ', '.join(str(r.expr) for r in routine.result_variables)
        wrapped.args=', '.join(str(a.name) for a in routine.arguments if not isinstance(a,OutputArgument))
        wrapped.returns=', '.join(str(r.result_var) if isinstance(r,OutputArgument) else 'nameless' for r in routine.result_variables)
        return wrapped


class F2PyCodeWrapper(CodeWrapper):
    def _compile(self,routine,directory):
        import numpy as np
        flags=list(self.flags)
        # Modern f2py supports Meson, avoiding numpy.distutils' incompatibility
        # with current setuptools. Preserve explicit flags/backend selections
        # and the legacy path for older NumPy or environments without Meson.
        if not flags and np.lib.NumpyVersion(np.__version__) >= '1.26.0' and importlib.util.find_spec('mesonbuild') is not None:
            flags=['--backend','meson']
        environment=os.environ.copy()
        environment['PATH']=sysconfig.get_path('scripts')+os.pathsep+environment.get('PATH','')
        _run([sys.executable,'-m','numpy.f2py','-c','-m',self.module_name,self.filename+'.f90',*flags],directory,self.verbose,env=environment)
        module=_load(self.module_name,directory)
        return getattr(module,routine.name.lower())


class CythonCodeWrapper(CodeWrapper):
    def _compile(self,routine,directory):
        if importlib.util.find_spec('Cython') is None:raise ImportError('Cython backend requires Cython')
        params=[];call=[];declarations=[];outputs=[]
        for arg in routine.arguments:
            if arg.dimensions:raise NotImplementedError('Cython wrapper currently requires scalar arguments')
            if isinstance(arg,OutputArgument):
                if isinstance(arg,InOutArgument):params.append('double '+str(arg.name))
                else:declarations.append('    cdef double '+str(arg.name))
                call.append('&'+str(arg.name));outputs.append(str(arg.name))
            else:params.append('double '+str(arg.name));call.append(str(arg.name))
        if routine.results:
            declarations.append('    cdef double result');outputs.insert(0,'result')
        prototype=self.generator._prototype(routine).replace(routine.name,'native "'+routine.name+'"',1)
        pyx='cdef extern from "'+self.filename+'.h":\n    '+prototype+'\n\ndef '+routine.name+'('+', '.join(params)+'):\n'
        pyx+='\n'.join(declarations)+'\n    '+('result = ' if routine.results else '')+'native('+', '.join(call)+')\n    return '+(', '.join(outputs) or 'None')+'\n'
        (directory/(self.module_name+'.pyx')).write_text(pyx)
        setup='from setuptools import setup, Extension\nfrom Cython.Build import cythonize\nsetup(ext_modules=cythonize([Extension('+repr(self.module_name)+', ['+repr(self.module_name+'.pyx')+', '+repr(self.filename+'.c')+'], libraries=["m"])], language_level=3))\n'
        (directory/'setup.py').write_text(setup)
        _run([sys.executable,'setup.py','build_ext','--inplace',*self.flags],directory,self.verbose)
        return getattr(_load(self.module_name,directory),routine.name)


class UfuncifyCodeWrapper(CodeWrapper):
    def wrap_code(self,routine,helpers=None):
        # Helpers are (name, expression, arguments) triples compiled as
        # static C functions that the loop body calls by name.
        self.helpers=list(helpers or ())
        directory=Path(self.filepath or tempfile.mkdtemp(prefix='symox-ufunc-')).resolve()
        directory.mkdir(parents=True,exist_ok=True)
        return self._compile(routine,directory)
    def _compile(self,routine,directory):
        import numpy as np
        from ..printing.c import C99CodePrinter
        expressions=[r.expr for r in routine.results]
        arguments=[a for a in routine.arguments if not isinstance(a,OutputArgument)]
        if len(arguments)!=len(routine.arguments):raise NotImplementedError('ufunc outputs must be scalar expressions')
        if not expressions:raise ValueError('a ufunc requires at least one output')
        nin,nout=len(arguments),len(expressions);ntypes=nin+nout
        if ntypes>32:raise ValueError('NumPy ufuncs support at most 32 arguments and outputs')
        helper_names={str(name):str(name) for name,_,_ in getattr(self,'helpers',())}
        printer=C99CodePrinter({'user_functions':helper_names});lines=[]
        helper_source=''
        for name,value,hargs in getattr(self,'helpers',()):
            helper_source+='static double '+str(name)+'('+', '.join('double '+str(a) for a in hargs)+') {\n return '+printer.doprint(value)+';\n}\n'
        for i,arg in enumerate(arguments):lines.append('double '+str(arg.name)+' = *(double *)(args['+str(i)+'] + j*steps['+str(i)+']);')
        for i,expr in enumerate(expressions,nin):lines.append('*(double *)(args['+str(i)+'] + j*steps['+str(i)+']) = '+printer.doprint(expr)+';')
        name=self.module_name
        source='#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION\n#include <Python.h>\n#include <numpy/arrayobject.h>\n#include <numpy/ufuncobject.h>\n#include <math.h>\n'
        source+=helper_source
        source+='static void loop(char **args, const npy_intp *dims, const npy_intp *steps, void *data) {\n for(npy_intp j=0;j<dims[0];++j){\n'+ '\n'.join(lines)+'\n }}\n'
        source+='static PyUFuncGenericFunction loops[1]={loop};\nstatic void *data[1]={NULL};\nstatic char types['+str(ntypes)+']={'+','.join(['NPY_DOUBLE']*ntypes)+'};\n'
        source+='static struct PyModuleDef def={PyModuleDef_HEAD_INIT,"'+name+'",NULL,-1,NULL};\n'
        source+='PyMODINIT_FUNC PyInit_'+name+'(void){import_array();import_umath();PyObject *m=PyModule_Create(&def);if(!m)return NULL;PyObject *u=PyUFunc_FromFuncAndData(loops,data,types,1,'+str(nin)+','+str(nout)+',PyUFunc_None,"'+routine.name+'","Generated by Symox",0);if(!u){Py_DECREF(m);return NULL;}if(PyModule_AddObject(m,"ufunc",u)<0){Py_DECREF(u);Py_DECREF(m);return NULL;}return m;}\n'
        (directory/(name+'.c')).write_text(source)
        _run(['cc','-shared','-fPIC','-O3','-std=c99','-I'+np.get_include(),'-I'+sysconfig.get_path('include'),name+'.c','-o',name+sysconfig.get_config_var('EXT_SUFFIX'),'-lm',*self.flags],directory,self.verbose)
        return _load(name,directory).ufunc


def autowrap(expr,language=None,backend='f2py',tempdir=None,args=None,flags=None,verbose=False,helpers=None,**kwargs):
    backend=backend.lower()
    wrappers={'f2py':F2PyCodeWrapper,'cython':CythonCodeWrapper,'dummy':DummyWrapper}
    if backend not in wrappers:raise ValueError('unsupported autowrap backend '+backend)
    language=(language or ('F95' if backend=='f2py' else 'C')).upper()
    if backend=='f2py' and language not in ('F95','F90','FORTRAN'):raise ValueError('f2py requires Fortran')
    if backend=='cython' and language not in ('C','C99','C89'):raise ValueError('Cython wrapper requires C')
    generator=FCodeGen() if language in ('F95','F90','FORTRAN') else CCodeGen()
    routine=generator.routine('autofunc',expr,argument_sequence=args)
    helper_routines=[generator.routine(name,value,arguments) for name,value,arguments in helpers or ()]
    return wrappers[backend](generator,tempdir,flags,verbose,**kwargs).wrap_code(routine,helper_routines)


def binary_function(symfunc,expr,**kwargs):return implemented_function(symfunc,autowrap(expr,**kwargs))


def ufuncify(args,expr,language=None,backend='numpy',tempdir=None,flags=None,verbose=False,helpers=None,**kwargs):
    if backend.lower()!='numpy':
        # The reference wraps an elementwise loop over its array arguments;
        # the compiled scalar routine is applied elementwise the same way.
        import numpy as np
        scalar=autowrap(expr,language,backend,tempdir,args,flags,verbose,helpers,**kwargs)
        return np.vectorize(scalar,otypes=[float])
    arguments=list(args) if isinstance(args,(list,tuple)) else [args]
    generator=CCodeGen();routine=generator.routine('autofunc',expr,argument_sequence=arguments)
    return UfuncifyCodeWrapper(generator,tempdir,flags,verbose,**kwargs).wrap_code(routine,helpers)


__all__=['CodeWrapError','CodeWrapper','CythonCodeWrapper','DummyWrapper','F2PyCodeWrapper','UfuncifyCodeWrapper','autowrap','binary_function','ufuncify']
