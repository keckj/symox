"""Failure-path CUDA ownership checks; no CUDA installation required."""
import ctypes as ct
import threading

import numpy as np
import pytest
import symox as s
from symox.codegen import cuda


class FakeDriver:
    def __init__(self):
        self.events = []
        self.failures = {}
        self.local = threading.local()
        self.address = 1000

    def call(self, name, *arguments):
        self.events.append(name)
        failure = self.failures.get(name)
        if failure:
            if isinstance(failure, list):
                failure = failure.pop(0) if failure else None
            if failure:
                raise cuda.CudaError(failure)
        if name in ('cuDeviceGet', 'cuDeviceGetAttribute'):
            ct.cast(arguments[0],ct.POINTER(ct.c_int))[0] = 8 if name.endswith('Attribute') and arguments[1] == 75 else 0
        elif name in ('cuDevicePrimaryCtxRetain','cuModuleLoadData','cuModuleGetFunction'):
            ct.cast(arguments[0],ct.POINTER(ct.c_void_p))[0] = 123
        elif name == 'cuMemAlloc_v2':
            self.address += 16
            ct.cast(arguments[0],ct.POINTER(ct.c_uint64))[0] = self.address
        elif name == 'cuCtxPushCurrent_v2':
            stack = getattr(self.local,'stack',[])
            stack.append(arguments[0].value)
            self.local.stack = stack
        elif name == 'cuCtxPopCurrent_v2':
            self.local.stack.pop()


@pytest.fixture
def mocked_runtime(monkeypatch):
    driver = FakeDriver()
    monkeypatch.setattr(cuda,'_driver',lambda: driver)
    monkeypatch.setattr(cuda,'_compile',lambda *args: (b'.version 8.0\0','test-key',False))
    return driver


def test_close_releases_module_after_async_error_and_preserves_failure(mocked_runtime):
    driver=mocked_runtime
    kernel=cuda.compile_cuda('extern "C" __global__ void test() {}','test')
    driver.failures.update(cuCtxSynchronize='original async launch error',cuModuleUnload='unload error',cuDevicePrimaryCtxRelease_v2='release error')
    with pytest.raises(cuda.CudaError,match='original async launch error'):
        kernel.close()
    assert driver.events.count('cuModuleUnload') == 1
    assert driver.events.count('cuDevicePrimaryCtxRelease_v2') == 1
    assert driver.local.stack == []
    kernel.close()
    assert driver.events.count('cuDevicePrimaryCtxRelease_v2') == 1


def test_load_failure_releases_retained_primary_context(mocked_runtime):
    driver=mocked_runtime
    driver.failures.update(cuModuleGetFunction='missing function',cuCtxSynchronize='async error')
    with pytest.raises(cuda.CudaError,match='missing function'):
        cuda.compile_cuda('extern "C" __global__ void test() {}','test')
    assert driver.events.count('cuModuleUnload') == 1
    assert driver.events.count('cuDevicePrimaryCtxRetain') == driver.events.count('cuDevicePrimaryCtxRelease_v2') == 1
    assert driver.local.stack == []


def test_context_stack_restored_after_body_error(mocked_runtime):
    driver=mocked_runtime
    with cuda.compile_cuda('extern "C" __global__ void test() {}','test') as kernel:
        driver.local.stack=[999]
        with pytest.raises(ValueError,match='body error'):
            with kernel._current():
                assert driver.local.stack == [999,123]
                raise ValueError('body error')
        assert driver.local.stack == [999]
        driver.local.stack=[]


def test_context_manager_preserves_body_error_during_cleanup(mocked_runtime):
    driver = mocked_runtime
    with pytest.raises(ValueError, match='original user error'):
        with cuda.compile_cuda('extern "C" __global__ void test() {}', 'test'):
            driver.failures['cuCtxSynchronize'] = 'cleanup error'
            raise ValueError('original user error')
    assert driver.events.count('cuModuleUnload') == 1
    assert driver.events.count('cuDevicePrimaryCtxRelease_v2') == 1


def test_stream_pointer_cannot_wrap(mocked_runtime):
    with cuda.compile_cuda('extern "C" __global__ void test() {}','test') as kernel:
        with pytest.raises(ValueError,match='stream'):
            kernel.launch([],stream=1 << (ct.sizeof(ct.c_void_p)*8))
        assert 'cuLaunchKernel' not in mocked_runtime.events


def test_host_call_attempts_all_buffer_releases_and_preserves_original_error(mocked_runtime):
    x=s.Symbol('x')
    f=cuda.cuda_lambdify(x,(x,x*x),cache_dir=False)
    driver=mocked_runtime
    driver.failures.update(cuMemcpyHtoD_v2='original upload error',cuMemFree_v2=['first free error',None,None])
    with pytest.raises(cuda.CudaError,match='original upload error'):
        f(np.arange(3.))
    assert driver.events.count('cuMemFree_v2') == 3
    f.close()


def test_host_call_reports_cleanup_error_after_releasing_all_buffers(mocked_runtime):
    x=s.Symbol('x')
    f=cuda.cuda_lambdify(x,(x,x*x),cache_dir=False)
    driver=mocked_runtime
    driver.failures['cuMemFree_v2']=['first free error',None,None]
    with pytest.raises(cuda.CudaError,match='first free error'):
        f(np.arange(3.))
    assert driver.events.count('cuMemFree_v2') == 3
    f.close()


def test_nvrtc_library_discovery_does_not_load_builtins(monkeypatch,tmp_path):
    directory=tmp_path/'lib64'
    directory.mkdir()
    builtins=directory/'libnvrtc-builtins.so.12'
    compiler=directory/'libnvrtc.so.12'
    builtins.touch(); compiler.touch()
    monkeypatch.delenv('SYMOX_NVRTC_LIBRARY',raising=False)
    monkeypatch.setenv('CUDA_HOME',str(tmp_path))
    monkeypatch.delenv('CUDA_PATH',raising=False)
    monkeypatch.setattr(cuda.ctypes.util,'find_library',lambda name: None)
    attempted=[]
    def load(name):
        attempted.append(name)
        if name in (str(compiler),str(builtins)): return name
        raise OSError('not installed')
    monkeypatch.setattr(cuda.ct,'CDLL',load)
    assert cuda._library('nvrtc') == str(compiler)
    assert str(builtins) not in attempted
