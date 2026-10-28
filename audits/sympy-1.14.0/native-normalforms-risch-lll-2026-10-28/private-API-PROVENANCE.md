The optional getter layer resolves public CPython exports at runtime. It never reads a private structure or field offset. Its verified admission is deliberately restricted to CPython3.10.12 and3.10.20 on Linuxx86_64. Other interpreters/platforms/versions or missing exports return None, keeping the complete historical arithmetic path.

The Python3.10 function C-API documentation specifies PyFunction_Type and the borrowed-return PyFunction_GetCode, PyFunction_GetGlobals, PyFunction_GetDefaults and PyFunction_GetClosure APIs:
https://docs.python.org/3.10/c-api/function.html

PyFunction_GetKwDefaults is public PyAPI_FUNC in the exactv3.10.12 andv3.10.20 headers, although absent from that documentation page. The header puts these exports outside Py_LIMITED_API:
https://raw.githubusercontent.com/python/cpython/v3.10.12/Include/funcobject.h
https://raw.githubusercontent.com/python/cpython/v3.10.20/Include/funcobject.h

The exactv3.10.12 implementation returns the corresponding borrowed field for these five getters, with function-type validation and no audit event:
https://raw.githubusercontent.com/python/cpython/v3.10.12/Objects/funcobject.c
https://raw.githubusercontent.com/python/cpython/v3.10.20/Objects/funcobject.c

The DSO retains its abi3-py39 build configuration; this optional, restricted optimization is not a claim that these symbols belong to the LimitedAPI or work on every abi3 interpreter. dlsym resolves the public function-type address and getters only after the verified version/platform check. No lambda is compiled/evaluated to obtain FunctionType for this new layer. Only immutable API addresses/type identity are retained; each admission reads the function's current code/globals/defaults/keyword-defaults/closure afresh, so mutations still trigger fallback.

Tests verify five exact identities, changed code/default/keyword-default and closure contents, foreign objects without attribute callbacks, refusal of missing/replaced metadata functions, three foreign objects injected before helper import, native HNF/LU/solve admission under rejecting metadata audit hooks, and complete old-body results on refusal. This restricted optional route may leave acceleration unavailable elsewhere; mathematical domains, exactness, errors and completeness use the existing fallback.
