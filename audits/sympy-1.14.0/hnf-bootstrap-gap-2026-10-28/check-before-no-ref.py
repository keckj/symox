"""Supplement the frozen failing pytest cases with checked import/audit state."""
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

stage=Path(__file__).parent
tree=ast.parse((stage/'tests/test_native_hnf_function_bootstrap.py').read_text())
script=ast.literal_eval(next(node.value for node in tree.body if isinstance(node,ast.Assign)
                             and any(isinstance(target,ast.Name) and target.id=='SCRIPT' for target in node.targets)))
wrapper='''
import json,sys,hashlib
from pathlib import Path
namespace={}
program,observed_target,observed_raises=sys.argv[1:4]
sys.argv=[sys.argv[0],observed_target,observed_raises]
try:
    exec(compile(program, '<isolated-before>', 'exec'),namespace)
except AssertionError as error:
    assert namespace.get('actual_events') == []
    assert namespace.get('expected_events')
    assert namespace.get('audit_events') == []
    reference=any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
    assert reference is False
    binary=Path(namespace['sx']._core.__file__)
    print(json.dumps({'target':observed_target,'raises':observed_raises=='True',
                      'actual':namespace['actual'],'expected':namespace['expected'],
                      'actual_events':namespace['actual_events'],'expected_events':namespace['expected_events'],
                      'reference_imported':reference,'audit_events':namespace['audit_events'],
                      'native_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
                      'expected_regression_failure':str(error)}))
else:
    raise AssertionError('baseline unexpectedly fixed')
'''
records=[]
for target in ('_gcdex','add_columns'):
    for raises in (False,True):
        result=subprocess.run([sys.executable,'-c',wrapper,script,target,str(raises)],
                              env=os.environ.copy(),text=True,capture_output=True)
        assert result.returncode==0,result.stdout+result.stderr
        records.append(json.loads(result.stdout))
Path(sys.argv[1]).write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps({'cases':len(records),'checked_no_reference_import':True,'audit_events':[]}))
