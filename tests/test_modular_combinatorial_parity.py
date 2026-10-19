"""Large held combinatorial expressions must use exact modular algorithms."""
import json
import subprocess
import sys

import pytest
from sympy.ntheory.residue_ntheory import binomial_mod as reference_binomial_mod
from symox.ntheory.residue_ntheory import binomial_mod


@pytest.mark.parametrize('n,m,k',[
    (10**18,10**12,100003**2),
    (156675,4433,100003),
    (10**30,10**15,2**32),
    (10**20+99,10**9+1,3**20),
    (10**20+999,10**10,5**12),
])
def test_arbitrary_integer_sizes_and_large_prime_powers(n,m,k):
    assert binomial_mod(n=n,m=m,k=k)==reference_binomial_mod(n=n,m=m,k=k)


@pytest.mark.parametrize('args',[(3.0,1,2),(3.5,1,2),('3',1,2),(True,1,2),
    (3,1,2.0),(3,1,True),(3.5,1,1),(3,1,1.0),(-1,1,2.0)])
def test_integer_input_protocol_does_not_truncate(args):
    try:expected=reference_binomial_mod(*args)
    except (TypeError,ValueError) as error:
        with pytest.raises(type(error)):binomial_mod(*args)
    else:assert binomial_mod(*args)==expected


def test_huge_held_documentation_examples_finish_without_reference_imports():
    source='''
import importlib.abc,json,sys
class BlockReference(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference imports are forbidden')
sys.meta_path.insert(0,BlockReference())
import symox as s
values=[s.Mod(s.Pow(2,10**16,evaluate=False),97),
        s.Mod(s.factorial(10**9,evaluate=False),10**9+9),
        s.Mod(s.binomial(10**18,10**12,evaluate=False),100003**2)]
huge=s.binomial(2**70,2**69)
assert huge.func.__name__=='binomial'
assert s.Mod(huge,4)==2
print(json.dumps([int(value) for value in values]))
'''
    result=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=20)
    assert result.returncode==0,result.stderr
    assert json.loads(result.stdout)==[61,712524808,3744312326]
