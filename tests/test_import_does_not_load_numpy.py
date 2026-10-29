"""Importing symox, like importing SymPy 1.14, does not load NumPy."""
import os
import subprocess
import sys


def test_import_symox_does_not_load_numpy():
    code = 'import sys, symox; assert "numpy" not in sys.modules'
    subprocess.run([sys.executable, '-c', code], check=True, timeout=60,
                   env=dict(os.environ))
