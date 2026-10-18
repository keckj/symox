"""Send presentation MathML to the optional GtkMathView executable."""
import shutil
import subprocess
import tempfile
from .mathml import mathml

def print_gtk(x,start_viewer=True):
    document='<math xmlns="http://www.w3.org/1998/Math/MathML">'+mathml(x,printer='presentation')+'</math>'
    with tempfile.NamedTemporaryFile('w',suffix='.mml',encoding='utf-8') as stream:
        stream.write(document);stream.flush()
        if start_viewer:
            viewer=shutil.which('mathmlviewer')
            if viewer is None:raise ImportError('GtkMathView requires the mathmlviewer executable')
            subprocess.check_call([viewer,stream.name])
