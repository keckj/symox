"""Native printing helpers adapted from SymPy 1.14; see LICENSE-SYMPY."""
import os
from os.path import join
import shutil
import tempfile
from pathlib import Path
try:
    from subprocess import STDOUT, CalledProcessError, check_output
except ImportError:
    pass

def debug(*args):
    pass
from .latex import latex
__doctest_requires__ = {('preview',): ['pyglet']}

def _check_output_no_window(*args, **kwargs):
    if os.name == 'nt':
        creation_flag = 134217728
    else:
        creation_flag = 0
    return check_output(*args, creationflags=creation_flag, **kwargs)

def system_default_viewer(fname, fmt):
    """Open fname with the default system viewer."""
    with tempfile.NamedTemporaryFile(prefix='sympy-preview-', suffix=os.path.splitext(fname)[1], delete=False) as temp_f:
        with open(fname, 'rb') as f:
            shutil.copyfileobj(f, temp_f)
    import platform
    if platform.system() == 'Darwin':
        import subprocess
        subprocess.call(('open', temp_f.name))
    elif platform.system() == 'Windows':
        os.startfile(temp_f.name)
    else:
        import subprocess
        subprocess.call(('xdg-open', temp_f.name))

def pyglet_viewer(fname, fmt):
    try:
        from pyglet import window, image, gl
        from pyglet.window import key
        from pyglet.image.codecs import ImageDecodeException
    except ImportError:
        raise ImportError('pyglet is required for preview.\n visit https://pyglet.org/')
    try:
        img = image.load(fname)
    except ImageDecodeException:
        raise ValueError("pyglet preview does not work for '{}' files.".format(fmt))
    offset = 25
    config = gl.Config(double_buffer=False)
    win = window.Window(width=img.width + 2 * offset, height=img.height + 2 * offset, caption='SymPy', resizable=False, config=config)
    win.set_vsync(False)
    try:

        def on_close():
            win.has_exit = True
        win.on_close = on_close

        def on_key_press(symbol, modifiers):
            if symbol in [key.Q, key.ESCAPE]:
                on_close()
        win.on_key_press = on_key_press

        def on_expose():
            gl.glClearColor(1.0, 1.0, 1.0, 1.0)
            gl.glClear(gl.GL_COLOR_BUFFER_BIT)
            img.blit((win.width - img.width) / 2, (win.height - img.height) / 2)
        win.on_expose = on_expose
        while not win.has_exit:
            win.dispatch_events()
            win.flip()
    except KeyboardInterrupt:
        pass
    win.close()

def _get_latex_main(expr, *, preamble=None, packages=(), extra_preamble=None, euler=True, fontsize=None, **latex_settings):
    """Generate string of a LaTeX document rendering ``expr``."""
    if preamble is None:
        actual_packages = packages + ('amsmath', 'amsfonts')
        if euler:
            actual_packages += ('euler',)
        package_includes = '\n' + '\n'.join(['\\usepackage{%s}' % p for p in actual_packages])
        if extra_preamble:
            package_includes += extra_preamble
        if not fontsize:
            fontsize = '12pt'
        elif isinstance(fontsize, int):
            fontsize = '{}pt'.format(fontsize)
        preamble = '\\documentclass[varwidth,%s]{standalone}\n%s\n\n\\begin{document}\n' % (fontsize, package_includes)
    elif packages or extra_preamble:
        raise ValueError('The "packages" or "extra_preamble" keywordsmust not be set if a custom LaTeX preamble was specified')
    if isinstance(expr, str):
        latex_string = expr
    else:
        latex_string = '$\\displaystyle ' + latex(expr, mode='plain', **latex_settings) + '$'
    return preamble + '\n' + latex_string + '\n\n' + '\\end{document}'

def preview(expr, output='png', viewer=None, euler=True, packages=(), filename=None, outputbuffer=None, preamble=None, dvioptions=None, outputTexFile=None, extra_preamble=None, fontsize=None, **latex_settings):
    """View expression or LaTeX markup in PNG, DVI, PostScript or PDF form."""
    if viewer is None and output == 'png':
        try:
            import pyglet
        except ImportError:
            pass
        else:
            viewer = pyglet_viewer
    if viewer is None:
        candidates = {'dvi': ['evince', 'okular', 'kdvi', 'xdvi'], 'ps': ['evince', 'okular', 'gsview', 'gv'], 'pdf': ['evince', 'okular', 'kpdf', 'acroread', 'xpdf', 'gv']}
        for candidate in candidates.get(output, []):
            path = shutil.which(candidate)
            if path is not None:
                viewer = path
                break
    if viewer is None:
        viewer = system_default_viewer
    if viewer == 'file':
        if filename is None:
            raise ValueError('filename has to be specified if viewer="file"')
    elif viewer == 'BytesIO':
        if outputbuffer is None:
            raise ValueError('outputbuffer has to be a BytesIO compatible object if viewer="BytesIO"')
    elif not callable(viewer) and (not shutil.which(viewer)):
        raise OSError('Unrecognized viewer: %s' % viewer)
    latex_main = _get_latex_main(expr, preamble=preamble, packages=packages, euler=euler, extra_preamble=extra_preamble, fontsize=fontsize, **latex_settings)
    debug('Latex code:')
    debug(latex_main)
    with tempfile.TemporaryDirectory() as workdir:
        Path(join(workdir, 'texput.tex')).write_text(latex_main, encoding='utf-8')
        if outputTexFile is not None:
            shutil.copyfile(join(workdir, 'texput.tex'), outputTexFile)
        if not shutil.which('latex'):
            raise RuntimeError('latex program is not installed')
        try:
            _check_output_no_window(['latex', '-halt-on-error', '-interaction=nonstopmode', 'texput.tex'], cwd=workdir, stderr=STDOUT)
        except CalledProcessError as e:
            raise RuntimeError("'latex' exited abnormally with the following output:\n%s" % e.output)
        src = 'texput.%s' % output
        if output != 'dvi':
            commandnames = {'ps': ['dvips'], 'pdf': ['dvipdfmx', 'dvipdfm', 'dvipdf'], 'png': ['dvipng'], 'svg': ['dvisvgm']}
            try:
                cmd_variants = commandnames[output]
            except KeyError:
                raise ValueError('Invalid output format: %s' % output) from None
            for cmd_variant in cmd_variants:
                cmd_path = shutil.which(cmd_variant)
                if cmd_path:
                    cmd = [cmd_path]
                    break
            else:
                if len(cmd_variants) > 1:
                    raise RuntimeError('None of %s are installed' % ', '.join(cmd_variants))
                else:
                    raise RuntimeError('%s is not installed' % cmd_variants[0])
            defaultoptions = {'dvipng': ['-T', 'tight', '-z', '9', '--truecolor'], 'dvisvgm': ['--no-fonts']}
            commandend = {'dvips': ['-o', src, 'texput.dvi'], 'dvipdf': ['texput.dvi', src], 'dvipdfm': ['-o', src, 'texput.dvi'], 'dvipdfmx': ['-o', src, 'texput.dvi'], 'dvipng': ['-o', src, 'texput.dvi'], 'dvisvgm': ['-o', src, 'texput.dvi']}
            if dvioptions is not None:
                cmd.extend(dvioptions)
            else:
                cmd.extend(defaultoptions.get(cmd_variant, []))
            cmd.extend(commandend[cmd_variant])
            try:
                _check_output_no_window(cmd, cwd=workdir, stderr=STDOUT)
            except CalledProcessError as e:
                raise RuntimeError("'%s' exited abnormally with the following output:\n%s" % (' '.join(cmd), e.output))
        if viewer == 'file':
            shutil.move(join(workdir, src), filename)
        elif viewer == 'BytesIO':
            s = Path(join(workdir, src)).read_bytes()
            outputbuffer.write(s)
        elif callable(viewer):
            viewer(join(workdir, src), fmt=output)
        else:
            try:
                _check_output_no_window([viewer, src], cwd=workdir, stderr=STDOUT)
            except CalledProcessError as e:
                raise RuntimeError("'%s %s' exited abnormally with the following output:\n%s" % (viewer, src, e.output))
