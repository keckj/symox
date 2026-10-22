"""Native samples and actual headless renders for every plotting entry point."""
import builtins
import io
import numpy as np
import pytest
import symox as s
from symox.plotting import (plot,plot_parametric,plot3d,plot_contour,
    plot3d_parametric_line,plot3d_parametric_surface,plot_implicit,
    PlotGrid,PygletPlot,textplot,plot_backends)

x,y,t=s.symbols('x y t')


def test_scalar_lines_constants_and_composition():
    p=plot(x*x,3,(x,-2,2),n=9,show=False)
    xx,yy=p[0].get_data()
    np.testing.assert_allclose(yy,xx*xx)
    np.testing.assert_allclose(p[1].get_data()[1],3)
    q=plot((s.sin(x),(x,-1,1),'wave'),show=False,n=7)
    p.extend(q)
    assert p[2].label=='wave'
    p[1]=q[0]
    del p[0]
    assert len(p)==2
    with pytest.raises(TypeError):p.append(x)


def test_callable_and_parameter_substitution():
    p=plot(lambda z:z*z,(x,0,2),show=False,n=4)
    xx,yy=p[0].get_data()
    np.testing.assert_allclose(yy,xx*xx)
    p=plot(y*x,(x,0,2),params={y:2},show=False,n=4)
    xx,yy=p[0].get_data()
    np.testing.assert_allclose(yy,2*xx)


def test_parametric_circle():
    p=plot_parametric(s.cos(t),s.sin(t),(t,0,2*s.pi),n=17,show=False)
    xx,yy,tt=p[0].get_data()
    np.testing.assert_allclose(xx*xx+yy*yy,1,atol=1e-14)
    np.testing.assert_allclose(tt,np.linspace(0,2*np.pi,17))


def test_surface_orientation_and_parametric_coordinates():
    p=plot3d(x+2*y,(x,-1,1),(y,-2,2),n=(5,7),show=False)
    xx,yy,zz=p[0].get_data()
    assert zz.shape==(7,5)
    np.testing.assert_allclose(zz,xx+2*yy)
    p=plot3d_parametric_surface(x,y,x*y,(x,-1,1),(y,0,1),n=(4,6),show=False)
    xx,yy,zz=p[0].get_data()
    np.testing.assert_allclose(zz,xx*yy)


@pytest.mark.parametrize('relation,expected',[
    (s.Eq(x*x+y*y,1),lambda z:None),
    (x*x+y*y<1,lambda z:z<0),
    (x*x+y*y<=1,lambda z:z<=0),
    (x*x+y*y>1,lambda z:z>0),
    (x*x+y*y>=1,lambda z:z>=0),
])
def test_implicit_samples(relation,expected):
    p=plot_implicit(relation,(x,-2,2),(y,-2,2),n=(5,7),show=False)
    xx,yy,zz,mode=p[0].get_data()
    value=xx*xx+yy*yy-1
    if mode=='contour':np.testing.assert_allclose(zz,value)
    else:np.testing.assert_array_equal(zz,expected(value))


def test_nonreal_and_singular_samples_are_masked():
    p=plot(1/x,(x,-1,1),n=5,show=False)
    assert np.isnan(p[0].get_data()[1][2])
    p=plot(s.sqrt(x),(x,-1,1),n=5,show=False)
    assert np.all(np.isnan(p[0].get_data()[1][:2]))
    p=plot(1/x,(x,-1,1),n=100,detect_poles=True,show=False)
    assert np.any(np.isnan(p[0].get_data()[1]))


@pytest.mark.parametrize('factory,args',[
    (plot,(s.sin(x),(x,-2,2))),
    (plot_parametric,(s.cos(t),s.sin(t),(t,0,2*s.pi))),
    (plot3d,(x*y,(x,-1,1),(y,-1,1))),
    (plot_contour,(x*x+y*y,(x,-1,1),(y,-1,1))),
    (plot3d_parametric_line,(s.cos(t),s.sin(t),t,(t,0,2*s.pi))),
    (plot3d_parametric_surface,(x,y,x*y,(x,-1,1),(y,-1,1))),
    (plot_implicit,(s.Eq(x*x+y*y,1),(x,-2,2),(y,-2,2))),
    (plot_implicit,(x*x+y*y<1,(x,-2,2),(y,-2,2))),
    (PygletPlot,(x*x,(x,-2,2))),
])
def test_png_render(factory,args,tmp_path):
    mpl=pytest.importorskip('matplotlib')
    mpl.use('Agg',force=True)
    import matplotlib.pyplot as plt
    before=set(plt.get_fignums())
    p=factory(*args,n=15,show=False)
    path=tmp_path/'plot.png'
    p.save(path)
    assert path.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert path.stat().st_size>1000
    p.close()
    assert set(plt.get_fignums())==before


def test_plot_grid_mixed_dimensions(tmp_path):
    mpl=pytest.importorskip('matplotlib');mpl.use('Agg',force=True)
    p=plot(x,(x,0,1),n=5,show=False)
    q=plot3d(x*y,(x,0,1),(y,0,1),n=5,show=False)
    grid=PlotGrid(1,2,p,q,show=False,size=(8,3))
    grid.save(tmp_path/'grid.svg')
    assert len(grid.fig.axes)==2
    assert grid.fig.axes[1].name=='3d'
    grid.close()
    with pytest.raises(ValueError):PlotGrid(1,1,p,q,show=False)


def test_text_backend_and_no_symbolic_dependency(monkeypatch,capsys,tmp_path):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('plotting imported SymPy')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    textplot(x*x,-1,1,W=17,H=9)
    lines=capsys.readouterr().out.splitlines()
    assert len(lines)==9 and all(len(line)==17 for line in lines)
    assert '*' in ''.join(lines)
    p=plot(x,(x,-1,1),backend='text',show=False)
    p.save(tmp_path/'plot.txt')
    assert '*' in (tmp_path/'plot.txt').read_text()
    assert set(plot_backends)=={'matplotlib','text','default'}


@pytest.mark.parametrize('kwargs',[{'n':1},{'n':2.5},{'n':(5,6)}])
def test_invalid_sample_counts(kwargs):
    with pytest.raises(ValueError):plot(x,(x,0,1),show=False,**kwargs)


def test_range_validation():
    with pytest.raises(ValueError):plot(x,(x,0,s.oo),show=False)
    # The reference samples a degenerate range (Truss.draw relies on it).
    assert list(plot(x,(x,1,1),show=False)[0].get_data()[0][:3])==[1.0,1.0,1.0]
    with pytest.raises(ValueError):plot(x+y,(x,0,1),show=False)
    with pytest.raises(ValueError):plot3d(x,(x,0,1),(x,0,2),show=False)
    # Inequality regions combine as boolean masks; an equality has no area.
    assert plot_implicit(s.And(x>0,y>0),show=False)[0].get_data()[3]=='contourf'
    with pytest.raises(NotImplementedError):plot_implicit(s.And(s.Eq(x,0),y>0),show=False)


def test_samples_do_not_require_matplotlib(monkeypatch):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='matplotlib' or name.startswith('matplotlib.'):
            raise ImportError('deliberately unavailable')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    p=plot(x*x,(x,-1,1),show=False,n=5)
    assert len(p[0].get_data()[0])==5
    with pytest.raises(ImportError,match=r'symox\[plotting\]'):
        p.save(io.BytesIO())


def test_failed_render_releases_figure():
    mpl=pytest.importorskip('matplotlib');mpl.use('Agg',force=True)
    import matplotlib.pyplot as plt
    before=set(plt.get_fignums())
    p=plot(x,(x,-1,1),show=False,rendering_kw={'invalid_option_for_test':True})
    with pytest.raises(AttributeError):p.save(io.BytesIO())
    assert set(plt.get_fignums())==before
