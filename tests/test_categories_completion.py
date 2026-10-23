"""Composition laws, diagram properties and actual Xy-pic compilation."""
import shutil
import subprocess
import pytest
from symox.categories import *
from symox.categories.diagram_drawing import ArrowStringDescription


def arrows():
    A,B,C,D = map(Object,'ABCD')
    return (A,B,C,D),(NamedMorphism(A,B,'f'),NamedMorphism(B,C,'g'),NamedMorphism(C,D,'h'))


def test_associativity_units_and_normalization():
    (A,B,C,D),(f,g,h) = arrows()
    assert (h*g)*f == h*(g*f) == CompositeMorphism(f,g,h)
    assert CompositeMorphism([f,g]) == g*f
    assert IdentityMorphism(B)*f == f == f*IdentityMorphism(A)
    assert CompositeMorphism(IdentityMorphism(A),IdentityMorphism(A)) == IdentityMorphism(A)
    assert (g*f).components == (f,g)
    assert (g*f).domain == A and (g*f).codomain == C
    assert (g*f).flatten('gf') == NamedMorphism(A,C,'gf')
    with pytest.raises(NotImplementedError):
        Morphism(A,B)
    assert len({f,NamedMorphism(A,B,'f')}) == 1
    with pytest.raises(ValueError):
        f*g
    with pytest.raises(TypeError):
        CompositeMorphism()
    with pytest.raises(AttributeError):
        f.args = ()


def test_finite_diagram_closure_and_properties():
    (A,B,C,D),(f,g,h) = arrows()
    diagram = Diagram({h:['mono'],f:['mono','epi'],g:['mono']},{h*g*f:'unique'})
    assert all(IdentityMorphism(obj) in diagram.premises for obj in (A,B,C,D))
    assert h*g*f in diagram.premises
    assert diagram.premises[h*g*f] == frozenset({'mono'})
    assert diagram.conclusions[h*g*f] == frozenset({'unique'})
    premise,conclusion = diagram.hom(A,D)
    assert h*g*f in premise and h*g*f in conclusion
    smaller = Diagram({f:['mono','epi']})
    assert diagram.is_subdiagram(smaller)
    assert not smaller.is_subdiagram(diagram)
    assert diagram.subdiagram_from_objects([A,B]) == smaller
    endpoints = diagram.subdiagram_from_objects([A,D])
    assert endpoints.objects == frozenset((A,D)) and h*g*f in endpoints.premises
    with pytest.raises(ValueError):
        diagram.subdiagram_from_objects([Object('missing')])
    with pytest.raises(ValueError):
        Diagram({IdentityMorphism(A):'mono'})
    category = Category('K',[A,B,C,D],[diagram])
    assert category.name == 'K' and diagram in category.commutative_diagrams
    assert A in category.objects
    with pytest.raises(NotImplementedError):
        category.hom(A,B)


def test_loops_are_finite_and_explicit_words_retained():
    A,B = map(Object,'AB')
    f,g = NamedMorphism(A,B,'f'),NamedMorphism(B,A,'g')
    loop = NamedMorphism(A,A,'u')
    diagram = Diagram([f,g,loop,g*f])
    assert g*f in diagram.premises and loop in diagram.premises
    assert len(diagram.premises) < 20
    # Reference SymPy 1.14 output.
    assert xypic_draw_diagram(diagram) == (
        '\\xymatrix{\nA \\ar[r]^{f} \\ar@(u,l)[]^{u} & B \\ar@/^3mm/[l]^{g} \n}\n')


def test_grids_rendering_masks_and_formatters():
    (A,B,C,D),(f,g,h) = arrows()
    diagram = Diagram([f,g],{g*f:'unique'})
    grid = DiagramGrid(diagram)
    # Reference SymPy 1.14 generic layout: a triangle on a 2x2 grid.
    assert grid.width == 2 and grid.height == 2
    assert [[grid[i,j] for j in range(2)] for i in range(2)] == [[A,B],[C,None]]
    transposed = DiagramGrid(diagram,transpose=True)
    assert transposed.width == 2 and transposed.height == 2
    drawer = XypicDiagramDrawer()
    drawer.arrow_formatters['unique'] = lambda arrow:setattr(arrow,'arrow_style','{-->}')
    text = drawer.draw(diagram,grid)
    assert text == '\\xymatrix{\nA \\ar@{-->}[d]_{g\\circ f} \\ar[r]^{f} & B \\ar[ld]^{g} \\\\\nC & \n}\n'
    assert drawer.draw(diagram,grid,masked=[g*f]) == (
        '\\xymatrix{\nA \\ar[r]^{f} & B \\ar[ld]^{g} \\\\\nC & \n}\n')
    assert xypic_draw_diagram(diagram,diagram_format='@R=1cm').startswith(r'\xymatrix@R=1cm')
    # SymPy 1.14 fails on an empty diagram and on groups missing an object.
    with pytest.raises(IndexError):
        xypic_draw_diagram(Diagram())
    with pytest.raises(KeyError):
        DiagramGrid(diagram,groups=[[B,A]])
    grouped = DiagramGrid(diagram,groups=[[B,A],C])
    assert {grouped[i,j] for i in range(grouped.height) for j in range(grouped.width)} >= {A,B,C}


@pytest.mark.skipif(not all(shutil.which(tool) for tool in ('latex','dvipng','kpsewhich')),reason='TeX rendering tools unavailable')
def test_actual_diagram_preview(tmp_path):
    probe = subprocess.run(['kpsewhich','xypic.sty'],capture_output=True,text=True)
    if probe.returncode:
        pytest.skip('Xy-pic unavailable')
    (A,B,C,D),(f,g,h) = arrows()
    parallel = NamedMorphism(A,B,'k')
    diagram = Diagram([f,g,parallel,NamedMorphism(C,C,'u')],{g*f:'unique'})
    destination = tmp_path/'diagram.png'
    preview_diagram(diagram,viewer='file',filename=str(destination),euler=False,transpose=True)
    assert destination.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
