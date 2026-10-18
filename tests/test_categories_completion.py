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
    assert Morphism(A,B).domain == A
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
    source = xypic_draw_diagram(diagram)
    assert '@(r,u)' in source


def test_grids_rendering_masks_and_formatters():
    (A,B,C,D),(f,g,h) = arrows()
    diagram = Diagram([f,g],{g*f:'unique'})
    grid = DiagramGrid(diagram)
    assert grid.width == 3 and grid.height == 1
    assert [grid[0,i] for i in range(3)] == [A,B,C]
    transposed = DiagramGrid(diagram,transpose=True)
    assert transposed.width == 1 and transposed.height == 3
    drawer = XypicDiagramDrawer()
    drawer.arrow_formatters['unique'] = lambda arrow:setattr(arrow,'arrow_style','{-->}')
    text = drawer.draw(diagram,grid)
    assert r'\xymatrix' in text and r'\circ' in text and r'@{-->}' in text
    assert r'\circ' not in drawer.draw(diagram,grid,masked=[g*f])
    assert xypic_draw_diagram(diagram,diagram_format='@R=1cm').startswith(r'\xymatrix@R=1cm')
    assert xypic_draw_diagram(Diagram()) == '\\xymatrix{\n\n}'
    assert r'\ar[r]^{}' in xypic_draw_diagram(Diagram([Morphism(A,B)]))
    grouped = DiagramGrid(diagram,groups=[[B,A]])
    assert {grouped[i,j] for i in range(grouped.height) for j in range(grouped.width)} >= {A,B,C}
    with pytest.raises(NotImplementedError):
        DiagramGrid(diagram,unsupported=True)


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
