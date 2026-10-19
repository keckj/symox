"""Reference keyword spellings must execute the corresponding operation."""
import pytest
import sympy as sp
import symox as sx

GEOMETRY = [
 ('Circle((0,0),1)','contains','other','Point(1,0)'),
 ('Circle((0,0),1)','encloses_point','p','Point(0,0)'),
 ('Circle((0,0),1)','equals','o','Circle((0,0),1)'),
 ('Circle((0,0),1)','intersection','o','Line((-2,0),(2,0))'),
 ('Circle((0,0),1)','is_tangent','o','Line((-1,1),(1,1))'),
 ('Circle((0,0),1)','tangent_lines','p','Point(2,0)'),
 ('Ellipse((0,0),2,1)','intersection','o','Line((-3,0),(3,0))'),
 ('Line((0,0),(1,0))','parallel_line','p','Point(1,1)'),
 ('Line((0,0),(1,0))','perpendicular_line','p','Point(1,1)'),
 ('Line((0,0),(1,0))','perpendicular_segment','p','Point(2,2)'),
 ('Segment((0,0),(2,0))','perpendicular_bisector','p','Point(1,2)'),
 ('Parabola(Point(0,1),Line((-1,-1),(1,-1)))','intersection','o','Line((-1,0),(1,0))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','angle_between','o','Plane((0,0,0),normal_vector=(1,0,0))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','distance','o','Point(0,0,3)'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','equals','o','Plane((1,1,0),normal_vector=(0,0,1))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','intersection','o','Plane((0,0,0),normal_vector=(1,0,0))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','is_parallel','l','Plane((0,0,1),normal_vector=(0,0,1))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','is_perpendicular','l','Plane((0,0,0),normal_vector=(1,0,0))'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','perpendicular_line','pt','Point(1,2,3)'),
 ('Plane((0,0,0),normal_vector=(0,0,1))','projection','pt','Point(1,2,3)'),
 ('Point(1,2)','dot','p','Point(2,3)'),
 ('Point(1,2)','midpoint','p','Point(2,3)'),
 ('Point(1,2)','taxicab_distance','p','Point(2,3)'),
 ('Polygon((0,0),(4,0),(4,4),(0,4))','distance','o','Point(5,1)'),
 ('Polygon((0,0),(4,0),(4,4),(0,4))','encloses','o','Point(1,1)'),
 ('Polygon((0,0),(4,0),(4,4),(0,4))','encloses_point','p','Point(1,1)'),
 ('Polygon((0,0),(4,0),(4,4),(0,4))','equals','o','Polygon((0,0),(4,0),(4,4),(0,4))'),
 ('Polygon((0,0),(4,0),(4,4),(0,4))','intersection','o','Line((-1,1),(5,1))'),
]


def converted(value):
 if isinstance(value,(list,tuple)):
  return type(value)(converted(v) for v in value)
 return value.to_sympy() if hasattr(value,'to_sympy') else value


def same_geometry(actual, expected):
 # Lines may use different defining points, segments opposite endpoint order,
 # and intersection lists different enumeration order. Compare their geometry.
 if isinstance(actual,list) and isinstance(expected,list):
  remaining=list(expected)
  for item in actual:
   for index,other in enumerate(remaining):
    if same_geometry(item,other):
     remaining.pop(index)
     break
   else:return False
  return not remaining
 return actual==expected or bool(getattr(actual,'equals',lambda other:False)(expected))


@pytest.mark.parametrize('owner,method,keyword,value', GEOMETRY)
def test_geometry_keyword_operations_match_reference(owner, method, keyword, value):
 actual=getattr(eval(owner,vars(sx)),method)(**{keyword:eval(value,vars(sx))})
 expected=getattr(eval(owner,vars(sp)),method)(**{keyword:eval(value,vars(sp))})
 assert same_geometry(converted(actual),expected)


@pytest.mark.parametrize('keyword', ['row','col'])
def test_sparse_deletion_keyword_mutates_the_requested_axis(keyword):
 method=keyword+'_del'
 actual=sx.MutableSparseMatrix([[1,0],[2,3]])
 expected=sp.MutableSparseMatrix([[1,0],[2,3]])
 getattr(actual,method)(**{keyword:0})
 getattr(expected,method)(**{keyword:0})
 assert actual.to_sympy()==expected


def test_array_and_sparse_callbacks_and_replacements():
 x=sx.Symbol('x');X=sp.Symbol('x')
 for cls,ref in [(sx.Array,sp.Array),(sx.ImmutableSparseMatrix,sp.ImmutableSparseMatrix)]:
  actual,expected=cls([[x,0],[2,3]]),ref([[X,0],[2,3]])
  if cls is sx.Array:
   assert actual.applyfunc(f=lambda a:a+1).to_sympy()==expected.applyfunc(f=lambda a:a+1)
  else:
   # SymPy 1.14's immutable sparse applyfunc skips implicit zero entries.
   # Verify the stated entrywise operation independently instead of copying it.
   assert actual.applyfunc(f=lambda a:a+1).to_sympy()==sp.ImmutableSparseMatrix(
       expected.rows,expected.cols,[a+1 for a in expected])
  assert actual.xreplace(rule={x:4}).to_sympy()==expected.xreplace(rule={X:4})
 matrix=sx.ImmutableSparseMatrix([[1,2],[3,4]])
 assert matrix.extract(rowsList=[1],colsList=[0])==sx.ImmutableSparseMatrix([[3]])


def test_ordered_set_reference_keyword_names():
 from symox.core.containers import OrderedSet
 values=OrderedSet([1,2])
 values.add(key=3);values.discard(key=2)
 assert list(values)==[1,3]


def test_parabola_contains_accepts_reference_keyword_and_checks_equation():
 parabola=sx.Parabola(sx.Point(0,1),sx.Line((-1,-1),(1,-1)))
 assert parabola.contains(other=sx.Point(0,0))
 assert not parabola.contains(other=sx.Point(0,1))
