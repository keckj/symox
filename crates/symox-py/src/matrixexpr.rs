//! Python interface to the native typed MatrixExpr DAG.
use pyo3::prelude::*;
use pyo3::types::{PyDict,PyList,PySet,PyTuple,PySlice};
use pyo3::exceptions::{PyTypeError,PyNotImplementedError,PyValueError,PyIndexError};
use rustc_hash::FxHashMap;
use symox::expr::Expr;
use symox::matrixexpr::{MatrixExpr,MatrixKind,MatrixOperand,BlockGrid,SliceSpec};
use crate::{wrap,to_expr,PyExpr};

#[pyclass(name="_MatrixExpr",module="symox._core",frozen)]
#[derive(Clone)]
pub(crate) struct PyMatrixExpr {pub m:MatrixExpr}
#[pyclass(name="_SliceSpec",module="symox._core",frozen)]
#[derive(Clone)]
pub(crate) struct PySliceSpec{pub spec:SliceSpec}
#[pymethods]
impl PySliceSpec{
    #[getter] fn args<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{PyTuple::new(py,self.spec.args().into_iter().map(|e|e.map(|e|wrap(e).into_pyobject(py).unwrap().into_any().unbind()).unwrap_or_else(||py.None())))}
    fn __reduce__<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyAny>,Bound<'py,PyTuple>)>{Ok((py.import("symox._core")?.getattr("_normalized_slice")?,(wrap(self.spec.start.clone()),wrap(self.spec.stop.clone()),wrap(self.spec.step.clone()),wrap(self.spec.length.clone()),self.spec.stop_is_omitted).into_pyobject(py)?))}
}
fn slice_object<'py>(py:Python<'py>,spec:&SliceSpec)->PyResult<Bound<'py,PyAny>>{py.import("symox._matrixexpr")?.getattr("_SliceTuple")?.call1((PySliceSpec{spec:spec.clone()},))}
fn slice_index(value:&Bound<'_,PyAny>)->PyResult<Expr>{
    if value.is_instance_of::<PyExpr>()||value.get_type().module().is_ok_and(|m|m.to_string().starts_with("sympy")){return to_expr(value);}
    to_expr(&value.py().import("operator")?.getattr("index")?.call1((value,))?)
}
fn slice_spec(value:&Bound<'_,PyAny>,size:&Expr)->PyResult<SliceSpec>{
    if let Ok(spec)=value.getattr("_symox_slice_spec"){if let Ok(spec)=spec.downcast::<PySliceSpec>(){return Ok(spec.get().spec.clone());}}
    if !value.is_instance_of::<PySlice>()&&value.downcast::<PyTuple>().is_err()&&value.downcast::<PyList>().is_err(){return SliceSpec::index(slice_index(value)?,size).map_err(|e|error(value.py(),e));}
    SliceSpec::from_raw(slice_argument(value,size)?,size).map_err(|e|error(value.py(),e))
}
#[pyfunction]
fn _normalized_slice(start:&Bound<'_,PyAny>,stop:&Bound<'_,PyAny>,step:&Bound<'_,PyAny>,length:&Bound<'_,PyAny>,omitted:bool)->PyResult<PySliceSpec>{Ok(PySliceSpec{spec:SliceSpec::normalized(to_expr(start)?,to_expr(stop)?,to_expr(step)?,to_expr(length)?,omitted).map_err(|e|error(start.py(),e))?})}
#[pyclass(name="_BlockGrid",module="symox._core",frozen)]
#[derive(Clone)]
pub(crate) struct PyBlockGrid{pub(crate) grid:BlockGrid}
#[pymethods]
impl PyBlockGrid{
    #[new]
    fn new(rows:usize,cols:usize,values:&Bound<'_,PyAny>)->PyResult<Self>{
        let entries=values.try_iter()?.map(|v|to_matrixexpr(&v?)).collect::<PyResult<Vec<_>>>()?;
        Ok(Self{grid:BlockGrid::new(rows,cols,entries).map_err(|e|error(values.py(),e))?})
    }
    #[getter] fn _symox_blockgrid(&self)->bool{true}
    #[getter] fn is_Matrix(&self)->bool{true}
    #[getter] fn is_Immutable(&self)->bool{true}
    #[getter] fn shape(&self)->(usize,usize){self.grid.shape()}
    #[getter] fn rows(&self)->usize{self.grid.shape().0}
    #[getter] fn cols(&self)->usize{self.grid.shape().1}
    fn __len__(&self)->usize{self.grid.data().len()}
    fn __getitem__<'py>(&self,key:&Bound<'py,PyAny>)->PyResult<Bound<'py,PyAny>>{
        let position=|value:&Bound<'_,PyAny>,size:usize|->PyResult<usize>{let index:i64=value.extract()?;let index=if index<0{size as i64+index}else{index};usize::try_from(index).ok().filter(|i|*i<size).ok_or_else(||PyIndexError::new_err("Block grid index out of range"))};
        let entry=if let Ok(pair)=key.downcast::<PyTuple>(){if pair.len()!=2{return Err(PyIndexError::new_err("Block grid requires two indices"));}self.grid.get(position(&pair.get_item(0)?,self.rows())?,position(&pair.get_item(1)?,self.cols())?)}else{&self.grid.data()[position(key,self.__len__())?]};
        matrix_object(key.py(),entry)
    }
    fn __iter__<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{PyList::new(py,self.grid.data().iter().map(|m|matrix_object(py,m)).collect::<PyResult<Vec<_>>>()?)?.call_method0("__iter__")}
    fn tolist<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyList>>{
        let rows=(0..self.rows()).map(|i|PyList::new(py,(0..self.cols()).map(|j|matrix_object(py,self.grid.get(i,j))).collect::<PyResult<Vec<_>>>()?)).collect::<PyResult<Vec<_>>>()?;
        PyList::new(py,rows)
    }
    fn transpose(&self)->Self{Self{grid:self.grid.transpose()}}
    #[getter] fn T(&self)->Self{self.transpose()}
    #[getter] fn free_symbols<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PySet>>{
        let result=PySet::empty(py)?;for matrix in self.grid.data(){for symbol in matrix.matrix_symbols(){result.add(wrap_matrixexpr(symbol))?;}for symbol in matrix.scalar_free_symbols(){result.add(crate::expr_object(py,symbol)?)?;}}Ok(result)
    }
    fn __add__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{let other=other.downcast::<Self>()?;Ok(Self{grid:self.grid.add(&other.get().grid).map_err(|e|error(other.py(),e))?})}
    fn __neg__(&self,py:Python<'_>)->PyResult<Self>{Ok(Self{grid:self.grid.scale(&Expr::minus_one()).map_err(|e|error(py,e))?})}
    fn __sub__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{let other=other.downcast::<Self>()?;let negative=other.get().grid.scale(&Expr::minus_one()).map_err(|e|error(other.py(),e))?;Ok(Self{grid:self.grid.add(&negative).map_err(|e|error(other.py(),e))?})}
    fn __mul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{
        let grid=if let Ok(other)=other.downcast::<Self>(){self.grid.matmul(&other.get().grid)}else{if is_matrix(other){return Err(PyTypeError::new_err("Multiply block grids with another block grid or a scalar"));}self.grid.scale(&to_expr(other)?)};
        Ok(Self{grid:grid.map_err(|e|error(other.py(),e))?})
    }
    fn __matmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{let other=other.downcast::<Self>()?;Ok(Self{grid:self.grid.matmul(&other.get().grid).map_err(|e|error(other.py(),e))?})}
    fn __rmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if is_matrix(other){return Err(PyTypeError::new_err("Multiply block grids with another block grid or a scalar"));}Ok(Self{grid:self.grid.scale(&to_expr(other)?).map_err(|e|error(other.py(),e))?})}
    fn __richcmp__(&self,other:&Bound<'_,PyAny>,op:pyo3::basic::CompareOp)->PyObject{
        let equal=other.downcast::<Self>().is_ok_and(|v|v.get().grid==self.grid);let py=other.py();match op{pyo3::basic::CompareOp::Eq=>equal.into_pyobject(py).unwrap().to_owned().into_any().unbind(),pyo3::basic::CompareOp::Ne=>(!equal).into_pyobject(py).unwrap().to_owned().into_any().unbind(),_=>py.NotImplemented()}
    }
    fn __hash__(&self)->isize{use std::hash::{Hash,Hasher};let mut hash=std::collections::hash_map::DefaultHasher::new();self.grid.hash(&mut hash);let result=hash.finish() as isize;if result== -1{-2}else{result}}
    fn __repr__(&self)->String{format!("BlockGrid({:?}, [{}])",self.shape(),self.grid.data().iter().map(MatrixExpr::to_str).collect::<Vec<_>>().join(", "))}
    fn to_sympy<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{
        let entries=PyList::new(py,self.grid.data().iter().map(|m|matrixexpr_to_sympy(py,m)).collect::<PyResult<Vec<_>>>()?)?;
        let options=PyDict::new(py);options.set_item("evaluate",false)?;
        py.import("sympy")?.getattr("ImmutableMatrix")?.call((self.rows(),self.cols(),entries),Some(&options))
    }
    fn _sympy_<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{self.to_sympy(py)}
    fn __reduce__<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyAny>,Bound<'py,PyTuple>)>{
        crate::matrix_dag::reduce_grid(py,&self.grid)
    }
}
pub(crate) fn wrap_matrixexpr(m:MatrixExpr)->PyMatrixExpr{PyMatrixExpr{m}}
fn error(py:Python<'_>,e:symox::matrix::MatrixError)->PyErr{crate::matrices::err(py,e)}
pub(crate) fn matrix_object<'py>(py:Python<'py>,m:&MatrixExpr)->PyResult<Bound<'py,PyAny>>{
    if let MatrixKind::Explicit(value)=m.kind(){return Ok(crate::matrices::wrapm(value.clone()).into_pyobject(py)?.call_method0("as_immutable")?);}
    Ok(wrap_matrixexpr(m.clone()).into_pyobject(py)?.into_any())
}
fn is_matrix(value:&Bound<'_,PyAny>)->bool{
    value.downcast::<PyMatrixExpr>().is_ok() || value.downcast::<crate::matrices::PyMatrix>().is_ok()
        || value.getattr("is_Matrix").is_ok_and(|v|v.is_truthy().unwrap_or(false))
        || value.getattr("_symox_sparse").is_ok_and(|v|v.is_truthy().unwrap_or(false))
}
fn slice_argument(value:&Bound<'_,PyAny>,size:&Expr)->PyResult<[Option<Expr>;3]>{
    if value.is_instance_of::<PySlice>(){
        let field=|name:&str|->PyResult<Option<Expr>>{let item=value.getattr(name)?;if item.is_none(){Ok(None)}else{slice_index(&item).map(Some)}};
        return Ok([field("start")?,field("stop")?,field("step")?]);
    }
    if value.downcast::<PyTuple>().is_ok()||value.downcast::<PyList>().is_ok(){
        let parts=value.try_iter()?.collect::<PyResult<Vec<_>>>()?;
        if parts.len()!=2&&parts.len()!=3{return Err(PyIndexError::new_err("A matrix slice needs two or three bounds"));}
        let mut result=[None,None,None];for (i,item) in parts.iter().enumerate(){if !item.is_none(){result[i]=Some(slice_index(item)?);}}
        return Ok(result);
    }
    let mut index=slice_index(value)?;
    MatrixExpr::validate_index(&index,size).map_err(|e|error(value.py(),e))?;
    if symox::assumptions::is_negative(&index)==Some(true){index=symox::ops::add2(size,&index);}
    Ok([Some(index.clone()),Some(symox::ops::add2(&index,&Expr::one())),Some(Expr::one())])
}
fn to_block_grid(value:&Bound<'_,PyAny>)->PyResult<BlockGrid>{
    if let Ok(grid)=value.downcast::<PyBlockGrid>(){return Ok(grid.get().grid.clone());}
    let value=if value.hasattr("tolist")?{value.call_method0("tolist")?}else{value.clone()};
    if is_matrix(&value){return Err(PyTypeError::new_err("BlockMatrix requires a grid of matrix-valued cells"));}
    let items=value.try_iter()?.collect::<PyResult<Vec<_>>>()?;
    if items.is_empty(){return Err(PyValueError::new_err("BlockMatrix requires a nonempty grid"));}
    if items.iter().all(is_matrix){return BlockGrid::new(1,items.len(),items.iter().map(to_matrixexpr).collect::<PyResult<_>>()?).map_err(|e|error(value.py(),e));}
    let mut cols=None;let mut cells=Vec::new();for row in &items{
        let row=row.try_iter()?.collect::<PyResult<Vec<_>>>()?;
        if row.is_empty()||cols.is_some_and(|n|n!=row.len()){return Err(PyValueError::new_err("Block matrix rows must be nonempty and have equal lengths"));}
        cols=Some(row.len());for cell in row{cells.push(to_matrixexpr(&cell)?);}
    }BlockGrid::new(items.len(),cols.unwrap(),cells).map_err(|e|error(value.py(),e))
}
pub(crate) fn to_matrixexpr(value:&Bound<'_,PyAny>)->PyResult<MatrixExpr>{
    if value.is_instance_of::<PyBlockGrid>(){return Err(PyTypeError::new_err("A block grid stores matrix-valued cells; use a block matrix constructor"));}
    if let Ok(v)=value.downcast::<PyMatrixExpr>(){return Ok(v.get().m.clone());}
    if value.downcast::<crate::matrices::PyMatrix>().is_ok() || value.getattr("_symox_sparse").is_ok_and(|v|v.is_truthy().unwrap_or(false)){
        return Ok(MatrixExpr::explicit(crate::matrices::to_matrix(value)?));
    }
    if value.getattr("is_MatrixExpr").is_ok_and(|v|v.is_truthy().unwrap_or(false)){
        let head:String=value.getattr("func")?.getattr("__name__")?.extract()?;
        let args=value.getattr("args")?;
        let args=args.downcast::<PyTuple>()?;
        if head=="MatrixSlice"{
            let parent=to_matrixexpr(&args.get_item(0)?)?;
            let normalized=|index:usize,name:&str|->PyResult<SliceSpec>{let values=args.get_item(index)?;let a=to_expr(&values.get_item(0)?)?;let h=to_expr(&values.get_item(2)?)?;let length=to_expr(&value.getattr(name)?)?;let b=symox::ops::add2(&a,&symox::ops::mul2(&h,&length));let omitted=symox::assumptions::is_negative(&h)==Some(true)&&symox::assumptions::is_negative(&b)==Some(true);SliceSpec::normalized(a,b,h,length,omitted).map_err(|e|error(value.py(),e))};
            return parent.slice_normalized(normalized(1,"rows")?,normalized(2,"cols")?).map_err(|e|error(value.py(),e));
        }
        if head=="MatrixSymbol"{
            let name:String=value.getattr("name")?.extract()?;
            return MatrixExpr::symbol(&name,to_expr(&value.getattr("rows")?)?,to_expr(&value.getattr("cols")?)?).map_err(|e|error(value.py(),e));
        }
        return construct(&head,args,false).and_then(|result|{
            result.bind(value.py()).downcast::<PyMatrixExpr>().map(|v|v.get().m.clone())
                .map_err(|_|PyTypeError::new_err("Expected a matrix expression"))
        });
    }
    if is_matrix(value){return Ok(MatrixExpr::explicit(crate::matrices::to_matrix(value)?));}
    Err(PyTypeError::new_err("Expected a matrix expression or explicit matrix"))
}
pub(crate) fn matrixexpr_to_sympy<'py>(py:Python<'py>,m:&MatrixExpr)->PyResult<Bound<'py,PyAny>>{
    let sympy=py.import("sympy")?;
    if let MatrixKind::Explicit(value)=m.kind(){return matrix_object(py,&MatrixExpr::explicit(value.clone()))?.call_method0("to_sympy");}
    if let MatrixKind::Conjugate(parent)=m.kind(){return matrixexpr_to_sympy(py,parent)?.call_method0("conjugate");}
    let args=matrix_args(py,m,true)?;
    if matches!(m.kind(),MatrixKind::ApplyFunction(..)){
        let head=py.import("sympy.matrices.expressions.applyfunc")?.getattr("ElementwiseApplyFunction")?;
        let mut held=vec![head];held.extend(args.iter());
        return sympy.getattr("Basic")?.getattr("__new__")?.call1(PyTuple::new(py,held)?);
    }
    let kwargs=PyDict::new(py);
    if matches!(m.kind(),MatrixKind::Add(_)|MatrixKind::Mul(..)|MatrixKind::Pow(..)){kwargs.set_item("evaluate",false)?;}
    if matches!(m.kind(),MatrixKind::Slice(..)){
        // SymPy's MatrixSlice.as_explicit infers shape from nested lists and
        // collapses empty axes. Empty explicit matrices preserve both axes.
        if let(Some(rows),Some(cols))=(m.rows().as_i64(),m.cols().as_i64()){
            if rows==0||cols==0{return sympy.getattr("ImmutableMatrix")?.call1((rows,cols,PyList::empty(py)));}
        }
        let mut raw=vec![sympy.getattr("MatrixSlice")?];raw.extend(args.iter());
        return sympy.getattr("Basic")?.getattr("__new__")?.call(PyTuple::new(py,raw)?,None);
    }
    sympy.getattr(m.kind_name())?.call(args,Some(&kwargs))
}
fn permutation_object<'py>(py:Python<'py>,images:&[usize],sympy:bool)->PyResult<Bound<'py,PyAny>>{
    py.import(if sympy{"sympy.combinatorics"}else{"symox.combinatorics"})?
        .getattr("Permutation")?.call1((images.to_vec(),))
}
fn permutation_images(value:&Bound<'_,PyAny>)->PyResult<Vec<usize>>{
    if let Ok(matrix)=value.downcast::<PyMatrixExpr>() {
        if let MatrixKind::Permutation(images)=matrix.get().m.kind(){return Ok(images.clone());}
    }
    value.getattr("array_form").and_then(|v|v.extract::<Vec<usize>>())
        .map_err(|_|PyValueError::new_err("Expected a finite Permutation or PermutationMatrix"))
}
fn matrix_args<'py>(py:Python<'py>,m:&MatrixExpr,sympy:bool)->PyResult<Bound<'py,PyTuple>>{
    let scalar=|e:&Expr|->PyResult<Bound<'py,PyAny>>{if sympy{crate::expr_to_sympy(py,e)}else{Ok(wrap(e.clone()).into_pyobject(py)?.into_any())}};
    let matrix=|e:&MatrixExpr|if sympy{matrixexpr_to_sympy(py,e)}else{matrix_object(py,e)};
    let mut args=Vec::new();
    match m.kind(){
        MatrixKind::Permutation(p)=>args.push(permutation_object(py,p,sympy)?),
        MatrixKind::Permute(parent,p,axis)=>{args.push(matrix(parent)?);args.push(permutation_object(py,p,sympy)?);args.push(scalar(&Expr::int(*axis as i64))?);},
        MatrixKind::Symbol(name)=>{args.push(name.into_pyobject(py)?.into_any());args.push(scalar(m.rows())?);args.push(scalar(m.cols())?);},
        MatrixKind::ApplyFunction(function,parent)=>{args.push(scalar(function)?);args.push(matrix(parent)?);},
        MatrixKind::Function(lambda)=>{args.push(scalar(m.rows())?);args.push(scalar(m.cols())?);args.push(scalar(lambda)?);},
        MatrixKind::Identity=>args.push(scalar(m.rows())?),MatrixKind::Zero|MatrixKind::One=>{args.push(scalar(m.rows())?);args.push(scalar(m.cols())?);},
        MatrixKind::Add(a)|MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)|MatrixKind::BlockDiag(a)=>for e in a{args.push(matrix(e)?);},
        MatrixKind::Block(grid)=>{let value=PyBlockGrid{grid:grid.clone()};args.push(if sympy{value.to_sympy(py)?}else{value.into_pyobject(py)?.into_any()});},
        MatrixKind::Slice(parent,rows,cols)=>{args.push(matrix(parent)?);for spec in [rows,cols]{args.push(if sympy{PyTuple::new(py,[scalar(&spec.start)?,scalar(&spec.stop)?,scalar(&spec.step)?])?.into_any()}else{slice_object(py,spec)?});}},
        MatrixKind::HadamardPower(a,b)=>for operand in [a,b]{args.push(match operand{MatrixOperand::Scalar(e)=>scalar(e)?,MatrixOperand::Matrix(m)=>matrix(m)?});},
        MatrixKind::Mul(c,a)=>{if !c.is_one(){args.push(scalar(c)?);}for e in a{args.push(matrix(e)?);}},
        MatrixKind::Pow(b,e)=>{args.push(matrix(b)?);args.push(scalar(e)?);},
        MatrixKind::Inverse(b)=>{args.push(matrix(b)?);args.push(scalar(&Expr::minus_one())?);},
        MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Conjugate(b)|MatrixKind::Diag(b)|MatrixKind::DiagonalOf(b)|MatrixKind::Diagonal(b)=>args.push(matrix(b)?),
        MatrixKind::Explicit(value)=>{args.push(value.rows().into_pyobject(py)?.into_any());args.push(value.cols().into_pyobject(py)?.into_any());args.push(PyList::new(py,value.data().iter().map(|e|wrap(e.clone())))?.into_any());},
    }
    PyTuple::new(py,args)
}
fn construct(kind:&str,args:&Bound<'_,PyTuple>,evaluate:bool)->PyResult<PyObject>{
    let py=args.py();let fail=||PyTypeError::new_err(format!("Invalid {kind} arguments"));
    let matrix=|i|to_matrixexpr(&args.get_item(i)?);
    let scalar=|i|to_expr(&args.get_item(i)?);
    let operand=|i|->PyResult<MatrixOperand>{let value=args.get_item(i)?;if is_matrix(&value){Ok(MatrixOperand::Matrix(to_matrixexpr(&value)?))}else{Ok(MatrixOperand::Scalar(to_expr(&value)?))}};
    let result=match kind{
        "PermutationMatrix" if args.len()==1=>MatrixExpr::permutation(permutation_images(&args.get_item(0)?)?,false),
        "MatrixPermute" if args.len()==3=>{let axis=scalar(2)?.as_i64().filter(|i|*i==0||*i==1).ok_or_else(||PyValueError::new_err("Permutation axis must be zero or one"))?;matrix(0)?.permute(permutation_images(&args.get_item(1)?)?,axis as u8,evaluate)},
        "MatrixSymbol" if args.len()==3=>{let name:String=args.get_item(0)?.extract()?;MatrixExpr::symbol(&name,scalar(1)?,scalar(2)?)},
        "ElementwiseApplyFunction" if args.len()==2=>MatrixExpr::apply_function(scalar(0)?,&matrix(1)?,false),
        "FunctionMatrix" if args.len()==3=>MatrixExpr::function(scalar(0)?,scalar(1)?,scalar(2)?),
        "Identity" if args.len()==1=>MatrixExpr::identity(scalar(0)?),"ZeroMatrix" if args.len()==2=>MatrixExpr::zero(scalar(0)?,scalar(1)?),
        "OneMatrix" if args.len()==2=>MatrixExpr::one(scalar(0)?,scalar(1)?),
        "DiagMatrix" if args.len()==1=>MatrixExpr::diag(&matrix(0)?),
        "DiagonalOf" if args.len()==1=>MatrixExpr::diagonal_of(&matrix(0)?),
        "DiagonalMatrix" if args.len()==1=>MatrixExpr::diagonal(&matrix(0)?),
        "BlockMatrix" if args.len()==1=>MatrixExpr::block(to_block_grid(&args.get_item(0)?)?),
        "BlockDiagMatrix"=>MatrixExpr::block_diag(&args.iter().map(|v|to_matrixexpr(&v)).collect::<PyResult<Vec<_>>>()?),
        "MatrixSlice" if args.len()==3=>{let parent=matrix(0)?;parent.slice_normalized(slice_spec(&args.get_item(1)?,parent.rows())?,slice_spec(&args.get_item(2)?,parent.cols())?)},
        "HadamardProduct"=>MatrixExpr::hadamard_product(&args.iter().map(|v|to_matrixexpr(&v)).collect::<PyResult<Vec<_>>>()?,evaluate),
        "KroneckerProduct"=>MatrixExpr::kronecker(&args.iter().map(|v|to_matrixexpr(&v)).collect::<PyResult<Vec<_>>>()?,evaluate),
        "HadamardPower" if args.len()==2=>{let a=operand(0)?;let b=operand(1)?;if let (MatrixOperand::Scalar(a),MatrixOperand::Scalar(b))=(&a,&b){return Ok(wrap(symox::ops::pow(a,b)).into_pyobject(py)?.into_any().unbind());}MatrixExpr::hadamard_power(a,b,evaluate)},
        "MatAdd"=>{if args.is_empty(){return Ok(wrap(Expr::zero()).into_pyobject(py)?.into_any().unbind());}let terms=args.iter().map(|v|to_matrixexpr(&v)).collect::<PyResult<Vec<_>>>()?;MatrixExpr::add(&terms,evaluate)},
        "MatMul"=>{let mut coefficient=Expr::one();let mut factors=Vec::new();for arg in args.iter(){if is_matrix(&arg){factors.push(to_matrixexpr(&arg)?);}else{coefficient=symox::ops::mul2(&coefficient,&to_expr(&arg)?);}}if factors.is_empty(){return Ok(wrap(coefficient).into_pyobject(py)?.into_any().unbind());}MatrixExpr::mul(coefficient,&factors,evaluate)},
        "MatPow" if args.len()==2=>matrix(0)?.pow(scalar(1)?,evaluate),"Transpose" if args.len()==1=>matrix(0)?.transpose(evaluate),"Adjoint" if args.len()==1=>matrix(0)?.adjoint(evaluate),"Inverse" if args.len()==1 || (args.len()==2 && scalar(1)?==Expr::minus_one())=>matrix(0)?.inverse(evaluate),"MatrixConjugate" if args.len()==1=>matrix(0)?.conjugate(evaluate),
        _=>return Err(fail()),
    }.map_err(|e|error(py,e))?;
    Ok(matrix_object(py,&result)?.unbind())
}
#[pyfunction]
#[pyo3(signature=(kind,args,evaluate=false))]
fn _matrixexpr_node(kind:&str,args:&Bound<'_,PyTuple>,evaluate:bool)->PyResult<PyObject>{construct(kind,args,evaluate)}
#[pyfunction]
fn _matrix_element(parent:&Bound<'_,PyAny>,i:&Bound<'_,PyAny>,j:&Bound<'_,PyAny>)->PyResult<PyExpr>{
    let m=to_matrixexpr(parent)?;Ok(wrap(m.element(&to_expr(i)?,&to_expr(j)?).map_err(|e|error(parent.py(),e))?))
}
#[pyfunction]
fn _as_matrixexpr(value:&Bound<'_,PyAny>)->PyResult<PyMatrixExpr>{Ok(wrap_matrixexpr(to_matrixexpr(value)?))}
#[pyfunction]
#[pyo3(signature=(parent,evaluate=false))]
fn _matrix_trace(parent:&Bound<'_,PyAny>,evaluate:bool)->PyResult<PyExpr>{
    let matrix=to_matrixexpr(parent)?;
    Ok(wrap((if evaluate{matrix.trace()}else{matrix.trace_unevaluated()}).map_err(|e|error(parent.py(),e))?))
}
#[pyfunction]
#[pyo3(signature=(parent,evaluate=false))]
fn _matrix_determinant(parent:&Bound<'_,PyAny>,evaluate:bool)->PyResult<PyExpr>{
    let matrix=to_matrixexpr(parent)?;
    Ok(wrap((if evaluate{matrix.determinant()}else{matrix.determinant_unevaluated()}).map_err(|e|error(parent.py(),e))?))
}
#[pyfunction]
fn _matrix_dot_product(a:&Bound<'_,PyAny>,b:&Bound<'_,PyAny>)->PyResult<PyExpr>{
    Ok(wrap(symox::matrixcontraction::dot_product(&to_matrixexpr(a)?,&to_matrixexpr(b)?,false).map_err(|e|error(a.py(),e))?))
}
#[pymethods]
impl PyMatrixExpr{
    #[getter] fn _dag_id(&self)->usize{self.m.dag_id()}
    #[getter] fn _symox_matrixexpr(&self)->bool{true}
    #[getter] fn is_Matrix(&self)->bool{true}
    #[getter] fn is_MatrixExpr(&self)->bool{true}
    #[getter] fn is_commutative(&self)->bool{false}
    #[getter] fn is_Symbol(&self)->bool{false}
    #[getter] fn is_Atom(&self)->bool{false}
    #[getter] fn is_Number(&self)->bool{false}
    #[getter] fn is_number(&self)->bool{false}
    #[getter] fn is_MatrixSymbol(&self)->bool{matches!(self.m.kind(),MatrixKind::Symbol(_))}
    #[getter] fn function(&self)->PyResult<PyExpr>{if let MatrixKind::ApplyFunction(function,_)=self.m.kind(){Ok(wrap(function.clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only ElementwiseApplyFunction has a function"))}}
    #[getter] fn expr(&self)->PyResult<Self>{if let MatrixKind::ApplyFunction(_,parent)=self.m.kind(){Ok(wrap_matrixexpr(parent.clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only ElementwiseApplyFunction has an expr"))}}
    #[getter] fn is_ElementwiseApplyFunction(&self)->bool{matches!(self.m.kind(),MatrixKind::ApplyFunction(..))}
    #[getter] fn is_FunctionMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::Function(_))}
    #[getter] fn lamda(&self)->PyResult<PyExpr>{self.m.function_lambda().cloned().map(wrap).ok_or_else(||pyo3::exceptions::PyAttributeError::new_err("Only FunctionMatrix has a lamda"))}
    #[getter] fn is_ZeroMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::Zero)}
    #[getter] fn is_Identity(&self)->bool{matches!(self.m.kind(),MatrixKind::Identity)||matches!(self.m.kind(),MatrixKind::Permutation(p) if p.iter().enumerate().all(|(i,j)|i==*j))}
    #[getter] fn is_PermutationMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::Permutation(_))}
    #[getter] fn is_MatrixPermute(&self)->bool{matches!(self.m.kind(),MatrixKind::Permute(..))}
    #[getter] fn is_OneMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::One)}
    #[getter] fn is_HadamardProduct(&self)->bool{matches!(self.m.kind(),MatrixKind::HadamardProduct(_))}
    #[getter] fn is_HadamardPower(&self)->bool{matches!(self.m.kind(),MatrixKind::HadamardPower(..))}
    #[getter] fn is_KroneckerProduct(&self)->bool{matches!(self.m.kind(),MatrixKind::Kronecker(_))}
    #[getter] fn is_BlockDiagMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::BlockDiag(_))}
    #[getter] fn is_MatrixSlice(&self)->bool{matches!(self.m.kind(),MatrixKind::Slice(..))}
    #[getter] fn is_MatAdd(&self)->bool{matches!(self.m.kind(),MatrixKind::Add(_))}
    #[getter] fn is_MatMul(&self)->bool{matches!(self.m.kind(),MatrixKind::Mul(..))}
    #[getter] fn is_MatPow(&self)->bool{matches!(self.m.kind(),MatrixKind::Pow(..)|MatrixKind::Inverse(_))}
    #[getter] fn kind(&self)->&'static str{self.m.kind_name()}
    #[getter] fn rows(&self)->PyExpr{wrap(self.m.rows().clone())}
    #[getter] fn cols(&self)->PyExpr{wrap(self.m.cols().clone())}
    #[getter] fn shape(&self)->(PyExpr,PyExpr){(self.rows(),self.cols())}
    #[getter] fn is_square(&self)->Option<bool>{if self.m.rows()==self.m.cols(){Some(true)}else{symox::matrix::ratfunc::is_zero_assume(&symox::ops::sub(self.m.rows(),self.m.cols()))}}
    #[getter] fn name(&self)->PyResult<String>{if let MatrixKind::Symbol(name)=self.m.kind(){Ok(name.clone())}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only MatrixSymbol has a name"))}}
    #[getter] fn base(&self,py:Python<'_>)->PyResult<PyObject>{match self.m.kind(){MatrixKind::Pow(b,_)|MatrixKind::Inverse(b)|MatrixKind::HadamardPower(MatrixOperand::Matrix(b),_)=>Ok(matrix_object(py,b)?.unbind()),MatrixKind::HadamardPower(MatrixOperand::Scalar(e),_)=>Ok(wrap(e.clone()).into_pyobject(py)?.into_any().unbind()),_=>Err(pyo3::exceptions::PyAttributeError::new_err("Only matrix powers have a base"))}}
    #[getter] fn exp(&self,py:Python<'_>)->PyResult<PyObject>{match self.m.kind(){MatrixKind::Pow(_,e)|MatrixKind::HadamardPower(_,MatrixOperand::Scalar(e))=>Ok(wrap(e.clone()).into_pyobject(py)?.into_any().unbind()),MatrixKind::Inverse(_)=>Ok(wrap(Expr::minus_one()).into_pyobject(py)?.into_any().unbind()),MatrixKind::HadamardPower(_,MatrixOperand::Matrix(m))=>Ok(matrix_object(py,m)?.unbind()),_=>Err(pyo3::exceptions::PyAttributeError::new_err("Only matrix powers have an exponent"))}}
    #[getter] fn arg(&self)->PyResult<PyMatrixExpr>{let c=self.m.children();if c.len()==1{Ok(wrap_matrixexpr(c[0].clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("Matrix expression does not have one matrix argument"))}}
    #[getter] fn parent(&self)->PyResult<PyMatrixExpr>{if let MatrixKind::Slice(parent,..)=self.m.kind(){Ok(wrap_matrixexpr(parent.clone()))}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only a matrix slice has a parent"))}}
    #[getter] fn rowslice<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{if let MatrixKind::Slice(_,row,_)=self.m.kind(){slice_object(py,row)}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only a matrix slice has rowslice"))}}
    #[getter] fn colslice<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{if let MatrixKind::Slice(_,_,col)=self.m.kind(){slice_object(py,col)}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only a matrix slice has colslice"))}}
    #[getter] fn on_diag(&self)->PyResult<bool>{if let MatrixKind::Slice(_,row,col)=self.m.kind(){Ok(row==col)}else{Err(pyo3::exceptions::PyAttributeError::new_err("Only a matrix slice has on_diag"))}}
    #[getter] fn blocks(&self,py:Python<'_>)->PyResult<PyBlockGrid>{Ok(PyBlockGrid{grid:self.m.block_grid().map_err(|e|error(py,e))?})}
    #[getter] fn blockshape(&self,py:Python<'_>)->PyResult<(usize,usize)>{Ok(self.m.block_grid().map_err(|e|error(py,e))?.shape())}
    #[getter] fn diag<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{if let MatrixKind::BlockDiag(blocks)=self.m.kind(){PyTuple::new(py,blocks.iter().map(|m|matrix_object(py,m)).collect::<PyResult<Vec<_>>>()?)}else{Err(pyo3::exceptions::PyAttributeError::new_err("This matrix has no diagonal blocks"))}}
    fn get_diag_blocks<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{self.diag(py)}
    #[getter] fn rowblocksizes(&self,py:Python<'_>)->PyResult<Vec<PyExpr>>{let grid=self.m.block_grid().map_err(|e|error(py,e))?;Ok((0..grid.shape().0).map(|i|wrap(grid.get(i,0).rows().clone())).collect())}
    #[getter] fn colblocksizes(&self,py:Python<'_>)->PyResult<Vec<PyExpr>>{let grid=self.m.block_grid().map_err(|e|error(py,e))?;Ok((0..grid.shape().1).map(|i|wrap(grid.get(0,i).cols().clone())).collect())}
    #[getter] fn is_BlockMatrix(&self)->bool{matches!(self.m.kind(),MatrixKind::Block(_)|MatrixKind::BlockDiag(_))}
    #[getter] fn is_structurally_symmetric(&self,py:Python<'_>)->PyResult<bool>{Ok(self.m.block_grid().map_err(|e|error(py,e))?.is_structurally_symmetric())}
    fn structurally_equal(&self,other:&Bound<'_,PyAny>)->PyResult<bool>{let first=self.m.block_grid().map_err(|e|error(other.py(),e))?;let second=to_matrixexpr(other)?.block_grid().map_err(|e|error(other.py(),e))?;Ok(first.same_partition(&second))}
    #[pyo3(signature=(mat="A",generalized=false))]
    fn schur(&self,py:Python<'_>,mat:&str,generalized:bool)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.schur(mat,generalized).map_err(|e|error(py,e))?))}
    fn LDUdecomposition(&self,py:Python<'_>)->PyResult<(Self,Self,Self)>{let(a,b,c)=self.m.ldu_decomposition().map_err(|e|error(py,e))?;Ok((wrap_matrixexpr(a),wrap_matrixexpr(b),wrap_matrixexpr(c)))}
    fn UDLdecomposition(&self,py:Python<'_>)->PyResult<(Self,Self,Self)>{let(a,b,c)=self.m.udl_decomposition().map_err(|e|error(py,e))?;Ok((wrap_matrixexpr(a),wrap_matrixexpr(b),wrap_matrixexpr(c)))}
    #[getter] fn args<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{matrix_args(py,&self.m,false)}
    fn class_key(&self,py:Python<'_>)->PyResult<PyObject>{
        crate::presentation_key_object(py,&symox::sortkey::matrix_class_key(&self.m))
    }
    #[pyo3(signature=(order=None))]
    fn sort_key(&self,py:Python<'_>,order:Option<&str>)->PyResult<PyObject>{
        crate::validate_sort_order(order)?;
        crate::presentation_key_object(py,&symox::sortkey::matrix_sort_key(&self.m))
    }
    #[getter] fn func<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{py.import("symox")?.getattr(self.m.kind_name())}
    #[getter] fn free_symbols<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PySet>>{
        let result=PySet::empty(py)?;for m in self.m.matrix_symbols(){result.add(wrap_matrixexpr(m))?;}
        for e in self.m.scalar_free_symbols(){result.add(crate::expr_object(py,e)?)?;}Ok(result)
    }
    fn __str__(&self)->String{self.m.to_str()}
    fn __repr__(&self)->String{self.m.to_str()}
    fn srepr(&self)->String{self.m.to_repr()}
    fn __hash__(&self)->isize{
        if let MatrixKind::Explicit(m)=self.m.kind(){
            let entries=m.data().iter().enumerate().filter(|(_,e)|!e.is_zero()).map(|(i,e)|((i/m.cols(),i%m.cols()),e.clone())).collect();
            return crate::matrices::matrix_value_hash(m.rows(),m.cols(),entries);
        }
        use std::hash::{Hash,Hasher};let mut h=std::collections::hash_map::DefaultHasher::new();self.m.hash(&mut h);let v=h.finish() as isize;if v== -1{-2}else{v}
    }
    fn __richcmp__(&self,other:&Bound<'_,PyAny>,op:pyo3::basic::CompareOp)->PyObject{
        use pyo3::basic::CompareOp::*;let py=other.py();let same=to_matrixexpr(other).is_ok_and(|m|m==self.m);
        match op{Eq=>same.into_pyobject(py).unwrap().to_owned().into_any().unbind(),Ne=>(!same).into_pyobject(py).unwrap().to_owned().into_any().unbind(),_=>py.NotImplemented()}
    }
    fn __add__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,false,false)}
    fn __radd__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,false,true)}
    fn __sub__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,true,false)}
    fn __rsub__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_add(other,true,true)}
    fn __mul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_mul(other,false)}
    fn __rmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{self.binary_mul(other,true)}
    fn __matmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if !is_matrix(other){return Err(PyTypeError::new_err("Matrix @ requires matrix operands"));}self.binary_mul(other,false)}
    fn __rmatmul__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if !is_matrix(other){return Err(PyTypeError::new_err("Matrix @ requires matrix operands"));}self.binary_mul(other,true)}
    fn __truediv__(&self,other:&Bound<'_,PyAny>)->PyResult<Self>{if is_matrix(other){return Err(PyTypeError::new_err("Matrix division requires a scalar divisor"));}let c=symox::ops::pow(&to_expr(other)?,&Expr::minus_one());Ok(wrap_matrixexpr(MatrixExpr::mul(c,&[self.m.clone()],true).map_err(|e|error(other.py(),e))?))}
    fn __pow__(&self,exponent:&Bound<'_,PyAny>,modulo:Option<&Bound<'_,PyAny>>)->PyResult<Self>{if modulo.is_some_and(|m|!m.is_none()){return Err(PyTypeError::new_err("Modular matrix powers are unsupported"));}if is_matrix(exponent){return Err(PyTypeError::new_err("Matrix exponent must be scalar"));}Ok(wrap_matrixexpr(self.m.pow(to_expr(exponent)?,true).map_err(|e|error(exponent.py(),e))?))}
    fn __neg__(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(MatrixExpr::mul(Expr::minus_one(),&[self.m.clone()],true).map_err(|e|error(py,e))?))}
    fn __pos__(&self)->Self{self.clone()}
    fn applyfunc(&self,function:&Bound<'_,PyAny>)->PyResult<PyObject>{
        let py=function.py();
        Ok(py.import("symox._matrixexpr")?.getattr("ElementwiseApplyFunction")?.call1((function,matrix_object(py,&self.m)?))?.unbind())
    }
    fn transpose(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.transpose(true).map_err(|e|error(py,e))?))}
    fn adjoint(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.adjoint(true).map_err(|e|error(py,e))?))}
    fn conjugate(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.conjugate(true).map_err(|e|error(py,e))?))}
    fn inverse(&self,py:Python<'_>)->PyResult<Self>{Ok(wrap_matrixexpr(self.m.inverse(true).map_err(|e|error(py,e))?))}
    fn trace(&self,py:Python<'_>)->PyResult<PyExpr>{Ok(wrap(self.m.trace().map_err(|e|error(py,e))?))}
    fn det(&self,py:Python<'_>)->PyResult<PyExpr>{Ok(wrap(self.m.determinant().map_err(|e|error(py,e))?))}
    fn trace_explicit(&self,py:Python<'_>)->PyResult<PyExpr>{Ok(wrap(self.m.trace_explicit().map_err(|e|error(py,e))?))}
    fn determinant_explicit(&self,py:Python<'_>)->PyResult<PyExpr>{Ok(wrap(self.m.determinant_explicit().map_err(|e|error(py,e))?))}
    fn inv(&self,py:Python<'_>)->PyResult<Self>{self.inverse(py)}
    #[getter] fn I(&self,py:Python<'_>)->PyResult<Self>{self.inverse(py)}
    // These are public Expr views, not the internal power/coefficient views
    // used to normalize matrix multiplication. In particular MatPow inherits
    // Expr.as_base_exp, and MatMul does not extract a scalar in as_coeff_Mul.
    fn as_base_exp<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        PyTuple::new(py,[matrix_object(py,&self.m)?,wrap(Expr::one()).into_pyobject(py)?.into_any()])
    }
    #[pyo3(signature=(rational=false))]
    fn as_coeff_Mul<'py>(&self,py:Python<'py>,rational:bool)->PyResult<Bound<'py,PyTuple>>{
        let _=rational;
        PyTuple::new(py,[wrap(Expr::one()).into_pyobject(py)?.into_any(),matrix_object(py,&self.m)?])
    }
    #[pyo3(signature=(rational=false))]
    fn as_coeff_Add<'py>(&self,py:Python<'py>,rational:bool)->PyResult<Bound<'py,PyTuple>>{
        let _=rational;
        PyTuple::new(py,[wrap(Expr::zero()).into_pyobject(py)?.into_any(),matrix_object(py,&self.m)?])
    }
    fn as_coeff_matrices<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        let MatrixKind::Mul(coefficient,factors)=self.m.kind() else{
            return Err(pyo3::exceptions::PyAttributeError::new_err("Only MatMul has as_coeff_matrices"));
        };
        if symox::decompose::commutativity(coefficient)==Some(false){
            return Err(PyNotImplementedError::new_err("Noncommutative scalars in MatMul are unsupported"));
        }
        PyTuple::new(py,[crate::expr_object(py,coefficient.clone())?.into_bound(py),PyList::new(py,factors.iter().map(|m|matrix_object(py,m)).collect::<PyResult<Vec<_>>>()?)?.into_any()])
    }
    fn as_coeff_mmul<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        let (coefficient,factors)=match self.m.kind(){
            MatrixKind::Mul(c,factors)=>(c.clone(),factors.clone()),
            _=>(Expr::one(),vec![self.m.clone()]),
        };
        if symox::decompose::commutativity(&coefficient)==Some(false){
            return Err(PyNotImplementedError::new_err("Noncommutative scalars in MatMul are unsupported"));
        }
        let product=MatrixExpr::mul(Expr::one(),&factors,false).map_err(|e|error(py,e))?;
        PyTuple::new(py,[crate::expr_object(py,coefficient)?.into_bound(py),matrix_object(py,&product)?])
    }
    #[pyo3(signature=(*deps,rational=true,**kwargs))]
    fn as_coeff_mul<'py>(&self,py:Python<'py>,deps:&Bound<'_,PyTuple>,rational:bool,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Bound<'py,PyTuple>>{
        if kwargs.is_some_and(|d|!d.is_empty()){
            return Err(PyNotImplementedError::new_err("Unsupported MatrixExpr.as_coeff_mul options"));
        }
        if !deps.is_empty(){
            let dependencies=deps.iter().map(|d|if is_matrix(&d){Ok(MatrixOperand::Matrix(to_matrixexpr(&d)?))}else{Ok(MatrixOperand::Scalar(to_expr(&d)?))}).collect::<PyResult<Vec<_>>>()?;
            let dependent=|value:&MatrixOperand|dependencies.iter().any(|target|match (value,target){
                (MatrixOperand::Matrix(value),MatrixOperand::Matrix(target))=>value.has_matrix(target),
                (MatrixOperand::Matrix(value),MatrixOperand::Scalar(target))=>value.has_scalar(target),
                (MatrixOperand::Scalar(value),MatrixOperand::Matrix(target))=>symox::matrixexpr::expression_has_matrix(value,target),
                (MatrixOperand::Scalar(value),MatrixOperand::Scalar(target))=>symox::matrixexpr::expression_has_scalar(value,target),
            });
            let args=if let MatrixKind::Mul(coefficient,factors)=self.m.kind(){
                let mut args=Vec::new();if !coefficient.is_one(){args.push(MatrixOperand::Scalar(coefficient.clone()));}
                args.extend(factors.iter().cloned().map(MatrixOperand::Matrix));args
            }else{vec![MatrixOperand::Matrix(self.m.clone())]};
            let mut independent_matrices=Vec::new();let mut independent_scalars=Vec::new();let mut tail=Vec::new();
            for argument in args{
                if dependent(&argument){tail.push(match argument{MatrixOperand::Scalar(e)=>crate::expr_object(py,e)?.into_bound(py),MatrixOperand::Matrix(m)=>matrix_object(py,&m)?});}
                else{match argument{MatrixOperand::Scalar(e)=>independent_scalars.push(e),MatrixOperand::Matrix(m)=>independent_matrices.push(m)}}
            }
            let scalar=symox::ops::mul(&independent_scalars);
            let head=if independent_matrices.is_empty(){
                if independent_scalars.is_empty()&&matches!(self.m.kind(),MatrixKind::Mul(..)){
                    return Err(PyNotImplementedError::new_err("An empty independent MatMul requires a dimensionless GenericIdentity"));
                }
                crate::expr_object(py,scalar)?.into_bound(py)
            }else{
                // Match _new_rawargs: singleton tails collapse, but multiple
                // factors keep their original matrix order and evaluation.
                if independent_matrices.len()==1&&scalar.is_one(){matrix_object(py,&independent_matrices[0])?}
                else{matrix_object(py,&MatrixExpr::mul(scalar,&independent_matrices,false).map_err(|e|error(py,e))?)?}
            };
            return PyTuple::new(py,[head,PyTuple::new(py,tail)?.into_any()]);
        }
        let mut tail=Vec::new();let mut head=Expr::one();
        if let MatrixKind::Mul(coefficient,factors)=self.m.kind(){
            let number=coefficient.is_number()||matches!(coefficient.kind(),symox::expr::Kind::Const(symox::expr::Constant::Infinity|symox::expr::Constant::NegInfinity|symox::expr::Constant::NaN));
            if number&&(!rational||coefficient.is_rational_number()){
                head=coefficient.clone();
            }else if number&&(coefficient.as_number().is_some_and(symox::number::Number::is_negative)||coefficient.is_constant(symox::expr::Constant::NegInfinity)){
                head=Expr::minus_one();tail.push(crate::expr_object(py,symox::ops::neg(coefficient))?.into_bound(py));
            }else if !coefficient.is_one(){tail.push(crate::expr_object(py,coefficient.clone())?.into_bound(py));}
            for matrix in factors{tail.push(matrix_object(py,matrix)?);}
        }else{tail.push(matrix_object(py,&self.m)?);}
        PyTuple::new(py,[crate::expr_object(py,head)?.into_bound(py),PyTuple::new(py,tail)?.into_any()])
    }
    #[pyo3(signature=(*patterns))]
    fn has(&self,patterns:&Bound<'_,PyTuple>)->PyResult<bool>{
        for pattern in patterns.iter(){
            if is_matrix(&pattern){if self.m.has_matrix(&to_matrixexpr(&pattern)?){return Ok(true);}}
            else if self.m.has_scalar(&to_expr(&pattern)?){return Ok(true);}
        }
        Ok(false)
    }
    fn as_two_terms<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyTuple>>{
        let args=matrix_args(py,&self.m,false)?;
        match self.m.kind(){
            MatrixKind::Add(terms)=>{
                if terms.len()==1{
                    return Err(PyNotImplementedError::new_err("A singleton MatAdd tail requires a dimensionless GenericZeroMatrix"));
                }
                if terms.len()==2{return Ok(args);}
                let tail=MatrixExpr::add(&terms[1..],false).map_err(|e|error(py,e))?;
                PyTuple::new(py,[matrix_object(py,&terms[0])?,matrix_object(py,&tail)?])
            },
            MatrixKind::Mul(coefficient,factors)=>{
                if args.len()==1{return PyTuple::new(py,[crate::expr_object(py,Expr::one())?.into_bound(py),matrix_object(py,&self.m)?]);}
                if args.len()==2{return Ok(args);}
                let remaining=if coefficient.is_one(){&factors[1..]}else{&factors[..]};
                let tail=MatrixExpr::mul(Expr::one(),remaining,false).map_err(|e|error(py,e))?;
                PyTuple::new(py,[args.get_item(0)?,matrix_object(py,&tail)?])
            },
            _=>Err(pyo3::exceptions::PyAttributeError::new_err("Only MatAdd and MatMul have as_two_terms")),
        }
    }
    #[getter] fn T(&self,py:Python<'_>)->PyResult<Self>{self.transpose(py)}
    #[getter] fn H(&self,py:Python<'_>)->PyResult<Self>{self.adjoint(py)}
    #[getter] fn C(&self,py:Python<'_>)->PyResult<Self>{self.conjugate(py)}
    fn as_explicit(&self,py:Python<'_>)->PyResult<crate::matrices::PyMatrix>{Ok(crate::matrices::wrapm(self.m.as_explicit().map_err(|e|error(py,e))?).into_pyobject(py)?.call_method0("as_immutable")?.extract()?)}
    fn as_mutable(&self,py:Python<'_>)->PyResult<crate::matrices::PyMatrix>{Ok(crate::matrices::wrapm(self.m.as_explicit().map_err(|e|error(py,e))?))}
    #[pyo3(signature=(deep=true,**kwargs))]
    fn doit(&self,py:Python<'_>,deep:bool,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<PyObject>{if kwargs.is_some_and(|d|!d.is_empty()){return Err(PyNotImplementedError::new_err("MatrixExpr doit options are unsupported"));}Ok(matrix_object(py,&self.m.doit(deep).map_err(|e|error(py,e))?)?.unbind())}
    #[pyo3(signature=(deep=true,**hints))]
    fn expand(&self,py:Python<'_>,deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<Self>{
        let (options,mut budget)=crate::parse_expand_options(deep,hints)?;
        Ok(wrap_matrixexpr(self.m.expand_with_options(&options,&mut budget).map_err(|e|error(py,e))?))
    }
    #[pyo3(signature=(*args,**kwargs))]
    fn diff(&self,args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Self>{
        if let Some(options)=kwargs {for (key,value) in options.iter(){
            if key.extract::<String>()? != "evaluate" || !value.is_truthy()? {return Err(PyNotImplementedError::new_err("MatrixExpr.diff supports evaluated scalar derivatives only"));}
        }}
        self.differentiate(args)
    }
    fn __getitem__(&self,key:&Bound<'_,PyAny>)->PyResult<PyObject>{
        let py=key.py();
        if key.is_instance_of::<PySlice>(){return Ok(wrap_matrixexpr(self.m.slice_normalized(slice_spec(key,self.m.rows())?,SliceSpec::from_raw([None,None,None],self.m.cols()).map_err(|e|error(py,e))?).map_err(|e|error(py,e))?).into_pyobject(py)?.into_any().unbind());}
        if let Ok(pair)=key.downcast::<PyTuple>(){
            if pair.len()!=2{return Err(PyTypeError::new_err("MatrixExpr indexing requires two indices"));}
            let row=pair.get_item(0)?;let col=pair.get_item(1)?;
            if row.is_instance_of::<PySlice>()||col.is_instance_of::<PySlice>(){return Ok(wrap_matrixexpr(self.m.slice_normalized(slice_spec(&row,self.m.rows())?,slice_spec(&col,self.m.cols())?).map_err(|e|error(py,e))?).into_pyobject(py)?.into_any().unbind());}
            return Ok(wrap(self.m.element(&to_expr(&row)?,&to_expr(&col)?).map_err(|e|error(py,e))?).into_pyobject(py)?.into_any().unbind());
        }
        let mut index:i64=key.extract().map_err(|_|PyIndexError::new_err("Single matrix indexing requires a concrete integer"))?;
        let rows=self.m.rows().as_i64();let cols=self.m.cols().as_i64().filter(|n|*n>0).ok_or_else(||PyNotImplementedError::new_err("Single indexing requires a concrete nonempty column dimension"))?;
        if index<0{index=rows.and_then(|r|r.checked_mul(cols)).and_then(|size|size.checked_add(index)).ok_or_else(||PyIndexError::new_err("Negative single indexing requires concrete finite dimensions"))?;}
        Ok(wrap(self.m.element(&Expr::int(index.div_euclid(cols)),&Expr::int(index.rem_euclid(cols))).map_err(|e|error(py,e))?).into_pyobject(py)?.into_any().unbind())
    }
    fn xreplace(&self,rule:&Bound<'_,PyDict>)->PyResult<Self>{let mut matrices=Vec::new();let mut scalars=FxHashMap::default();for (a,b) in rule.iter(){if is_matrix(&a){matrices.push((to_matrixexpr(&a)?,to_matrixexpr(&b)?));}else{scalars.insert(to_expr(&a)?,to_expr(&b)?);}}Ok(wrap_matrixexpr(self.m.xreplace(&matrices,&scalars).map_err(|e|error(rule.py(),e))?))}
    #[pyo3(signature=(*args,**kwargs))]
    fn subs(&self,args:&Bound<'_,PyTuple>,kwargs:Option<&Bound<'_,PyDict>>)->PyResult<Self>{
        if kwargs.is_some_and(|d|!d.is_empty()){return Err(PyNotImplementedError::new_err("MatrixExpr subs options are unsupported"));}
        let mut pairs=Vec::new();if args.len()==2{pairs.push((args.get_item(0)?,args.get_item(1)?));}else if args.len()==1{let rule=args.get_item(0)?;if let Ok(d)=rule.downcast::<PyDict>(){pairs.extend(d.iter());}else{for pair in rule.try_iter()?{let pair=pair?;let a=pair.get_item(0)?;let b=pair.get_item(1)?;pairs.push((a,b));}}}else{return Err(PyTypeError::new_err("subs expects a mapping, pairs or (old,new)"));}
        let mut result=self.m.clone();for (a,b) in pairs{result=if is_matrix(&a){result.xreplace(&[(to_matrixexpr(&a)?,to_matrixexpr(&b)?)],&FxHashMap::default())}else{let old=to_expr(&a)?;let new=to_expr(&b)?;result.scalar_map(&mut |e|symox::subs::subs(e,&[(old.clone(),new.clone())]))}.map_err(|e|error(args.py(),e))?;}Ok(wrap_matrixexpr(result))
    }
    fn to_sympy<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{matrixexpr_to_sympy(py,&self.m)}
    fn _sympy_<'py>(&self,py:Python<'py>)->PyResult<Bound<'py,PyAny>>{self.to_sympy(py)}
    fn __reduce__<'py>(&self,py:Python<'py>)->PyResult<(Bound<'py,PyAny>,Bound<'py,PyTuple>)>{crate::matrix_dag::reduce_matrix(py,&self.m)}
}
impl PyMatrixExpr{
    pub(crate) fn differentiate(&self,args:&Bound<'_,PyTuple>)->PyResult<Self>{
        let mut variables=crate::diff_vars(args)?;
        if variables.is_empty(){
            let scalar=self.m.scalar_free_symbols();let matrices=self.m.matrix_symbols();
            if scalar.len()==1&&matrices.is_empty(){variables.push((scalar[0].clone(),1));}
            else if scalar.is_empty()&&matrices.len()==1{return Err(PyNotImplementedError::new_err("Matrix-symbol derivatives require tensor expressions; supply a scalar differentiation variable"));}
            else{return Err(PyValueError::new_err("Supply the scalar variable(s) of differentiation"));}
        }
        Ok(wrap_matrixexpr(self.m.diff_many(&variables).map_err(|e|error(args.py(),e))?))
    }
    fn binary_add(&self,other:&Bound<'_,PyAny>,subtract:bool,reverse:bool)->PyResult<Self>{let py=other.py();let mut a=self.m.clone();let mut b=to_matrixexpr(other)?;if reverse{std::mem::swap(&mut a,&mut b);}if subtract{b=MatrixExpr::mul(Expr::minus_one(),&[b],true).map_err(|e|error(py,e))?;}Ok(wrap_matrixexpr(MatrixExpr::add(&[a,b],true).map_err(|e|error(py,e))?))}
    fn binary_mul(&self,other:&Bound<'_,PyAny>,reverse:bool)->PyResult<Self>{let py=other.py();let result=if is_matrix(other){let peer=to_matrixexpr(other)?;let terms=if reverse{vec![peer,self.m.clone()]}else{vec![self.m.clone(),peer]};MatrixExpr::mul(Expr::one(),&terms,true)}else{MatrixExpr::mul(to_expr(other)?,&[self.m.clone()],true)};Ok(wrap_matrixexpr(result.map_err(|e|error(py,e))?))}
}
#[pyfunction]
fn _block_collapse(value:&Bound<'_,PyAny>)->PyResult<PyObject>{
    let py=value.py();
    if let Ok(matrix)=value.downcast::<PyMatrixExpr>(){return Ok(wrap_matrixexpr(matrix.get().m.block_collapse().map_err(|e|error(py,e))?).into_pyobject(py)?.into_any().unbind());}
    if let Ok(expression)=value.downcast::<PyExpr>(){return crate::expr_object(py,symox::matrixexpr::block_collapse_scalar(&expression.get().e).map_err(|e|error(py,e))?);}
    if value.downcast::<crate::matrices::PyMatrix>().is_ok(){return Ok(value.clone().unbind());}
    if value.is_instance_of::<PyBlockGrid>(){return Err(PyTypeError::new_err("Block grids require a BlockMatrix wrapper before collapse"));}
    let expression=to_expr(value)?;crate::expr_object(py,symox::matrixexpr::block_collapse_scalar(&expression).map_err(|e|error(py,e))?)
}
pub(crate) fn register(m:&Bound<'_,PyModule>)->PyResult<()>{m.add_class::<PyMatrixExpr>()?;m.add_class::<PyBlockGrid>()?;m.add_class::<PySliceSpec>()?;m.add_function(wrap_pyfunction!(_normalized_slice,m)?)?;m.add_function(wrap_pyfunction!(_matrixexpr_node,m)?)?;m.add_function(wrap_pyfunction!(_block_collapse,m)?)?;m.add_function(wrap_pyfunction!(_matrix_element,m)?)?;m.add_function(wrap_pyfunction!(_as_matrixexpr,m)?)?;m.add_function(wrap_pyfunction!(_matrix_trace,m)?)?;m.add_function(wrap_pyfunction!(_matrix_determinant,m)?)?;m.add_function(wrap_pyfunction!(_matrix_dot_product,m)?)?;Ok(())}
