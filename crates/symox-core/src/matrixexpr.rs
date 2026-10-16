//! Typed, noncommutative symbolic matrix expressions over scalar Expr dimensions.
use std::cmp::Ordering;
use std::sync::{Arc,OnceLock};
use std::hash::{Hash,Hasher};
use rustc_hash::{FxHashMap, FxHashSet};
use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::matrix::{Matrix, MatrixError, MResult};
use crate::matrix::ratfunc::is_zero_assume;
use crate::ops::{add, mul2, neg, pow, sub};

#[derive(Clone, Debug)]
pub struct MatrixExpr(Arc<MatrixNode>);
#[derive(Debug)]
struct MatrixNode { kind: MatrixKind, rows: Expr, cols: Expr, hash:OnceLock<u64>, scalar_symbols:OnceLock<Vec<Expr>> }
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MatrixKind {
    Symbol(String), Explicit(Matrix), Identity, Zero, One,
    Add(Vec<MatrixExpr>), Mul(Expr, Vec<MatrixExpr>), Pow(MatrixExpr, Expr),
    Transpose(MatrixExpr), Adjoint(MatrixExpr), Inverse(MatrixExpr), Conjugate(MatrixExpr),
    HadamardProduct(Vec<MatrixExpr>), HadamardPower(MatrixOperand,MatrixOperand),
    Kronecker(Vec<MatrixExpr>),
    Diag(MatrixExpr), DiagonalOf(MatrixExpr), Diagonal(MatrixExpr),
    BlockDiag(Vec<MatrixExpr>), Block(BlockGrid), Slice(MatrixExpr,SliceSpec,SliceSpec),
}
/// A normalized affine range. `stop` is exclusive and may overshoot the
/// parent by less than a stride; length, not truncating division, is authoritative.
#[derive(Clone,Debug,PartialEq,Eq,Hash,PartialOrd,Ord)]
pub struct SliceSpec{pub start:Expr,pub stop:Expr,pub step:Expr,pub length:Expr,pub stop_is_omitted:bool}
impl SliceSpec{
    fn sign(step:&Expr)->MResult<bool>{
        if step.is_zero(){return Err(MatrixError::Value("Slice step cannot be zero".into()));}
        if crate::assumptions::is_integer(step)==Some(false){return Err(MatrixError::Type("Slice step must be an integer".into()));}
        if crate::assumptions::is_positive(step)==Some(true){Ok(true)}else if crate::assumptions::is_negative(step)==Some(true){Ok(false)}else{Err(MatrixError::NotImplemented("Symbolic slice step needs a known positive or negative sign".into()))}
    }
    fn parts(start:Expr,step:Expr,length:Expr)->MResult<Self>{
        let positive=Self::sign(&step)?;
        if crate::assumptions::is_integer(&start)==Some(false)||crate::assumptions::is_integer(&length)==Some(false)||crate::assumptions::is_negative(&length)==Some(true){return Err(MatrixError::Type("Normalized slice indices and length must be nonnegative-length integers".into()));}
        if length.is_zero(){return Ok(Self{start:Expr::zero(),stop:Expr::zero(),step:Expr::one(),length,stop_is_omitted:false});}
        let stop=crate::ops::add2(&start,&mul2(&step,&length));let stop_is_omitted=!positive&&crate::assumptions::is_negative(&stop)==Some(true);
        Ok(Self{start,stop,step,length,stop_is_omitted})
    }
    pub fn normalized(start:Expr,stop:Expr,step:Expr,length:Expr,omitted:bool)->MResult<Self>{
        let result=Self::parts(start.clone(),step.clone(),length.clone())?;
        if start!=result.start||stop!=result.stop||step!=result.step||length!=result.length||omitted!=result.stop_is_omitted{return Err(MatrixError::Value("Inconsistent normalized slice metadata".into()));}
        Ok(result)
    }
    pub fn args(&self)->[Option<Expr>;3]{[Some(self.start.clone()),if self.stop_is_omitted{None}else{Some(self.stop.clone())},Some(self.step.clone())]}
    pub fn fields(&self)->[Expr;4]{[self.start.clone(),self.stop.clone(),self.step.clone(),self.length.clone()]}
    pub fn index(index:Expr,size:&Expr)->MResult<Self>{
        MatrixExpr::validate_index(&index,size)?;
        let start=if crate::assumptions::is_negative(&index)==Some(true){crate::ops::add2(size,&index)}else if crate::assumptions::is_nonnegative(&index)==Some(true){index}else{crate::relational::piecewise(vec![(crate::ops::add2(size,&index),crate::relational::lt(&index,&Expr::zero())),(index,Expr::boolean(true))])};
        Self::parts(start,Expr::one(),Expr::one())
    }
    pub fn from_raw(raw:[Option<Expr>;3],size:&Expr)->MResult<Self>{
        let[start,stop,step]=raw;let step=step.unwrap_or_else(Expr::one);let positive=Self::sign(&step)?;
        let last=sub(size,&Expr::one());let(lo,hi)=if positive{(Expr::zero(),size.clone())}else{(Expr::minus_one(),last.clone())};
        let bound=|value:Expr|->MResult<Expr>{
            if crate::assumptions::is_integer(&value)==Some(false){return Err(MatrixError::Type("Slice bounds must be integers".into()));}
            let shifted=if crate::assumptions::is_negative(&value)==Some(true){crate::ops::add2(size,&value)}else if crate::assumptions::is_nonnegative(&value)==Some(true){value}else{crate::relational::piecewise(vec![(crate::ops::add2(size,&value),crate::relational::lt(&value,&Expr::zero())),(value,Expr::boolean(true))])};
            Ok(crate::functions::max(vec![lo.clone(),crate::functions::min(vec![shifted,hi.clone()])]))
        };
        let start=match start{Some(v)=>bound(v)?,None=>if positive{Expr::zero()}else{last}};
        let stop=match stop{Some(v)=>bound(v)?,None=>if positive{size.clone()}else{Expr::minus_one()}};
        let length=crate::functions::max(vec![Expr::zero(),func(Func::Ceiling,vec![crate::ops::div(&sub(&stop,&start),&step)])]);
        Self::parts(start,step,length)
    }
    fn same_presentation(&self,other:&Self)->bool{self.stop_is_omitted==other.stop_is_omitted&&self.fields().iter().zip(other.fields()).all(|(a,b)|a.same_presentation(&b))}
    fn mapped(&self,f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        let a=f(&self.start);let b=f(&self.stop);let h=f(&self.step);let length=f(&self.length);
        if a.same_presentation(&self.start)&&b.same_presentation(&self.stop)&&h.same_presentation(&self.step)&&length.same_presentation(&self.length){return Ok(self.clone());}
        Self::sign(&h)?;let expected=crate::ops::add2(&a,&mul2(&h,&length));
        let length=if expected==b{length}else{crate::functions::max(vec![Expr::zero(),func(Func::Ceiling,vec![crate::ops::div(&sub(&b,&a),&h)])])};
        Self::parts(a,h,length)
    }
    fn compose(&self,next:&Self)->MResult<Self>{Self::parts(crate::ops::add2(&self.start,&mul2(&self.step,&next.start)),mul2(&self.step,&next.step),next.length.clone())}
    fn validate_parent(&self,size:&Expr)->MResult<()>{
        if self.length.is_zero(){return Ok(());}
        if crate::assumptions::is_positive(&self.length)==Some(true){let last=crate::ops::add2(&self.start,&mul2(&self.step,&sub(&self.length,&Expr::one())));for index in [&self.start,&last]{if crate::assumptions::is_negative(index)==Some(true)||crate::assumptions::is_nonnegative(&sub(index,size))==Some(true){return Err(MatrixError::Index("Normalized slice selects an index outside its parent".into()));}}}
        Ok(())
    }
    fn repr(&self)->String{format!("NormalizedSlice({}, {}, {}, {}, {})",crate::printer::srepr(&self.start),crate::printer::srepr(&self.stop),crate::printer::srepr(&self.step),crate::printer::srepr(&self.length),if self.stop_is_omitted{"True"}else{"False"})}
}

/// Immutable, matrix-valued grid. Its cells use noncommutative matrix
/// arithmetic rather than the commutative scalar expression engine.
#[derive(Clone,Debug,PartialEq,Eq,Hash)]
pub struct BlockGrid{rows:usize,cols:usize,data:Vec<MatrixExpr>}
impl BlockGrid{
    pub fn new(rows:usize,cols:usize,data:Vec<MatrixExpr>)->MResult<Self>{
        if rows.checked_mul(cols)!=Some(data.len()){return Err(MatrixError::Shape("Block grid entries do not match its shape".into()));}
        Ok(Self{rows,cols,data})
    }
    pub fn shape(&self)->(usize,usize){(self.rows,self.cols)}
    pub fn data(&self)->&[MatrixExpr]{&self.data}
    pub fn get(&self,row:usize,col:usize)->&MatrixExpr{&self.data[row*self.cols+col]}
    pub fn same_partition(&self,other:&Self)->bool{
        self.shape()==other.shape()&&(self.rows==0||self.cols==0||(0..self.rows).all(|i|MatrixExpr::same_dimension(self.get(i,0).rows(),other.get(i,0).rows()))&&(0..self.cols).all(|j|MatrixExpr::same_dimension(self.get(0,j).cols(),other.get(0,j).cols())))
    }
    pub fn is_structurally_symmetric(&self)->bool{self.rows==self.cols&&(0..self.rows).all(|i|MatrixExpr::same_dimension(self.get(i,0).rows(),self.get(0,i).cols()))}
    fn can_multiply(&self,other:&Self)->bool{self.cols==other.rows&&(0..self.cols).all(|i|MatrixExpr::same_dimension(self.get(0,i).cols(),other.get(i,0).rows()))}
    fn diagonal(&self)->bool{self.is_structurally_symmetric()&&(0..self.rows).all(|i|(0..self.cols).all(|j|i==j||matches!(self.get(i,j).kind(),MatrixKind::Zero)))}
    pub fn transpose(&self)->Self{Self{rows:self.cols,cols:self.rows,data:(0..self.cols).flat_map(|i|(0..self.rows).map(move|j|self.get(j,i).clone())).collect()}}
    pub fn add(&self,other:&Self)->MResult<Self>{
        if self.shape()!=other.shape(){return Err(MatrixError::Shape("Block grid sizes do not match".into()));}
        Self::new(self.rows,self.cols,self.data.iter().zip(&other.data).map(|(a,b)|MatrixExpr::add(&[a.clone(),b.clone()],true)).collect::<MResult<Vec<_>>>()?)
    }
    pub fn scale(&self,scalar:&Expr)->MResult<Self>{Self::new(self.rows,self.cols,self.data.iter().map(|m|MatrixExpr::mul(scalar.clone(),&[m.clone()],true)).collect::<MResult<Vec<_>>>()?)}
    pub fn matmul(&self,other:&Self)->MResult<Self>{
        if self.cols!=other.rows{return Err(MatrixError::Shape("Block grid multiplication sizes do not match".into()));}
        let mut data=Vec::new();for i in 0..self.rows{for j in 0..other.cols{
            let terms=(0..self.cols).map(|k|MatrixExpr::mul(Expr::one(),&[self.get(i,k).clone(),other.get(k,j).clone()],true)).collect::<MResult<Vec<_>>>()?;
            if terms.is_empty(){return Err(MatrixError::NotImplemented("The shape of an empty block sum is unknown".into()));}
            data.push(MatrixExpr::add(&terms,true)?);
        }}Self::new(self.rows,other.cols,data)
    }
}
#[derive(Clone,Debug,PartialEq,Eq,Hash,PartialOrd,Ord)]
pub enum MatrixOperand {Scalar(Expr),Matrix(MatrixExpr)}
impl MatrixOperand {
    pub fn matrix(&self)->Option<&MatrixExpr>{if let Self::Matrix(m)=self{Some(m)}else{None}}
    fn same_presentation(&self,other:&Self)->bool{match(self,other){(Self::Scalar(a),Self::Scalar(b))=>a.same_presentation(b),(Self::Matrix(a),Self::Matrix(b))=>a.dag_id()==b.dag_id(),_=>false}}
    fn scalar(&self)->Option<&Expr>{if let Self::Scalar(e)=self{Some(e)}else{None}}
    fn entry(&self,i:&Expr,j:&Expr)->MResult<Expr>{match self{Self::Scalar(e)=>Ok(e.clone()),Self::Matrix(m)=>m.element(i,j)}}
    fn map(&self,matrices:&mut dyn FnMut(&MatrixExpr)->MatrixExpr,scalars:&mut dyn FnMut(&Expr)->Expr)->Self{match self{Self::Scalar(e)=>Self::Scalar(scalars(e)),Self::Matrix(m)=>Self::Matrix(matrices(m))}}
}
impl Hash for MatrixExpr {
    fn hash<H:Hasher>(&self,state:&mut H){
        let hash=self.0.hash.get_or_init(||{
            let mut hasher=std::collections::hash_map::DefaultHasher::new();
            self.kind().hash(&mut hasher);self.rows().hash(&mut hasher);self.cols().hash(&mut hasher);
            hasher.finish()
        });
        state.write_u64(*hash);
    }
}
impl PartialEq for MatrixExpr {
    fn eq(&self,other:&Self)->bool{
        if Arc::ptr_eq(&self.0,&other.0){return true;}
        if self.0.hash.get().zip(other.0.hash.get()).is_some_and(|(a,b)|a!=b){return false;}
        self.cmp(other)==Ordering::Equal
    }
}
impl Eq for MatrixExpr {}
impl Ord for MatrixExpr {
    fn cmp(&self,other:&Self)->Ordering{self.cmp_cached(other,&mut FxHashMap::default())}
}
impl MatrixExpr {
    fn cmp_cached(&self,other:&Self,memo:&mut FxHashMap<(usize,usize),Ordering>)->Ordering{
        if Arc::ptr_eq(&self.0,&other.0){return Ordering::Equal;}
        let pair=(self.dag_id(),other.dag_id());
        if let Some(result)=memo.get(&pair){return *result;}
        let kind=self.kind_name().cmp(other.kind_name());if kind!=Ordering::Equal{return kind;}
        let result=match (self.kind(),other.kind()){
            (MatrixKind::Symbol(a),MatrixKind::Symbol(b))=>a.cmp(b).then_with(||self.rows().cmp(other.rows())).then_with(||self.cols().cmp(other.cols())),
            (MatrixKind::Explicit(a),MatrixKind::Explicit(b))=>a.shape().cmp(&b.shape()).then_with(||a.data().cmp(b.data())),
            (MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One,MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One)=>self.rows().cmp(other.rows()).then_with(||self.cols().cmp(other.cols())),
            (MatrixKind::Add(a),MatrixKind::Add(b))|(MatrixKind::HadamardProduct(a),MatrixKind::HadamardProduct(b))|(MatrixKind::Kronecker(a),MatrixKind::Kronecker(b))=>Self::cmp_children(a,b,memo),
            (MatrixKind::BlockDiag(a),MatrixKind::BlockDiag(b))=>Self::cmp_children(a,b,memo),
            (MatrixKind::Block(a),MatrixKind::Block(b))=>a.shape().cmp(&b.shape()).then_with(||Self::cmp_children(a.data(),b.data(),memo)),
            (MatrixKind::Slice(a,ar,ac),MatrixKind::Slice(b,br,bc))=>a.cmp_cached(b,memo).then_with(||ar.cmp(br)).then_with(||ac.cmp(bc)),
            (MatrixKind::HadamardPower(a,aa),MatrixKind::HadamardPower(b,bb))=>Self::cmp_operand(a,b,memo).then_with(||Self::cmp_operand(aa,bb,memo)),
            (MatrixKind::Mul(a,aa),MatrixKind::Mul(b,bb))=>a.cmp(b).then_with(||Self::cmp_children(aa,bb,memo)),
            (MatrixKind::Pow(a,aa),MatrixKind::Pow(b,bb))=>a.cmp_cached(b,memo).then_with(||aa.cmp(bb)),
            (MatrixKind::Transpose(a),MatrixKind::Transpose(b))|(MatrixKind::Adjoint(a),MatrixKind::Adjoint(b))|(MatrixKind::Inverse(a),MatrixKind::Inverse(b))|(MatrixKind::Conjugate(a),MatrixKind::Conjugate(b))|(MatrixKind::Diag(a),MatrixKind::Diag(b))|(MatrixKind::DiagonalOf(a),MatrixKind::DiagonalOf(b))|(MatrixKind::Diagonal(a),MatrixKind::Diagonal(b))=>a.cmp_cached(b,memo),
            _=>unreachable!("distinct matrix kinds have distinct names"),
        };
        memo.insert(pair,result);result
    }
    fn cmp_children(a:&[Self],b:&[Self],memo:&mut FxHashMap<(usize,usize),Ordering>)->Ordering{
        for (a,b) in a.iter().zip(b){let order=a.cmp_cached(b,memo);if order!=Ordering::Equal{return order;}}
        a.len().cmp(&b.len())
    }
    fn cmp_operand(a:&MatrixOperand,b:&MatrixOperand,memo:&mut FxHashMap<(usize,usize),Ordering>)->Ordering{
        match(a,b){(MatrixOperand::Scalar(a),MatrixOperand::Scalar(b))=>a.cmp(b),(MatrixOperand::Matrix(a),MatrixOperand::Matrix(b))=>a.cmp_cached(b,memo),(MatrixOperand::Scalar(_),_)=>Ordering::Less,(_,MatrixOperand::Scalar(_))=>Ordering::Greater}
    }
}

impl PartialOrd for MatrixExpr { fn partial_cmp(&self, other:&Self)->Option<Ordering>{Some(self.cmp(other))} }

fn dimension(e:&Expr)->MResult<()> {
    if crate::assumptions::is_integer(e) == Some(false) || crate::assumptions::is_negative(e) == Some(true) {
        return Err(MatrixError::Value("Matrix dimensions must be nonnegative integers".into()));
    }
    Ok(())
}
fn compatible(a:&Expr,b:&Expr)->MResult<()> {
    if a != b && is_zero_assume(&sub(a,b)) == Some(false) {
        return Err(MatrixError::Shape(format!("Incompatible matrix dimensions: {a} and {b}")));
    }
    Ok(())
}
fn finite_dim(e:&Expr)->MResult<usize> {
    e.as_i64().and_then(|n| usize::try_from(n).ok()).ok_or_else(||MatrixError::NotImplemented("Explicit matrices require finite integer dimensions".into()))
}
impl MatrixExpr {
    fn new(kind:MatrixKind,rows:Expr,cols:Expr)->Self {Self(Arc::new(MatrixNode{kind,rows,cols,hash:OnceLock::new(),scalar_symbols:OnceLock::new()}))}
    /// Allocation identity for a scoped DAG memo. Keep a clone alive while
    /// using this identifier; it is not a persistent symbolic identifier.
    pub fn dag_id(&self)->usize{Arc::as_ptr(&self.0) as usize}
    fn visit_nodes(&self,f:&mut dyn FnMut(&Self)){
        let mut seen=FxHashSet::default();let mut pending=vec![self.clone()];
        while let Some(node)=pending.pop(){if !seen.insert(node.dag_id()){continue;}f(&node);pending.extend(node.children());}
    }
    pub fn kind(&self)->&MatrixKind {&self.0.kind}
    pub fn rows(&self)->&Expr {&self.0.rows}
    pub fn cols(&self)->&Expr {&self.0.cols}
    pub fn shape(&self)->(Expr,Expr) {(self.rows().clone(),self.cols().clone())}
    pub fn symbol(name:&str,rows:Expr,cols:Expr)->MResult<Self> {
        dimension(&rows)?;dimension(&cols)?;
        Ok(Self::new(MatrixKind::Symbol(name.to_string()),rows,cols))
    }
    pub fn identity(n:Expr)->MResult<Self> { dimension(&n)?;Ok(Self::new(MatrixKind::Identity,n.clone(),n)) }
    pub fn zero(rows:Expr,cols:Expr)->MResult<Self> {dimension(&rows)?;dimension(&cols)?;Ok(Self::new(MatrixKind::Zero,rows,cols))}
    pub fn one(rows:Expr,cols:Expr)->MResult<Self> {dimension(&rows)?;dimension(&cols)?;Ok(Self::new(MatrixKind::One,rows,cols))}
    /// Construct a diagonal matrix from a proven row/column vector (or 0x0).
    pub fn diag(vector:&Self)->MResult<Self>{
        let size=if vector.cols().is_one(){vector.rows().clone()}else if vector.rows().is_one(){vector.cols().clone()}else if vector.rows().is_zero()&&vector.cols().is_zero(){Expr::zero()}else{return Err(MatrixError::Shape("DiagMatrix requires a row or column vector".into()));};
        Ok(Self::new(MatrixKind::Diag(vector.clone()),size.clone(),size))
    }
    pub fn diagonal_of(parent:&Self)->MResult<Self>{
        let size=crate::functions::min(vec![parent.rows().clone(),parent.cols().clone()]);
        Ok(Self::new(MatrixKind::DiagonalOf(parent.clone()),size,Expr::one()))
    }
    pub fn diagonal(parent:&Self)->MResult<Self>{Ok(Self::new(MatrixKind::Diagonal(parent.clone()),parent.rows().clone(),parent.cols().clone()))}
    pub fn block_diag(blocks:&[Self])->MResult<Self>{
        let rows=add(&blocks.iter().map(|m|m.rows().clone()).collect::<Vec<_>>());
        let cols=add(&blocks.iter().map(|m|m.cols().clone()).collect::<Vec<_>>());
        Ok(Self::new(MatrixKind::BlockDiag(blocks.to_vec()),rows,cols))
    }
    fn same_dimension(a:&Expr,b:&Expr)->bool{a==b||is_zero_assume(&sub(a,b))==Some(true)}
    pub fn block(grid:BlockGrid)->MResult<Self>{
        let(rows,cols)=grid.shape();if rows==0||cols==0{return Err(MatrixError::Value("BlockMatrix requires one or more nonempty rows of matrices".into()));}
        for i in 0..rows{for j in 0..cols{
            if !Self::same_dimension(grid.get(i,j).rows(),grid.get(i,0).rows())||!Self::same_dimension(grid.get(i,j).cols(),grid.get(0,j).cols()){
                return Err(MatrixError::Shape("Block matrix cells must have matching row heights and column widths".into()));
            }
        }}
        let height=add(&(0..rows).map(|i|grid.get(i,0).rows().clone()).collect::<Vec<_>>());
        let width=add(&(0..cols).map(|j|grid.get(0,j).cols().clone()).collect::<Vec<_>>());
        Ok(Self::new(MatrixKind::Block(grid),height,width))
    }
    pub fn block_grid(&self)->MResult<BlockGrid>{
        if let MatrixKind::Block(grid)=self.kind(){return Ok(grid.clone());}
        let MatrixKind::BlockDiag(blocks)=self.kind() else{return Err(MatrixError::Type("This matrix has no block grid".into()));};
        let mut entries=Vec::new();for (i,a) in blocks.iter().enumerate(){for (j,b) in blocks.iter().enumerate(){entries.push(if i==j{a.clone()}else{Self::zero(a.rows().clone(),b.cols().clone())?});}}
        BlockGrid::new(blocks.len(),blocks.len(),entries)
    }
    pub fn slice(&self,rows:[Option<Expr>;3],cols:[Option<Expr>;3])->MResult<Self>{
        self.slice_normalized(SliceSpec::from_raw(rows,self.rows())?,SliceSpec::from_raw(cols,self.cols())?)
    }
    pub fn slice_normalized(&self,rows:SliceSpec,cols:SliceSpec)->MResult<Self>{
        rows.validate_parent(self.rows())?;cols.validate_parent(self.cols())?;
        if let MatrixKind::Slice(parent,old_rows,old_cols)=self.kind(){return parent.slice_normalized(old_rows.compose(&rows)?,old_cols.compose(&cols)?);}
        Ok(Self::new(MatrixKind::Slice(self.clone(),rows.clone(),cols.clone()),rows.length,cols.length))
    }
    pub fn explicit(m:Matrix)->Self {Self::new(MatrixKind::Explicit(m.clone()),Expr::int(m.rows() as i64),Expr::int(m.cols() as i64))}
    fn scalar_nonfinite(e:&Expr)->bool{
        match e.kind(){Kind::Const(crate::expr::Constant::Infinity|crate::expr::Constant::NegInfinity|crate::expr::Constant::ComplexInfinity|crate::expr::Constant::NaN)=>true,Kind::Add(..)|Kind::Mul(..)=>e.args().iter().any(Self::scalar_nonfinite),_=>false}
    }
    fn has_nonfinite(&self)->bool{
        let mut scalars=Vec::new();self.scalar_fields(&mut scalars);
        scalars.iter().any(Self::scalar_nonfinite)
    }
    fn entry_product(a:&Expr,b:&Expr)->Expr{
        if (a.is_zero()&&Self::scalar_nonfinite(b))||(b.is_zero()&&Self::scalar_nonfinite(a)){Expr::nan()}else{mul2(a,b)}
    }
    pub fn hadamard_product(factors:&[Self],evaluate:bool)->MResult<Self>{
        let first=factors.first().ok_or_else(||MatrixError::Value("HadamardProduct needs at least one argument".into()))?;
        let (rows,cols)=first.shape();
        for factor in factors{compatible(&rows,factor.rows())?;compatible(&cols,factor.cols())?;}
        if !evaluate{return Ok(Self::new(MatrixKind::HadamardProduct(factors.to_vec()),rows,cols));}
        let mut flat=Vec::new();
        for factor in factors{if let MatrixKind::HadamardProduct(args)=factor.kind(){flat.extend(args.iter().cloned());}else{flat.push(factor.clone());}}
        if flat.iter().any(|f|matches!(f.kind(),MatrixKind::Zero))&&!flat.iter().any(Self::has_nonfinite){return Self::zero(rows,cols);}
        let mut explicit:Option<Matrix>=None;let mut rest=Vec::new();
        for factor in flat{match factor.kind(){
            MatrixKind::One=>{},
            MatrixKind::Explicit(matrix)=>{explicit=Some(if let Some(previous)=explicit{Matrix::new(matrix.rows(),matrix.cols(),previous.data().iter().zip(matrix.data()).map(|(a,b)|Self::entry_product(a,b)).collect())?}else{matrix.clone()});},
            _=>rest.push(factor),
        }}
        if let Some(matrix)=explicit{rest.push(Self::explicit(matrix));}
        rest.sort();
        let mut grouped=Vec::new();let mut i=0;
        while i<rest.len(){let mut end=i+1;while end<rest.len()&&rest[end]==rest[i]{end+=1;}
            grouped.push(if end-i==1{rest[i].clone()}else{Self::hadamard_power(MatrixOperand::Matrix(rest[i].clone()),MatrixOperand::Scalar(Expr::int((end-i) as i64)),true)?});i=end;
        }
        grouped.sort();
        Ok(match grouped.len(){0=>Self::one(rows,cols)?,1=>grouped.remove(0),_=>Self::new(MatrixKind::HadamardProduct(grouped),rows,cols)})
    }
    pub fn hadamard_power(base:MatrixOperand,exponent:MatrixOperand,evaluate:bool)->MResult<Self>{
        let matrix=base.matrix().or_else(||exponent.matrix()).ok_or_else(||MatrixError::Type("HadamardPower requires a matrix operand".into()))?;
        let (rows,cols)=matrix.shape();
        if let (Some(a),Some(b))=(base.matrix(),exponent.matrix()){compatible(a.rows(),b.rows())?;compatible(a.cols(),b.cols())?;}
        if evaluate{
            if exponent.scalar().is_some_and(Expr::is_zero){return Self::one(rows,cols);}
            if exponent.scalar().is_some_and(Expr::is_one){return Ok(base.matrix().unwrap().clone());}
            let explicit=|operand:&MatrixOperand|match operand{MatrixOperand::Scalar(_)=>true,MatrixOperand::Matrix(m)=>matches!(m.kind(),MatrixKind::Explicit(_))};
            if explicit(&base)&&explicit(&exponent){
                let mut entries=Vec::new();for i in 0..finite_dim(&rows)?{for j in 0..finite_dim(&cols)?{let(i,j)=(Expr::int(i as i64),Expr::int(j as i64));entries.push(pow(&base.entry(&i,&j)?,&exponent.entry(&i,&j)?));}}
                return Ok(Self::explicit(Matrix::new(finite_dim(&rows)?,finite_dim(&cols)?,entries)?));
            }
        }
        Ok(Self::new(MatrixKind::HadamardPower(base,exponent),rows,cols))
    }
    pub fn kronecker(factors:&[Self],evaluate:bool)->MResult<Self>{
        if factors.is_empty(){return Ok(Self::explicit(Matrix::eye(1)));}
        let rows=crate::ops::mul(&factors.iter().map(|m|m.rows().clone()).collect::<Vec<_>>());
        let cols=crate::ops::mul(&factors.iter().map(|m|m.cols().clone()).collect::<Vec<_>>());
        if factors.iter().all(|m|matches!(m.kind(),MatrixKind::Identity)){return Self::identity(rows);}
        if !evaluate{return Ok(Self::new(MatrixKind::Kronecker(factors.to_vec()),rows,cols));}
        let mut flat=Vec::new();let mut coefficient=Expr::one();
        for factor in factors{match factor.kind(){
            MatrixKind::Kronecker(args)=>flat.extend(args.iter().cloned()),
            MatrixKind::Mul(c,args)=>{coefficient=Self::entry_product(&coefficient,c);flat.push(Self::mul(Expr::one(),args,true)?);},
            _=>flat.push(factor.clone()),
        }}
        if flat.iter().any(|m|matches!(m.kind(),MatrixKind::Zero))&&!flat.iter().any(Self::has_nonfinite){return Self::zero(rows,cols);}
        if flat.len()==1{return Self::mul(coefficient,&flat,true);}
        if flat.iter().all(|m|matches!(m.kind(),MatrixKind::Explicit(_))){
            let node=Self::new(MatrixKind::Kronecker(flat),rows,cols);
            return Ok(Self::explicit(node.as_explicit()?.scale(&coefficient)));
        }
        let node=Self::new(MatrixKind::Kronecker(flat),rows,cols);
        if coefficient.is_one(){Ok(node)}else{Self::mul(coefficient,&[node],true)}
    }
    pub fn add(terms:&[Self],evaluate:bool)->MResult<Self> {
        let first=terms.first().ok_or_else(||MatrixError::Type("Matrix addition requires a matrix".into()))?;
        let (rows,cols)=first.shape();
        for term in terms {compatible(&rows,term.rows())?;compatible(&cols,term.cols())?;}
        if !evaluate {return Ok(Self::new(MatrixKind::Add(terms.to_vec()),rows,cols));}
        let mut flat=Vec::new();
        for term in terms {if let MatrixKind::Add(args)=term.kind(){flat.extend(args.iter().cloned());}else{flat.push(term.clone());}}
        let mut groups:Vec<(Self,Expr)>=Vec::new();
        let mut explicit:Option<Matrix>=None;
        for term in flat {
            if matches!(term.kind(),MatrixKind::Zero){continue;}
            if let MatrixKind::Explicit(m)=term.kind(){explicit=Some(if let Some(prev)=explicit{prev.add(m)?}else{m.clone()});continue;}
            let (base,coefficient)=match term.kind(){
                MatrixKind::Mul(c,factors)=>(Self::mul(Expr::one(),factors,true)?,c.clone()),
                _=>(term.clone(),Expr::one()),
            };
            if let Some((_,c))=groups.iter_mut().find(|(b,_)|b==&base){*c=crate::ops::add2(c,&coefficient);}else{groups.push((base,coefficient));}
        }
        let mut result=Vec::new();
        if let Some(m)=explicit {if !m.is_zero_matrix(){result.push(Self::explicit(m));}}
        for (base,c) in groups {if !c.is_zero(){result.push(Self::mul(c,&[base],true)?);}}
        result.sort();
        Ok(match result.len(){0=>Self::zero(rows,cols)?,1=>result.remove(0),_=>Self::new(MatrixKind::Add(result),rows,cols)})
    }
    pub fn mul(coefficient:Expr,factors:&[Self],evaluate:bool)->MResult<Self> {
        let first=factors.first().ok_or_else(||MatrixError::Type("Matrix multiplication requires a matrix".into()))?;
        let rows=first.rows().clone(); let cols=factors.last().unwrap().cols().clone();
        for pair in factors.windows(2){compatible(pair[0].cols(),pair[1].rows())?;}
        if !evaluate{return Ok(Self::new(MatrixKind::Mul(coefficient,factors.to_vec()),rows,cols));}
        let mut scalar=coefficient;
        let mut flat=Vec::new();
        for factor in factors {if let MatrixKind::Mul(c,fs)=factor.kind(){scalar=mul2(&scalar,c);flat.extend(fs.iter().cloned());}else{flat.push(factor.clone());}}
        if scalar.is_zero() || flat.iter().any(|f|matches!(f.kind(),MatrixKind::Zero)) {return Self::zero(rows,cols);}
        let mut result:Vec<Self>=Vec::new();
        for factor in flat {
            if matches!(factor.kind(),MatrixKind::Identity){continue;}
            if let Some(prev)=result.last().cloned(){
                if let (MatrixKind::Explicit(a),MatrixKind::Explicit(b))=(prev.kind(),factor.kind()){
                    result.pop();result.push(Self::explicit(a.matmul(b)?));continue;
                }
                let (a,ae)=prev.as_base_exp();let(b,be)=factor.as_base_exp();
                if a==b {
                    result.pop();let merged=a.pow(crate::ops::add2(&ae,&be),true)?;
                    if !matches!(merged.kind(),MatrixKind::Identity){result.push(merged);}continue;
                }
            }
            result.push(factor);
        }
        if result.is_empty(){return Ok(if scalar.is_one(){Self::identity(rows)?}else{Self::new(MatrixKind::Mul(scalar,vec![Self::identity(rows.clone())?]),rows,cols)});}
        if result.len()==1 {
            if scalar.is_one(){return Ok(result.remove(0));}
            if let MatrixKind::Explicit(m)=result[0].kind(){return Ok(Self::explicit(m.scale(&scalar)));}
            if let MatrixKind::Add(args)=result[0].kind(){return Self::add(&args.iter().map(|a|Self::mul(scalar.clone(),&[a.clone()],true)).collect::<MResult<Vec<_>>>()?,true);}
        }
        Ok(Self::new(MatrixKind::Mul(scalar,result),rows,cols))
    }
    pub fn as_base_exp(&self)->(Self,Expr){match self.kind(){MatrixKind::Pow(b,e)=>(b.clone(),e.clone()),MatrixKind::Inverse(b)=>(b.clone(),Expr::minus_one()),_=>(self.clone(),Expr::one())}}
    pub fn pow(&self, exponent:Expr, evaluate:bool)->MResult<Self>{
        compatible(self.rows(),self.cols())?;
        if evaluate {
            if exponent.is_zero(){return Self::identity(self.rows().clone());}
            if exponent.is_one(){return Ok(self.clone());}
            if exponent==Expr::minus_one(){return self.inverse(true);}
            if matches!(self.kind(),MatrixKind::Identity){return Ok(self.clone());}
            if let Some(n)=exponent.as_i64(){
                if let MatrixKind::Explicit(m)=self.kind(){return Ok(Self::explicit(m.powi(n)?));}
                if matches!(self.kind(),MatrixKind::Zero){if n<0{return Err(MatrixError::NonInvertible("Zero matrix is not invertible".into()));}return Ok(self.clone());}
                if let MatrixKind::Pow(b,e)=self.kind(){if e.is_integer_number(){return b.pow(mul2(e,&exponent),true);}}
                if let MatrixKind::Inverse(b)=self.kind(){return b.pow(neg(&exponent),true);}
            }
        }
        Ok(Self::new(MatrixKind::Pow(self.clone(),exponent),self.rows().clone(),self.cols().clone()))
    }
    pub fn inverse(&self,evaluate:bool)->MResult<Self>{
        compatible(self.rows(),self.cols())?;
        if evaluate {match self.kind(){
            MatrixKind::Identity=>return Ok(self.clone()),MatrixKind::Zero=>return Err(MatrixError::NonInvertible("Zero matrix is not invertible".into())),
            MatrixKind::Inverse(b)=>return Ok(b.clone()),
            MatrixKind::Explicit(m)=>return Ok(Self::explicit(m.inv()?)),
            MatrixKind::One=>{if let Some(n)=self.rows().as_i64(){return if n==0{Self::identity(Expr::zero())}else if n==1{Ok(self.clone())}else{Err(MatrixError::NonInvertible("An all-ones matrix of size greater than one is singular".into()))};}},
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>return Self::kronecker(&factors.iter().map(|m|m.inverse(true)).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::BlockDiag(blocks) if blocks.iter().all(|m|m.rows()==m.cols())=>return Self::block_diag(&blocks.iter().map(|m|m.inverse(true)).collect::<MResult<Vec<_>>>()?),
            MatrixKind::BlockDiag(blocks) if blocks.iter().any(|m|is_zero_assume(&sub(m.rows(),m.cols()))==Some(false))=>return Err(MatrixError::NonInvertible("A square block diagonal matrix with rectangular diagonal blocks is singular".into())),
            MatrixKind::Mul(c,fs) if fs.iter().all(|f|f.rows()==f.cols())=>return Self::mul(pow(c,&Expr::minus_one()),&fs.iter().rev().map(|f|f.inverse(true)).collect::<MResult<Vec<_>>>()?,true),
            _=>{},
        }}
        Ok(Self::new(MatrixKind::Inverse(self.clone()),self.rows().clone(),self.cols().clone()))
    }
    pub fn transpose(&self,evaluate:bool)->MResult<Self>{self.unary(0,evaluate)}
    pub fn adjoint(&self,evaluate:bool)->MResult<Self>{self.unary(1,evaluate)}
    pub fn conjugate(&self,evaluate:bool)->MResult<Self>{self.unary(2,evaluate)}
    fn require_square(&self)->MResult<()>{
        compatible(self.rows(),self.cols()).map_err(|_|MatrixError::NonSquare("Trace and determinant require a square matrix".into()))
    }
    pub fn trace_unevaluated(&self)->MResult<Expr>{
        self.require_square()?;
        Ok(Expr::from_kind(Kind::Func(Func::MatrixTrace(Arc::new(self.clone())),vec![])))
    }
    pub fn determinant_unevaluated(&self)->MResult<Expr>{
        self.require_square()?;
        Ok(Expr::from_kind(Kind::Func(Func::MatrixDeterminant(Arc::new(self.clone())),vec![])))
    }
    /// Evaluate certified trace identities, retaining a typed scalar trace
    /// when the matrix is abstract. No implicit expansion into scalar entries.
    pub fn trace(&self)->MResult<Expr>{
        self.require_square()?;
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),MatrixKind::Identity|MatrixKind::One=>Ok(self.rows().clone()),
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>Ok(crate::ops::mul(&factors.iter().map(Self::trace).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Block(grid) if grid.is_structurally_symmetric()=>Ok(add(&(0..grid.shape().0).map(|i|grid.get(i,i).trace()).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::BlockDiag(blocks) if blocks.iter().all(|m|m.rows()==m.cols())=>Ok(add(&blocks.iter().map(Self::trace).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Explicit(matrix)=>matrix.trace(),
            MatrixKind::Add(terms)=>Ok(add(&terms.iter().map(Self::trace).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Mul(coefficient,factors)=>{
                if factors.len()==1{return Ok(mul2(coefficient,&factors[0].trace()?));}
                // The trace is invariant under cyclic rotations, including
                // rectangular compatible factors forming a closed chain.
                let transposed=factors.iter().rev().map(|factor|factor.transpose(true)).collect::<MResult<Vec<_>>>()?;
                let mut canonical=factors.clone();
                for orientation in [factors,&transposed]{
                    for i in 0..orientation.len(){
                        let rotation=orientation[i..].iter().chain(&orientation[..i]).cloned().collect::<Vec<_>>();
                        if rotation<canonical{canonical=rotation;}
                    }
                }
                let parent=Self::mul(Expr::one(),&canonical,false)?;
                Ok(mul2(coefficient,&parent.trace_unevaluated()?))
            },
            MatrixKind::Transpose(base)=>base.trace(),
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>Ok(func(Func::Conjugate,vec![base.trace()?])),
            _=>self.trace_unevaluated(),
        }
    }
    pub fn determinant(&self)->MResult<Expr>{
        self.require_square()?;
        if self.rows().is_zero(){return Ok(Expr::one());}
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),MatrixKind::Identity=>Ok(Expr::one()),
            MatrixKind::One=>match self.rows().as_i64(){Some(1)=>Ok(Expr::one()),Some(n) if n>1=>Ok(Expr::zero()),_=>self.determinant_unevaluated()},
            MatrixKind::BlockDiag(blocks) if blocks.iter().all(|m|m.rows()==m.cols())=>Ok(crate::ops::mul(&blocks.iter().map(Self::determinant).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::BlockDiag(blocks) if blocks.iter().any(|m|is_zero_assume(&sub(m.rows(),m.cols()))==Some(false))=>Ok(Expr::zero()),
            MatrixKind::Kronecker(factors) if factors.iter().all(|m|m.rows()==m.cols())=>{
                let mut terms=Vec::new();for (index,factor) in factors.iter().enumerate(){let exponent=crate::ops::mul(&factors.iter().enumerate().filter(|(i,_)|*i!=index).map(|(_,m)|m.rows().clone()).collect::<Vec<_>>());terms.push(pow(&factor.determinant()?,&exponent));}Ok(crate::ops::mul(&terms))
            },
            MatrixKind::Explicit(matrix)=>matrix.det(),
            MatrixKind::Mul(coefficient,factors) if factors.iter().all(|factor|factor.rows()==factor.cols())=>{
                let mut values=vec![pow(coefficient,self.rows())];
                values.extend(factors.iter().map(Self::determinant).collect::<MResult<Vec<_>>>()?);
                Ok(crate::ops::mul(&values))
            },
            MatrixKind::Transpose(base)=>base.determinant(),
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>Ok(func(Func::Conjugate,vec![base.determinant()?])),
            MatrixKind::Inverse(base)=>{
                let determinant=base.determinant()?;
                if determinant.is_zero(){return Err(MatrixError::NonInvertible("A singular matrix has no inverse".into()));}
                Ok(pow(&determinant,&Expr::minus_one()))
            },
            MatrixKind::Pow(base,exponent) if crate::assumptions::is_integer(exponent)==Some(true)=>{
                let determinant=base.determinant()?;
                if determinant.is_zero()&&crate::assumptions::is_negative(exponent)==Some(true){return Err(MatrixError::NonInvertible("A singular matrix has no negative power".into()));}
                Ok(pow(&determinant,exponent))
            },
            _=>self.determinant_unevaluated(),
        }
    }
    pub fn trace_explicit(&self)->MResult<Expr>{
        self.require_square()?;
        let n=finite_dim(self.rows())?;
        Ok(add(&(0..n).map(|i|self.element(&Expr::int(i as i64),&Expr::int(i as i64))).collect::<MResult<Vec<_>>>()?))
    }
    pub fn determinant_explicit(&self)->MResult<Expr>{self.require_square()?;self.as_explicit()?.det()}
    /// Jacobi's adjugate formula is valid even at singular matrices. In
    /// contrast, det(A)*trace(A^-1*A') loses those cases through division.
    pub fn determinant_diff(&self,variable:&Expr)->MResult<Expr>{
        self.require_square()?;
        let derivative=self.diff(variable)?;
        if matches!(derivative.kind(),MatrixKind::Zero){return Ok(Expr::zero());}
        self.as_explicit()?.adjugate()?.matmul(&derivative.as_explicit()?)?.trace()
    }
    /// Derivative with respect to a scalar variable. Matrix factors retain
    /// their order; dimensions describe the result and are not scalar entries.
    pub fn diff(&self,variable:&Expr)->MResult<Self>{
        if !matches!(variable.kind(),Kind::Sym(_)|Kind::Func(Func::Undefined(_),_)){
            return Err(MatrixError::NotImplemented("MatrixExpr differentiation requires a scalar symbol or scalar function; matrix-variable derivatives require tensor expressions".into()));
        }
        let mut cache=FxHashMap::default();
        self.diff_cached(variable,&mut cache)
    }
    pub fn diff_many(&self,variables:&[(Expr,u32)])->MResult<Self>{
        let mut result=self.clone();
        for (variable,order) in variables {for _ in 0..*order {result=result.diff(variable)?;if matches!(result.kind(),MatrixKind::Zero){break;}}}
        Ok(result)
    }
    fn scalar_diff(e:&Expr,variable:&Expr)->MResult<Expr>{
        if crate::assumptions::is_real(variable)!=Some(true){
            let mut nonholomorphic=false;
            crate::traverse::preorder(e,&mut |node|{
                if let Kind::Func(Func::Conjugate|Func::Re|Func::Im,args)=node.kind(){
                    if args.iter().any(|a|!crate::diff::diff(a,variable).is_zero()){nonholomorphic=true;}
                }
            });
            if nonholomorphic{return Err(MatrixError::NotImplemented("Differentiating conjugate, real or imaginary matrix coefficients requires a real differentiation variable".into()));}
        }
        Ok(crate::diff::diff(e,variable))
    }
    fn diff_cached(&self,variable:&Expr,cache:&mut FxHashMap<usize,(Self,Self)>)->MResult<Self>{
        if let Some((_,result))=cache.get(&self.dag_id()){return Ok(result.clone());}
        let zero=||Self::zero(self.rows().clone(),self.cols().clone());
        let result=match self.kind(){
            MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>zero()?,
            MatrixKind::Explicit(matrix)=>{
                let derived=matrix.try_applyfunc::<MatrixError>(|e|Self::scalar_diff(e,variable))?;
                if derived.is_zero_matrix(){zero()?}else{Self::explicit(derived)}
            },
            MatrixKind::Add(terms)=>Self::add(&terms.iter().map(|a|a.diff_cached(variable,cache)).collect::<MResult<Vec<_>>>()?,true)?,
            MatrixKind::Block(grid)=>{let cells=grid.data().iter().map(|m|m.diff_cached(variable,cache)).collect::<MResult<Vec<_>>>()?;if cells.iter().all(|m|matches!(m.kind(),MatrixKind::Zero)){zero()?}else{Self::block(BlockGrid::new(grid.shape().0,grid.shape().1,cells)?)?}},
            MatrixKind::BlockDiag(blocks)=>{let derivatives=blocks.iter().map(|m|m.diff_cached(variable,cache)).collect::<MResult<Vec<_>>>()?;if derivatives.iter().all(|m|matches!(m.kind(),MatrixKind::Zero)){zero()?}else{Self::block_diag(&derivatives)?}},
            MatrixKind::Diag(parent)|MatrixKind::DiagonalOf(parent)|MatrixKind::Diagonal(parent)=>{
                let derivative=parent.diff_cached(variable,cache)?;if matches!(derivative.kind(),MatrixKind::Zero){zero()?}else{match self.kind(){MatrixKind::Diag(_)=>Self::diag(&derivative)?,MatrixKind::DiagonalOf(_)=>Self::diagonal_of(&derivative)?,_=>Self::diagonal(&derivative)?}}
            },
            MatrixKind::Slice(parent,rows,cols)=>{
                if !Self::scalar_diff(&rows.start,variable)?.is_zero()||!Self::scalar_diff(&cols.start,variable)?.is_zero()||!Self::scalar_diff(&rows.step,variable)?.is_zero()||!Self::scalar_diff(&cols.step,variable)?.is_zero(){return Err(MatrixError::NotImplemented("Differentiation of a slice with a moving start index or stride is unsupported".into()));}
                let derivative=parent.diff_cached(variable,cache)?;if matches!(derivative.kind(),MatrixKind::Zero){zero()?}else{derivative.slice_normalized(rows.clone(),cols.clone())?}
            },
            MatrixKind::HadamardProduct(factors)|MatrixKind::Kronecker(factors)=>{
                let mut terms=Vec::new();
                for (i,factor) in factors.iter().enumerate(){let derivative=factor.diff_cached(variable,cache)?;if matches!(derivative.kind(),MatrixKind::Zero){continue;}let mut product=factors.clone();product[i]=derivative;
                    terms.push(if matches!(self.kind(),MatrixKind::Kronecker(_)){Self::kronecker(&product,true)?}else{Self::hadamard_product(&product,true)?});}
                if terms.is_empty(){zero()?}else{Self::add(&terms,true)?}
            },
            MatrixKind::HadamardPower(base,exponent)=>{
                if exponent.scalar().is_some_and(Expr::is_zero){return zero();}
                let mut derive=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>Self::scalar_diff(e,variable).map(MatrixOperand::Scalar),MatrixOperand::Matrix(m)=>m.diff_cached(variable,cache).map(MatrixOperand::Matrix)};
                let db=derive(base)?;let de=derive(exponent)?;
                let is_zero=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>e.is_zero(),MatrixOperand::Matrix(m)=>matches!(m.kind(),MatrixKind::Zero)};
                if is_zero(&db)&&is_zero(&de){zero()?}
                else if is_zero(&de){
                    let previous=match exponent{MatrixOperand::Scalar(e)=>MatrixOperand::Scalar(sub(e,&Expr::one())),MatrixOperand::Matrix(m)=>MatrixOperand::Matrix(Self::add(&[m.clone(),Self::mul(Expr::minus_one(),&[Self::one(m.rows().clone(),m.cols().clone())?],true)?],true)?)};
                    let mut factors=vec![Self::hadamard_power(base.clone(),previous,true)?];let mut coefficient=Expr::one();
                    for operand in [exponent,&db]{match operand{MatrixOperand::Scalar(e)=>coefficient=mul2(&coefficient,e),MatrixOperand::Matrix(m)=>factors.push(m.clone())}}
                    Self::mul(coefficient,&[Self::hadamard_product(&factors,true)?],true)?
                }else if let MatrixOperand::Scalar(a)=base{
                    let MatrixOperand::Matrix(exponent)=exponent else{unreachable!()};let MatrixOperand::Matrix(de)=de else{unreachable!()};let MatrixOperand::Scalar(da)=db else{unreachable!()};
                    let mut terms=vec![Self::mul(func(Func::Log,vec![a.clone()]),&[Self::hadamard_product(&[self.clone(),de],true)?],true)?];
                    if !da.is_zero(){terms.push(Self::mul(crate::ops::div(&da,a),&[Self::hadamard_product(&[self.clone(),exponent.clone()],true)?],true)?);}
                    Self::add(&terms,true)?
                }else{
                    let explicit=self.as_explicit().map_err(|_|MatrixError::NotImplemented("Variable Hadamard exponents with a matrix base require finite dimensions for exact entrywise differentiation".into()))?;
                    Self::explicit(explicit.try_applyfunc::<MatrixError>(|e|Self::scalar_diff(e,variable))?)
                }
            },
            MatrixKind::Mul(coefficient,factors)=>{
                if coefficient.is_zero(){return zero();}
                let mut terms=Vec::new();
                let dc=Self::scalar_diff(coefficient,variable)?;
                if !dc.is_zero(){terms.push(Self::mul(dc,factors,true)?);}
                for (index,factor) in factors.iter().enumerate(){
                    let derivative=factor.diff_cached(variable,cache)?;
                    if matches!(derivative.kind(),MatrixKind::Zero){continue;}
                    let mut product=factors.clone();product[index]=derivative;
                    terms.push(Self::mul(coefficient.clone(),&product,true)?);
                }
                if terms.is_empty(){zero()?}else{Self::add(&terms,true)?}
            },
            MatrixKind::Pow(base,exponent)=>{
                if exponent.is_zero(){return zero();}
                let derivative=base.diff_cached(variable,cache)?;
                let de=Self::scalar_diff(exponent,variable)?;
                if matches!(derivative.kind(),MatrixKind::Zero)&&de.is_zero(){zero()?}
                else{
                    let n=exponent.as_i64().ok_or_else(||MatrixError::NotImplemented("Derivatives of matrix powers require a constant integer exponent".into()))?;
                    if n<0{
                        let inverse=base.inverse(true)?;
                        let dinverse=Self::mul(Expr::minus_one(),&[inverse.clone(),derivative,inverse.clone()],true)?;
                        Self::diff_power(&inverse,&dinverse,n.unsigned_abs(),&mut FxHashMap::default())?
                    }else{Self::diff_power(base,&derivative,n as u64,&mut FxHashMap::default())?}
                }
            },
            MatrixKind::Inverse(base)=>{
                let derivative=base.diff_cached(variable,cache)?;
                if matches!(derivative.kind(),MatrixKind::Zero){zero()?}
                else{let inverse=base.inverse(true)?;Self::mul(Expr::minus_one(),&[inverse.clone(),derivative,inverse],true)?}
            },
            MatrixKind::Transpose(base)=>base.diff_cached(variable,cache)?.transpose(true)?,
            MatrixKind::Adjoint(base)|MatrixKind::Conjugate(base)=>{
                let derivative=base.diff_cached(variable,cache)?;
                if matches!(derivative.kind(),MatrixKind::Zero){zero()?}
                else{
                    if crate::assumptions::is_real(variable)!=Some(true){return Err(MatrixError::NotImplemented("Differentiating a dependent adjoint or conjugate requires a real scalar variable".into()));}
                    if matches!(self.kind(),MatrixKind::Adjoint(_)){derivative.adjoint(true)?}else{derivative.conjugate(true)?}
                }
            },
        };
        cache.insert(self.dag_id(),(self.clone(),result.clone()));Ok(result)
    }
    /// A balanced product rule keeps very large integer powers compact.
    fn diff_power(base:&Self,derivative:&Self,n:u64,cache:&mut FxHashMap<u64,Self>)->MResult<Self>{
        if n==0{return Self::zero(base.rows().clone(),base.cols().clone());}
        if n==1{return Ok(derivative.clone());}
        if let Some(result)=cache.get(&n){return Ok(result.clone());}
        let left=n/2;let right=n-left;
        let exponent=|n:u64|Expr::num(crate::number::Number::from_bigint(num_bigint::BigInt::from(n)));
        let dl=Self::diff_power(base,derivative,left,cache)?;
        let dr=Self::diff_power(base,derivative,right,cache)?;
        let result=Self::add(&[
            Self::mul(Expr::one(),&[dl,base.pow(exponent(right),true)?],true)?,
            Self::mul(Expr::one(),&[base.pow(exponent(left),true)?,dr],true)?,
        ],true)?;
        cache.insert(n,result.clone());Ok(result)
    }
    fn unary(&self,op:u8,evaluate:bool)->MResult<Self>{
        let apply=|m:&Self|m.unary(op,true);
        if evaluate {match self.kind(){
            MatrixKind::Identity=>return Ok(self.clone()),MatrixKind::One=>return if op==2{Ok(self.clone())}else{Self::one(self.cols().clone(),self.rows().clone())},MatrixKind::Zero=>return if op==2{Ok(self.clone())}else{Self::zero(self.cols().clone(),self.rows().clone())},
            MatrixKind::HadamardProduct(factors)=>return Self::hadamard_product(&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Kronecker(factors)=>return Self::kronecker(&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Diag(parent)=>return if op==0{Ok(self.clone())}else{Self::diag(&parent.conjugate(true)?)},
            MatrixKind::Diagonal(parent)=>return Self::diagonal(&apply(parent)?),
            MatrixKind::DiagonalOf(parent) if op==2=>return Self::diagonal_of(&parent.conjugate(true)?),
            MatrixKind::Block(grid)=>{let grid=if op==2{grid.clone()}else{grid.transpose()};return Self::block(BlockGrid::new(grid.shape().0,grid.shape().1,grid.data().iter().map(apply).collect::<MResult<Vec<_>>>()?)?);},
            MatrixKind::BlockDiag(blocks)=>return Self::block_diag(&blocks.iter().map(apply).collect::<MResult<Vec<_>>>()?),
            MatrixKind::Slice(parent,rows,cols)=>return if op==2{apply(parent)?.slice_normalized(rows.clone(),cols.clone())}else{apply(parent)?.slice_normalized(cols.clone(),rows.clone())},
            MatrixKind::HadamardPower(base,exponent) if op==0 || exponent.scalar().is_some_and(|e|crate::assumptions::is_integer(e)==Some(true))=>{
                let transform=|o:&MatrixOperand|->MResult<MatrixOperand>{Ok(match o{MatrixOperand::Matrix(m)=>MatrixOperand::Matrix(apply(m)?),MatrixOperand::Scalar(e)=>MatrixOperand::Scalar(if op==0{e.clone()}else{func(Func::Conjugate,vec![e.clone()])})})};
                return Self::hadamard_power(transform(base)?,transform(exponent)?,true);
            },
            MatrixKind::Explicit(m)=>return Ok(Self::explicit(match op{0=>m.transpose(),1=>m.adjoint(),_=>m.conjugate()})),
            MatrixKind::Transpose(b) if op==0=>return Ok(b.clone()),MatrixKind::Adjoint(b) if op==1=>return Ok(b.clone()),MatrixKind::Conjugate(b) if op==2=>return Ok(b.clone()),
            MatrixKind::Transpose(b) if op==1=>return b.conjugate(true),MatrixKind::Adjoint(b) if op==0=>return b.conjugate(true),
            MatrixKind::Conjugate(b) if op==0=>return b.adjoint(true),MatrixKind::Conjugate(b) if op==1=>return b.transpose(true),
            MatrixKind::Add(args)=>return Self::add(&args.iter().map(apply).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Mul(c,fs)=>{let mut factors=fs.clone();if op!=2{factors.reverse();}let c=if op==0{c.clone()}else{func(Func::Conjugate,vec![c.clone()])};return Self::mul(c,&factors.iter().map(apply).collect::<MResult<Vec<_>>>()?,true);},
            MatrixKind::Pow(b,e)=>return apply(b)?.pow(if op==1{func(Func::Conjugate,vec![e.clone()])}else{e.clone()},true),
            MatrixKind::Inverse(b)=>return apply(b)?.inverse(true),_=>{},
        }}
        let (r,c)=if op==2{self.shape()}else{(self.cols().clone(),self.rows().clone())};
        Ok(Self::new(match op{0=>MatrixKind::Transpose(self.clone()),1=>MatrixKind::Adjoint(self.clone()),_=>MatrixKind::Conjugate(self.clone())},r,c))
    }
    fn same_children(a:&[Self],b:&[Self])->bool{a.len()==b.len()&&a.iter().zip(b).all(|(a,b)|a.dag_id()==b.dag_id())}
    /// Collapse compatible block arithmetic without assuming any Schur pivot invertible.
    pub fn block_collapse(&self)->MResult<Self>{self.collapse_cached(&mut FxHashMap::default())}
    fn collapse_cached(&self,memo:&mut FxHashMap<usize,Self>)->MResult<Self>{
        if let Some(node)=memo.get(&self.dag_id()){return Ok(node.clone());}
        let old=self.children();let children=old.iter().map(|c|c.collapse_cached(memo)).collect::<MResult<Vec<_>>>()?;
        let node=if Self::same_children(&old,&children){self.clone()}else{let mut mapped=children.iter();self.map_children(&mut |_|mapped.next().unwrap().clone(),&mut |e|e.clone())?};
        let result=match node.kind(){
            MatrixKind::Block(grid) if grid.shape()==(1,1)=>grid.get(0,0).clone(),
            MatrixKind::Add(terms)=>{
                let grids=terms.iter().map(|m|m.block_grid()).collect::<MResult<Vec<_>>>();
                if let Ok(grids)=grids{if grids.iter().all(|g|g.same_partition(&grids[0])){let mut grid=grids[0].clone();for other in &grids[1..]{grid=grid.add(other)?;}Self::block(grid)?}else{node.clone()}}else{node.clone()}
            },
            MatrixKind::Mul(coefficient,factors)=>{
                let grids=factors.iter().map(|m|m.block_grid()).collect::<MResult<Vec<_>>>();
                if let Ok(grids)=grids{if !grids.is_empty()&&grids.windows(2).all(|p|p[0].can_multiply(&p[1])){let mut grid=grids[0].clone();for other in &grids[1..]{grid=grid.matmul(other)?;}grid=grid.scale(coefficient)?;if factors.iter().all(|m|matches!(m.kind(),MatrixKind::BlockDiag(_)))&&grid.diagonal(){Self::block_diag(&(0..grid.shape().0).map(|i|grid.get(i,i).clone()).collect::<Vec<_>>())?}else{Self::block(grid)?}}else{node.clone()}}else{node.clone()}
            },
            MatrixKind::Pow(base,exponent)=>{
                if let(Some(mut n),Ok(mut square))=(exponent.as_i64().filter(|n|*n>=0),base.block_grid()){
                    if n==0{Self::identity(base.rows().clone())?}else if n==1{base.clone()}else if square.is_structurally_symmetric(){let mut answer:Option<BlockGrid>=None;while n>0{if n&1==1{answer=Some(if let Some(a)=answer{a.matmul(&square)?}else{square.clone()});}n>>=1;if n>0{square=square.matmul(&square)?;}}let answer=answer.unwrap();if answer.shape()==(0,0){Self::block_diag(&[])?}else{Self::block(answer)?}}else{node.clone()}
                }else{node.clone()}
            },
            MatrixKind::Inverse(base)=>{
                if let Ok(grid)=base.block_grid(){if grid.diagonal(){Self::block_diag(&(0..grid.shape().0).map(|i|grid.get(i,i).inverse(true)).collect::<MResult<Vec<_>>>()?)?}else{node.clone()}}else{node.clone()}
            },
            MatrixKind::Transpose(base)=>base.transpose(true)?,MatrixKind::Adjoint(base)=>base.adjoint(true)?,MatrixKind::Conjugate(base)=>base.conjugate(true)?,
            _=>node.clone(),
        };
        let result=if let MatrixKind::Block(grid)=result.kind(){if grid.shape()==(1,1){grid.get(0,0).clone()}else{result}}else{result};
        memo.insert(self.dag_id(),result.clone());Ok(result)
    }
    pub fn children(&self)->Vec<Self>{match self.kind(){MatrixKind::Add(a)|MatrixKind::Mul(_,a)|MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)|MatrixKind::BlockDiag(a)=>a.clone(),MatrixKind::Block(grid)=>grid.data().to_vec(),MatrixKind::Slice(parent,..)=>vec![parent.clone()],MatrixKind::HadamardPower(a,b)=>[a,b].into_iter().filter_map(|o|o.matrix().cloned()).collect(),MatrixKind::Pow(b,_)|MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)|MatrixKind::Diag(b)|MatrixKind::DiagonalOf(b)|MatrixKind::Diagonal(b)=>vec![b.clone()],_=>Vec::new()}}
    pub fn map_children(&self, matrices:&mut dyn FnMut(&Self)->Self, scalars:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        match self.kind(){
            MatrixKind::Symbol(name)=>{let r=scalars(self.rows());let c=scalars(self.cols());if r.same_presentation(self.rows())&&c.same_presentation(self.cols()){Ok(self.clone())}else{Self::symbol(name,r,c)}},
            MatrixKind::Identity=>{let n=scalars(self.rows());if n.same_presentation(self.rows()){Ok(self.clone())}else{Self::identity(n)}},
            MatrixKind::Zero|MatrixKind::One=>{let r=scalars(self.rows());let c=scalars(self.cols());if r.same_presentation(self.rows())&&c.same_presentation(self.cols()){Ok(self.clone())}else if matches!(self.kind(),MatrixKind::One){Self::one(r,c)}else{Self::zero(r,c)}},
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)=>{let b:Vec<_>=a.iter().map(matrices).collect();if Self::same_children(&b,a){Ok(self.clone())}else if matches!(self.kind(),MatrixKind::Kronecker(_)){Self::kronecker(&b,true)}else{Self::hadamard_product(&b,true)}},
            MatrixKind::Block(grid)=>{let cells:Vec<_>=grid.data().iter().map(matrices).collect();if Self::same_children(&cells,grid.data()){Ok(self.clone())}else{Self::block(BlockGrid::new(grid.shape().0,grid.shape().1,cells)?)}},
            MatrixKind::BlockDiag(a)=>{let b:Vec<_>=a.iter().map(matrices).collect();if Self::same_children(&b,a){Ok(self.clone())}else{Self::block_diag(&b)}},
            MatrixKind::Slice(parent,rows,cols)=>{let mapped=matrices(parent);let r=rows.mapped(scalars)?;let c=cols.mapped(scalars)?;if mapped.dag_id()==parent.dag_id()&&r.same_presentation(rows)&&c.same_presentation(cols){Ok(self.clone())}else{mapped.slice_normalized(r,c)}},
            MatrixKind::HadamardPower(a,b)=>{let aa=a.map(matrices,scalars);let bb=b.map(matrices,scalars);if aa.same_presentation(a)&&bb.same_presentation(b){Ok(self.clone())}else{Self::hadamard_power(aa,bb,true)}},
            MatrixKind::Explicit(m)=>{let entries:Vec<_>=m.data().iter().map(scalars).collect();if entries.iter().zip(m.data()).all(|(a,b)|a.same_presentation(b)){Ok(self.clone())}else{Ok(Self::explicit(Matrix::new(m.rows(),m.cols(),entries)?))}},
            MatrixKind::Add(a)=>{let b:Vec<_>=a.iter().map(matrices).collect();if Self::same_children(&b,a){Ok(self.clone())}else{Self::add(&b,true)}},
            MatrixKind::Mul(c,a)=>{let coefficient=scalars(c);let b:Vec<_>=a.iter().map(matrices).collect();if Self::same_children(&b,a)&&coefficient.same_presentation(c){Ok(self.clone())}else{Self::mul(coefficient,&b,true)}},
            MatrixKind::Pow(b,e)=>{let base=matrices(b);let exponent=scalars(e);if base.dag_id()==b.dag_id()&&exponent.same_presentation(e){Ok(self.clone())}else{base.pow(exponent,true)}},
            MatrixKind::Diag(parent)|MatrixKind::DiagonalOf(parent)|MatrixKind::Diagonal(parent)=>{let mapped=matrices(parent);if mapped.dag_id()==parent.dag_id(){Ok(self.clone())}else{match self.kind(){MatrixKind::Diag(_)=>Self::diag(&mapped),MatrixKind::DiagonalOf(_)=>Self::diagonal_of(&mapped),_=>Self::diagonal(&mapped)}}},
            MatrixKind::Transpose(b)|MatrixKind::Adjoint(b)|MatrixKind::Inverse(b)|MatrixKind::Conjugate(b)=>{let base=matrices(b);if base.dag_id()==b.dag_id(){return Ok(self.clone());}match self.kind(){MatrixKind::Transpose(_)=>base.transpose(true),MatrixKind::Adjoint(_)=>base.adjoint(true),MatrixKind::Inverse(_)=>base.inverse(true),_=>base.conjugate(true)}},
        }
    }
    pub fn matrix_symbols(&self)->Vec<Self>{
        let mut out=FxHashSet::default();
        self.visit_nodes(&mut |node|{if matches!(node.kind(),MatrixKind::Symbol(_)){out.insert(node.clone());}});
        let mut out:Vec<_>=out.into_iter().collect();out.sort();out
    }
    pub fn scalar_free_symbols(&self)->Vec<Expr>{
        self.0.scalar_symbols.get_or_init(||{
            let mut out=FxHashSet::default();self.visit_nodes(&mut |node|node.collect_scalars(&mut out));
            let mut out:Vec<_>=out.into_iter().collect();out.sort();out
        }).clone()
    }
    fn collect_scalars(&self,out:&mut FxHashSet<Expr>){
        let mut push=|e:&Expr|out.extend(crate::traverse::free_symbols(e));
        match self.kind(){
            MatrixKind::Slice(_,rows,cols)=>{for e in rows.fields().iter().chain(cols.fields().iter()){push(e);}},
            MatrixKind::Symbol(_)=>{}, MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>{push(self.rows());push(self.cols());},
            MatrixKind::Explicit(m)=>{for e in m.data(){push(e);}},MatrixKind::Mul(c,_)=>push(c),MatrixKind::Pow(_,e)=>push(e),MatrixKind::HadamardPower(a,b)=>{for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{push(e);}}},_=>{},
        }
    }
    pub fn matrix_map(&self,f:&mut dyn FnMut(&Self)->Self)->MResult<Self>{self.matrix_map_cached(f,&mut FxHashMap::default())}
    fn matrix_map_cached(&self,f:&mut dyn FnMut(&Self)->Self,memo:&mut FxHashMap<usize,Self>)->MResult<Self>{
        if let Some(result)=memo.get(&self.dag_id()){return Ok(result.clone());}
        let replacement=f(self);if replacement!=*self{memo.insert(self.dag_id(),replacement.clone());return Ok(replacement);}
        let old=self.children();let children=old.iter().map(|c|c.matrix_map_cached(f,memo)).collect::<MResult<Vec<_>>>()?;
        let result=if Self::same_children(&children,&old){self.clone()}else{self.map_children(&mut |child|old.iter().position(|c|c==child).map(|i|children[i].clone()).unwrap_or_else(||child.clone()),&mut |e|e.clone())?};
        memo.insert(self.dag_id(),result.clone());Ok(result)
    }

    pub(crate) fn scalar_fields(&self,out:&mut Vec<Expr>){
        self.visit_nodes(&mut |node|match node.kind(){MatrixKind::Slice(_,rows,cols)=>out.extend(rows.fields().iter().chain(cols.fields().iter()).cloned()),MatrixKind::Symbol(_)|MatrixKind::Identity|MatrixKind::Zero|MatrixKind::One=>{out.push(node.rows().clone());out.push(node.cols().clone());},MatrixKind::Explicit(m)=>out.extend(m.data().iter().cloned()),MatrixKind::Mul(c,_)=>out.push(c.clone()),MatrixKind::Pow(_,e)=>out.push(e.clone()),MatrixKind::HadamardPower(a,b)=>{for operand in [a,b]{if let MatrixOperand::Scalar(e)=operand{out.push(e.clone());}}},_=>{}});
    }
    pub fn scalar_map(&self,f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        let mut fields=Vec::new();self.scalar_fields(&mut fields);
        let values:FxHashMap<_,_>=fields.into_iter().map(|e|{let mapped=f(&e);(e.dag_id(),(e,mapped))}).collect();
        if values.values().all(|(a,b)|a.same_presentation(b)){return Ok(self.clone());}
        self.remap_internal(&[],&mut |e|values.get(&e.dag_id()).map(|(_,mapped)|mapped.clone()).unwrap_or_else(||e.clone()))
    }
    pub fn map_all(&self,replacements:&[(Self,Self)],f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{
        self.remap_internal(replacements,f)
    }
    fn remap_internal(&self,replacements:&[(Self,Self)],f:&mut dyn FnMut(&Expr)->Expr)->MResult<Self>{self.remap_cached(replacements,f,&mut FxHashMap::default())}
    fn remap_cached(&self,replacements:&[(Self,Self)],f:&mut dyn FnMut(&Expr)->Expr,memo:&mut FxHashMap<usize,Self>)->MResult<Self>{
        if let Some(result)=memo.get(&self.dag_id()){return Ok(result.clone());}
        if let Some((_,value))=replacements.iter().find(|(key,_)|key==self){memo.insert(self.dag_id(),value.clone());return Ok(value.clone());}
        let children=self.children();
        let mapped=children.iter().map(|c|c.remap_cached(replacements,f,memo)).collect::<MResult<Vec<_>>>()?;
        let result=self.map_children(&mut |child|children.iter().position(|c|c==child).map(|i|mapped[i].clone()).unwrap_or_else(||child.clone()),f)?;
        memo.insert(self.dag_id(),result.clone());Ok(result)
    }

    pub fn xreplace(&self,matrices:&[(Self,Self)],scalars:&FxHashMap<Expr,Expr>)->MResult<Self>{
        let mut rules=scalars.clone();
        for (old,new) in matrices{
            let proxy=|m:&Self|Expr::from_kind(Kind::Func(Func::MatrixProxy(Arc::new(m.clone())),vec![]));
            rules.insert(proxy(old),proxy(new));
        }
        self.map_all(matrices,&mut |e|crate::subs::xreplace(e,&rules))
    }
    pub fn doit(&self,deep:bool)->MResult<Self>{self.doit_cached(deep,&mut FxHashMap::default())}
    fn doit_cached(&self,deep:bool,memo:&mut FxHashMap<usize,Self>)->MResult<Self>{
        if let Some(result)=memo.get(&self.dag_id()){return Ok(result.clone());}
        let children=if deep{self.children().iter().map(|c|c.doit_cached(true,memo)).collect::<MResult<Vec<_>>>()?}else{self.children()};
        let result=match self.kind(){MatrixKind::Diag(_)=>Self::diag(&children[0]),MatrixKind::DiagonalOf(_)=>Self::diagonal_of(&children[0]),MatrixKind::Diagonal(_)=>Self::diagonal(&children[0]),MatrixKind::Block(grid)=>Self::block(BlockGrid::new(grid.shape().0,grid.shape().1,children.clone())?),MatrixKind::BlockDiag(_)=>Self::block_diag(&children),MatrixKind::Slice(_,rows,cols)=>children[0].slice_normalized(rows.clone(),cols.clone()),MatrixKind::HadamardProduct(_)=>Self::hadamard_product(&children,true),MatrixKind::Kronecker(_)=>Self::kronecker(&children,true),MatrixKind::HadamardPower(a,b)=>{let mut index=0;let mut map=|_:&Self|{let child=children[index].clone();index+=1;child};Self::hadamard_power(a.map(&mut map,&mut |e|e.clone()),b.map(&mut map,&mut |e|e.clone()),true)},MatrixKind::Add(_)=>Self::add(&children,true),MatrixKind::Mul(c,_)=>Self::mul(c.clone(),&children,true),MatrixKind::Pow(_,e)=>children[0].pow(e.clone(),true),MatrixKind::Inverse(_)=>children[0].inverse(true),MatrixKind::Transpose(_)=>children[0].transpose(true),MatrixKind::Adjoint(_)=>children[0].adjoint(true),MatrixKind::Conjugate(_)=>children[0].conjugate(true),_=>Ok(self.clone())}?;
        let result=if result==*self{self.clone()}else{result};
        memo.insert(self.dag_id(),result.clone());Ok(result)
    }
    pub fn expand(&self)->MResult<Self>{
        let node=self.doit(true)?;
        match node.kind(){
            MatrixKind::Add(args)=>Self::add(&args.iter().map(|a|a.expand()).collect::<MResult<Vec<_>>>()?,true),
            MatrixKind::Mul(c,fs)=>{let mut terms=vec![Vec::new()];for factor in fs{let expanded=factor.expand()?;let choices=if let MatrixKind::Add(a)=expanded.kind(){a.clone()}else{vec![expanded]};let mut next=Vec::new();for prefix in terms{for term in &choices{let mut chain=prefix.clone();chain.push(term.clone());next.push(chain);}}terms=next;}Self::add(&terms.iter().map(|chain|Self::mul(c.clone(),chain,true)).collect::<MResult<Vec<_>>>()?,true)},_=>Ok(node)
        }
    }
    pub fn validate_index(index:&Expr,dim:&Expr)->MResult<()> {
        if crate::assumptions::is_integer(index)==Some(false){return Err(MatrixError::Index("Matrix indices must be integers".into()));}
        if crate::assumptions::is_nonnegative(&sub(index,dim))==Some(true)||crate::assumptions::is_negative(&crate::ops::add2(index,dim))==Some(true){return Err(MatrixError::Index("Matrix index out of range".into()));}
        Ok(())
    }
    fn kronecker_element(&self,factors:&[Self],i:&Expr,j:&Expr)->MResult<Expr>{
        if factors.iter().any(|m|crate::assumptions::is_positive(m.rows())!=Some(true)||crate::assumptions::is_positive(m.cols())!=Some(true)){return Ok(self.opaque_element(i,j));}
        let split=|index:&Expr,size:&Expr|{
            if let (Some(i),Some(n))=(index.as_i64(),size.as_i64()){(Expr::int(i.div_euclid(n)),Expr::int(i.rem_euclid(n)))}
            else{let quotient=func(Func::Floor,vec![crate::ops::div(index,size)]);let remainder=sub(index,&mul2(size,&quotient));(quotient,remainder)}
        };
        let mut i=i.clone();let mut j=j.clone();let mut value=Expr::one();
        for factor in factors.iter().rev(){let(quotient_i,row)=split(&i,factor.rows());let(quotient_j,col)=split(&j,factor.cols());value=Self::entry_product(&value,&factor.element(&row,&col)?);i=quotient_i;j=quotient_j;}
        Ok(value)
    }
    fn block_element(&self,grid:&BlockGrid,i:&Expr,j:&Expr)->MResult<Expr>{
        let locate=|index:&Expr,total:&Expr,sizes:Vec<&Expr>|{
            let mut local=if crate::assumptions::is_negative(index)==Some(true){crate::ops::add2(index,total)}else{index.clone()};
            for(position,size)in sizes.iter().enumerate(){let beyond=sub(&local,size);if crate::assumptions::is_nonnegative(&beyond)==Some(true){local=beyond;continue;}if crate::assumptions::is_negative(&beyond)==Some(true)||position+1==sizes.len(){return Some((position,local));}return None;}None
        };
        let row=locate(i,self.rows(),(0..grid.shape().0).map(|i|grid.get(i,0).rows()).collect());
        let col=locate(j,self.cols(),(0..grid.shape().1).map(|j|grid.get(0,j).cols()).collect());
        match(row,col){(Some((r,i)),Some((c,j)))=>grid.get(r,c).element(&i,&j),_=>Ok(self.opaque_element(i,j))}
    }
    fn block_diag_element(&self,blocks:&[Self],i:&Expr,j:&Expr)->MResult<Expr>{
        fn locate(blocks:&[MatrixExpr],index:&Expr,total:&Expr,rows:bool)->Option<(usize,Expr)>{
            let mut local=if crate::assumptions::is_negative(index)==Some(true){crate::ops::add2(index,total)}else{index.clone()};
            for (position,block) in blocks.iter().enumerate(){
                let size=if rows{block.rows()}else{block.cols()};let beyond=sub(&local,size);
                if crate::assumptions::is_nonnegative(&beyond)==Some(true){local=beyond;continue;}
                if crate::assumptions::is_negative(&beyond)==Some(true)||position+1==blocks.len(){return Some((position,local));}
                return None;
            }None
        }
        match(locate(blocks,i,self.rows(),true),locate(blocks,j,self.cols(),false)){
            (Some((row,i)),Some((col,j)))=>if row==col{blocks[row].element(&i,&j)}else{Ok(Expr::zero())},
            _=>Ok(self.opaque_element(i,j)),
        }
    }
    fn normalize_local_index(index:&Expr,size:&Expr)->Expr{if crate::assumptions::is_negative(index)==Some(true){crate::ops::add2(index,size)}else{index.clone()}}
    fn diagonal_element(&self,parent:&Self,i:&Expr,j:&Expr,vector:bool)->MResult<Expr>{
        let i=Self::normalize_local_index(i,self.rows());let j=Self::normalize_local_index(j,self.cols());
        let condition=crate::relational::eq(&i,&j);
        if condition==Expr::boolean(false){return Ok(Expr::zero());}
        let value=if vector{if parent.rows().is_one(){parent.element(&Expr::zero(),&i)?}else{parent.element(&i,&Expr::zero())?}}else{parent.element(&i,&j)?};
        Ok(crate::relational::piecewise(vec![(value,condition),(Expr::zero(),Expr::boolean(true))]))
    }
    pub fn element(&self,i:&Expr,j:&Expr)->MResult<Expr>{
        Self::validate_index(i,self.rows())?;Self::validate_index(j,self.cols())?;
        match self.kind(){
            MatrixKind::Zero=>Ok(Expr::zero()),
            MatrixKind::One=>Ok(Expr::one()),
            MatrixKind::Diag(parent)=>self.diagonal_element(parent,i,j,true),
            MatrixKind::Diagonal(parent)=>self.diagonal_element(parent,i,j,false),
            MatrixKind::DiagonalOf(parent)=>{let i=Self::normalize_local_index(i,self.rows());parent.element(&i,&i)},
            MatrixKind::Block(grid)=>self.block_element(grid,i,j),
            MatrixKind::BlockDiag(blocks)=>self.block_diag_element(blocks,i,j),
            MatrixKind::Slice(parent,rows,cols)=>{let normalize=|index:&Expr,size:&Expr|if crate::assumptions::is_negative(index)==Some(true){crate::ops::add2(index,size)}else{index.clone()};parent.element(&crate::ops::add2(&rows.start,&mul2(&rows.step,&normalize(i,self.rows()))),&crate::ops::add2(&cols.start,&mul2(&cols.step,&normalize(j,self.cols()))))},
            MatrixKind::HadamardProduct(factors)=>{let mut value=Expr::one();for factor in factors{value=Self::entry_product(&value,&factor.element(i,j)?);}Ok(value)},
            MatrixKind::HadamardPower(base,exponent)=>Ok(pow(&base.entry(i,j)?,&exponent.entry(i,j)?)),
            MatrixKind::Kronecker(factors)=>self.kronecker_element(factors,i,j),
            MatrixKind::Identity if i==j=>Ok(Expr::one()),
            MatrixKind::Identity if i.as_i64().is_some()&&j.as_i64().is_some()=>Ok(Expr::zero()),
            MatrixKind::Explicit(m)=>{if let(Some(i),Some(j))=(i.as_i64(),j.as_i64()){return Ok(m.get(i.rem_euclid(m.rows() as i64) as usize,j.rem_euclid(m.cols() as i64) as usize).clone());}Ok(self.opaque_element(i,j))},
            MatrixKind::Add(a)=>Ok(add(&a.iter().map(|e|e.element(i,j)).collect::<MResult<Vec<_>>>()?)),
            MatrixKind::Mul(c,fs)=>{if fs.len()==1{return Ok(mul2(c,&fs[0].element(i,j)?));}let first=&fs[0];if let Ok(n)=finite_dim(first.cols()){let rest=Self::mul(Expr::one(),&fs[1..],true)?;let mut terms=Vec::with_capacity(n);for k in 0..n{let k=Expr::int(k as i64);terms.push(mul2(&first.element(i,&k)?,&rest.element(&k,j)?));}Ok(mul2(c,&add(&terms)))}else{Ok(self.opaque_element(i,j))}},
            MatrixKind::Transpose(b)=>b.element(j,i),MatrixKind::Adjoint(b)=>Ok(func(Func::Conjugate,vec![b.element(j,i)?])),MatrixKind::Conjugate(b)=>Ok(func(Func::Conjugate,vec![b.element(i,j)?])),
            MatrixKind::Pow(_,e) if e.as_i64().is_none()=>Ok(self.opaque_element(i,j)),
            MatrixKind::Pow(..)|MatrixKind::Inverse(_)=>{if finite_dim(self.rows()).is_ok(){let m=self.as_explicit()?;if let(Some(i),Some(j))=(i.as_i64(),j.as_i64()){return Ok(m.get(i.rem_euclid(m.rows() as i64) as usize,j.rem_euclid(m.cols() as i64) as usize).clone());}}Ok(self.opaque_element(i,j))},
            _=>Ok(self.opaque_element(i,j)),
        }
    }
    /// Validate and construct a held scalar entry without materializing its parent.
    pub fn element_unevaluated(&self,i:&Expr,j:&Expr)->MResult<Expr>{
        Self::validate_index(i,self.rows())?;Self::validate_index(j,self.cols())?;
        Ok(self.opaque_element(i,j))
    }
    fn opaque_element(&self,i:&Expr,j:&Expr)->Expr{Expr::from_kind(crate::expr::Kind::Func(Func::MatrixElement(Arc::new(self.clone())),vec![i.clone(),j.clone()]))}
    pub fn as_explicit(&self)->MResult<Matrix>{
        let (rows,cols)=(finite_dim(self.rows())?,finite_dim(self.cols())?);
        match self.kind(){
            MatrixKind::Explicit(m)=>Ok(m.clone()),MatrixKind::Identity=>Ok(Matrix::eye(rows)),MatrixKind::Zero=>Ok(Matrix::zeros(rows,cols)),
            MatrixKind::One=>Ok(Matrix::from_fn(rows,cols,|_,_|Expr::one())),
            MatrixKind::HadamardProduct(_)|MatrixKind::HadamardPower(..)|MatrixKind::Kronecker(_)|MatrixKind::BlockDiag(_)|MatrixKind::Block(_)|MatrixKind::Diag(_)|MatrixKind::DiagonalOf(_)|MatrixKind::Diagonal(_)|MatrixKind::Slice(..)=>{let mut entries=Vec::with_capacity(rows*cols);for i in 0..rows{for j in 0..cols{entries.push(self.element(&Expr::int(i as i64),&Expr::int(j as i64))?);}}Matrix::new(rows,cols,entries)},
            MatrixKind::Add(a)=>{let mut m=Matrix::zeros(rows,cols);for e in a{m=m.add(&e.as_explicit()?)?;}Ok(m)},
            MatrixKind::Mul(c,fs)=>{let mut m=fs[0].as_explicit()?;for e in &fs[1..]{m=m.matmul(&e.as_explicit()?)?;}Ok(m.scale(c))},
            MatrixKind::Pow(b,e)=>{let n=e.as_i64().ok_or_else(||MatrixError::NotImplemented("Explicit matrix powers require an integer exponent".into()))?;b.as_explicit()?.powi(n)},
            MatrixKind::Inverse(b)=>b.as_explicit()?.inv(),MatrixKind::Transpose(b)=>Ok(b.as_explicit()?.transpose()),MatrixKind::Adjoint(b)=>Ok(b.as_explicit()?.adjoint()),MatrixKind::Conjugate(b)=>Ok(b.as_explicit()?.conjugate()),
            MatrixKind::Symbol(_)=>{let mut entries=Vec::with_capacity(rows*cols);for i in 0..rows{for j in 0..cols{entries.push(self.opaque_element(&Expr::int(i as i64),&Expr::int(j as i64)));}}Matrix::new(rows,cols,entries)},
        }
    }
    pub fn kind_name(&self)->&'static str{match self.kind(){MatrixKind::Symbol(_)=>"MatrixSymbol",MatrixKind::Explicit(_)=>"ImmutableMatrix",MatrixKind::Identity=>"Identity",MatrixKind::Zero=>"ZeroMatrix",MatrixKind::One=>"OneMatrix",MatrixKind::HadamardProduct(_)=>"HadamardProduct",MatrixKind::HadamardPower(..)=>"HadamardPower",MatrixKind::Kronecker(_)=>"KroneckerProduct",MatrixKind::BlockDiag(_)=>"BlockDiagMatrix",MatrixKind::Block(_)=>"BlockMatrix",MatrixKind::Diag(_)=>"DiagMatrix",MatrixKind::DiagonalOf(_)=>"DiagonalOf",MatrixKind::Diagonal(_)=>"DiagonalMatrix",MatrixKind::Slice(..)=>"MatrixSlice",MatrixKind::Add(_)=>"MatAdd",MatrixKind::Mul(..)=>"MatMul",MatrixKind::Pow(..)=>"MatPow",MatrixKind::Transpose(_)=>"Transpose",MatrixKind::Adjoint(_)=>"Adjoint",MatrixKind::Inverse(_)=>"Inverse",MatrixKind::Conjugate(_)=>"MatrixConjugate"}}
    pub fn to_repr(&self)->String{
        let s=crate::printer::srepr;let me=|e:&Self|e.to_repr();
        let args=match self.kind(){
            MatrixKind::Symbol(name)=>format!("{name:?}, {}, {}",s(self.rows()),s(self.cols())),MatrixKind::Identity=>s(self.rows()),MatrixKind::Zero|MatrixKind::One=>format!("{}, {}",s(self.rows()),s(self.cols())),
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)|MatrixKind::BlockDiag(a)=>a.iter().map(me).collect::<Vec<_>>().join(", "),
            MatrixKind::Block(grid)=>format!("Tuple({})",(0..grid.shape().0).map(|i|format!("Tuple({})",(0..grid.shape().1).map(|j|me(grid.get(i,j))).collect::<Vec<_>>().join(", "))).collect::<Vec<_>>().join(", ")),
            MatrixKind::Slice(parent,rows,cols)=>format!("{}, {}, {}",parent.to_repr(),rows.repr(),cols.repr()),
            MatrixKind::HadamardPower(a,b)=>{let operand=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>s(e),MatrixOperand::Matrix(m)=>m.to_repr()};format!("{}, {}",operand(a),operand(b))},
            MatrixKind::Explicit(m)=>format!("{}, {}, [{}]",m.rows(),m.cols(),m.data().iter().map(s).collect::<Vec<_>>().join(", ")),
            MatrixKind::Add(a)=>a.iter().map(me).collect::<Vec<_>>().join(", "),MatrixKind::Mul(c,a)=>{let mut a=a.iter().map(me).collect::<Vec<_>>();if !c.is_one(){a.insert(0,s(c));}a.join(", ")},
            MatrixKind::Pow(b,e)=>format!("{}, {}",b.to_repr(),s(e)),_=>self.children()[0].to_repr(),
        };format!("{}({args})",self.kind_name())
    }
    pub fn to_str(&self)->String{
        let child=|e:&Self|if matches!(e.kind(),MatrixKind::Add(_)){format!("({})",e.to_str())}else{e.to_str()};
        match self.kind(){MatrixKind::Symbol(n)=>n.clone(),MatrixKind::Explicit(m)=>m.to_str(),MatrixKind::Identity=>"I".into(),MatrixKind::Zero=>"0".into(),MatrixKind::One=>"1".into(),
            MatrixKind::HadamardProduct(a)|MatrixKind::Kronecker(a)|MatrixKind::BlockDiag(a)=>format!("{}({})",self.kind_name(),a.iter().map(Self::to_str).collect::<Vec<_>>().join(", ")),
            MatrixKind::Block(grid)=>format!("BlockMatrix([{}])",(0..grid.shape().0).map(|i|format!("[{}]",(0..grid.shape().1).map(|j|grid.get(i,j).to_str()).collect::<Vec<_>>().join(", "))).collect::<Vec<_>>().join(", ")),
            MatrixKind::Diag(parent)|MatrixKind::DiagonalOf(parent)|MatrixKind::Diagonal(parent)=>format!("{}({})",self.kind_name(),parent.to_str()),
            MatrixKind::Slice(parent,rows,cols)=>format!("{}[{}:{}:{}, {}:{}:{}]",parent.to_str(),rows.start,if rows.stop_is_omitted{String::new()}else{rows.stop.to_string()},rows.step,cols.start,if cols.stop_is_omitted{String::new()}else{cols.stop.to_string()},cols.step),
            MatrixKind::HadamardPower(a,b)=>{let operand=|o:&MatrixOperand|match o{MatrixOperand::Scalar(e)=>e.to_string(),MatrixOperand::Matrix(m)=>m.to_str()};format!("HadamardPower({}, {})",operand(a),operand(b))},
            MatrixKind::Add(a)=>{let mut out=String::new();for e in a{let v=e.to_str();if !out.is_empty(){if let Some(positive)=v.strip_prefix('-'){out.push_str(" - ");out.push_str(positive);continue;}out.push_str(" + ");}out.push_str(&v);}out},
            MatrixKind::Mul(c,a)=>{let mut terms:Vec<_>=a.iter().map(child).collect();if !c.is_one(){let scalar=if c.is_add(){format!("({c})")}else{c.to_string()};if c==&Expr::minus_one(){return format!("-{}",terms.join("*"));}terms.insert(0,scalar);}terms.join("*")},
            MatrixKind::Pow(b,e)=>format!("{}**{}",child(b),if e.is_atom()&&!e.as_number().is_some_and(|n|n.is_negative()){e.to_string()}else{format!("({e})")}),MatrixKind::Inverse(b)=>format!("{}**(-1)",child(b)),MatrixKind::Transpose(b)=>format!("{}.T",child(b)),MatrixKind::Adjoint(b)=>format!("Adjoint({})",b.to_str()),MatrixKind::Conjugate(b)=>format!("Adjoint({}.T)",child(b)),
        }
    }
}

impl std::fmt::Display for MatrixExpr { fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result{f.write_str(&self.to_str())} }

/// Block simplification inside scalar expressions, including metadata parents.
pub fn block_collapse_scalar(expression:&Expr)->MResult<Expr>{
    fn rec(e:&Expr,cache:&mut FxHashMap<usize,(Expr,Expr)>)->MResult<Expr>{
        if let Some((_,result))=cache.get(&e.dag_id()){return Ok(result.clone());}
        let result=match e.kind(){
            Kind::Func(Func::MatrixProxy(parent),_)=>Expr::from_kind(Kind::Func(Func::MatrixProxy(Arc::new(parent.block_collapse()?)),vec![])),
            Kind::Func(Func::MatrixTrace(parent),_)=>parent.block_collapse()?.trace()?,
            Kind::Func(Func::MatrixDeterminant(parent),_)=>parent.block_collapse()?.determinant()?,
            Kind::Func(Func::MatrixElement(parent),args)=>parent.block_collapse()?.element_unevaluated(&rec(&args[0],cache)?,&rec(&args[1],cache)?)?,
            _=>{let mut error=None;let result=crate::traverse::map_children(e,&mut |child|match rec(child,cache){Ok(result)=>result,Err(e)=>{error=Some(e);child.clone()}});if let Some(error)=error{return Err(error);}result},
        };cache.insert(e.dag_id(),(e.clone(),result.clone()));Ok(result)
    }rec(expression,&mut FxHashMap::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn symbol(name:&str,rows:i64,cols:i64)->MatrixExpr {MatrixExpr::symbol(name,Expr::int(rows),Expr::int(cols)).unwrap()}
    #[test]
    fn noncommutative_canonicalization() {
        let (a,b)=(symbol("A",2,2),symbol("B",2,2));
        assert_ne!(MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap(),MatrixExpr::mul(Expr::one(),&[b.clone(),a.clone()],true).unwrap());
        assert_eq!(MatrixExpr::add(&[a.clone(),a.clone()],true).unwrap(),MatrixExpr::mul(Expr::int(2),&[a.clone()],true).unwrap());
        assert_eq!(MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone(),b.inverse(true).unwrap()],true).unwrap(),a);
        let p=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        assert_eq!(p.transpose(true).unwrap(),MatrixExpr::mul(Expr::one(),&[b.transpose(true).unwrap(),a.transpose(true).unwrap()],true).unwrap());
        assert!(MatrixExpr::mul(Expr::one(),&[a.clone(),symbol("C",3,2)],true).is_err());
    }
    #[test]
    fn finite_entries_and_explicit_matrices() {
        let a=symbol("A",2,3);let b=symbol("B",3,2);
        let product=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        let value=product.element(&Expr::zero(),&Expr::one()).unwrap();
        let expected=add(&(0..3).map(|k|mul2(&a.element(&Expr::zero(),&Expr::int(k)).unwrap(),&b.element(&Expr::int(k),&Expr::one()).unwrap())).collect::<Vec<_>>());
        assert_eq!(value,expected);
        assert_eq!(product.as_explicit().unwrap(),a.as_explicit().unwrap().matmul(&b.as_explicit().unwrap()).unwrap());
        let vars=crate::traverse::free_symbols(&value);assert_eq!(vars.len(),2);
        assert!(vars.iter().all(|e|matches!(e.kind(),crate::expr::Kind::Func(Func::MatrixProxy(_),_))));
    }
    #[test]
    fn dimension_substitutions_and_single_pass_mapping() {
        let n=Expr::symbol("n");let m=Expr::symbol("m");let p=Expr::symbol("p");
        let a=MatrixExpr::symbol("A",n.clone(),n.clone()).unwrap();
        assert!(a.scalar_free_symbols().is_empty());assert_eq!(a.matrix_symbols(),vec![a.clone()]);
        let map=FxHashMap::from_iter([(n.clone(),m.clone()),(m.clone(),p)]);
        let actual=a.xreplace(&[],&map).unwrap();assert_eq!(actual.rows(),&m);
        assert!(a.as_explicit().is_err());
        let raw=MatrixExpr::mul(Expr::one(),&[a.clone(),MatrixExpr::identity(n).unwrap()],false).unwrap();
        assert_eq!(raw.scalar_map(&mut |e|e.clone()).unwrap(),raw);
        assert_eq!(raw.map_children(&mut |e|e.clone(),&mut |e|e.clone()).unwrap(),raw);
    }
    #[test]
    fn ordered_scalar_derivatives_and_inverse() {
        let x=Expr::symbol("x");let a=symbol("A",2,2);let b=symbol("B",2,2);
        let xa=MatrixExpr::mul(x.clone(),&[a.clone()],true).unwrap();
        let base=MatrixExpr::add(&[xa.clone(),b.clone()],true).unwrap();
        let derivative=base.pow(Expr::int(3),true).unwrap().diff(&x).unwrap().expand().unwrap();
        let expected=MatrixExpr::add(&[
            MatrixExpr::mul(Expr::one(),&[a.clone(),base.pow(Expr::int(2),true).unwrap()],true).unwrap(),
            MatrixExpr::mul(Expr::one(),&[base.clone(),a.clone(),base.clone()],true).unwrap(),
            MatrixExpr::mul(Expr::one(),&[base.pow(Expr::int(2),true).unwrap(),a.clone()],true).unwrap(),
        ],true).unwrap().expand().unwrap();
        let actual=derivative.as_explicit().unwrap();let wanted=expected.as_explicit().unwrap();
        for (a,b) in actual.data().iter().zip(wanted.data()){assert!(crate::expand::expand(&sub(a,b)).is_zero());}
        let inverse=base.inverse(true).unwrap();
        assert_eq!(inverse.diff(&x).unwrap(),MatrixExpr::mul(Expr::minus_one(),&[inverse.clone(),a.clone(),inverse],true).unwrap());
        assert!(matches!(a.diff(&x).unwrap().kind(),MatrixKind::Zero));
        assert!(base.pow(x.clone(),true).unwrap().diff(&x).is_err());
        assert!(base.adjoint(false).unwrap().diff(&x).is_err());
        assert!(xa.adjoint(true).unwrap().diff(&x).is_err());
    }
    #[test]
    fn abstract_trace_determinant_and_explicit_certificates(){
        let a=symbol("A",2,3);let b=symbol("B",3,2);
        let ab=MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap();
        let ba=MatrixExpr::mul(Expr::one(),&[b,a],true).unwrap();
        assert_eq!(ab.trace().unwrap(),ba.trace().unwrap());
        assert!(matches!(ab.determinant().unwrap().kind(),Kind::Func(Func::MatrixDeterminant(_),_)));
        assert_eq!(MatrixExpr::identity(Expr::int(0)).unwrap().determinant().unwrap(),Expr::one());
        assert_eq!(MatrixExpr::zero(Expr::int(0),Expr::int(0)).unwrap().determinant().unwrap(),Expr::one());
        let a=symbol("C",2,2);let b=symbol("D",2,2);
        let product=MatrixExpr::mul(Expr::int(3),&[a.clone(),b.clone()],true).unwrap();
        assert_eq!(product.determinant().unwrap(),crate::ops::mul(&[Expr::int(9),a.determinant().unwrap(),b.determinant().unwrap()]));
        assert_eq!(a.trace_explicit().unwrap(),a.as_explicit().unwrap().trace().unwrap());
        assert_eq!(a.determinant_explicit().unwrap(),a.as_explicit().unwrap().det().unwrap());
    }
    #[test]
    fn large_integer_power_derivatives_share_subexpressions(){
        let x=Expr::symbol("x");let a=symbol("A",2,2);let b=symbol("B",2,2);
        let base=MatrixExpr::add(&[MatrixExpr::mul(x.clone(),&[a],true).unwrap(),b],true).unwrap();
        let result=base.pow(Expr::int(1_000_000_000),true).unwrap().diff(&x).unwrap();
        assert_eq!(result.shape(),(Expr::int(2),Expr::int(2)));
    }
    #[test]
    fn hadamard_operands_canonicalization_and_entries(){
        let a=symbol("A",2,3);let b=symbol("B",2,3);
        let product=MatrixExpr::hadamard_product(&[b.clone(),a.clone(),a.clone()],true).unwrap();
        let reversed=MatrixExpr::hadamard_product(&[a.clone(),b.clone(),a.clone()],true).unwrap();assert_eq!(product,reversed);
        assert_eq!(product.element(&Expr::zero(),&Expr::one()).unwrap(),mul2(&a.element(&Expr::zero(),&Expr::one()).unwrap().powi(2),&b.element(&Expr::zero(),&Expr::one()).unwrap()));
        assert_eq!(MatrixExpr::hadamard_product(&[a.clone(),MatrixExpr::one(Expr::int(2),Expr::int(3)).unwrap()],true).unwrap(),a);
        let power=MatrixExpr::hadamard_power(MatrixOperand::Scalar(Expr::int(2)),MatrixOperand::Matrix(b),false).unwrap();assert_eq!(power.shape(),(Expr::int(2),Expr::int(3)));
        assert_eq!(power.scalar_map(&mut |e|e.clone()).unwrap(),power);
        assert!(MatrixExpr::hadamard_product(&[a,symbol("C",3,2)],false).is_err());
    }
    #[test]
    fn rectangular_kronecker_and_scalar_derivative(){
        let a=symbol("A",2,3);let b=symbol("B",3,2);let x=Expr::symbol("x");
        let product=MatrixExpr::kronecker(&[a.clone(),b.clone()],true).unwrap();assert_eq!(product.shape(),(Expr::int(6),Expr::int(6)));
        assert_eq!(product.element(&Expr::zero(),&Expr::int(3)).unwrap(),mul2(&a.element(&Expr::zero(),&Expr::one()).unwrap(),&b.element(&Expr::zero(),&Expr::one()).unwrap()));
        let scaled=MatrixExpr::kronecker(&[MatrixExpr::mul(x.clone(),&[a],true).unwrap(),b],false).unwrap();assert_eq!(scaled.diff(&x).unwrap(),product);
        assert!(matches!(product.determinant().unwrap().kind(),Kind::Func(Func::MatrixDeterminant(_),_)));
        assert!(product.as_explicit().is_ok());
    }
    #[test]
    fn elementwise_nonfinite_products_keep_explicit_zero_terms(){
        let zero=MatrixExpr::zero(Expr::one(),Expr::one()).unwrap();let infinite=MatrixExpr::explicit(Matrix::new(1,1,vec![Expr::oo()]).unwrap());
        for factors in [[zero.clone(),infinite.clone()],[infinite,zero]]{
            assert_eq!(MatrixExpr::hadamard_product(&factors,true).unwrap().as_explicit().unwrap().get(0,0),&Expr::nan());
            assert_eq!(MatrixExpr::kronecker(&factors,true).unwrap().as_explicit().unwrap().get(0,0),&Expr::nan());
        }
    }
    #[test]
    fn cached_matrix_dag_hash_order_symbols_and_rewrites(){
        use std::collections::hash_map::DefaultHasher;
        fn hash(value:&MatrixExpr)->u64{let mut h=DefaultHasher::new();value.hash(&mut h);h.finish()}
        fn derivative()->MatrixExpr{
            let x=Expr::symbol("x");let a=symbol("A",2,2);let b=symbol("B",2,2);
            MatrixExpr::add(&[MatrixExpr::mul(x.clone(),&[a],true).unwrap(),b],true).unwrap().pow(Expr::int(1_000_000_000),true).unwrap().diff(&x).unwrap()
        }
        let a=derivative();let b=derivative();
        assert_eq!(hash(&a),hash(&b));assert!(a==b);assert_eq!(a.cmp(&b),Ordering::Equal);
        assert_eq!(a.matrix_symbols().len(),2);assert_eq!(a.scalar_free_symbols(),vec![Expr::symbol("x")]);
        assert_eq!(a.scalar_free_symbols(),a.scalar_free_symbols());
        assert_eq!(a.doit(true).unwrap().dag_id(),a.dag_id());
        assert_eq!(a.scalar_map(&mut |e|e.clone()).unwrap().dag_id(),a.dag_id());
        let mut visits=0;a.matrix_map(&mut |m|{visits+=1;m.clone()}).unwrap();assert!(visits<1000);
        let replacements=FxHashMap::from_iter([(Expr::symbol("x"),Expr::symbol("y"))]);
        assert_eq!(a.xreplace(&[],&replacements).unwrap().scalar_free_symbols(),vec![Expr::symbol("y")]);
        let cloned=a.clone();let warm=hash(&a);
        assert_eq!(std::thread::spawn(move||{assert_eq!(cloned.scalar_free_symbols(),vec![Expr::symbol("x")]);hash(&cloned)}).join().unwrap(),warm);
        let n=Expr::symbol("n");let matrix=MatrixExpr::symbol("N",n.clone(),n.clone()).unwrap();
        assert!(matrix.xreplace(&[],&FxHashMap::from_iter([(n,Expr::minus_one())])).is_err());
    }
    #[test]
    fn block_diagonal_rectangles_and_slice_boundaries(){
        let a=symbol("A",2,3);let b=symbol("B",3,2);
        let block=MatrixExpr::block_diag(&[a.clone(),b.clone()]).unwrap();
        assert_eq!(block.shape(),(Expr::int(5),Expr::int(5)));
        assert_eq!(block.element(&Expr::int(1),&Expr::int(3)).unwrap(),Expr::zero());
        assert_eq!(block.element(&Expr::int(3),&Expr::int(4)).unwrap(),b.element(&Expr::one(),&Expr::one()).unwrap());
        assert_eq!(block.determinant().unwrap(),Expr::zero());assert!(block.inverse(true).is_err());
        assert_eq!(block.block_grid().unwrap().get(0,1).shape(),(Expr::int(2),Expr::int(2)));
        let parent=symbol("M",5,6);let slice=parent.slice([Some(Expr::int(2)),Some(Expr::int(4)),None],[Some(Expr::one()),Some(Expr::int(5)),None]).unwrap();
        assert_eq!(slice.shape(),(Expr::int(2),Expr::int(4)));
        assert_eq!(slice.element(&Expr::minus_one(),&Expr::minus_one()).unwrap(),parent.element(&Expr::int(3),&Expr::int(4)).unwrap());
        assert_eq!(parent.slice([Some(Expr::zero()),Some(Expr::zero()),None],[None,None,None]).unwrap().shape(),(Expr::zero(),Expr::int(6)));
        assert_eq!(parent.slice([None,None,Some(Expr::int(2))],[None,None,None]).unwrap().shape(),(Expr::int(3),Expr::int(6)));
        let nested=slice.slice([Some(Expr::one()),None,None],[Some(Expr::one()),None,None]).unwrap();
        assert_eq!(nested.element(&Expr::zero(),&Expr::zero()).unwrap(),parent.element(&Expr::int(3),&Expr::int(2)).unwrap());
    }
    #[test]
    fn cold_shared_matrix_caches_are_thread_safe(){
        use std::hash::{Hash,Hasher};
        let x=Expr::symbol("x");let a=MatrixExpr::mul(x,&[symbol("A",3,3)],true).unwrap().pow(Expr::int(1_000_000_000),true).unwrap();
        let workers:Vec<_>=(0..8).map(|_|{let a=a.clone();std::thread::spawn(move||{let mut h=std::collections::hash_map::DefaultHasher::new();a.hash(&mut h);(h.finish(),a.scalar_free_symbols())})}).collect();
        let answers:Vec<_>=workers.into_iter().map(|w|w.join().unwrap()).collect();assert!(answers.windows(2).all(|w|w[0]==w[1]));
    }

    #[test]
    fn true_block_matrix_and_safe_partition_collapse(){
        let a=symbol("A",2,2);let b=symbol("B",2,1);let c=symbol("C",1,2);let d=symbol("D",1,1);
        let matrix=MatrixExpr::block(BlockGrid::new(2,2,vec![a.clone(),b.clone(),c.clone(),d.clone()]).unwrap()).unwrap();
        assert_eq!(matrix.shape(),(Expr::int(3),Expr::int(3)));
        assert_eq!(matrix.element(&Expr::minus_one(),&Expr::minus_one()).unwrap(),d.element(&Expr::zero(),&Expr::zero()).unwrap());
        assert_eq!(matrix.trace().unwrap(),add(&[a.trace().unwrap(),d.trace().unwrap()]));
        let product=MatrixExpr::mul(Expr::one(),&[matrix.clone(),matrix.clone()],true).unwrap().block_collapse().unwrap();
        let grid=product.block_grid().unwrap();assert_eq!(grid.get(0,1),&MatrixExpr::add(&[MatrixExpr::mul(Expr::one(),&[a.clone(),b.clone()],true).unwrap(),MatrixExpr::mul(Expr::one(),&[b,d.clone()],true).unwrap()],true).unwrap());
        assert_eq!(product.as_explicit().unwrap(),matrix.as_explicit().unwrap().matmul(&matrix.as_explicit().unwrap()).unwrap());
        assert!(MatrixExpr::block(BlockGrid::new(2,2,vec![a.clone(),a.clone(),symbol("E",1,1),symbol("F",1,3)]).unwrap()).is_err());
        let nested=MatrixExpr::block(BlockGrid::new(1,1,vec![a.clone()]).unwrap()).unwrap();assert_eq!(nested.block_collapse().unwrap(),a);
        let power=matrix.pow(Expr::int(8),true).unwrap().block_collapse().unwrap();assert_eq!(power.shape(),matrix.shape());
        let diagonal=MatrixExpr::block_diag(&[symbol("U",1,1),symbol("V",2,2)]).unwrap();
        let unaligned=MatrixExpr::mul(Expr::one(),&[matrix,diagonal],true).unwrap();assert_eq!(unaligned.block_collapse().unwrap(),unaligned);
    }

    #[test]
    fn diagonal_maps_preserve_shapes_indices_and_derivatives(){
        let vector=symbol("v",3,1);let row=symbol("r",1,3);
        for vector in [vector,row]{let diagonal=MatrixExpr::diag(&vector).unwrap();assert_eq!(diagonal.shape(),(Expr::int(3),Expr::int(3)));assert_eq!(diagonal.element(&Expr::one(),&Expr::zero()).unwrap(),Expr::zero());assert_eq!(diagonal.transpose(true).unwrap(),diagonal);assert!(diagonal.as_explicit().is_ok());}
        let a=symbol("A",2,3);let projected=MatrixExpr::diagonal(&a).unwrap();assert_eq!(projected.shape(),a.shape());assert_eq!(projected.element(&Expr::minus_one(),&Expr::minus_one()).unwrap(),Expr::zero());
        let extracted=MatrixExpr::diagonal_of(&a).unwrap();assert_eq!(extracted.shape(),(Expr::int(2),Expr::one()));assert_eq!(extracted.element(&Expr::minus_one(),&Expr::zero()).unwrap(),a.element(&Expr::one(),&Expr::one()).unwrap());
        assert!(MatrixExpr::diag(&a).is_err());
        for(rows,cols)in[(0,0),(0,3),(3,0)]{let a=symbol("Z",rows,cols);assert_eq!(MatrixExpr::diagonal(&a).unwrap().as_explicit().unwrap().shape(),(rows as usize,cols as usize));assert_eq!(MatrixExpr::diagonal_of(&a).unwrap().as_explicit().unwrap().shape(),(0,1));}
        let n=Expr::symbol_with("n",crate::symbol::Assumptions::INTEGER.with(crate::symbol::Assumptions::NONNEGATIVE));let m=Expr::symbol_with("m",crate::symbol::Assumptions::INTEGER.with(crate::symbol::Assumptions::NONNEGATIVE));
        let symbolic=MatrixExpr::diagonal_of(&MatrixExpr::symbol("S",n.clone(),m.clone()).unwrap()).unwrap();assert_eq!(crate::assumptions::is_integer(symbolic.rows()),Some(true));
        let x=Expr::symbol("x");let scaled=MatrixExpr::mul(x.clone(),&[a.clone()],true).unwrap();assert_eq!(MatrixExpr::diagonal(&scaled).unwrap().diff(&x).unwrap(),MatrixExpr::diagonal(&a).unwrap());
    }

    #[test]
    fn slices_match_python_integer_ranges_and_composition(){
        // Independent imperative list-slice oracle, without expression math.
        fn selected(n:i64,a:Option<i64>,b:Option<i64>,h:i64)->Vec<i64>{
            let bound=|v:i64|{let v=if v<0{v+n}else{v};if h>0{v.clamp(0,n)}else{v.clamp(-1,n-1)}};
            let mut i=a.map(bound).unwrap_or(if h>0{0}else{n-1});
            let stop=b.map(bound).unwrap_or(if h>0{n}else{-1});let mut out=Vec::new();
            while if h>0{i<stop}else{i>stop}{out.push(i);i+=h;}out
        }
        let huge=crate::parser::parse_str("10000000000000000000000000000000000000000000").unwrap();
        assert!(SliceSpec::index(huge.clone(),&Expr::int(5)).is_err());
        assert!(SliceSpec::index(neg(&huge),&Expr::int(5)).is_err());
        let bounds=[None,Some(-10),Some(-2),Some(-1),Some(0),Some(1),Some(2),Some(10)];
        for n in 0..7{let a=symbol("A",n,2);for start in bounds{for stop in bounds{for step in [-3,-2,-1,1,2,3]{
            let rows=[start.map(Expr::int),stop.map(Expr::int),Some(Expr::int(step))];
            let view=a.slice(rows,[None,None,None]).unwrap();let expected=selected(n,start,stop,step);
            assert_eq!(view.rows().as_i64(),Some(expected.len()as i64),"n={n},start={start:?},stop={stop:?},step={step}");
            for(i,k)in expected.iter().enumerate(){assert_eq!(view.element(&Expr::int(i as i64),&Expr::zero()).unwrap(),a.element(&Expr::int(*k),&Expr::zero()).unwrap());}
            if let MatrixKind::Slice(_,r,c)=view.kind(){assert_eq!(view,a.slice_normalized(SliceSpec::normalized(r.start.clone(),r.stop.clone(),r.step.clone(),r.length.clone(),r.stop_is_omitted).unwrap(),c.clone()).unwrap());}
            let reversed=view.slice([None,None,Some(Expr::minus_one())],[None,None,None]).unwrap();
            assert_eq!(reversed.rows().as_i64(),Some(expected.len()as i64));
            for(i,k)in expected.iter().rev().enumerate(){assert_eq!(reversed.element(&Expr::int(i as i64),&Expr::zero()).unwrap(),a.element(&Expr::int(*k),&Expr::zero()).unwrap());}
        }}}}
    }
    #[test]
    fn symbolic_slice_lengths_and_integer_signs(){
        use crate::symbol::Assumptions as A;
        let n=Expr::symbol_with("n",A::INTEGER.with(A::NONNEGATIVE));
        let k=Expr::symbol_with("k",A::INTEGER);let h=Expr::symbol_with("h",A::INTEGER.with(A::POSITIVE));
        let a=MatrixExpr::symbol("A",n.clone(),Expr::int(4)).unwrap();
        for step in [h.clone(),neg(&h)]{let view=a.slice([Some(k.clone()),None,Some(step)],[None,None,None]).unwrap();
            assert_eq!(crate::assumptions::is_integer(view.rows()),Some(true));
            assert_eq!(crate::assumptions::is_nonnegative(view.rows()),Some(true));
            let concrete=view.xreplace(&[],&FxHashMap::from_iter([(n.clone(),Expr::int(5)),(k.clone(),Expr::int(-2)),(h.clone(),Expr::int(2))])).unwrap();
            assert!(concrete.as_explicit().is_ok());
        }
        assert!(a.slice([None,None,Some(k)],[None,None,None]).is_err());
    }

}
