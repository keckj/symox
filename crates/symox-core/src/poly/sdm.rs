//! Sparse exact algorithms: division-free Berkowitz and field Gauss--Jordan.
//! No n-by-n storage is constructed. Isolated zero rows AND columns can be
//! removed for Berkowitz: a simultaneous basis permutation splits off a zero direct summand,
//! multiplying the characteristic polynomial by a power of lambda. Its sparse
//! coefficient indices, counted from the leading monic term, do not change.
use std::collections::{BTreeMap, BTreeSet};
use super::Coeff;

pub type SparseVector<C> = BTreeMap<usize,C>;
pub type SparseRows<C> = BTreeMap<usize,SparseVector<C>>;

fn accumulate<C:Coeff>(values:&mut SparseVector<C>,index:usize,value:C) {
    if value.is_zero(){return;}
    let value=match values.remove(&index){Some(old)=>old.add(&value),None=>value};
    if !value.is_zero(){values.insert(index,value);}
}

fn dot<C:Coeff,E>(row:&SparseVector<C>,vector:&SparseVector<C>,start:usize,
    check:&mut impl FnMut()->Result<(),E>)->Result<C,E> {
    let mut result=C::zero();
    for (step,(&j,a)) in row.range(start..).enumerate() {
        if step%64==0{check()?;}
        if let Some(b)=vector.get(&j){result.add_mul(a,b);}
    }
    Ok(result)
}

fn matvec<C:Coeff,E>(rows:&[SparseVector<C>],vector:&SparseVector<C>,start:usize,
    check:&mut impl FnMut()->Result<(),E>)->Result<SparseVector<C>,E> {
    let mut result=SparseVector::new();
    for (i,row) in rows.iter().enumerate().skip(start) {
        check()?;
        let value=dot(row,vector,start,check)?;
        if !value.is_zero(){result.insert(i,value);}
    }
    Ok(result)
}

/// Multiply an already computed sparse coefficient vector by a lower
/// triangular Toeplitz column. The same routine can finish a legacy prefix:
/// q and t are supplied results, never reconstructed from the input matrix.
pub fn toeplitz_product_checked<C:Coeff,E>(q:&SparseVector<C>,t:&SparseVector<C>,degree:usize,
    check:&mut impl FnMut()->Result<(),E>)->Result<SparseVector<C>,E> {
    check()?;
    let mut next=SparseVector::new();
    for (&i,a) in q {
        check()?;
        for (step,(&j,b)) in t.iter().enumerate() {
            if step%64==0{check()?;}
            if i<=degree && j<=degree-i{accumulate(&mut next,i+j,a.mul(b));}
        }
    }
    check()?;
    Ok(next)
}

/// Sparse leading-coefficient map of det(lambda*I-M). Missing entries are zero.
/// None declines malformed coordinates; an interruption is a distinct Err.
/// The historical dimension-zero case returns one without reading entries.
/// Input rows/coefficients are immutable and are never published partially.
pub fn berkowitz_checked<C:Coeff,E>(matrix:&SparseRows<C>,n:usize,
    check:&mut impl FnMut()->Result<(),E>)->Result<Option<SparseVector<C>>,E> {
    check()?;
    if n==0{return Ok(Some(BTreeMap::from([(0,C::one())])));}
    let mut active=BTreeSet::new();
    for (&i,row) in matrix {
        check()?;
        if i>=n{return Ok(None);}
        for (step,(&j,value)) in row.iter().enumerate() {
            if step%64==0{check()?;}
            if j>=n{return Ok(None);}
            if !value.is_zero(){active.insert(i);active.insert(j);}
        }
    }
    let indices:BTreeMap<_,_>=active.into_iter().enumerate().map(|(new,old)|(old,new)).collect();
    let size=indices.len();
    // O(active vertices), not O(the declared dimension). At most two active
    // vertices occur per stored nonzero. Column views also stay sparse.
    let mut rows=vec![SparseVector::new();size];
    let mut columns=vec![SparseVector::new();size];
    for (&i,row) in matrix {
        check()?;
        for (step,(&j,value)) in row.iter().enumerate() {
            if step%64==0{check()?;}
            if !value.is_zero(){
                let (i,j)=(indices[&i],indices[&j]);
                rows[i].insert(j,value.clone());columns[j].insert(i,value.clone());
            }
        }
    }
    let mut q=BTreeMap::from([(0,C::one())]);
    // M=[[a,R],[C,A]]. The Toeplitz first column is
    // (1,-a,-R*C,-R*A*C,...). Multiplying it by charpoly(A) gives charpoly(M).
    // The suffix is a view of the same sparse matrix, rather than a copied A.
    for pivot in (0..size).rev() {
        check()?;
        let tail=size-pivot-1;
        let mut t=BTreeMap::from([(0,C::one())]);
        if let Some(a)=rows[pivot].get(&pivot){t.insert(1,a.neg());}
        let r=&rows[pivot];
        let mut vector:SparseVector<C>=columns[pivot].range((pivot+1)..)
            .map(|(&i,value)|(i,value.clone())).collect();
        if r.range((pivot+1)..).next().is_some() {
            for power in 0..tail {
                check()?;
                if vector.is_empty(){break;}
                let value=dot(r,&vector,pivot+1,check)?;
                if !value.is_zero(){t.insert(power+2,value.neg());}
                if power+1<tail{vector=matvec(&rows,&vector,pivot+1,check)?;}
            }
        }
        // Sparse Toeplitz convolution, without constructing the matrix T.
        q=toeplitz_product_checked(&q,&t,tail+1,check)?;
    }
    check()?;
    Ok(Some(q))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{Int,Rat};
    fn sparse<C:Coeff>(a:&[Vec<C>])->SparseRows<C> {
        a.iter().enumerate().filter_map(|(i,row)|{
            let row=row.iter().enumerate().filter(|(_,v)|!v.is_zero()).map(|(j,v)|(j,v.clone())).collect::<BTreeMap<_,_>>();
            (!row.is_empty()).then_some((i,row))
        }).collect()
    }
    fn rational(n:i64,d:i64)->Rat{Rat::new(Int::from(n),Int::from(d))}
    #[test]
    fn complete_exact_coefficients_match_the_existing_dense_kernel() {
        use crate::matrix::{domain::RatMat,int::{Int as MI,Rat as MR}};
        for n in 0..8 {
            for seed in 0..12 {
                let mut a=vec![vec![Rat::zero();n];n];
                for i in 0..n {for j in 0..n {
                    if (i*17+j*11+seed)%4==0{a[i][j]=rational(((i*13+j*7+seed)%19) as i64-9,((i+j+seed)%7) as i64+1);}
                }}
                let flat=a.iter().flatten().map(|v|MR{num:MI::from_big(v.n.to_big()),den:MI::from_big(v.d.to_big())}).collect();
                let expected=RatMat{rows:n,cols:n,data:flat}.charpoly();
                let actual=berkowitz_checked(&sparse(&a),n,&mut ||Ok::<_,()>(())).unwrap().unwrap();
                for (i,v) in expected.into_iter().enumerate(){
                    assert_eq!(actual.get(&i).cloned().unwrap_or_else(Rat::zero),Rat::new(Int::from(v.num.to_big()),Int::from(v.den.to_big())));
                }
            }
        }
    }
    #[test]
    fn huge_declared_dimension_does_not_allocate_dense_storage() {
        let n=usize::MAX;
        let a=BTreeMap::from([
            (n-9,BTreeMap::from([(n-9,Int::from(2)),(n-2,Int::from(3))])),
            (n-2,BTreeMap::from([(n-9,Int::from(5)),(n-2,Int::from(7))])),
        ]);
        assert_eq!(berkowitz_checked(&a,n,&mut ||Ok::<_,()>(())).unwrap().unwrap(),
            BTreeMap::from([(0,Int::one()),(1,Int::from(-9)),(2,Int::from(-1))]));
        assert_eq!(berkowitz_checked::<Int,_>(&BTreeMap::new(),n,&mut ||Ok::<_,()>(())).unwrap().unwrap(),BTreeMap::from([(0,Int::one())]));
    }
    #[test]
    fn zero_cancellations_nilpotence_alias_free_inputs_and_interruptions() {
        let a=BTreeMap::from([(0,BTreeMap::from([(1,Int::from(2))])),(1,BTreeMap::from([(2,Int::from(3))]))]);
        let original=a.clone();
        assert_eq!(berkowitz_checked(&a,3,&mut ||Ok::<_,()>(())).unwrap().unwrap(),BTreeMap::from([(0,Int::one())]));
        let mut calls=0;
        let result=berkowitz_checked(&a,3,&mut ||{calls+=1;if calls==13{Err("interrupted")}else{Ok(())}});
        assert_eq!(result,Err("interrupted"));assert_eq!(a,original);
        let invalid=BTreeMap::from([(2,BTreeMap::from([(0,Int::one())]))]);
        assert!(berkowitz_checked(&invalid,2,&mut ||Ok::<_,()>(())).unwrap().is_none());
    }
    #[test]
    fn sparse_prefix_convolution_preserves_gaps_truncation_and_overflow_safety() {
        let q=BTreeMap::from([(0,Int::one()),(2,Int::from(-3))]);
        let t=BTreeMap::from([(0,Int::one()),(1,Int::from(2)),(3,Int::from(5))]);
        assert_eq!(toeplitz_product_checked(&q,&t,4,&mut ||Ok::<_,()>(())).unwrap(),
            BTreeMap::from([(0,Int::one()),(1,Int::from(2)),(2,Int::from(-3)),(3,Int::from(-1))]));
        let edge=BTreeMap::from([(usize::MAX,Int::one())]);
        assert_eq!(toeplitz_product_checked(&edge,&t,usize::MAX,&mut ||Ok::<_,()>(())).unwrap(),edge);
        assert!(toeplitz_product_checked::<Int,_>(&BTreeMap::new(),&t,4,&mut ||Ok::<_,()>(())).unwrap().is_empty());
        let mut calls=0;
        assert_eq!(toeplitz_product_checked(&q,&t,4,&mut ||{calls+=1;if calls==3{Err("stop")}else{Ok(())}}),Err("stop"));
    }
}


// Field Gauss--Jordan is deliberately separate from the division-free
// Berkowitz above and the fraction-free elimination in fraction_free.rs.
// Opaque origins preserve representative key objects without algebra in Python.
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub struct IrrefKey { pub index:usize,pub origin:usize }
pub type IrrefRow<C> = super::fraction_free::OrderedMap<(usize,C)>;
#[derive(Clone,Copy,Debug)]
pub enum IrrefSetOp { Difference,Intersection }

/// Representation only. A host can maintain its original persistent sets;
/// every returned membership/order is checked against independent Rust sets.
pub trait IrrefSets<E> {
    fn create(&mut self,keys:&[IrrefKey],from_dictionary:bool)->Result<usize,E>;
    fn release(&mut self,_set:usize)->Result<(),E>{Ok(())}
    fn insert(&mut self,set:usize,key:IrrefKey)->Result<(),E>;
    fn remove(&mut self,set:usize,key:IrrefKey)->Result<(),E>;
    fn members(&mut self,set:usize)->Result<Vec<IrrefKey>,E>;
    fn binary(&mut self,op:IrrefSetOp,left:usize,right:usize)->Result<Vec<IrrefKey>,E>;
}
struct TrackedSet { id:usize,keys:BTreeSet<usize> }
fn tracked<E>(sets:&mut impl IrrefSets<E>,keys:&[IrrefKey],from_dictionary:bool)->Result<TrackedSet,E> {
    Ok(TrackedSet{id:sets.create(keys,from_dictionary)?,keys:keys.iter().map(|k|k.index).collect()})
}
fn checked_members(keys:Vec<IrrefKey>,expected:&BTreeSet<usize>)->Option<Vec<IrrefKey>> {
    let actual=keys.iter().map(|k|k.index).collect::<BTreeSet<_>>();
    (actual.len()==keys.len()&&&actual==expected).then_some(keys)
}
fn set_members<E>(sets:&mut impl IrrefSets<E>,set:&TrackedSet)->Result<Option<Vec<IrrefKey>>,E> {
    Ok(checked_members(sets.members(set.id)?,&set.keys))
}
fn set_binary<E>(sets:&mut impl IrrefSets<E>,op:IrrefSetOp,left:&TrackedSet,right:&TrackedSet)->Result<Option<Vec<IrrefKey>>,E> {
    let expected=match op {
        IrrefSetOp::Difference=>left.keys.difference(&right.keys).copied().collect(),
        IrrefSetOp::Intersection=>left.keys.intersection(&right.keys).copied().collect(),
    };
    Ok(checked_members(sets.binary(op,left.id,right.id)?,&expected))
}
fn set_insert<E>(sets:&mut impl IrrefSets<E>,set:&mut TrackedSet,key:IrrefKey)->Result<(),E> {
    sets.insert(set.id,key)?;set.keys.insert(key.index);Ok(())
}
fn set_remove<E>(sets:&mut impl IrrefSets<E>,set:&mut TrackedSet,key:IrrefKey)->Result<bool,E> {
    if !set.keys.remove(&key.index){return Ok(false);}
    sets.remove(set.id,key)?;Ok(true)
}
fn irref_keys<C>(row:&IrrefRow<C>)->Vec<IrrefKey> {
    row.iter().map(|(index,(origin,_))|IrrefKey{index,origin:*origin}).collect()
}
fn irref_store<C>(row:&mut IrrefRow<C>,key:IrrefKey,value:C) {
    let origin=row.get(key.index).map_or(key.origin,|(origin,_)|*origin);
    row.insert(key.index,(origin,value));
}
fn column_set<E>(sets:&mut impl IrrefSets<E>,columns:&mut super::fraction_free::OrderedMap<(usize,usize)>,
    storage:&mut Vec<TrackedSet>,key:IrrefKey)->Result<usize,E> {
    if let Some((_,index))=columns.get(key.index){return Ok(*index);}
    let index=storage.len();storage.push(tracked(sets,&[],false)?);columns.insert(key.index,(key.origin,index));Ok(index)
}
pub struct IrrefResult<C> {
    pub rows:Vec<IrrefRow<C>>,
    pub pivots:Vec<IrrefKey>,
    pub nonzero_columns:Vec<(IrrefKey,Vec<usize>)>,
}

/// Complete sparse field Gauss--Jordan, mirroring sdm_irref's actual sequence:
/// stable first-column sort then pop, previous-row cancellation, inverse/pivot
/// normalization and cancellation into every previous nonzero pivot row.
/// The third return retains insertion order and even historical empty sets.
/// None means a representation adapter failed validation, never rank zero.
pub fn irref<C:super::Field,E>(mut source:Vec<IrrefRow<C>>,sets:&mut impl IrrefSets<E>,
    check:&mut impl FnMut()->Result<(),E>)->Result<Option<IrrefResult<C>>,E> {
    check()?;
    if source.iter().any(|r|r.is_empty()||r.iter().any(|(_,(_,v))|v.is_zero())){return Ok(None);}
    source.sort_by_key(|row|row.iter().map(|(j,_)|j).min().unwrap());
    let mut pivot_rows:BTreeMap<usize,(IrrefKey,IrrefRow<C>)>=BTreeMap::new();
    let mut reduced=tracked(sets,&[],false)?;let mut unreduced=tracked(sets,&[],false)?;
    let mut columns:super::fraction_free::OrderedMap<(usize,usize)>=Default::default();
    let mut storage=Vec::new();
    while let Some(source)=source.pop() {
        check()?;
        let mut a:IrrefRow<C>=Default::default();
        for (j,(origin,value)) in source.iter() {
            if !reduced.keys.contains(&j){a.insert(j,(*origin,value.clone()));}
        }
        let snapshot=tracked(sets,&irref_keys(&a),true)?;
        let Some(cancel)=set_binary(sets,IrrefSetOp::Intersection,&unreduced,&snapshot)? else{return Ok(None)};
        sets.release(snapshot.id)?;drop(snapshot);
        for j in cancel {
            check()?;
            let Some((_,coefficient))=a.get(j.index).cloned() else{return Ok(None)};
            let Some((_,previous))=pivot_rows.get(&j.index) else{return Ok(None)};
            let mut an=tracked(sets,&irref_keys(&a),true)?;
            let pn=tracked(sets,&irref_keys(previous),true)?;
            let Some(new)=set_binary(sets,IrrefSetOp::Difference,&pn,&an)? else{return Ok(None)};
            for k in new {
                check()?;
                let Some((_,value))=previous.get(k.index) else{return Ok(None)};
                irref_store(&mut a,k,coefficient.neg().mul(value));
            }
            a.remove(j.index);
            if !set_remove(sets,&mut an,j)?{return Ok(None);}
            let Some(shared)=set_binary(sets,IrrefSetOp::Intersection,&pn,&an)? else{return Ok(None)};
            for k in shared {
                check()?;
                let (Some((_,old)),Some((_,previous)))=(a.get(k.index),previous.get(k.index)) else{return Ok(None)};
                let value=old.sub(&coefficient.mul(previous));
                if value.is_zero(){a.remove(k.index);}else{irref_store(&mut a,k,value);}
            }
            sets.release(an.id)?;sets.release(pn.id)?;
        }
        if a.is_empty(){continue;}
        let (j,(origin,pivot))=a.iter().min_by_key(|(j,_)|*j).unwrap();
        let j=IrrefKey{index:j,origin:*origin};let inverse=pivot.inv();
        if pivot_rows.contains_key(&j.index){return Ok(None);}
        let an=tracked(sets,&irref_keys(&a),true)?;
        for key in irref_keys(&a) {
            check()?;let value=a.get(key.index).unwrap().1.mul(&inverse);irref_store(&mut a,key,value);
        }
        if let Some((_,index))=columns.remove(j.index) {
            let Some(affected)=set_members(sets,&storage[index])? else{return Ok(None)};
            sets.release(storage[index].id)?;storage[index].keys.clear();
            for k in affected {
                check()?;
                let Some((_,previous))=pivot_rows.get_mut(&k.index) else{return Ok(None)};
                let Some((_,coefficient))=previous.get(j.index).cloned() else{return Ok(None)};
                let mut pn=tracked(sets,&irref_keys(previous),true)?;
                let Some(new)=set_binary(sets,IrrefSetOp::Difference,&an,&pn)? else{return Ok(None)};
                for l in new {
                    check()?;let Some((_,value))=a.get(l.index) else{return Ok(None)};
                    irref_store(previous,l,coefficient.neg().mul(value));
                    let column=column_set(sets,&mut columns,&mut storage,l)?;
                    set_insert(sets,&mut storage[column],k)?;
                }
                previous.remove(j.index);
                if !set_remove(sets,&mut pn,j)?{return Ok(None);}
                let Some(shared)=set_binary(sets,IrrefSetOp::Intersection,&an,&pn)? else{return Ok(None)};
                for l in shared {
                    check()?;
                    let (Some((_,old)),Some((_,value)))=(previous.get(l.index),a.get(l.index)) else{return Ok(None)};
                    let value=old.sub(&coefficient.mul(value));
                    if !value.is_zero(){irref_store(previous,l,value);}
                    else {
                        previous.remove(l.index);
                        if l.index!=j.index {
                            let column=column_set(sets,&mut columns,&mut storage,l)?;
                            if !set_remove(sets,&mut storage[column],k)?{return Ok(None);}
                        }
                    }
                }
                sets.release(pn.id)?;
                if previous.len()==1 {
                    set_insert(sets,&mut reduced,k)?;
                    if !set_remove(sets,&mut unreduced,k)?{return Ok(None);}
                }
            }
        }
        sets.release(an.id)?;
        if a.len()==1 {set_insert(sets,&mut reduced,j)?;}
        else {
            set_insert(sets,&mut unreduced,j)?;
            for l in irref_keys(&a) {
                if l.index!=j.index {
                    let column=column_set(sets,&mut columns,&mut storage,l)?;
                    set_insert(sets,&mut storage[column],j)?;
                }
            }
        }
        pivot_rows.insert(j.index,(j,a));
    }
    let union=reduced.keys.union(&unreduced.keys).copied().collect::<Vec<_>>();
    if union!=pivot_rows.keys().copied().collect::<Vec<_>>(){return Ok(None);}
    let pivots=pivot_rows.values().map(|(key,_)|*key).collect::<Vec<_>>();
    let numbering=pivots.iter().enumerate().map(|(i,key)|(key.index,i)).collect::<BTreeMap<_,_>>();
    let mut nonzero_columns=Vec::new();
    for (column,(origin,index)) in columns.iter() {
        check()?;
        let Some(members)=set_members(sets,&storage[*index])? else{return Ok(None)};
        let Some(indices)=members.iter().map(|key|numbering.get(&key.index).copied()).collect::<Option<Vec<_>>>() else{return Ok(None)};
        nonzero_columns.push((IrrefKey{index:column,origin:*origin},indices));
    }
    check()?;
    Ok(Some(IrrefResult{rows:pivot_rows.into_values().map(|(_,row)|row).collect(),pivots,nonzero_columns}))
}

#[cfg(test)]
mod irref_tests {
    use super::*;
    use crate::poly::{Int,Rat};
    #[derive(Default)]
    struct Repr(Vec<BTreeMap<usize,usize>>);
    impl IrrefSets<()> for Repr {
        fn create(&mut self,keys:&[IrrefKey],_:bool)->Result<usize,()> {
            let index=self.0.len();self.0.push(keys.iter().map(|k|(k.index,k.origin)).collect());Ok(index)
        }
        fn insert(&mut self,set:usize,key:IrrefKey)->Result<(),()> {self.0[set].entry(key.index).or_insert(key.origin);Ok(())}
        fn remove(&mut self,set:usize,key:IrrefKey)->Result<(),()> {self.0[set].remove(&key.index);Ok(())}
        fn members(&mut self,set:usize)->Result<Vec<IrrefKey>,()> {
            Ok(self.0[set].iter().map(|(&index,&origin)|IrrefKey{index,origin}).collect())
        }
        fn binary(&mut self,op:IrrefSetOp,left:usize,right:usize)->Result<Vec<IrrefKey>,()> {
            Ok(self.0[left].iter().filter(|(key,_)|match op {
                IrrefSetOp::Difference=>!self.0[right].contains_key(key),
                IrrefSetOp::Intersection=>self.0[right].contains_key(key),
            }).map(|(&index,&origin)|IrrefKey{index,origin}).collect())
        }
    }
    fn q(n:i64)->Rat{Rat::from(Int::from(n))}
    fn row(values:&[(usize,usize,i64)])->IrrefRow<Rat> {
        let mut result=IrrefRow::default();for &(index,origin,n)in values{result.insert(index,(origin,q(n)));}result
    }
    #[test]
    fn exact_reduction_matches_independent_dense_elimination() {
        let mut state=3719u64;
        for m in 0..6 {for n in 0..7 {for _ in 0..5 {
            let mut dense=vec![vec![Rat::zero();n];m];let mut sparse=Vec::new();
            for (i,values)in dense.iter_mut().enumerate(){
                let mut row=IrrefRow::default();
                for (j,value)in values.iter_mut().enumerate(){
                    state=state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let v=((state>>32)%9)as i64-4;*value=q(v);
                    if v!=0{row.insert(j,(i*128+j,value.clone()));}
                }
                if !row.is_empty(){sparse.push(row);}
            }
            let result=irref(sparse,&mut Repr::default(),&mut ||Ok::<_,()>(())).unwrap().unwrap();
            let pivots=crate::poly::ddm::irref(&mut dense,false,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap();
            assert_eq!(result.pivots.iter().map(|k|k.index).collect::<Vec<_>>(),pivots);
            for (i,row)in result.rows.iter().enumerate(){for j in 0..n{
                assert_eq!(row.get(j).map(|(_,v)|v.clone()).unwrap_or_else(Rat::zero),dense[i][j]);
            }}
            for (key,indices)in result.nonzero_columns{
                let expected=result.rows.iter().enumerate().filter_map(|(i,row)|row.get(key.index).map(|_|i)).collect::<BTreeSet<_>>();
                assert_eq!(indices.into_iter().collect::<BTreeSet<_>>(),expected);
            }
        }}}
    }
    #[test]
    fn stable_sort_then_pop_preserves_representative_key_origins() {
        let rows=vec![row(&[(1000,10,1),(2000,11,2)]),row(&[(1000,20,1),(2000,21,2)])];
        let result=irref(rows,&mut Repr::default(),&mut ||Ok::<_,()>(())).unwrap().unwrap();
        assert_eq!(result.pivots,vec![IrrefKey{index:1000,origin:20}]);
        assert_eq!(irref_keys(&result.rows[0]),vec![IrrefKey{index:1000,origin:20},IrrefKey{index:2000,origin:21}]);
        assert_eq!(result.nonzero_columns,vec![(IrrefKey{index:2000,origin:21},vec![0])]);
    }
    #[test]
    fn empty_rows_explicit_zeros_and_interruption_do_not_mean_rank_zero() {
        assert!(irref::<Rat,_>(vec![IrrefRow::default()],&mut Repr::default(),&mut ||Ok::<_,()>(())).unwrap().is_none());
        assert!(irref(vec![row(&[(0,0,0)])],&mut Repr::default(),&mut ||Ok::<_,()>(())).unwrap().is_none());
        assert!(irref(vec![row(&[(0,0,1)])],&mut Repr::default(),&mut ||Err(())).is_err());
        let result=irref::<Rat,_>(vec![],&mut Repr::default(),&mut ||Ok::<_,()>(())).unwrap().unwrap();
        assert!(result.rows.is_empty()&&result.pivots.is_empty()&&result.nonzero_columns.is_empty());
    }
    #[test]
    fn representation_oracle_cannot_invent_a_pivot() {
        struct Invalid(Repr);
        impl IrrefSets<()> for Invalid {
            fn create(&mut self,k:&[IrrefKey],d:bool)->Result<usize,()>{self.0.create(k,d)}
            fn insert(&mut self,s:usize,k:IrrefKey)->Result<(),()>{self.0.insert(s,k)}
            fn remove(&mut self,s:usize,k:IrrefKey)->Result<(),()>{self.0.remove(s,k)}
            fn members(&mut self,s:usize)->Result<Vec<IrrefKey>,()>{self.0.members(s)}
            fn binary(&mut self,_:IrrefSetOp,_:usize,_:usize)->Result<Vec<IrrefKey>,()>{Ok(vec![IrrefKey{index:99,origin:99}])}
        }
        assert!(irref(vec![row(&[(0,0,1)])],&mut Invalid(Repr::default()),&mut ||Ok::<_,()>(())).unwrap().is_none());
    }
}

#[cfg(test)]
mod irref_storage_tests {
    use super::*;
    use crate::poly::{Int,Rat};
    #[derive(Default)]
    struct Sets { data:Vec<BTreeMap<usize,usize>>,live:BTreeSet<usize>,peak:usize }
    impl IrrefSets<()> for Sets {
        fn create(&mut self,keys:&[IrrefKey],_:bool)->Result<usize,()> {
            let id=self.data.len();self.data.push(keys.iter().map(|k|(k.index,k.origin)).collect());
            self.live.insert(id);self.peak=self.peak.max(self.live.len());Ok(id)
        }
        fn release(&mut self,id:usize)->Result<(),()> {
            assert!(self.live.remove(&id));self.data[id].clear();Ok(())
        }
        fn insert(&mut self,id:usize,key:IrrefKey)->Result<(),()> {
            assert!(self.live.contains(&id));self.data[id].entry(key.index).or_insert(key.origin);Ok(())
        }
        fn remove(&mut self,id:usize,key:IrrefKey)->Result<(),()> {
            assert!(self.live.contains(&id));assert!(self.data[id].remove(&key.index).is_some());Ok(())
        }
        fn members(&mut self,id:usize)->Result<Vec<IrrefKey>,()> {
            assert!(self.live.contains(&id));Ok(self.data[id].iter().map(|(&index,&origin)|IrrefKey{index,origin}).collect())
        }
        fn binary(&mut self,op:IrrefSetOp,left:usize,right:usize)->Result<Vec<IrrefKey>,()> {
            assert!(self.live.contains(&left)&&self.live.contains(&right));
            Ok(self.data[left].iter().filter(|(key,_)|match op {
                IrrefSetOp::Difference=>!self.data[right].contains_key(key),
                IrrefSetOp::Intersection=>self.data[right].contains_key(key),
            }).map(|(&index,&origin)|IrrefKey{index,origin}).collect())
        }
    }
    #[test]
    fn temporary_set_storage_is_bounded_by_live_columns_not_operation_count() {
        let (m,n)=(28,35);let mut rows=Vec::new();
        for i in 0..m {
            let mut row=IrrefRow::default();
            for j in 0..n {
                let value=((i*31+j*17+i*j*3)%29+1)as i64;
                row.insert(j,(i*n+j,Rat::from(Int::from(value))));
            }
            rows.push(row);
        }
        let mut sets=Sets::default();
        let result=irref(rows,&mut sets,&mut ||Ok::<_,()>(())).unwrap().unwrap();
        assert!(sets.data.len()>n*3);
        assert!(sets.peak<=n+6,"peak {}",sets.peak);
        assert_eq!(sets.live.len(),2+result.nonzero_columns.len());
    }
}
