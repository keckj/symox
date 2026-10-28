//! Fraction-free reductions with their distinct dense/sparse denominator rules.
use super::Coeff;

#[derive(Debug)]
pub enum Failure<E,C> { Interrupted(E), NotExact(C,C) }

fn quotient<C:Coeff,E>(value:C,divisor:Option<&C>)->Result<C,Failure<E,C>> {
    let Some(divisor)=divisor else{return Ok(value)};
    value.div_exact(divisor).ok_or_else(||Failure::NotExact(value,divisor.clone()))
}

pub enum DenseEvent<'a,C> { Swap(usize,usize), Set(usize,usize,&'a C,Option<usize>) }

/// Faithful dense fraction-free Gauss-Jordan over an integral domain. Shared pivot assignment tokens
/// preserve the caller's object aliasing; all published writes remain ordered.
pub fn dense<C:Coeff,E>(a:&mut [Vec<C>],check:&mut impl FnMut()->Result<(),E>,
    publish:&mut impl FnMut(DenseEvent<'_,C>)->Result<(),E>)->Result<Option<(C,Vec<usize>)>,Failure<E,C>> {
    check().map_err(Failure::Interrupted)?;
    let m=a.len();let Some(first)=a.first() else{return Ok(Some((C::one(),vec![])))};
    let n=first.len();
    if a.iter().any(|row|row.len()!=n){return Ok(None);}
    let mut divisor=None;let mut pivots=Vec::new();let mut no_pivots=Vec::new();let mut i=0;
    for j in 0..n {
        check().map_err(Failure::Interrupted)?;
        if a[i][j].is_zero() {
            let mut found=None;
            for (ip,row) in a.iter().enumerate().skip(i+1) {
                check().map_err(Failure::Interrupted)?;
                if !row[j].is_zero(){found=Some(ip);break;}
            }
            let Some(ip)=found else{no_pivots.push(j);continue};
            publish(DenseEvent::Swap(i,ip)).map_err(Failure::Interrupted)?;a.swap(i,ip);
        }
        let pivot=a[i][j].clone();
        if let Some(&first)=pivots.first() {
            let value=quotient(pivot.mul(&a[0][first]),divisor.as_ref())?;
            for (ip,&jp) in pivots.iter().enumerate() {
                check().map_err(Failure::Interrupted)?;
                publish(DenseEvent::Set(ip,jp,&value,Some(i))).map_err(Failure::Interrupted)?;a[ip][jp]=value.clone();
            }
        }
        for &jnp in &no_pivots {
            for (ip,row) in a.iter_mut().enumerate().take(i) {
                check().map_err(Failure::Interrupted)?;
                if !row[jnp].is_zero() {
                    let value=quotient(row[jnp].mul(&pivot),divisor.as_ref())?;
                    publish(DenseEvent::Set(ip,jnp,&value,None)).map_err(Failure::Interrupted)?;row[jnp]=value;
                }
            }
        }
        for ip in 0..m {
            check().map_err(Failure::Interrupted)?;
            if ip==i{continue;}
            for kp in (j+1)..n {
                check().map_err(Failure::Interrupted)?;
                let value=quotient(pivot.mul(&a[ip][kp]).sub(&a[ip][j].mul(&a[i][kp])),divisor.as_ref())?;
                publish(DenseEvent::Set(ip,kp,&value,None)).map_err(Failure::Interrupted)?;a[ip][kp]=value;
            }
            let zero=C::zero();
            publish(DenseEvent::Set(ip,j,&zero,None)).map_err(Failure::Interrupted)?;a[ip][j]=zero;
        }
        pivots.push(j);i+=1;if i>=m{break;}
        divisor=if pivot.is_one(){None}else{Some(pivot)};
    }
    check().map_err(Failure::Interrupted)?;
    let denominator=pivots.first().map(|&j|a[0][j].clone()).unwrap_or_else(C::one);
    Ok(Some((denominator,pivots)))
}

#[derive(Clone,Debug)]
pub struct OrderedMap<V> { entries:Vec<Option<(usize,V)>>, positions:rustc_hash::FxHashMap<usize,usize> }
impl<V> Default for OrderedMap<V> {fn default()->Self{Self{entries:vec![],positions:Default::default()}}}
impl<V> OrderedMap<V> {
    pub fn get(&self,key:usize)->Option<&V>{self.positions.get(&key).map(|&i|&self.entries[i].as_ref().unwrap().1)}
    pub fn contains(&self,key:usize)->bool{self.positions.contains_key(&key)}
    pub fn insert(&mut self,key:usize,value:V){if let Some(&i)=self.positions.get(&key){self.entries[i]=Some((key,value));}else{self.positions.insert(key,self.entries.len());self.entries.push(Some((key,value)));}}
    pub fn remove(&mut self,key:usize)->Option<V>{let i=self.positions.remove(&key)?;Some(self.entries[i].take().unwrap().1)}
    pub fn iter(&self)->impl Iterator<Item=(usize,&V)>{self.entries.iter().filter_map(|entry|entry.as_ref().map(|(k,v)|(*k,v)))}
    pub fn keys(&self)->Vec<usize>{self.iter().map(|(k,_)|k).collect()}
    pub fn len(&self)->usize{self.positions.len()}
    pub fn is_empty(&self)->bool{self.positions.is_empty()}
}

#[derive(Clone,Copy)]
pub enum SetOrder { ViewIntersection, Difference, Intersection }
fn ordered<E,C:Coeff>(operation:SetOrder,left:&[usize],right:&[usize],
    order:&mut impl FnMut(SetOrder,&[usize],&[usize])->Result<Vec<usize>,E>)->Result<Option<Vec<usize>>,Failure<E,C>> {
    let l=left.iter().copied().collect::<std::collections::BTreeSet<_>>();
    let r=right.iter().copied().collect::<std::collections::BTreeSet<_>>();
    let expected=match operation{SetOrder::Difference=>l.difference(&r).copied().collect::<Vec<_>>(),_=>l.intersection(&r).copied().collect()};
    let result=order(operation,left,right).map_err(Failure::Interrupted)?;
    let mut check=result.clone();check.sort_unstable();
    // An ordering adapter cannot discard or invent a mathematical term.
    if check!=expected{return Ok(None);}
    Ok(Some(result))
}

/// Incremental sparse fraction-free Gauss-Jordan over an integral domain.
/// Input rows contain nonzero coefficients; no densification is performed.
/// The order callback serves representation compatibility only; every returned
/// key permutation is checked against the independently computed set operation.
pub fn sparse<C:Coeff,E>(mut rows:Vec<OrderedMap<C>>,check:&mut impl FnMut()->Result<(),E>,
    order:&mut impl FnMut(SetOrder,&[usize],&[usize])->Result<Vec<usize>,E>)->Result<Option<(Vec<OrderedMap<C>>,C,Vec<usize>)>,Failure<E,C>> {
    check().map_err(Failure::Interrupted)?;
    if rows.is_empty(){return Ok(Some((vec![],C::one(),vec![])));}
    if rows.iter().any(OrderedMap::is_empty){return Ok(None);}
    if rows.len()==1 {let row=rows.pop().unwrap();let j=*row.keys().iter().min().unwrap();let den=row.get(j).unwrap().clone();return Ok(Some((vec![row],den,vec![j])));}
    rows.sort_by_key(|row|row.keys().into_iter().min().unwrap());
    let mut reduced:OrderedMap<usize>=Default::default();let mut unreduced:OrderedMap<usize>=Default::default();
    let mut result:Vec<OrderedMap<C>>=Vec::new();let mut denom:Option<C>=None;let mut divisor=None;
    for source in rows {
        check().map_err(Failure::Interrupted)?;
        let mut a:OrderedMap<C>=Default::default();
        for (j,value) in source.iter(){if !reduced.contains(j){a.insert(j,value.clone());}}
        let mut cancel:OrderedMap<C>=Default::default();
        let Some(keys)=ordered(SetOrder::ViewIntersection,&unreduced.keys(),&a.keys(),order)? else{return Ok(None)};
        for j in keys {
            let coefficient=a.remove(j).unwrap();let previous=&result[*unreduced.get(j).unwrap()];
            for (k,value) in previous.iter() {
                check().map_err(Failure::Interrupted)?;
                let product=coefficient.mul(value);
                if let Some(old)=cancel.get(k){let value=old.add(&product);if !value.is_zero(){cancel.insert(k,value);}else{cancel.remove(k);}}
                else{cancel.insert(k,product);}
            }
        }
        let ak=a.keys();let ck=cancel.keys();let scale=denom.clone().unwrap_or_else(C::one);
        let Some(keys)=ordered(SetOrder::Difference,&ck,&ak,order)? else{return Ok(None)};
        for k in keys{check().map_err(Failure::Interrupted)?;a.insert(k,cancel.get(k).unwrap().neg());}
        let Some(keys)=ordered(SetOrder::Difference,&ak,&ck,order)? else{return Ok(None)};
        for k in keys{check().map_err(Failure::Interrupted)?;a.insert(k,a.get(k).unwrap().mul(&scale));}
        let Some(keys)=ordered(SetOrder::Intersection,&ck,&ak,order)? else{return Ok(None)};
        for k in keys{check().map_err(Failure::Interrupted)?;let value=a.get(k).unwrap().mul(&scale).sub(cancel.get(k).unwrap());if value.is_zero(){a.remove(k);}else{a.insert(k,value);}}
        if a.is_empty(){continue;}
        let j=a.keys().into_iter().min().unwrap();let pivot=a.remove(j).unwrap();
        let previous=unreduced.iter().map(|(p,&r)|(p,r)).collect::<Vec<_>>();
        for (pk,k) in previous {
            check().map_err(Failure::Interrupted)?;
            let old=&mut result[k];
            let Some(factor)=old.remove(j) else{
                for l in old.keys(){check().map_err(Failure::Interrupted)?;let value=quotient(old.get(l).unwrap().mul(&pivot),divisor.as_ref())?;old.insert(l,value);}continue;
            };
            let ak=a.keys();let ok=old.keys();
            let Some(keys)=ordered(SetOrder::Difference,&ak,&ok,order)? else{return Ok(None)};
            for l in keys{check().map_err(Failure::Interrupted)?;let value=quotient(factor.neg().mul(a.get(l).unwrap()),divisor.as_ref())?;old.insert(l,value);}
            let Some(keys)=ordered(SetOrder::Difference,&ok,&ak,order)? else{return Ok(None)};
            for l in keys{check().map_err(Failure::Interrupted)?;let value=quotient(pivot.mul(old.get(l).unwrap()),divisor.as_ref())?;old.insert(l,value);}
            let Some(keys)=ordered(SetOrder::Intersection,&ak,&ok,order)? else{return Ok(None)};
            for l in keys{check().map_err(Failure::Interrupted)?;let value=pivot.mul(old.get(l).unwrap()).sub(&factor.mul(a.get(l).unwrap()));if value.is_zero(){old.remove(l);}else{old.insert(l,quotient(value,divisor.as_ref())?);}}
            if old.is_empty(){unreduced.remove(pk);reduced.insert(pk,k);}
        }
        let i=result.len();let empty=a.is_empty();result.push(a);
        if empty{reduced.insert(j,i);}else{unreduced.insert(j,i);}
        if !pivot.is_one(){denom=Some(denom.map(|d|d.mul(&pivot)).unwrap_or(pivot));}
        if let Some(divisor)=&divisor{denom=Some(quotient(denom.unwrap(),Some(divisor))?);}
        divisor=denom.clone();
    }
    let denominator=denom.unwrap_or_else(C::one);
    let mut indexed=reduced.iter().chain(unreduced.iter()).map(|(p,&i)|(p,i)).collect::<Vec<_>>();indexed.sort_by_key(|&(p,_)|p);
    let mut out=Vec::with_capacity(indexed.len());let mut pivots=Vec::with_capacity(indexed.len());
    for (pivot,row) in indexed{let mut row=std::mem::take(&mut result[row]);row.insert(pivot,denominator.clone());out.push(row);pivots.push(pivot);}
    check().map_err(Failure::Interrupted)?;
    Ok(Some((out,denominator,pivots)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{Int,Rat};
    fn i(n:i64)->Int{Int::from(n)}
    fn row<C:Coeff>(entries:&[(usize,i64)])->OrderedMap<C>{let mut row=OrderedMap::default();for &(j,n) in entries{row.insert(j,C::from_i64(n));}row}
    fn sorted(operation:SetOrder,left:&[usize],right:&[usize])->Result<Vec<usize>,()> {
        let left=left.iter().copied().collect::<std::collections::BTreeSet<_>>();let right=right.iter().copied().collect::<std::collections::BTreeSet<_>>();
        Ok(match operation{SetOrder::Difference=>left.difference(&right).copied().collect(),_=>left.intersection(&right).copied().collect()})
    }
    #[test]
    fn dense_denominator_and_missing_pivot_columns_are_exact() {
        let mut a=vec![vec![i(0),i(2),i(4),i(6)],vec![i(0),i(1),i(3),i(5)],vec![i(0),i(0),i(0),i(0)]];
        let (den,pivots)=dense(&mut a,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap();
        assert_eq!(den,i(2));assert_eq!(pivots,vec![1,2]);
        assert_eq!(a,vec![vec![i(0),i(2),i(0),i(-2)],vec![i(0),i(0),i(2),i(4)],vec![i(0),i(0),i(0),i(0)]]);
        let mut a=vec![vec![i(1),i(2)],vec![i(3),i(4)]];
        assert_eq!(dense(&mut a,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap(),(i(-2),vec![0,1]));
    }
    #[test]
    fn sparse_preserves_denominator_rank_and_checked_order() {
        let rows=vec![row::<Int>(&[(0,1),(1,2)]),row::<Int>(&[(0,3),(1,4)])];
        let (out,den,pivots)=sparse(rows,&mut ||Ok::<_,()>(()),&mut sorted).unwrap().unwrap();
        assert_eq!(den,i(-2));assert_eq!(pivots,vec![0,1]);assert_eq!(out[0].get(0),Some(&i(-2)));assert_eq!(out[1].get(1),Some(&i(-2)));
        let rows=vec![row::<Rat>(&[(0,1),(9,2)]),row::<Rat>(&[(0,2),(9,4)])];
        let (_,den,pivots)=sparse(rows.clone(),&mut ||Ok::<_,()>(()),&mut sorted).unwrap().unwrap();
        assert_eq!(den,Rat::one());assert_eq!(pivots,vec![0]);
        assert!(sparse(rows,&mut ||Ok::<_,()>(()),&mut |_,_,_|Ok(vec![999])).unwrap().is_none());
    }
    #[test]
    fn interruptions_are_not_empty_rank_and_dense_pivot_tokens_share() {
        let rows=vec![row::<Int>(&[(0,1)]),row::<Int>(&[(1,1)])];
        assert!(matches!(sparse(rows,&mut ||Err("stop"),&mut |_,_,_|Ok(vec![])),Err(Failure::Interrupted("stop"))));
        let mut a=vec![vec![i(2),i(0),i(0)],vec![i(0),i(3),i(0)],vec![i(0),i(0),i(4)]];let mut tokens=vec![];
        dense(&mut a,&mut ||Ok::<_,()>(()),&mut |event|{if let DenseEvent::Set(r,c,_,Some(token))=event{tokens.push((r,c,token));}Ok(())}).unwrap().unwrap();
        assert_eq!(tokens,vec![(0,0,1),(0,0,2),(1,1,2)]);
    }
    #[test]
    fn dense_interruption_preserves_the_completed_cell_prefix() {
        let mut a=vec![vec![i(2),i(3)],vec![i(4),i(5)]];let mut writes=0;
        let result=dense(&mut a,&mut ||Ok::<_,&str>(()),&mut |event| {
            if matches!(event,DenseEvent::Set(..)){writes+=1;if writes==2{return Err("stop");}}
            Ok(())
        });
        assert!(matches!(result,Err(Failure::Interrupted("stop"))));
        assert_eq!(a,vec![vec![i(2),i(3)],vec![i(4),i(-2)]]);
    }

}

/// Exact denominator preparation for QQ matrices. Rows contain the existing
/// ordered nonzero entries; no rectangular or dense intermediate is allocated.
/// In global mode the returned denominator vector has one entry, even for an
/// empty matrix. Rowwise mode retains the input row order.
pub fn clear_denominators<E>(rows: &[Vec<super::Rat>], rowwise: bool,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<(Vec<super::Int>, Vec<Vec<super::Int>>), E> {
    use super::Int;
    check()?;
    let mut denominators = if rowwise { Vec::with_capacity(rows.len()) } else { vec![Int::one()] };
    for row in rows {
        check()?;
        let mut denominator = if rowwise { Int::one() } else { denominators[0].clone() };
        for (index, value) in row.iter().enumerate() {
            if index % 64 == 0 { check()?; }
            denominator = denominator.lcm(&value.d);
        }
        if rowwise { denominators.push(denominator); } else { denominators[0] = denominator; }
    }
    let mut numerators = Vec::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        check()?;
        let denominator = &denominators[if rowwise { i } else { 0 }];
        let mut values = Vec::with_capacity(row.len());
        for (index, value) in row.iter().enumerate() {
            if index % 64 == 0 { check()?; }
            values.push(value.n.mul(&denominator.div_exact_unchecked(&value.d)));
        }
        numerators.push(values);
    }
    check()?;
    Ok((denominators, numerators))
}

#[cfg(test)]
mod clear_denominators_tests {
    use super::*;
    use crate::poly::{Int, Rat};
    #[test]
    fn exact_row_and_global_preparation() {
        let rows = vec![vec![Rat::new(Int::from(-3), Int::from(14)), Rat::new(Int::from(5), Int::from(21))], vec![],
                        vec![Rat::new(Int::from(1), Int::from(25))]];
        let (d,n) = clear_denominators(&rows, true, &mut || Ok::<_, ()>(())).unwrap();
        assert_eq!(d, vec![Int::from(42), Int::one(), Int::from(25)]);
        assert_eq!(n, vec![vec![Int::from(-9),Int::from(10)],vec![],vec![Int::one()]]);
        let (d,n) = clear_denominators(&rows, false, &mut || Ok::<_, ()>(())).unwrap();
        assert_eq!(d, vec![Int::from(1050)]);
        assert_eq!(n, vec![vec![Int::from(-225),Int::from(250)],vec![],vec![Int::from(42)]]);
        let (d,n) = clear_denominators::<()>(&[], false, &mut || Ok(())).unwrap();
        assert_eq!(d, vec![Int::one()]); assert!(n.is_empty());
    }
    #[test]
    fn preparation_propagates_interruptions() {
        assert_eq!(clear_denominators::<u8>(&[], false, &mut || Err(7)), Err(7));
    }
}


/// Stock ZZ/QQ AUTO decisions, preserving the Python heuristic's binary64
/// operation order. Statistics count literal sparse rows/entries, including
/// explicit zeros and empty rows; dense callers omit zero rows/entries.
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum AutoMethod { GJ, FF, CD }
impl AutoMethod {
    pub fn name(self)->&'static str {match self {Self::GJ=>"GJ",Self::FF=>"FF",Self::CD=>"CD"}}
}
#[derive(Clone,Copy,Debug)]
pub struct AutoStats { pub rows:usize,pub entries:u128,pub columns:usize }
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum AutoPrefix { Done(AutoMethod), IntegerBits, RationalDenominators }

fn integer_ratio(n:num_bigint::BigInt,d:num_bigint::BigInt)->f64 {
    // CPython int/int rounds the exact rational once, rather than converting
    // each integer to f64 first. Reuse the certified destination rounding.
    crate::Number::rational(n,d).expect("positive AUTO ratio denominator").to_f64()
}
impl AutoStats {
    fn density(self)->f64 {
        if self.rows==0 {0.0}
        else {integer_ratio(self.entries.into(),self.rows.into())}
    }
}
pub fn auto_prefix(stats:AutoStats,rational:bool)->AutoPrefix {
    let density=stats.density();
    if rational {
        if density<(5.0_f64).min(integer_ratio(stats.columns.into(),2.into())) {
            AutoPrefix::Done(AutoMethod::GJ)
        }else {AutoPrefix::RationalDenominators}
    }else if stats.rows<10 {
        AutoPrefix::Done(if density<integer_ratio(stats.columns.into(),2.into()) {AutoMethod::GJ}else{AutoMethod::FF})
    }else if density<5.0 {AutoPrefix::Done(AutoMethod::GJ)}
    else if density>5.0+integer_ratio(10000.into(),stats.rows.into()) {AutoPrefix::Done(AutoMethod::FF)}
    else {AutoPrefix::IntegerBits}
}
pub fn auto_integer_bits(stats:AutoStats,bits:u64)->Option<AutoMethod> {
    if stats.rows==0||bits==0 {return None;} // Preserve malformed-storage fallback.
    let wideness=(1.0_f64).max(((2.0_f64/3.0)*stats.columns as f64)/stats.rows as f64);
    let divisor=num_bigint::BigInt::from(stats.rows)*num_bigint::BigInt::from(bits).pow(2);
    let max_density=(5.0+integer_ratio(10000.into(),divisor))*wideness;
    Some(if stats.density()<max_density {AutoMethod::GJ}else{AutoMethod::FF})
}
pub fn auto_rational_denominators<E>(values:&[super::Rat],check:&mut impl FnMut()->Result<(),E>)->Result<AutoMethod,E> {
    let mut maximum=None;
    for value in values {check()?;maximum=Some(maximum.unwrap_or(0).max(value.n.bits()));}
    let bits=maximum.unwrap_or(1);
    let mut lcm=super::Int::one();
    for value in values {
        check()?;lcm=lcm.lcm(&value.d);
        if u128::from(lcm.bits())>5*u128::from(bits) {return Ok(AutoMethod::GJ);}
    }
    check()?;
    Ok(if lcm.bits()<50 {AutoMethod::CD}else{AutoMethod::FF})
}

#[cfg(test)]
mod auto_tests {
    use super::*;
    use super::super::{Int,Rat};
    #[test]
    fn literal_sparse_rows_and_zero_width_choices() {
        assert_eq!(auto_prefix(AutoStats{rows:0,entries:0,columns:0},false),AutoPrefix::Done(AutoMethod::FF));
        assert_eq!(auto_prefix(AutoStats{rows:0,entries:0,columns:0},true),AutoPrefix::RationalDenominators);
        assert_eq!(auto_prefix(AutoStats{rows:0,entries:0,columns:1},true),AutoPrefix::Done(AutoMethod::GJ));
        assert_eq!(auto_prefix(AutoStats{rows:10,entries:50,columns:5},false),AutoPrefix::IntegerBits);
        assert_eq!(auto_integer_bits(AutoStats{rows:10,entries:50,columns:5},0),None);
    }
    #[test]
    fn rational_lcm_thresholds_and_cancellation() {
        let mut check=||Ok::<(),()>(());
        assert_eq!(auto_rational_denominators(&[],&mut check),Ok(AutoMethod::CD));
        assert_eq!(auto_rational_denominators(&[Rat::zero()],&mut check),Ok(AutoMethod::GJ));
        assert_eq!(auto_rational_denominators(&[Rat::new(Int::from(1),Int::from(3))],&mut check),Ok(AutoMethod::CD));
        let big=Int::from(num_bigint::BigInt::from(1)<<50usize); // 51 bits, numerator large enough for FF.
        assert_eq!(auto_rational_denominators(&[Rat::new(big.add(&Int::one()),big)],&mut check),Ok(AutoMethod::FF));
        assert_eq!(auto_rational_denominators(&[Rat::one()],&mut ||Err(17)),Err(17));
    }
    #[test]
    fn rational_rounding_is_not_rounded_integer_operands() {
        let n=num_bigint::BigInt::from(9007199254740993u64);
        let d=num_bigint::BigInt::from(9007199254740994u64);
        assert_eq!(integer_ratio(n,d).to_bits(),1.0_f64.to_bits()-1);
    }
}


#[cfg(test)]
mod auto_rounding_fixtures {
    use super::*;
    #[test]
    fn python_binary64_statistics_and_threshold_boundaries() {
        // Independently generated with CPython int/int, then its exact
        // binary64 add/mul/div sequence; includes integers beyond 2^53.
        let cases=[
            (10u64,48u128,9u64,1u64,"GJ",4617090337980232499u64),
            (10u64,57u128,9u64,3u64,"GJ",4618103647896390861u64),
            (10u64,53u128,20u64,17u64,"GJ",4617653287933653811u64),
            (10u64,52u128,9u64,139u64,"FF",4617540697942969549u64),
            (10u64,53u128,67u64,211u64,"GJ",4617653287933653811u64),
            (13u64,66u128,18446744073709551615u64,1u64,"GJ",4617402125646742764u64),
            (13u64,64u128,18446744073709551615u64,3u64,"GJ",4617228910276459284u64),
            (13u64,65u128,18446744073709551615u64,17u64,"GJ",4617315517961601024u64),
            (13u64,71u128,9u64,139u64,"FF",4617835164072451466u64),
            (13u64,64u128,18446744073709551615u64,211u64,"GJ",4617228910276459284u64),
            (29u64,149u128,18446744073709551615u64,1u64,"GJ",4617470814500475869u64),
            (29u64,151u128,20u64,3u64,"GJ",4617548462769913291u64),
            (29u64,144u128,67u64,17u64,"GJ",4617276693826882313u64),
            (29u64,144u128,67u64,139u64,"GJ",4617276693826882313u64),
            (29u64,143u128,67u64,211u64,"GJ",4617237869692163602u64),
            (997u64,4987u128,67u64,1u64,"GJ",4617317776537141330u64),
            (997u64,4987u128,18446744073709551615u64,3u64,"GJ",4617317776537141330u64),
            (997u64,4987u128,9007199254740993u64,17u64,"GJ",4617317776537141330u64),
            (997u64,4982u128,67u64,139u64,"GJ",4617312130098290565u64),
            (997u64,4983u128,9u64,211u64,"GJ",4617313259386060718u64),
            (9007199254740993u64,45035996273704963u128,9007199254740993u64,1u64,"GJ",4617315517961601024u64),
            (9007199254740993u64,45035996273704968u128,20u64,3u64,"GJ",4617315517961601024u64),
            (9007199254740993u64,45035996273704965u128,9007199254740993u64,17u64,"GJ",4617315517961601024u64),
            (9007199254740993u64,45035996273704963u128,20u64,139u64,"FF",4617315517961601024u64),
            (9007199254740993u64,45035996273704969u128,9u64,211u64,"FF",4617315517961601024u64),
            (9007199254740994u64,45035996273704975u128,9007199254740993u64,1u64,"GJ",4617315517961601025u64),
            (9007199254740994u64,45035996273704977u128,9007199254740993u64,3u64,"GJ",4617315517961601025u64),
            (9007199254740994u64,45035996273704969u128,67u64,17u64,"GJ",4617315517961601024u64),
            (9007199254740994u64,45035996273704970u128,9007199254740993u64,139u64,"FF",4617315517961601024u64),
            (9007199254740994u64,45035996273704970u128,67u64,211u64,"FF",4617315517961601024u64),
            (18446744073709551615u64,92233720368547758076u128,9007199254740993u64,1u64,"GJ",4617315517961601024u64),
            (18446744073709551615u64,92233720368547758074u128,18446744073709551615u64,3u64,"FF",4617315517961601024u64),
            (18446744073709551615u64,92233720368547758073u128,67u64,17u64,"FF",4617315517961601024u64),
            (18446744073709551615u64,92233720368547758075u128,67u64,139u64,"FF",4617315517961601024u64),
            (18446744073709551615u64,92233720368547758073u128,20u64,211u64,"FF",4617315517961601024u64),
        ];
        for (rows,entries,columns,bits,expected,density_bits) in cases {
            let (Ok(rows),Ok(columns))=(usize::try_from(rows),usize::try_from(columns)) else{continue};
            let stats=AutoStats{rows,entries,columns};
            assert_eq!(stats.density().to_bits(),density_bits);
            let result=match auto_prefix(stats,false) {
                AutoPrefix::Done(method)=>method,
                AutoPrefix::IntegerBits=>auto_integer_bits(stats,bits).unwrap(),
                AutoPrefix::RationalDenominators=>unreachable!(),
            };
            assert_eq!(result.name(),expected,"{stats:?} bits={bits}");
        }
    }
}
