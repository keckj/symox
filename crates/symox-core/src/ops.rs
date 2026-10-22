//! Canonical constructors: `add`, `mul`, `pow` and friends.
//!
//! These implement SymPy's automatic simplification rules (flattening,
//! collecting like terms and powers, numeric folding, exact radicals,
//! distribution of a numeric coefficient over a sum, ...). Every expression
//! created through them is in canonical form, which is what makes
//! structural equality meaningful.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, ToPrimitive};

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::{exact_root, trial_factor, Number};

mod memo;

// ----------------------------------------------------------------------
// helpers on the internal representation
// ----------------------------------------------------------------------

/// `k * t` where `t` is a term of an `Add` (no numeric coefficient).
pub fn add_term(t: &Expr, k: &Number) -> Expr {
    if k.is_one() {
        return t.clone();
    }
    if let Kind::Mul(c,factors)=t.kind(){
        if !t.is_canonical(){
            let mut tail=Vec::with_capacity(factors.len()+1);
            if !c.is_one(){tail.push((Expr::num(c.clone()),Expr::one()));}
            tail.extend(factors.iter().cloned());
            return Expr::from_unevaluated_kind(Kind::Mul(k.clone(),tail));
        }
    }
    match t.kind() {
        Kind::Mul(c, facs) if c.is_one() => Expr::from_kind(Kind::Mul(k.clone(), facs.clone())),
        Kind::Pow(b, e) => Expr::from_kind(Kind::Mul(k.clone(), vec![(b.clone(), e.clone())])),
        _ => Expr::from_kind(Kind::Mul(k.clone(), vec![(t.clone(), Expr::one())])),
    }
}

/// `b ** e` where `(b, e)` is a factor of a `Mul`.
pub fn mul_factor(b: &Expr, e: &Expr) -> Expr {
    if e.is_one() {
        b.clone()
    } else {
        Expr::from_kind(Kind::Pow(b.clone(), e.clone()))
    }
}

/// Split an expression into `(numeric coefficient, rest)`.
pub fn as_coeff_mul(e: &Expr) -> (Number, Expr) {
    match e.kind() {
        Kind::Num(n) => (n.clone(), Expr::one()),
        Kind::Mul(c, facs) if !c.is_one() => (c.clone(), mul_from_factors_like(e,facs.clone())),
        _ => (Number::ONE, e.clone()),
    }
}

/// Split an expression into `(numeric constant term, rest)`.
pub fn as_coeff_add(e: &Expr) -> (Number, Expr) {
    match e.kind() {
        Kind::Num(n) => (n.clone(), Expr::zero()),
        Kind::Add(c, terms) if !c.is_zero() => {
            (c.clone(), Expr::from_kind(Kind::Add(Number::ZERO, terms.clone())).normalize_add())
        }
        _ => (Number::ZERO, e.clone()),
    }
}

/// `(base, exponent)` view of an expression.
pub fn as_base_exp(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        _ => (e.clone(), Expr::one()),
    }
}

/// Build a Mul with coefficient 1 from already canonical sorted factors.
fn mul_from_factors(facs: Vec<(Expr, Expr)>) -> Expr {
    if facs.len() == 1 {
        let (b, e) = &facs[0];
        return mul_factor(b, e);
    }
    Expr::from_kind(Kind::Mul(Number::ONE, facs))
}

/// Extracting a coefficient must not certify a raw product tail as canonical.
fn mul_from_factors_like(parent:&Expr,factors:Vec<(Expr,Expr)>)->Expr{
    if parent.is_canonical(){return mul_from_factors(factors);}
    match factors.len(){
        0=>Expr::one(),
        1=>mul_factor(&factors[0].0,&factors[0].1),
        _=>Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,factors)),
    }
}

impl Expr {
    /// Fix an `Add` node which may have become degenerate (0 or 1 term).
    fn normalize_add(self) -> Expr {
        match self.kind() {
            Kind::Add(c, terms) if terms.len() == 1 && c.is_zero() => add_term(&terms[0].0, &terms[0].1),
            Kind::Add(c, terms) if terms.is_empty() => Expr::num(c.clone()),
            _ => self,
        }
    }
}

// ----------------------------------------------------------------------
// Add
// ----------------------------------------------------------------------

#[derive(Default)]
struct InfFlags {
    pos: bool,
    neg: bool,
    complex: bool,
}

/// Canonical sum of `args`.
pub fn add(args: &[Expr]) -> Expr {
    add_with_cache(args,&mut rustc_hash::FxHashMap::default())
}

type ArithmeticMemo=rustc_hash::FxHashMap<usize,(Expr,Expr)>;

/// Flatten held sums once per physical node, reclassifying their literal
/// arguments rather than treating a raw payload as a canonical term table.
fn normalized_sum(expression:&Expr,cache:&mut ArithmeticMemo)->Expr{
    if let Some((_,value))=cache.get(&expression.dag_id()){return value.clone();}
    let values=expression.args();
    let value=if values.len()==1 && values[0].is_add() && values[0].has_held_arithmetic(){normalized_sum(&values[0],cache)}else{add_with_cache(&values,cache)};
    cache.insert(expression.dag_id(),(expression.clone(),value.clone()));value
}

/// Preserve exact arithmetic when a compact float's exponent cannot fit the
/// native counter. Ordinary expressions use only a cached feature-bit check.
fn wide_arithmetic(args:&[Expr],product:bool)->Option<Expr>{
    if !args.iter().any(Expr::has_wide_exponent){return None;}
    if args.len()==1{return Some(args[0].clone());}
    let mut coefficient=if product{Number::ONE}else{Number::ZERO};
    let complete=args.iter().all(|e|if let Some(number)=e.as_number(){
        let next=if product{coefficient.checked_mul(number)}else{coefficient.checked_add(number)};
        if let Some(next)=next{coefficient=next;true}else{false}
    }else{false});
    if complete{return Some(Expr::num(coefficient));}
    // Literal complex values also have a compact, checked coefficient path.
    // In particular a huge finite coefficient times I must remain usable as a
    // numerical-field value, rather than needing general symbolic expansion.
    let mut real=if product{Number::ONE}else{Number::ZERO};let mut imaginary=Number::ZERO;
    let literal=args.iter().all(|argument|{
        let Some((re,im))=crate::complex_parts::literal_parts(argument)else{return false;};
        let (re,im)=(re.as_number().unwrap(),im.as_number().unwrap());
        let next=if product{(||{
            let r=real.checked_mul(re)?.checked_add(&imaginary.checked_mul(im)?.neg())?;
            let i=real.checked_mul(im)?.checked_add(&imaginary.checked_mul(re)?)?;
            Some((r,i))
        })()}else{real.checked_add(re).zip(imaginary.checked_add(im))};
        if let Some((r,i))=next{real=r;imaginary=i;true}else{false}
    });
    if literal{
        return Some(if imaginary.is_zero(){Expr::num(real)}else if real.is_zero(){
            if imaginary.is_one(){Expr::i()}else{Expr::from_kind(Kind::Mul(imaginary,vec![(Expr::i(),Expr::one())]))}
        }else{Expr::from_kind(Kind::Add(real,vec![(Expr::i(),imaginary)]))});
    }
    Some(Expr::from_unevaluated_kind(if product{
        Kind::Mul(Number::ONE,args.iter().cloned().map(|e|(e,Expr::one())).collect())
    }else{
        Kind::Add(Number::ZERO,args.iter().cloned().map(|e|(e,Number::ONE)).collect())
    }))
}

fn add_with_cache(args:&[Expr],cache:&mut ArithmeticMemo)->Expr{
    if let Some(result)=wide_arithmetic(args,false){return result;}
    match args.len() {
        0 => return Expr::zero(),
        1 => return args[0].clone(),
        _ => {}
    }
    if args.iter().any(Expr::may_have_order) {
        if let Some(r) = crate::series::order::add_hook(args) {
            return r;
        }
    }
    let mut coeff = Number::ZERO;
    let mut terms: Vec<(Expr, Number)> = Vec::with_capacity(args.len() + 4);
    let mut inf = InfFlags::default();
    for argument in args {
        let normalized=if argument.is_add() && argument.has_held_arithmetic(){Some(normalized_sum(argument,cache))}else{None};
        let a=normalized.as_ref().unwrap_or(argument);
        match a.kind() {
            Kind::Num(n) => coeff = coeff.add(n),
            Kind::Add(c, ts) => {
                coeff = coeff.add(c);
                if a.has_nonfinite(){
                    for(t,k)in ts{
                        if k.is_nan(){return Expr::nan();}
                        if k.is_zero(){continue;}
                        match t.kind(){
                            Kind::Const(Constant::NaN)=>return Expr::nan(),
                            Kind::Const(Constant::ComplexInfinity)=>inf.complex=true,
                            Kind::Const(Constant::Infinity)=>if k.is_negative(){inf.neg=true}else{inf.pos=true},
                            Kind::Const(Constant::NegInfinity)=>if k.is_negative(){inf.pos=true}else{inf.neg=true},
                            _=>terms.push((t.clone(),k.clone())),
                        }
                    }
                }else{terms.extend(ts.iter().cloned());}
            }
            Kind::Mul(c, facs) if !c.is_one() => {
                // a term `c * oo` never exists: mul normalises it
                terms.push((mul_from_factors_like(a,facs.clone()), c.clone()));
            }
            Kind::Const(Constant::NaN) => return Expr::nan(),
            Kind::Const(Constant::Infinity) => inf.pos = true,
            Kind::Const(Constant::NegInfinity) => inf.neg = true,
            Kind::Const(Constant::ComplexInfinity) => inf.complex = true,
            _ => terms.push((a.clone(), Number::ONE)),
        }
    }
    finish_add(coeff, terms, inf)
}

/// Fast binary addition.
pub fn add2(a: &Expr, b: &Expr) -> Expr {
    if a.has_wide_exponent()||b.has_wide_exponent(){return wide_arithmetic(&[a.clone(),b.clone()],false).unwrap();}
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    if let (Kind::Num(x), Kind::Num(y)) = (a.kind(), b.kind()) {
        return Expr::num(x.add(y));
    }
    memo::apply(memo::Operation::Add,a,b,||add2_nontrivial(a,b))
}

fn add2_nontrivial(a:&Expr,b:&Expr)->Expr{
    // fast path: merge already sorted term lists in linear time
    if let Some(r) = add2_merge(a, b) {
        return r;
    }
    add(&[a.clone(), b.clone()])
}

/// View `e` as a sorted list of `(term, coeff)` plus a constant, when it is
/// a plain finite term (no infinities / nan involved).
fn as_terms(e: &Expr) -> Option<(Number, std::borrow::Cow<'_, [(Expr, Number)]>)> {
    use std::borrow::Cow;
    match e.kind() {
        Kind::Add(..) if !e.is_canonical() => None,
        Kind::Add(c, terms) => Some((c.clone(), Cow::Borrowed(terms.as_slice()))),
        Kind::Num(n) => Some((n.clone(), Cow::Owned(Vec::new()))),
        Kind::Const(_) => {
            if matches!(e.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN) {
                None
            } else {
                Some((Number::ZERO, Cow::Owned(vec![(e.clone(), Number::ONE)])))
            }
        }
        Kind::Mul(c, facs) if !c.is_one() => {
            Some((Number::ZERO, Cow::Owned(vec![(mul_from_factors_like(e,facs.clone()), c.clone())])))
        }
        _ => Some((Number::ZERO, Cow::Owned(vec![(e.clone(), Number::ONE)]))),
    }
}

fn add2_merge(a: &Expr, b: &Expr) -> Option<Expr> {
    if !(a.is_add() || b.is_add()) || a.may_have_order() || b.may_have_order() || a.has_nonfinite() || b.has_nonfinite() {
        return None;
    }
    let (ca, ta) = as_terms(a)?;
    let (cb, tb) = as_terms(b)?;
    let coeff = ca.add(&cb);
    if coeff.is_nan() {
        return None;
    }
    let mut out: Vec<(Expr, Number)> = Vec::with_capacity(ta.len() + tb.len());
    let (mut i, mut j) = (0, 0);
    while i < ta.len() && j < tb.len() {
        match ta[i].0.cmp(&tb[j].0) {
            std::cmp::Ordering::Less => {
                out.push(ta[i].clone());
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(tb[j].clone());
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let k = ta[i].1.add(&tb[j].1);
                if !k.is_zero() {
                    out.push((ta[i].0.clone(), k));
                }
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&ta[i..]);
    out.extend_from_slice(&tb[j..]);
    if out.is_empty() {
        return Some(Expr::num(coeff));
    }
    if out.len() == 1 && coeff.is_zero() {
        return Some(add_term(&out[0].0, &out[0].1));
    }
    Some(Expr::from_kind(Kind::Add(coeff, out)))
}

fn finish_add(coeff: Number, mut terms: Vec<(Expr, Number)>, inf: InfFlags) -> Expr {
    if coeff.is_nan() {
        return Expr::nan();
    }
    // sort + merge equal terms
    // stable sort: linear on the common case of concatenated sorted runs
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out: Vec<(Expr, Number)> = Vec::with_capacity(terms.len());
    for (t, k) in terms {
        match out.last_mut() {
            Some((lt, lk)) if *lt == t => *lk = lk.add(&k),
            _ => out.push((t, k)),
        }
    }
    out.retain(|(_, k)| !k.is_zero());
    let mut coeff = coeff;
    if inf.pos || inf.neg || inf.complex {
        if (inf.pos && inf.neg) || (inf.complex && (inf.pos || inf.neg)) {
            return Expr::nan();
        }
        crate::assumptions::exact_queries(||out.retain(|(term,_)|{
            if inf.complex{
                !(crate::assumptions::is_finite(term)==Some(true) && crate::assumptions::is_real(term).is_some())
            }else{
                crate::assumptions::is_real(term)!=Some(true)
                    && if inf.pos{crate::assumptions::is_nonnegative(term)!=Some(true)}else{crate::assumptions::is_nonpositive(term)!=Some(true)}
            }
        }));
        coeff = Number::ZERO;
        let c = if inf.complex {
            Constant::ComplexInfinity
        } else if inf.pos {
            Constant::Infinity
        } else {
            Constant::NegInfinity
        };
        if out.is_empty() {
            return Expr::constant(c);
        }
        out.push((Expr::constant(c), Number::ONE));
        out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    }
    if out.is_empty() {
        return Expr::num(coeff);
    }
    if out.len() == 1 && coeff.is_zero() {
        let (t, k) = &out[0];
        return add_term(t, k);
    }
    Expr::from_kind(Kind::Add(coeff, out))
}

/// Build a sum directly from a map `term -> coefficient` (terms must be
/// canonical "coefficient-free" terms, distinct). Used by fast paths.
pub fn add_from_terms(coeff: Number, terms: Vec<(Expr, Number)>) -> Expr {
    finish_add(coeff, terms, InfFlags::default())
}

// ----------------------------------------------------------------------
// Mul
// ----------------------------------------------------------------------

/// A compact active-product summary. Repeated DAG edges increase a factor's
/// multiplicity instead of materializing repeated occurrences. Finalizing a
/// child product here would distribute coefficients over Add too early.
#[derive(Clone)]
struct ProductParts {
    coefficient: Number,
    factors: Vec<(Expr, Expr, BigInt)>,
    positions: rustc_hash::FxHashMap<(usize, usize), usize>,
    infinite: bool,
    negative_infinity: bool,
    complex_infinity: bool,
    zero_times_infinite_sum: bool,
}
impl Default for ProductParts {
    fn default() -> Self {
        Self { coefficient: Number::ONE, factors: Vec::new(), positions: rustc_hash::FxHashMap::default(),
            infinite: false, negative_infinity: false, complex_infinity: false,
            zero_times_infinite_sum: false }
    }
}
impl ProductParts {
    fn factor(&mut self, base: &Expr, exponent: &Expr, count: &BigInt) {
        let key = (base.dag_id(), if exponent.is_one() { 0 } else { exponent.dag_id() });
        if let Some(&index) = self.positions.get(&key) { self.factors[index].2 += count; }
        else {
            self.positions.insert(key, self.factors.len());
            self.factors.push((base.clone(), exponent.clone(), count.clone()));
        }
    }
    fn merge(&mut self, other: &Self) {
        self.coefficient = self.coefficient.mul(&other.coefficient);
        self.infinite |= other.infinite;
        self.negative_infinity ^= other.negative_infinity;
        self.complex_infinity |= other.complex_infinity;
        self.zero_times_infinite_sum |= other.zero_times_infinite_sum;
        for (base, exponent, count) in &other.factors { self.factor(base, exponent, count); }
    }
    fn infinity(&mut self, constant: Constant) {
        self.infinite = true;
        self.negative_infinity ^= constant == Constant::NegInfinity;
        self.complex_infinity |= constant == Constant::ComplexInfinity;
    }
    fn absorb(&mut self, argument: &Expr, cache: &mut ProductMemo) {
        match argument.kind() {
            Kind::Num(number) => self.coefficient = self.coefficient.mul(number),
            Kind::Mul(..) => self.merge(&product_parts(argument, cache)),
            Kind::Pow(base, exponent) => {
                if let Some(value) = collection_scalar_power(base, exponent) {
                    self.absorb(&value, cache);
                } else { self.factor(base, exponent, &BigInt::one()); }
            }
            Kind::Const(Constant::NaN) => self.coefficient = Number::float(f64::NAN),
            Kind::Const(constant) if constant.is_infinite() => self.infinity(*constant),
            _ => {
                self.zero_times_infinite_sum |= zero_add_infinity(argument);
                self.factor(argument, &Expr::one(), &BigInt::one());
            }
        }
    }
}
#[derive(Default)]
struct ProductMemo {
    // Owners prevent pointer reuse when public arguments synthesize temporaries.
    parts: rustc_hash::FxHashMap<usize, (Expr, ProductParts)>,
}
fn product_parts(expression: &Expr, cache: &mut ProductMemo) -> ProductParts {
    if let Some((_, value)) = cache.parts.get(&expression.dag_id()) { return value.clone(); }
    let Kind::Mul(coefficient, factors) = expression.kind() else { unreachable!() };
    let mut parts = ProductParts { coefficient: coefficient.clone(), ..ProductParts::default() };
    for (base, exponent) in factors {
        if exponent.is_one() { parts.absorb(base, cache); }
        else { parts.factor(base, exponent, &BigInt::one()); }
    }
    cache.parts.insert(expression.dag_id(), (expression.clone(), parts.clone()));
    parts
}

/// Literal infinity in the public arguments of an Add/Mul. This intentionally
/// stops at functions, powers and nested products, as SymPy's Mul.flatten does.
fn direct_arithmetic_infinity(e: &Expr) -> bool {
    if !(e.is_add() || e.is_mul()) || !e.has_nonfinite() { return false; }
    e.args().iter().any(|argument| matches!(argument.kind(), Kind::Const(c) if c.is_infinite()))
}
fn zero_add_infinity(e: &Expr) -> bool {
    if !e.is_add() || !e.has_nonfinite() { return false; }
    e.args().iter().any(|term| {
        matches!(term.kind(), Kind::Const(c) if c.is_infinite())
            || (term.is_mul() && direct_arithmetic_infinity(term))
    })
}
fn needs_product_scan(e: &Expr) -> bool {
    e.is_mul() && (e.has_held_arithmetic() || e.has_nonfinite()
        || matches!(e.kind(), Kind::Mul(Number::Float(value), _) if !value.is_finite()))
}

/// Canonical product of `args`. Shared held products are normalized only along
/// their active multiplicative spine; powers and functions remain boundaries.
pub fn mul(args: &[Expr]) -> Expr {
    mul_with_cache(args, 0, &mut ProductMemo::default())
}
fn mul_with_cache(args: &[Expr], depth: u32, cache: &mut ProductMemo) -> Expr {
    if args.iter().any(Expr::has_noncommutative_host) {
        if let Some(result) = mul_noncommutative(args, depth, cache) { return result; }
    }
    mul_commutative(args, depth, cache)
}

/// SymPy's `Mul.flatten` over factors of host classes that may not commute:
/// a factor whose `is_commutative` is not True keeps its position after the
/// commutative part, and only adjacent powers of one base combine, unless
/// their exponents add to a sum. `None` when every factor commutes.
fn mul_noncommutative(args: &[Expr], depth: u32, cache: &mut ProductMemo) -> Option<Expr> {
    use crate::decompose::commutativity;
    if args.iter().all(|argument| !argument.has_noncommutative_host() || commutativity(argument) == Some(true)) {
        return None;
    }
    fn push(factor: Expr, commutative: &mut Vec<Expr>, ordered: &mut Vec<Expr>) {
        let mut pending = std::collections::VecDeque::from([factor]);
        while let Some(factor) = pending.pop_front() {
            let Some(last) = ordered.pop() else { ordered.push(factor); continue };
            let ((b1, e1), (b2, e2)) = (crate::decompose::as_base_exp(&last), crate::decompose::as_base_exp(&factor));
            let exponent = add2(&e1, &e2);
            if b1 == b2 && !exponent.is_add() {
                let combined = pow(&b1, &exponent);
                if crate::decompose::commutativity(&combined) == Some(true) { commutative.push(combined); }
                else { pending.push_front(combined); }
            } else { ordered.push(last); ordered.push(factor); }
        }
    }
    let (mut commutative, mut ordered) = (Vec::with_capacity(args.len()), Vec::new());
    for argument in args {
        if commutativity(argument) == Some(true) { commutative.push(argument.clone()); }
        else if argument.is_mul() {
            for factor in argument.args() {
                if commutativity(&factor) == Some(true) { commutative.push(factor); }
                else { push(factor, &mut commutative, &mut ordered); }
            }
        } else { push(argument.clone(), &mut commutative, &mut ordered); }
    }
    let scalar = mul_commutative(&commutative, depth, cache);
    if ordered.is_empty() || scalar.is_zero() || scalar.is_constant(Constant::NaN)
        || scalar.as_number().is_some_and(Number::is_nan) { return Some(scalar); }
    let (coefficient, mut factors) = match scalar.kind() {
        Kind::Num(number) => (number.clone(), Vec::new()),
        Kind::Mul(coefficient, factors) => (coefficient.clone(), factors.clone()),
        _ => (Number::ONE, vec![as_base_exp(&scalar)]),
    };
    if coefficient.is_one() && factors.is_empty() && ordered.len() == 1 { return ordered.pop(); }
    factors.extend(ordered.iter().map(as_base_exp));
    Some(Expr::from_kind(Kind::Mul(coefficient, factors)))
}

fn mul_commutative(args: &[Expr], depth: u32, cache: &mut ProductMemo) -> Expr {
    if let Some(result)=wide_arithmetic(args,true){return result;}
    // Removing exact identity factors must not evaluate the surviving operand.
    let mut nontrivial = args.iter().filter(|argument| !argument.is_one());
    let Some(first) = nontrivial.next() else { return Expr::one(); };
    if nontrivial.next().is_none() { return first.clone(); }
    if args.iter().any(Expr::may_have_order) {
        if let Some(result) = crate::series::order::mul_hook(args) { return result; }
    }
    let mut coefficient = Number::ONE;
    let mut factors = Vec::with_capacity(args.len() + 4);
    let mut infinities = Vec::new();
    let mut zero_infinite_sum = false;
    if args.iter().any(needs_product_scan) {
        let mut parts = ProductParts::default();
        for argument in args { parts.absorb(argument, cache); }
        coefficient = parts.coefficient;
        zero_infinite_sum = parts.zero_times_infinite_sum;
        if parts.infinite {
            infinities.push(if parts.complex_infinity { Constant::ComplexInfinity }
                else if parts.negative_infinity { Constant::NegInfinity } else { Constant::Infinity });
        }
        for (base, exponent, count) in parts.factors {
            let exponent = if count.is_one() { exponent }
                else { mul2(&Expr::num(Number::from_bigint(count)), &exponent) };
            factors.push((base, exponent));
        }
    } else {
        // Keep the usual canonical path allocation-light.
        for argument in args {
            match argument.kind() {
                Kind::Num(number) => coefficient = coefficient.mul(number),
                Kind::Mul(c, fs) => { coefficient = coefficient.mul(c); factors.extend(fs.iter().cloned()); }
                Kind::Pow(base, exponent) => {
                    if let Some(value) = collection_scalar_power(base, exponent) {
                        match value.kind() {
                            Kind::Num(number) => coefficient = coefficient.mul(number),
                            Kind::Const(Constant::NaN) => return Expr::nan(),
                            Kind::Const(constant) if constant.is_infinite() => infinities.push(*constant),
                            _ => factors.push((value, Expr::one())),
                        }
                    } else { factors.push((base.clone(), exponent.clone())); }
                }
                Kind::Const(Constant::NaN) => return Expr::nan(),
                Kind::Const(constant) if constant.is_infinite() => infinities.push(*constant),
                _ => {
                    zero_infinite_sum |= zero_add_infinity(argument);
                    factors.push((argument.clone(), Expr::one()));
                }
            }
        }
    }
    if coefficient.is_nan() { return Expr::nan(); }
    if zero_infinite_sum && (coefficient.is_zero() || factors.iter().any(|(base, exponent)| {
        exceptional_scalar_power(base, exponent).is_some_and(|value| value.is_zero())
    })) { return Expr::nan(); }
    if !infinities.is_empty() { return mul_infinite(coefficient, factors, &infinities, cache); }
    finish_mul(coefficient, factors, depth, cache)
}

/// Fast binary multiplication.
pub fn mul2(a: &Expr, b: &Expr) -> Expr {
    if a.has_wide_exponent()||b.has_wide_exponent(){return wide_arithmetic(&[a.clone(),b.clone()],true).unwrap();}
    if a.is_one() {
        return b.clone();
    }
    if b.is_one() {
        return a.clone();
    }
    match (a.kind(), b.kind()) {
        (Kind::Num(x), Kind::Num(y)) => return Expr::num(x.mul(y)),
        (Kind::Num(x), _) if x.is_zero() && !matches!(b.kind(), Kind::Const(_)) && !b.has_nonfinite() && !b.has_held_arithmetic() => {
            return Expr::zero()
        }
        (_, Kind::Num(y)) if y.is_zero() && !matches!(a.kind(), Kind::Const(_)) && !a.has_nonfinite() && !a.has_held_arithmetic() => {
            return Expr::zero()
        }
        _ => {}
    }
    if a.may_have_order() || b.may_have_order() || a.has_noncommutative_host() || b.has_noncommutative_host() {
        return mul(&[a.clone(), b.clone()]);
    }
    memo::apply(memo::Operation::Multiply,a,b,||mul2_nontrivial(a,b))
}

fn mul2_nontrivial(a:&Expr,b:&Expr)->Expr{
    if let Some(r) = mul2_fast(a, b).or_else(|| mul2_fast(b, a)) {
        return r;
    }
    mul(&[a.clone(), b.clone()])
}

/// A factor that can never merge with a factor of a different base.
#[inline]
fn plain_factor(e: &Expr) -> Option<(Expr, Expr)> {
    let plain_base = |b: &Expr| match b.kind() {
        Kind::Sym(_) | Kind::Add(..) => true,
        Kind::Func(f, _) => !matches!(f, Func::Exp | Func::ExpPolar),
        Kind::Const(c) => matches!(c, Constant::Pi | Constant::EulerGamma | Constant::I | Constant::Catalan | Constant::GoldenRatio | Constant::TribonacciConstant),
        _ => false,
    };
    match e.kind() {
        Kind::Pow(b, x) if plain_base(b) && !matches!(b.kind(), Kind::Const(_)) => {
            Some((b.clone(), x.clone()))
        }
        Kind::Add(..) => None,
        _ if plain_base(e) => Some((e.clone(), Expr::one())),
        _ => None,
    }
}

/// Allocation-light products for the most common shapes.
fn mul2_fast(a: &Expr, b: &Expr) -> Option<Expr> {
    if needs_product_scan(a) || needs_product_scan(b) { return None; }
    match (a.kind(), b.kind()) {
        // number * plain factor
        (Kind::Num(n), _) if !n.is_nan() && !n.is_zero() => {
            if let Kind::Mul(c, facs) = b.kind() {
                let k = n.mul(c);
                if k.is_zero() {
                    return None;
                }
                if k.is_one() {
                    return Some(mul_from_factors(facs.clone()));
                }
                return Some(Expr::from_kind(Kind::Mul(k, facs.clone())));
            }
            let f = plain_factor(b)?;
            Some(Expr::from_kind(Kind::Mul(n.clone(), vec![f])))
        }
        // plain * plain with different bases
        (_, _) if !a.is_mul() && !b.is_mul() => {
            let fa = plain_factor(a)?;
            let fb = plain_factor(b)?;
            if fa.0 == fb.0 {
                return None;
            }
            let v = if cmp_factor(&fa, &fb) == std::cmp::Ordering::Less { vec![fa, fb] } else { vec![fb, fa] };
            Some(Expr::from_kind(Kind::Mul(Number::ONE, v)))
        }
        // product * plain factor whose base is new
        (Kind::Mul(c, facs), _) if !b.is_mul() => {
            let f = plain_factor(b)?;
            let pos = facs.partition_point(|g| g.0 < f.0);
            if pos < facs.len() && facs[pos].0 == f.0 {
                return None;
            }
            let mut v = Vec::with_capacity(facs.len() + 1);
            v.extend_from_slice(&facs[..pos]);
            v.push(f);
            v.extend_from_slice(&facs[pos..]);
            Some(Expr::from_kind(Kind::Mul(c.clone(), v)))
        }
        _ => None,
    }
}

fn mul_infinite(coeff: Number, facs: Vec<(Expr, Expr)>, infs: &[Constant], cache: &mut ProductMemo) -> Expr {
    if coeff.is_zero() {
        return Expr::nan();
    }
    let mut sign_neg = coeff.is_negative();
    let mut complex = false;
    for c in infs {
        match c {
            Constant::NegInfinity => sign_neg = !sign_neg,
            Constant::ComplexInfinity => complex = true,
            _ => {}
        }
    }
    // Finish the commutative powers before filtering their signs. A zoo
    // exposed here is a factor, rather than an initial scalar coefficient:
    // oo*Pow(zoo,1,False) retains (oo,zoo), while an initial zoo absorbs it.
    // Sending the completed factors through mul again loses that phase.
    let rest = finish_mul(Number::ONE, facs, 0, cache);
    if rest.is_zero() || rest.is_constant(Constant::NaN) || rest.as_number().is_some_and(Number::is_nan) {
        return Expr::nan();
    }
    let (number, tail) = as_coeff_mul(&rest);
    if number.is_negative() { sign_neg = !sign_neg; }
    let factors = if tail.is_one() { Vec::new() } else {
        match tail.kind() {
            Kind::Mul(_, factors) => factors.clone(),
            _ => vec![as_base_exp(&tail)],
        }
    };
    let mut retained = Vec::with_capacity(factors.len() + 1);
    for (base, exponent) in factors {
        let factor = mul_factor(&base, &exponent);
        let (zero, positive, negative, real) = crate::assumptions::exact_queries(|| (
            crate::assumptions::is_zero(&factor),
            crate::assumptions::fact(&factor,"extended_positive").or_else(||crate::assumptions::is_positive(&factor)),
            crate::assumptions::fact(&factor,"extended_negative").or_else(||crate::assumptions::is_negative(&factor)),
            crate::assumptions::fact(&factor,"extended_real"),
        ));
        if zero == Some(true) { return Expr::nan(); }
        if complex {
            // Infinity atoms are exactly nonzero, even where the generic
            // assumption query conservatively leaves is_zero unknown.
            if (zero == Some(false) && real.is_some())
                || matches!(factor.kind(),Kind::Const(c) if c.is_infinite()) { continue; }
        } else if positive == Some(true) { continue; }
        else if negative == Some(true) { sign_neg = !sign_neg; continue; }
        retained.push((base, exponent));
    }
    let constant = if complex { Constant::ComplexInfinity }
        else if sign_neg { Constant::NegInfinity } else { Constant::Infinity };
    if retained.is_empty() { return Expr::constant(constant); }
    // The scalar coefficient has its own slot. Only the remaining public
    // factors use Basic.compare, including late zoo and nested held nodes.
    retained.sort_by(|(a,x),(b,y)| crate::traverse::args_cmp(&mul_factor(a,x),&mul_factor(b,y)));
    retained.insert(0, (Expr::constant(constant), Expr::one()));
    let held = retained.iter().any(|(base, exponent)| exponent.is_one()
        && (base.is_mul() || base.is_constant(Constant::ComplexInfinity)));
    let kind = Kind::Mul(Number::ONE, retained);
    if held { Expr::from_unevaluated_kind(kind) } else { Expr::from_kind(kind) }
}

/// Canonical order of `Mul` factors: by base, then by exponent.
#[inline]
pub(crate) fn cmp_factor(a: &(Expr, Expr), b: &(Expr, Expr)) -> std::cmp::Ordering {
    a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1))
}

fn is_positive_int_base(b: &Expr) -> bool {
    matches!(b.kind(), Kind::Num(n) if n.is_integer() && n.is_positive())
}

fn exceptional_scalar_power(base: &Expr, exponent: &Expr) -> Option<Expr> {
    let exceptional = matches!(base.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN)
        || base.as_number().is_some_and(|number| number.is_zero()
            || matches!(number, Number::Float(value) if !value.is_finite()));
    exceptional.then(|| pow(base, exponent))
}

/// Numeric-base exceptional powers contribute to the coefficient during
/// collection. ComplexInfinity is not a SymPy Number, so its held powers are
/// deliberately evaluated later, with the other commutative powers.
fn collection_scalar_power(base: &Expr, exponent: &Expr) -> Option<Expr> {
    if base.is_constant(Constant::ComplexInfinity) || !exponent.is_rational_number() { return None; }
    let value = exceptional_scalar_power(base, exponent)?;
    if value.is_pow() { None } else { Some(value) }
}

fn finish_mul(mut coeff: Number, mut facs: Vec<(Expr, Expr)>, depth: u32, cache: &mut ProductMemo) -> Expr {
    // Exceptional scalar powers must be checked before exponent cancellation
    // and before zero absorption. Ordinary finite powers need no evaluation
    // here, particularly when their (possibly enormous) value is discarded.
    let mut scalar_zero = coeff.is_zero();
    let mut scalar_infinite = false;
    for (base, exponent) in &facs {
        if let Some(value) = exceptional_scalar_power(base, exponent) {
            if value.is_constant(Constant::NaN) || value.as_number().is_some_and(Number::is_nan) { return Expr::nan(); }
            scalar_zero |= value.is_zero();
            scalar_infinite |= matches!(value.kind(),Kind::Const(c) if c.is_infinite());
        }
    }
    if scalar_zero && scalar_infinite { return Expr::nan(); }
    facs.sort_unstable_by(cmp_factor);
    // merge equal bases
    let mut merged: Vec<(Expr, Expr, bool)> = Vec::with_capacity(facs.len());
    {
        let mut pending_exps: Vec<Expr> = Vec::new();
        let mut it = facs.into_iter().peekable();
        while let Some((b, e)) = it.next() {
            pending_exps.clear();
            pending_exps.push(e);
            while let Some((nb, _)) = it.peek() {
                if *nb == b {
                    pending_exps.push(it.next().unwrap().1);
                } else {
                    break;
                }
            }
            if pending_exps.len() == 1 {
                merged.push((b, pending_exps.pop().unwrap(), false));
            } else {
                // SymPy only combines exponents sharing the same non-numeric
                // part: x**y*x**(2*y) = x**(3*y) but x*x**y stays.
                let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
                for e in pending_exps.drain(..) {
                    let (_, t) = as_coeff_mul(&e);
                    match groups.iter_mut().find(|(gt, _)| *gt == t) {
                        Some((_, v)) => v.push(e),
                        None => groups.push((t, vec![e])),
                    }
                }
                for (_, mut v) in groups {
                    if v.len() == 1 {
                        merged.push((b.clone(), v.pop().unwrap(), false));
                    } else {
                        merged.push((b.clone(), add(&v), true));
                    }
                }
            }
        }
    }

    if merged.iter().any(|(base, exponent, _)| exponent.is_zero() && direct_arithmetic_infinity(base)) {
        return Expr::nan();
    }
    if merged.iter().any(|(base, exponent, _)| base.is_constant(Constant::NaN) && !exponent.is_zero()) {
        return Expr::nan();
    }
    if coeff.is_zero() || scalar_zero {
        // Scalar powers can expose 0**-1 or oo**1 before zero absorption.
        // Composite bases stay held boundaries, e.g. 0*(x*oo)**1 == 0.
        for (base, exponent, _) in &merged {
            if exponent.is_one() && base.has_nonfinite()
                && crate::assumptions::exact_queries(||crate::assumptions::is_finite(base))==Some(false){return Expr::nan();}
            if let Some(value) = exceptional_scalar_power(base, exponent) {
                if matches!(value.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN)
                    || value.as_number().is_some_and(Number::is_nan) { return Expr::nan(); }
            }
        }
        return Expr::zero();
    }

    let mut out: Vec<(Expr, Expr)> = Vec::with_capacity(merged.len());
    let mut pending: Vec<Expr> = Vec::new();
    let mut radicals: Vec<(Expr, Expr)> = Vec::new();

    // exp_polar(a)*exp_polar(b): SymPy's flatten sees both as powers of
    // exp_polar(1) and gathers exponents sharing the same non-numeric term
    // (I*pi with I*pi/2, x with 2*x), so opposite windings cancel to 1.
    // exp_polar(0) is its own base and is never merged.
    let is_winding = |b: &Expr| match b.kind() {
        Kind::Func(Func::ExpPolar, a) => !a[0].is_zero(),
        _ => false,
    };
    if merged.iter().filter(|(b, _, _)| is_winding(b)).nth(1).is_some() {
        let mut groups: Vec<(Expr, Vec<Expr>, Vec<(Expr, Expr, bool)>)> = Vec::new();
        merged.retain(|item| {
            if !is_winding(&item.0) {
                return true;
            }
            let Kind::Func(_, a) = item.0.kind() else { return true };
            let (c, term) = as_coeff_mul(&mul2(&a[0], &item.1));
            match groups.iter_mut().find(|(t, _, _)| *t == term) {
                Some((_, cs, orig)) => { cs.push(Expr::num(c)); orig.push(item.clone()); }
                None => groups.push((term, vec![Expr::num(c)], vec![item.clone()])),
            }
            false
        });
        for (term, cs, mut orig) in groups {
            if orig.len() == 1 {
                merged.push(orig.pop().unwrap());
                continue;
            }
            let total = mul2(&add(&cs), &term);
            if !total.is_zero() {
                pending.push(crate::functions::func(Func::ExpPolar, vec![total]));
            }
        }
    }

    // exp(a)*exp(b)*E: SymPy combines the arguments, then splits the sum
    // into one exp per term again (so exp(x)*exp(-x) = 1).
    let n_exp = merged
        .iter()
        .filter(|(b, e, _)| {
            e.is_one() && (b.is_func(&Func::Exp) || b.is_constant(Constant::E))
        })
        .count();
    if n_exp >= 2 {
        let mut exp_args: Vec<Expr> = Vec::with_capacity(n_exp);
        merged.retain(|(b, e, _)| {
            if !e.is_one() {
                return true;
            }
            match b.kind() {
                Kind::Func(Func::Exp, a) => {
                    exp_args.push(a[0].clone());
                    false
                }
                Kind::Const(Constant::E) => {
                    exp_args.push(Expr::one());
                    false
                }
                _ => true,
            }
        });
        let s = add(&exp_args);
        let mut parts: Vec<Expr> = Vec::new();
        match s.kind() {
            Kind::Add(c, terms) => {
                for (t, k) in terms {
                    parts.push(add_term(t, k));
                }
                if !c.is_zero() {
                    parts.push(Expr::num(c.clone()));
                }
            }
            _ => parts.push(s.clone()),
        }
        for p in parts {
            let v = crate::functions::exp(&p);
            match v.kind() {
                Kind::Func(Func::Exp, _) | Kind::Const(Constant::E) => out.push((v, Expr::one())),
                Kind::Num(n) => coeff = coeff.mul(n),
                _ => pending.push(v),
            }
        }
    }

    for (b, e, was_merged) in merged {
        if e.is_zero() {
            continue;
        }
        if e.is_one() && !b.is_number() {
            if b.is_constant(Constant::Infinity) || b.is_constant(Constant::NegInfinity) {
                pending.push(b);
            } else { out.push((b, e)); }
            continue;
        }
        let special_base = matches!(b.kind(), Kind::Num(_) | Kind::Const(_));
        if !(was_merged || special_base) {
            out.push((b, e));
            continue;
        }
        let v = pow(&b, &e);
        match v.kind() {
            Kind::Pow(vb, ve) if *vb == b => {
                if b.is_number() && ve.is_rational_number() {
                    radicals.push((b, ve.clone()));
                } else {
                    out.push((b, ve.clone()));
                }
            }
            Kind::Num(n) => coeff = coeff.mul(n),
            Kind::Const(Constant::ComplexInfinity) if b.is_constant(Constant::ComplexInfinity) => {
                out.push((v.clone(), Expr::one()));
            }
            Kind::Mul(c, fs) if fs.iter().all(|(fb, _)| fb.is_number()) => {
                coeff = coeff.mul(c);
                radicals.extend(fs.iter().cloned());
            }
            _ => pending.push(v),
        }
    }

    // SymPy's Mul.flatten splits a negative numeric base off as a power of
    // -1 and lets the total power of -1 join a positive base with the same
    // exponent: (-1)**(1/6)*2**(1/6) = (-2)**(1/6).
    let minus_one_exponent = signed_radicals(&mut radicals, &mut out, &mut pending);
    let first_out = out.len();
    // sqrt(2)*sqrt(3) -> sqrt(6): merge positive integer bases with equal
    // rational exponents.
    if radicals.len() > 1 {
        radicals.sort_unstable_by(|a, b| a.1.cmp(&b.1));
        let mut i = 0;
        while i < radicals.len() {
            let mut j = i + 1;
            while j < radicals.len() && radicals[j].1 == radicals[i].1 {
                j += 1;
            }
            let group = &radicals[i..j];
            let ints: Vec<&(Expr, Expr)> =
                group.iter().filter(|(b, _)| is_positive_int_base(b)).collect();
            if ints.len() > 1 {
                let mut prod = Number::ONE;
                for (b, _) in &ints {
                    prod = prod.mul(b.as_number().unwrap());
                }
                let v = pow(&Expr::num(prod), &group[0].1);
                match v.kind() {
                    Kind::Num(n) => coeff = coeff.mul(n),
                    Kind::Mul(c, fs) => {
                        coeff = coeff.mul(c);
                        out.extend(fs.iter().cloned());
                    }
                    Kind::Pow(b, e) => out.push((b.clone(), e.clone())),
                    _ => pending.push(v),
                }
                for (b, e) in group.iter().filter(|(b, _)| !is_positive_int_base(b)) {
                    out.push((b.clone(), e.clone()));
                }
            } else {
                out.extend(group.iter().cloned());
            }
            i = j;
        }
    } else {
        out.extend(radicals);
    }
    if let Some(total) = minus_one_exponent {
        join_minus_one(total, &mut coeff, &mut out, first_out, &mut pending);
    }

    if !pending.is_empty() {
        if depth < 8 {
            pending.push(Expr::num(coeff));
            for (b, e) in out {
                pending.push(mul_factor(&b, &e));
            }
            return mul_depth(&pending, depth + 1, cache);
        }
        // give up simplifying further but never drop factors
        for p in pending {
            match p.kind() {
                Kind::Num(n) => coeff = coeff.mul(n),
                Kind::Mul(c, fs) => {
                    coeff = coeff.mul(c);
                    out.extend(fs.iter().cloned());
                }
                _ => out.push(as_base_exp(&p)),
            }
        }
    }

    if coeff.is_nan() { return Expr::nan(); }
    if coeff.is_zero() {
        if out.iter().any(|(base, exponent)| exponent.is_one()
            && matches!(base.kind(),Kind::Const(c) if c.is_infinite())) { return Expr::nan(); }
        return Expr::zero();
    }
    if out.is_empty() {
        return Expr::num(coeff);
    }
    out.sort_unstable_by(cmp_factor);
    if out.len() == 1 {
        let (b, e) = &out[0];
        if coeff.is_one() {
            return mul_factor(b, e);
        }
        // Number * (a + b) distributes
        if e.is_one() {
            if let Kind::Add(c, terms) = b.kind() {
                if !b.is_canonical(){
                    let normalized=normalized_sum(b,&mut ArithmeticMemo::default());
                    return mul2(&Expr::num(coeff),&normalized);
                }
                return distribute(&coeff, c, terms);
            }
        }
    }
    // A power can reveal a Mul or zoo only after the active flattening pass.
    // Those public trees need evaluate=False on reconstruction; they are not
    // fixed points of another multiplication, even though this operation is
    // complete. Keep their presentation metadata honest for interop/pickling.
    let held = out.iter().any(|(base, exponent)| exponent.is_one()
        && (base.is_mul() || base.is_constant(Constant::ComplexInfinity)));
    if held {
        out.sort_by(|(a,x),(b,y)| crate::traverse::args_cmp(&mul_factor(a,x),&mul_factor(b,y)));
    }
    let kind = Kind::Mul(coeff, out);
    if held { Expr::from_unevaluated_kind(kind) } else { Expr::from_kind(kind) }
}

/// Splits `(-n)**e` radicals into `n**e` and returns the total exponent of
/// -1, counting a `(-1)**t` factor and `I` as `(-1)**(1/2)` like SymPy's
/// `neg1e`, or `None` when no radical has a negative base, which leaves the
/// product as it was.
fn signed_radicals(radicals: &mut Vec<(Expr, Expr)>, out: &mut Vec<(Expr, Expr)>, pending: &mut Vec<Expr>) -> Option<Number> {
    let rational_exponent = |e: &Expr| e.as_number().filter(|n| n.is_rational() && n.is_positive()).cloned();
    let negative_base = |b: &Expr| b.as_number().is_some_and(|n| n.is_integer() && n.is_negative());
    if !radicals.iter().any(|(b, e)| negative_base(b) && rational_exponent(e).is_some()) { return None; }
    let mut total = Number::ZERO;
    let half = Number::rational_i64(1, 2).unwrap();
    if let Some(index) = pending.iter().position(|p| p.is_constant(Constant::I)) {
        pending.remove(index);
        total = half.clone();
    } else if let Some(index) = out.iter().position(|(b, e)| b.is_constant(Constant::I) && e.is_one()) {
        out.remove(index);
        total = half.clone();
    }
    let mut kept = Vec::with_capacity(radicals.len());
    for (b, e) in radicals.drain(..) {
        match rational_exponent(&e) {
            Some(n) if b.is_minus_one() => total = total.add(&n),
            Some(n) if negative_base(&b) => {
                total = total.add(&n);
                kept.push((Expr::num(b.as_number().unwrap().neg()), e));
            }
            _ => kept.push((b, e)),
        }
    }
    *radicals = kept;
    Some(total)
}

/// SymPy's handling of `neg1e`: odd integer parts flip the coefficient, a
/// half gives `I`, and any other remainder joins the first positive integer
/// base with that exponent or stays a separate power of -1.
fn join_minus_one(total: Number, coeff: &mut Number, out: &mut Vec<(Expr, Expr)>, first: usize, pending: &mut Vec<Expr>) {
    let two = Number::int(2);
    let mut rest = total.sub(&total.div(&two).unwrap().floor().mul(&two));
    if rest.cmp_value(&Number::ONE) != Some(std::cmp::Ordering::Less) {
        *coeff = coeff.neg();
        rest = rest.sub(&Number::ONE);
    }
    if rest.is_zero() { return; }
    if rest == Number::rational_i64(1, 2).unwrap() {
        pending.push(Expr::i());
        return;
    }
    let exponent = Expr::num(rest);
    if let Some(slot) = out[first..].iter_mut().find(|(b, e)| is_positive_int_base(b) && *e == exponent) {
        slot.0 = Expr::num(slot.0.as_number().unwrap().neg());
        return;
    }
    out.push((Expr::minus_one(), exponent));
}

fn mul_depth(args: &[Expr], depth: u32, cache: &mut ProductMemo) -> Expr {
    mul_with_cache(args, depth, cache)
}

fn distribute(k: &Number, c: &Number, terms: &[(Expr, Number)]) -> Expr {
    if terms.iter().any(|(t, _)| t.may_have_order()||t.has_nonfinite()) {
        let mut v: Vec<Expr> = terms.iter().map(|(t, n)| mul2(&Expr::num(n.mul(k)), t)).collect();
        v.push(Expr::num(c.mul(k)));
        return add(&v);
    }
    let terms: Vec<(Expr, Number)> = terms.iter().map(|(t, n)| (t.clone(), n.mul(k))).collect();
    let c = c.mul(k);
    // coefficients may vanish only through float underflow; be safe
    if terms.iter().any(|(_, n)| n.is_zero()) {
        return add_from_terms(c, terms);
    }
    Expr::from_kind(Kind::Add(c, terms))
}

// ----------------------------------------------------------------------
// Pow
// ----------------------------------------------------------------------

/// Canonical power `b ** e`.
pub fn pow(b: &Expr, e: &Expr) -> Expr {
    if !e.is_zero() {
        if let Kind::Func(Func::ExpPolar, args) = b.kind() {
            return crate::functions::func(Func::ExpPolar, vec![mul2(&args[0], e)]);
        }
    }
    if let Kind::Num(en) = e.kind() {
        if en.is_zero() && en.is_rational() {
            return Expr::one();
        }
        if en.is_one() {
            return b.clone();
        }
        if let Kind::Num(bn) = b.kind() {
            return num_pow(bn, en);
        }
    }
    if e.is_constant(Constant::NaN) || b.is_constant(Constant::NaN) {
        return Expr::nan();
    }
    match b.kind() {
        Kind::Num(n) if n.is_one() => {
            if matches!(e.kind(), Kind::Const(c) if c.is_infinite()) {
                return Expr::nan();
            }
            return Expr::one();
        }
        Kind::Num(n) if n.is_zero() => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::zero() } else { Expr::zoo() };
            }
            if e.is_constant(Constant::Infinity) {
                return Expr::zero();
            }
        }
        Kind::Const(Constant::E) => return crate::functions::exp(e),
        // Piecewise._eval_power: the power of every piece
        Kind::Piecewise(pieces) => {
            return crate::relational::piecewise(pieces.iter().map(|(v, c)| (pow(v, e), c.clone())).collect());
        }
        Kind::Func(Func::Order, _) => {
            if let Some(r) = crate::series::order::pow_hook(b, e) {
                return r;
            }
        }
        Kind::Const(Constant::I) => {
            if let Some(k) = e.as_number().and_then(|n| n.to_i64()) {
                return i_pow(k);
            }
        }
        Kind::Const(Constant::Infinity) => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::oo() } else { Expr::zero() };
            }
        }
        Kind::Const(Constant::NegInfinity) => {
            if let Some(k) = e.as_i64() {
                return if k > 0 {
                    if k % 2 == 0 {
                        Expr::oo()
                    } else {
                        Expr::neg_oo()
                    }
                } else {
                    Expr::zero()
                };
            }
        }
        Kind::Const(Constant::ComplexInfinity) => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::zoo() } else { Expr::zero() };
            }
        }
        Kind::Pow(bb, be) => {
            if e.is_integer_number() {
                return pow(bb, &mul2(be, e));
            }
            // Pow._eval_power: a polar base never needs a branch factor.
            if crate::special::polar::is_polar(bb) {
                return pow(bb, &mul2(be, e));
            }
            // Pow._eval_power: the branch factor is 1 for a real exponent
            // a with |a| < 1, so (b**a)**e = b**(a*e) for every b and e.
            if let Kind::Num(a) = be.kind() {
                if a.abs().cmp_value(&Number::ONE) == Some(std::cmp::Ordering::Less) {
                    return pow(bb, &mul2(be, e));
                }
            }
            // (b**a)**e = b**(a*e) when b >= 0 and a is real
            if crate::assumptions::is_nonnegative(bb) == Some(true)
                && crate::assumptions::is_real(be) == Some(true)
            {
                return pow(bb, &mul2(be, e));
            }
            // (r**(2k))**e = Abs(r)**(2k*e) for real r
            if e.is_number()
                && crate::assumptions::is_real(bb) == Some(true)
                && be.as_i64().is_some_and(|k| k % 2 == 0)
            {
                return pow(&crate::functions::abs(bb), &mul2(be, e));
            }
        }
        Kind::Mul(..) if b.has_noncommutative_host() && crate::decompose::commutativity(b) != Some(true) => {
            return noncommutative_product_power(b, e);
        }
        Kind::Mul(c, facs) => {
            if facs.len() == 1 && facs[0].0.is_constant(Constant::I) && facs[0].1.is_one() && c.is_rational() {
                if let Some(value) = imaginary_square_root_power(c, e) {
                    return value;
                }
            }
            if e.is_integer_number() {
                // A held positive numeric product can collect to a unit.
                // Keep its power as a multiplicative boundary: multiplying
                // another held product by it must perform active collection,
                // while multiplication by a literal 1 preserves that operand.
                if !b.is_canonical() && c.is_positive() && facs.iter().all(|(base, exponent)|
                    exponent.is_one() && base.as_number().is_some_and(Number::is_positive)) {
                    let coefficient = facs.iter().fold(c.clone(), |value, (base, _)|
                        value.mul(base.as_number().unwrap()));
                    if coefficient.is_one() {
                        return Expr::from_unevaluated_kind(Kind::Pow(Expr::one(), e.clone()));
                    }
                }
                let mut v = Vec::with_capacity(facs.len() + 1);
                v.push(num_pow(c, e.as_number().unwrap()));
                for (fb, fe) in facs {
                    v.push(pow(fb, &mul2(fe, e)));
                }
                return mul(&v);
            }
            // pull positive factors out of a radical: sqrt(p*x) = sqrt(p)*sqrt(x)
            if e.is_rational_number() {
                let (pos, other): (Vec<&(Expr, Expr)>, Vec<&(Expr, Expr)>) = facs
                    .iter()
                    .partition(|(fb, fe)| !fb.is_number() && (crate::special::polar::is_polar(fb)
                        || crate::assumptions::is_positive(&mul_factor(fb, fe)) == Some(true)));
                if !pos.is_empty() {
                    let mut v: Vec<Expr> = pos.iter().map(|(fb, fe)| pow(&mul_factor(fb, fe), e)).collect();
                    let mut rest = vec![Expr::num(c.clone())];
                    rest.extend(other.iter().map(|(fb, fe)| mul_factor(fb, fe)));
                    v.push(pow(&mul(&rest), e));
                    return mul(&v);
                }
            }
            if e.is_number() && !c.is_one() && (c.is_rational() || c.is_float()) {
                let rest = mul_from_factors(facs.clone());
                if c.is_positive() && facs.iter().all(|(fb, fe)| fb.as_number().is_some_and(|n| n.is_rational() && n.is_positive()) && fe.is_rational_number()) {
                    // SymPy's expand_power_base: every factor is positive,
                    // so sqrt(sqrt(2)/2) = sqrt(1/2)*2**(1/4) = 2**(3/4)/2.
                    let mut v = vec![num_pow(c, e.as_number().unwrap())];
                    v.extend(facs.iter().map(|(fb, fe)| pow(fb, &mul2(fe, e))));
                    return mul(&v);
                }
                if c.is_positive() {
                    return mul2(&num_pow(c, e.as_number().unwrap()), &pow_raw(&rest, e));
                } else if !c.is_minus_one() {
                    let neg_rest = Expr::from_kind(Kind::Mul(Number::MINUS_ONE, facs.clone()));
                    return mul2(&num_pow(&c.neg(), e.as_number().unwrap()), &pow_raw(&neg_rest, e));
                }
            }
        }
        Kind::Add(c, terms)
            if terms.len() == 1 && terms[0].0.is_constant(Constant::I) && e.is_rational_number() =>
        {
            if let Some(value) = pure_complex_power(c, &terms[0].1, e.as_number().unwrap()) {
                return value;
            }
        }
        Kind::Func(Func::Abs, args) => {
            // Abs(r)**(2k) = r**(2k) for real r
            if e.as_i64().is_some_and(|k| k % 2 == 0) && crate::assumptions::is_real(&args[0]) == Some(true) {
                return pow(&args[0], e);
            }
            // Abs(r)**k = r**(k-1)*Abs(r) for odd k != -1 (SymPy's Abs._eval_power)
            if let Some(k) = e.as_i64().filter(|k| k % 2 != 0 && *k != 1 && *k != -1) {
                if crate::assumptions::is_real(&args[0]) == Some(true) {
                    return mul2(&pow(&args[0], &Expr::int(k - 1)), b);
                }
            }
        }
        Kind::Num(n) if n.is_minus_one() => {
            // (-1)**(2k) = 1, (-1)**(2k+1) = -1
            if crate::assumptions::is_integer(e) == Some(true) {
                match crate::assumptions::parity_even(e) {
                    Some(true) => return Expr::one(),
                    Some(false) => return Expr::minus_one(),
                    None => {}
                }
            }
        }
        Kind::Func(Func::Special(crate::special::SpecialFn::KroneckerDelta), _) => {
            // KroneckerDelta._eval_power: idempotent for positive powers
            if crate::assumptions::is_positive(e) == Some(true) {
                return b.clone();
            }
            if crate::assumptions::is_negative(e) == Some(true) && !e.is_minus_one() {
                return pow(b, &Expr::minus_one());
            }
        }
        Kind::Func(Func::Exp, args) => {
            if let Some(value) = exp_power(&args[0], e) { return value; }
        }
        _ => {}
    }
    memo::apply(memo::Operation::Power,b,e,||pow_raw(b,e))
}

/// SymPy's `Add._eval_power` for a literal `r + i*I` with a rational
/// exponent: the reciprocal is the held product `(r - i*I)/(r**2 + i**2)`,
/// and a half-integer power of a Pythagorean base is its exact root.
fn pure_complex_power(r: &Number, i: &Number, e: &Number) -> Option<Expr> {
    let finite = |n: &Number| !n.is_nan() && (!n.is_float() || n.to_f64().is_finite());
    if !finite(r) || !finite(i) {
        return None;
    }
    let (p, q) = e.numer_denom()?;
    if q == num_bigint::BigInt::from(2) && r.is_rational() && i.is_rational() {
        let norm = r.mul(r).add(&i.mul(i));
        let (n, d) = norm.numer_denom()?;
        let (n, d) = (crate::number::exact_root(&n, 2)?, crate::number::exact_root(&d, 2)?);
        let modulus = Number::rational(n, d)?;
        let p = i64::try_from(&p).ok()?;
        let half = modulus.sub(r).div(&Number::int(2))?;
        let root = pow(&pow(&Expr::num(half), &Expr::num(Number::rational_i64(1, 2)?)), &Expr::int(p));
        let a = modulus.add(r).div(&i.abs())?;
        let b = if i.is_negative() { Number::MINUS_ONE } else { Number::ONE };
        let (u, v) = gaussian_rational_power(&a, &b, p.unsigned_abs());
        let expansion = if p < 0 {
            gaussian_reciprocal(&u, &v)?
        } else {
            add(&[Expr::num(u), mul2(&Expr::num(v), &Expr::i())])
        };
        return Some(mul2(&root, &expansion));
    }
    if e.is_minus_one() {
        return gaussian_reciprocal(r, i);
    }
    None
}

/// `1/(r + i*I)` as the reference's unevaluated `(r - i*I)*(1/(r**2+i**2))`.
fn gaussian_reciprocal(r: &Number, i: &Number) -> Option<Expr> {
    let norm = r.mul(r).add(&i.mul(i));
    let conjugate = add(&[Expr::num(r.clone()), mul2(&Expr::num(i.neg()), &Expr::i())]);
    Some(crate::simptools::unevaluated_mul(&[conjugate, Expr::num(norm.inv()?)]))
}

/// `(a + b*I)**k` with exact rational parts.
fn gaussian_rational_power(a: &Number, b: &Number, k: u64) -> (Number, Number) {
    let (mut u, mut v) = (Number::ONE, Number::ZERO);
    for _ in 0..k {
        let next = (u.mul(a).sub(&v.mul(b)), u.mul(b).add(&v.mul(a)));
        u = next.0;
        v = next.1;
    }
    (u, v)
}

/// SymPy's `Mul._eval_power` for `a*I` with `a` rational and a
/// half-integer exponent: `(2*a*I)**(1/2)` with a square `|a|/2` is exact.
fn imaginary_square_root_power(a: &Number, e: &Expr) -> Option<Expr> {
    let (p, q) = e.as_number()?.numer_denom()?;
    if q != num_bigint::BigInt::from(2) {
        return None;
    }
    let (n, d) = a.abs().div(&Number::int(2))?.numer_denom()?;
    let r = Number::rational(crate::number::exact_root(&n, 2)?, crate::number::exact_root(&d, 2)?)?;
    let p = Expr::int(i64::try_from(&p).ok()?);
    let unit = add(&[Expr::one(), if a.is_negative() { mul2(&Expr::minus_one(), &Expr::i()) } else { Expr::i() }]);
    Some(crate::simptools::unevaluated_mul(&[pow(&Expr::num(r), &p), pow(&unit, &p)]))
}

/// SymPy's exp._eval_power, i.e. Pow._eval_power of E**z: exp(z)**w is
/// exp(z*w) for an integer w or an extended real z; for a non-real z the
/// branch factor exp(2*pi*I*w*floor(1/2 - im(z)/(2*pi))) must be +-1.
fn exp_power(z: &Expr, w: &Expr) -> Option<Expr> {
    use crate::assumptions::{exact_queries, is_integer, is_real, real_signs};
    let product = || crate::functions::exp(&mul2(z, w));
    if w.is_integer_number() || exact_queries(|| is_integer(w)) == Some(true) { return Some(product()); }
    if real_signs(z).is_some() { return Some(product()); }
    if exact_queries(|| is_real(z)) != Some(false) { return None; }
    let imaginary = crate::functions::im(z);
    if !crate::functions::is_numeric_value(&imaginary) { return None; }
    let k = crate::functions::floor(&sub(&Expr::rational(1, 2), &div(&imaginary, &mul2(&Expr::int(2), &Expr::pi()))));
    if !k.is_integer_number() { return None; }
    if k.is_zero() { return Some(product()); }
    let turns = mul(&[Expr::int(2), w.clone(), k]);
    let (n, _) = turns.as_number().filter(|n| n.is_integer())?.numer_denom()?;
    Some(if num_integer::Integer::is_even(&n) { product() } else { neg(&product()) })
}

/// `b ** e` without simplification beyond trivial exponents.
/// SymPy's `Mul._eval_power` of a product that does not commute: an
/// integer power distributes over the commutative factors only, a rational
/// power over its nonnegative commutative factors (`expand_power_base`),
/// and the noncommutative product keeps its power.
fn noncommutative_product_power(b: &Expr, e: &Expr) -> Expr {
    let (commutative, ordered) = crate::decompose::args_cnc(b, false);
    let integer = e.is_integer_number();
    if !integer && !e.is_rational_number() { return Expr::from_kind(Kind::Pow(b.clone(), e.clone())); }
    let (mut outside, mut inside) = (Vec::new(), Vec::new());
    for factor in commutative {
        if integer || crate::assumptions::is_nonnegative(&factor) == Some(true) { outside.push(pow(&factor, e)); }
        else { inside.push(factor); }
    }
    inside.extend(ordered);
    outside.push(Expr::from_kind(Kind::Pow(mul(&inside), e.clone())));
    mul(&outside)
}

fn pow_raw(b: &Expr, e: &Expr) -> Expr {
    if e.is_one() {
        return b.clone();
    }
    if e.is_zero() {
        return Expr::one();
    }
    Expr::from_kind(Kind::Pow(b.clone(), e.clone()))
}

fn i_pow(k: i64) -> Expr {
    match k.rem_euclid(4) {
        0 => Expr::one(),
        1 => Expr::i(),
        2 => Expr::minus_one(),
        _ => Expr::from_kind(Kind::Mul(Number::MINUS_ONE, vec![(Expr::i(), Expr::one())])),
    }
}

/// Power of two numbers (exact radicals are simplified like SymPy).
pub fn num_pow(b: &Number, e: &Number) -> Expr {
    if b.is_rational() && (b.is_zero() || b.is_one()) {
        if e.is_zero() || b.is_one() { return Expr::one(); }
        return if e.is_positive() { Expr::zero() } else { Expr::zoo() };
    }
    if e.is_zero() && e.is_rational() {
        return Expr::one();
    }
    if e.is_one() {
        return Expr::num(b.clone());
    }
    if let Some(k) = e.to_i64() {
        return match b.pow_i64(k) {
            Some(v) => Expr::num(v),
            None if matches!(b, Number::MpFloat(..)) && !b.is_zero() => Expr::from_kind(Kind::Pow(Expr::num(b.clone()), Expr::num(e.clone()))),
            None => Expr::zoo(),
        };
    }
    if e.is_integer() {
        // huge integer exponent
        if b.is_one() {
            return Expr::one();
        }
        if b.is_zero() {
            return if e.is_positive() { Expr::zero() } else { Expr::zoo() };
        }
        if b.is_minus_one() {
            let odd = e.to_bigint().unwrap().is_odd();
            return Expr::int(if odd { -1 } else { 1 });
        }
        return Expr::from_kind(Kind::Pow(Expr::num(b.clone()), Expr::num(e.clone())));
    }
    if b.is_float() || e.is_float() {
        // SymPy's Float._eval_power: a negative float to 1/q with q odd keeps
        // the principal root of -1 symbolic.
        if b.is_float() && b.is_negative() && e.is_rational() && e.numer().is_one()
            && e.denom().to_bigint().is_some_and(|q| q.is_odd()) {
            return mul2(&Expr::from_kind(Kind::Pow(Expr::minus_one(), Expr::num(e.clone()))), &num_pow(&b.neg(), e));
        }
        if matches!(b, Number::MpFloat(..)) || matches!(e, Number::MpFloat(..)) || (b.is_negative() && !b.is_zero()) {
            // A complex result of a negative base is mpmath's mpc_pow at the
            // float precision; double precision trigonometry would leave a
            // spurious real part.
            let precision = b.float_precision().unwrap_or(0).max(e.float_precision().unwrap_or(0));
            let raw = Expr::from_kind(Kind::Pow(Expr::num(b.clone()), Expr::num(e.clone())));
            return crate::evalf_mp::eval_mp(&raw, precision)
                .map(|value| crate::evalf_mp::to_expr(value, precision)).unwrap_or(raw);
        }
        let (x, y) = (b.to_f64(), e.to_f64());
        if x >= 0.0 || y.fract() == 0.0 {
            return Expr::float(x.powf(y));
        }
        let r = (-x).powf(y);
        let t = std::f64::consts::PI * y;
        return add2(
            &Expr::float(r * t.cos()),
            &Expr::from_kind(Kind::Mul(Number::float(r * t.sin()), vec![(Expr::i(), Expr::one())])),
        );
    }
    // exact base, non-integer rational exponent
    if b.is_zero() {
        return if e.is_positive() { Expr::zero() } else { Expr::zoo() };
    }
    if b.is_one() {
        return Expr::one();
    }
    if b.is_negative() {
        if b.is_minus_one() {
            return minus_one_pow(e);
        }
        if e.is_positive() {
            // SymPy's Integer._eval_power: a negative base stays a power
            // ((-2)**(1/6)) unless sqrt or a root extraction changes the
            // positive part; Rational._eval_power keeps the sign on the
            // numerator: (-2/3)**(1/6) = (-2)**(1/6)*3**(5/6)/3.
            if !b.is_integer() {
                if let Number::Rational(r) = b {
                    let n = Number::from_bigint(r.numer().clone());
                    let d = Number::from_bigint(r.denom().clone());
                    return mul2(&num_pow(&n, e), &num_pow(&d, &e.neg()));
                }
            } else if *e != Number::rational_i64(1, 2).unwrap() {
                let positive = b.neg();
                let value = num_pow(&positive, e);
                if let Kind::Pow(base, exponent) = value.kind() {
                    if base.as_number() == Some(&positive) && exponent.as_number() == Some(e) {
                        return Expr::from_kind(Kind::Pow(Expr::num(b.clone()), Expr::num(e.clone())));
                    }
                }
                return mul2(&minus_one_pow(e), &value);
            }
        }
        return mul2(&minus_one_pow(e), &num_pow(&b.neg(), e));
    }
    match b {
        Number::Rational(r) => {
            let n = Number::from_bigint(r.numer().clone());
            let d = Number::from_bigint(r.denom().clone());
            mul2(&num_pow(&n, e), &num_pow(&d, &e.neg()))
        }
        _ => int_rational_pow(&b.to_bigint().unwrap(), e),
    }
}

/// `(-1) ** (p/q)` with the exponent reduced into (-1, 1].
fn minus_one_pow(e: &Number) -> Expr {
    // t = e - 2*round_half_down(e/2)  in (-1, 1]
    let two = Number::int(2);
    let half = e.div(&two).unwrap();
    let mut k = half.floor();
    let mut t = e.sub(&k.mul(&two));
    // t in [0, 2): bring into (-1, 1]
    if t.cmp_value(&Number::ONE) == Some(std::cmp::Ordering::Greater) {
        k = k.add(&Number::ONE);
        t = e.sub(&k.mul(&two));
    }
    if t.is_zero() {
        return Expr::one();
    }
    if t.is_one() {
        return Expr::minus_one();
    }
    let half = Number::rational_i64(1, 2).unwrap();
    if t == half {
        return Expr::i();
    }
    if t == half.neg() {
        return i_pow(3);
    }
    // SymPy keeps the exponent in (0, 1): (-1)**(-1/3) = -(-1)**(2/3)
    if t.is_negative() {
        let t2 = t.add(&Number::ONE);
        return Expr::from_kind(Kind::Mul(
            Number::MINUS_ONE,
            vec![(Expr::minus_one(), Expr::num(t2))],
        ));
    }
    Expr::from_kind(Kind::Pow(Expr::minus_one(), Expr::num(t)))
}

const TRIAL_LIMIT: u64 = 1 << 15;

/// `n ** (p/q)` for a positive integer `n` and a non-integer rational.
fn int_rational_pow(n: &BigInt, e: &Number) -> Expr {
    let (p, q) = e.numer_denom().unwrap();
    // 1. perfect power reduction: n = m**g  =>  n**e = m**(e*g)
    let (factors, cof) = trial_factor(n, TRIAL_LIMIT);
    let mut g: u64 = 0;
    for (_, m) in &factors {
        g = g.gcd(&(*m as u64));
    }
    if !cof.is_one() {
        // the unfactored cofactor: find its largest perfect-power exponent
        let bits = cof.bits();
        let mut k = (bits / 16).clamp(1, 64);
        let mut cof_root_exp = 1u64;
        while k >= 2 {
            if exact_root(&cof, k as u32).is_some() {
                cof_root_exp = k;
                break;
            }
            k -= 1;
        }
        g = if g == 0 { cof_root_exp } else { g.gcd(&cof_root_exp) };
    }
    if g > 1 {
        let m = exact_root(n, g as u32).unwrap();
        let e2 = e.mul(&Number::int(g as i64));
        return num_pow(&Number::from_bigint(m), &e2);
    }
    // 2. split integer part: e = k + r/q
    let k = e.floor();
    let r = p - k.to_bigint().unwrap() * &q;
    let mut outside = Number::from_bigint(n.clone()).pow_i64(k.to_i64().unwrap_or(0)).unwrap();
    // 3. extract q-th powers: n = A**q * M
    let q_u = q.to_u64().unwrap_or(u64::MAX);
    let mut a = BigInt::one();
    let mut m = BigInt::one();
    for (prime, mult) in &factors {
        let mult = *mult as u64;
        let out_e = mult / q_u;
        let in_e = mult % q_u;
        if out_e > 0 {
            a *= num_traits::pow::Pow::pow(prime, out_e);
        }
        if in_e > 0 {
            m *= num_traits::pow::Pow::pow(prime, in_e);
        }
    }
    if !cof.is_one() {
        match q.to_u32().and_then(|qq| exact_root(&cof, qq)) {
            Some(rt) => a *= rt,
            None => m *= &cof,
        }
    }
    if !a.is_one() {
        let ar = num_traits::pow::Pow::pow(&a, r.to_u64().unwrap());
        outside = outside.mul(&Number::from_bigint(ar));
    }
    if m.is_one() {
        return Expr::num(outside);
    }
    let rad_e = Number::rational(r, q).unwrap();
    let rad = (Expr::num(Number::from_bigint(m)), Expr::num(rad_e));
    if outside.is_one() {
        Expr::from_kind(Kind::Pow(rad.0, rad.1))
    } else {
        Expr::from_kind(Kind::Mul(outside, vec![rad]))
    }
}

// ----------------------------------------------------------------------
// derived operations
// ----------------------------------------------------------------------

pub fn neg(a: &Expr) -> Expr {
    if (a.has_held_arithmetic() && (a.is_add()||a.is_mul())) || (a.is_add()&&a.has_nonfinite()){
        return neg_held(a,&mut ArithmeticMemo::default(),&mut ArithmeticMemo::default());
    }
    neg_canonical(a)
}

fn neg_held(a:&Expr,negative:&mut ArithmeticMemo,sums:&mut ArithmeticMemo)->Expr{
    if let Some((_,value))=negative.get(&a.dag_id()){return value.clone();}
    let value=if a.is_add(){
        let values=a.args().iter().map(|term|neg_held(term,negative,sums)).collect::<Vec<_>>();
        add_with_cache(&values,sums)
    }else if a.is_mul() && !a.is_canonical(){
        // Keep the tuple tail: a single factor can itself be a held Mul.
        let(coefficient,mut factors)=crate::decompose::as_coeff_mul(a,&[],true);
        let coefficient=if factors.first().is_some_and(|factor|factor.is_constant(Constant::ComplexInfinity)){coefficient}else{neg_canonical(&coefficient)};
        if !coefficient.is_one(){
            if factors.first().is_some_and(|factor|matches!(factor.kind(),Kind::Num(_)|Kind::Const(Constant::Infinity|Constant::NegInfinity|Constant::NaN))){
                factors[0]=mul2(&coefficient,&factors[0]);
            }else{factors.insert(0,coefficient);}
        }
        match factors.len(){0=>Expr::one(),1=>factors.pop().unwrap(),_=>Expr::from_unevaluated_kind(Kind::Mul(Number::ONE,factors.into_iter().map(|factor|(factor,Expr::one())).collect()))}
    }else{neg_canonical(a)};
    negative.insert(a.dag_id(),(a.clone(),value.clone()));value
}

fn neg_canonical(a:&Expr)->Expr{
    if a.is_mul() && a.has_nonfinite() {
        return mul(&[Expr::minus_one(),a.clone()]);
    }
    if a.may_have_order() {
        return mul(&[Expr::minus_one(), a.clone()]);
    }
    match a.kind() {
        Kind::Num(n) => Expr::num(n.neg()),
        Kind::Mul(c, facs) => {
            let c = c.neg();
            if c.is_one() {
                mul_from_factors(facs.clone())
            } else {
                Expr::from_kind(Kind::Mul(c, facs.clone()))
            }
        }
        Kind::Add(c, terms) => distribute(&Number::MINUS_ONE, c, terms),
        _ => mul2(&Expr::minus_one(), a),
    }
}

pub fn sub(a: &Expr, b: &Expr) -> Expr {
    add2(a, &neg(b))
}

pub fn div(a: &Expr, b: &Expr) -> Expr {
    if let (Some(numerator),Some(denominator))=(a.as_number(),b.as_number()){
        if !denominator.is_zero(){if let Some(value)=numerator.div(denominator){return Expr::num(value);}}
    }
    mul2(a, &pow(b, &Expr::minus_one()))
}

pub fn inv(a: &Expr) -> Expr {
    pow(a, &Expr::minus_one())
}

pub fn sqrt(a: &Expr) -> Expr {
    pow(a, &Expr::rational(1, 2))
}

pub fn powi(a: &Expr, n: i64) -> Expr {
    pow(a, &Expr::int(n))
}

/// Is `e` (syntactically) "negative looking", i.e. has a negative leading
/// numeric coefficient (SymPy's `could_extract_minus_sign` heuristic)?
pub fn could_extract_minus_sign(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Const(Constant::NegInfinity) => true,
        Kind::Mul(c, _) => c.is_negative(),
        Kind::Add(c, terms) => {
            let neg = terms.iter().filter(|(_, k)| k.is_negative()).count();
            let pos = terms.len() - neg;
            if neg != pos {
                return neg > pos;
            }
            if !c.is_zero() {
                return c.is_negative();
            }
            // tie: decide on the structurally first term for determinism
            let first = terms.iter().min_by(|a, b| a.0.cmp_struct(&b.0)).unwrap();
            first.1.is_negative()
        }
        _ => false,
    }
}

// operator sugar ---------------------------------------------------------

macro_rules! impl_binop {
    ($tr:ident, $m:ident, $f:expr) => {
        impl std::ops::$tr<Expr> for Expr {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(&self, &o)
            }
        }
        impl std::ops::$tr<&Expr> for &Expr {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(self, o)
            }
        }
        impl std::ops::$tr<&Expr> for Expr {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(&self, o)
            }
        }
        impl std::ops::$tr<Expr> for &Expr {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(self, &o)
            }
        }
        impl std::ops::$tr<i64> for Expr {
            type Output = Expr;
            fn $m(self, o: i64) -> Expr {
                $f(&self, &Expr::int(o))
            }
        }
        impl std::ops::$tr<i64> for &Expr {
            type Output = Expr;
            fn $m(self, o: i64) -> Expr {
                $f(self, &Expr::int(o))
            }
        }
        impl std::ops::$tr<Expr> for i64 {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(&Expr::int(self), &o)
            }
        }
        impl std::ops::$tr<&Expr> for i64 {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(&Expr::int(self), o)
            }
        }
    };
}

impl_binop!(Add, add, add2);
impl_binop!(Sub, sub, sub);
impl_binop!(Mul, mul, mul2);
impl_binop!(Div, div, div);

impl std::ops::Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        neg(&self)
    }
}
impl std::ops::Neg for &Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        neg(self)
    }
}

impl Expr {
    pub fn pow(&self, e: &Expr) -> Expr {
        pow(self, e)
    }
    pub fn powi(&self, n: i64) -> Expr {
        powi(self, n)
    }
    pub fn sqrt(&self) -> Expr {
        sqrt(self)
    }
}

#[cfg(test)]
mod held_arithmetic_tests {
    use super::*;
    fn p(source:&str)->Expr{crate::parser::parse_str(source).unwrap()}

    #[test]
    fn raw_product_tails_survive_coefficient_extraction_and_negation(){
        for(source,expected)in[
            ("Mul(-1,x,1/2,evaluate=False)","Mul(x,1/2,evaluate=False)"),
            ("Mul(y,-2,I,evaluate=False)","Mul(-1,y,-2,I,evaluate=False)"),
            ("Mul(-2.5,x,2,evaluate=False)","Mul(2.5,x,2,evaluate=False)"),
        ]{
            let value=p(source);let wanted=p(expected);
            assert_eq!(crate::sortkey::sort_key(&neg(&value)),crate::sortkey::sort_key(&wanted),"{source}");
            assert!(value.has_held_arithmetic());
        }
        assert!(!p("x+y").has_held_arithmetic());
        assert!(p("sin(Mul(x,2,evaluate=False))").has_held_arithmetic());
        let value=p("Mul(-1,x,1/2,evaluate=False)");
        let(c,tail)=as_coeff_mul(&value);
        assert!(c.is_minus_one());assert!(!tail.is_canonical());
        assert_eq!(tail.args(),vec![p("x"),p("1/2")]);
        assert_eq!(add_term(&tail,&c).args(),value.args());
    }

    #[test]
    fn held_sum_flattening_reclassifies_constants_and_shares_deep_dags(){
        assert_eq!(add(&[p("Add(1,1,evaluate=False)"),Expr::int(2)]),Expr::int(4));
        assert_eq!(mul2(&Expr::int(2),&p("Add(1,1,evaluate=False)")),Expr::int(4));
        let mut value=p("x");
        for _ in 0..80{value=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![(value.clone(),Number::ONE),(value,Number::ONE)]));}
        let coefficient=Expr::num(Number::from_bigint(num_bigint::BigInt::from(1)<<80));
        let expected=mul2(&coefficient,&p("x"));
        assert_eq!(add(&[value.clone(),Expr::one()]),add(&[expected.clone(),Expr::one()]));
        assert_eq!(neg(&value),neg(&expected));
        assert_eq!(mul2(&Expr::int(2),&value),mul2(&Expr::int(2),&expected));
    }
}

#[cfg(test)]
mod nonfinite_add_tests{
    use super::*;
    #[test]
    fn infinities_are_reclassified_when_flattening_and_distributing(){
        let x=Expr::symbol("x");let real=Expr::symbol_with("r",crate::Assumptions::REAL);
        assert_eq!(add2(&Expr::oo(),&real),Expr::oo());
        assert_eq!(add2(&Expr::zoo(),&Expr::i()),Expr::zoo());
        let value=add2(&Expr::oo(),&x);
        assert_eq!(add2(&value,&Expr::neg_oo()),Expr::nan());
        assert_eq!(add2(&value,&Expr::int(2)),value);
        assert_eq!(neg(&value),add2(&Expr::neg_oo(),&neg(&x)));
        assert_eq!(mul2(&Expr::int(2),&value),add2(&Expr::oo(),&mul2(&Expr::int(2),&x)));
        let held=crate::parser::parse_str("Add(-2,Add(I,2,oo,evaluate=False),evaluate=False)").unwrap();
        assert_eq!(neg(&held),add2(&Expr::neg_oo(),&neg(&Expr::i())));
    }
}

#[cfg(test)]
mod active_product_tests {
    use super::*;
    fn p(source: &str) -> Expr { crate::parser::parse_str(source).unwrap() }
    fn raw_mul(arguments: Vec<Expr>) -> Expr {
        Expr::from_unevaluated_kind(Kind::Mul(Number::ONE, arguments.into_iter().map(|argument| (argument,Expr::one())).collect()))
    }
    #[test]
    fn active_products_reclassify_raw_numeric_and_nested_factors() {
        let x=p("x");let y=p("y");let z=p("z");
        let raw=raw_mul(vec![x.clone(),Expr::int(2)]);
        assert_eq!(mul2(&Expr::int(2),&raw),p("4*x"));
        assert_eq!(mul(&[raw.clone(),z.clone()]),p("2*x*z"));
        assert_eq!(mul2(&raw,&y),p("2*x*y"));
        assert_eq!(mul(&[raw_mul(vec![raw,y.clone()]),z]),p("2*x*y*z"));
        let wrong_order=raw_mul(vec![y,x.clone()]);
        assert_eq!(mul2(&wrong_order,&x),p("x**2*y"));
        // A child coefficient must not distribute across its Add prematurely.
        let raw=raw_mul(vec![Expr::int(2),p("x+y")]);
        assert_eq!(mul(&[raw,p("z")]),p("2*z*(x+y)"));
    }
    #[test]
    fn identity_and_unary_negation_preserve_the_raw_tail() {
        let raw=raw_mul(vec![p("x"),Expr::int(2)]);
        assert!(mul2(&Expr::one(),&raw).ptr_eq(&raw));
        assert!(mul(&[Expr::one(),raw.clone(),Expr::one()]).ptr_eq(&raw));
        let negative=neg(&raw);
        assert!(!negative.is_canonical());
        assert_eq!(negative.args(),vec![Expr::minus_one(),p("x"),Expr::int(2)]);
        assert_eq!(mul2(&Expr::minus_one(),&raw),p("-2*x"));
    }
    #[test]
    fn explicit_infinities_are_seen_before_zero_absorption() {
        for infinite in [Expr::oo(),Expr::neg_oo(),Expr::zoo(),Expr::nan()] {
            let raw=raw_mul(vec![p("x"),infinite.clone()]);
            assert_eq!(mul2(&Expr::zero(),&raw),Expr::nan());
            assert_eq!(mul(&[raw,Expr::zero()]),Expr::nan());
            assert_eq!(mul2(&Expr::zero(),&mul(&[p("x"),infinite])),Expr::nan());
        }
        let body=raw_mul(vec![Expr::zero(),Expr::oo()]);
        assert_eq!(mul2(&Expr::int(2),&body),Expr::nan());
        assert_eq!(mul2(&Expr::zero(),&p("f(oo)")),Expr::zero());
        assert_eq!(mul2(&Expr::zero(),&p("UnevaluatedExpr(oo)")),Expr::zero());
        let sum=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![(p("x"),Number::ONE),(Expr::oo(),Number::ONE)]));
        assert_eq!(mul2(&Expr::zero(),&sum),Expr::nan());
    }
    #[test]
    fn held_powers_are_boundaries_with_scalar_zero_exceptions() {
        let raw=raw_mul(vec![p("x"),Expr::oo()]);
        let power=|base:Expr,n:i64|Expr::from_unevaluated_kind(Kind::Pow(base,Expr::int(n)));
        assert_eq!(mul2(&Expr::zero(),&power(raw.clone(),1)),Expr::zero());
        let scaled=mul2(&Expr::int(2),&power(raw.clone(),1));
        assert_eq!(scaled.args(),vec![Expr::int(2),raw.clone()]);
        for coefficient in [Expr::zero(),Expr::int(2)] {
            assert_eq!(mul2(&coefficient,&power(raw.clone(),0)),Expr::nan());
        }
        assert_eq!(mul2(&Expr::zero(),&power(Expr::oo(),0)),Expr::zero());
        assert_eq!(mul2(&Expr::zero(),&power(Expr::oo(),1)),Expr::nan());
        assert_eq!(mul2(&Expr::zero(),&power(Expr::zero(),-1)),Expr::nan());
        assert_eq!(mul2(&Expr::int(2),&power(Expr::oo(),1)),Expr::oo());
        assert_eq!(mul(&[power(Expr::zero(),1),power(Expr::zoo(),1)]),Expr::nan());
        assert_eq!(mul(&[power(Expr::oo(),1),power(Expr::oo(),-1)]),Expr::nan());
        assert_eq!(mul(&[power(Expr::oo(),0),power(Expr::oo(),0)]),Expr::one());
        assert_eq!(mul2(&Expr::zero(),&power(Expr::int(2),1_000_000_000)),Expr::zero());
        assert_eq!(mul(&[power(Expr::zero(),1),power(Expr::int(2),1_000_000_000)]),Expr::zero());
    }
    #[test]
    fn nonfinite_collection_distinguishes_coefficients_from_late_factors() {
        let power=|base:Expr,n:i64|Expr::from_unevaluated_kind(Kind::Pow(base,Expr::int(n)));
        let late=power(Expr::zoo(),1);
        for (infinity,negative) in [(Expr::oo(),false),(Expr::neg_oo(),true)] {
            for reversed in [false,true] {
                let args=if reversed {vec![infinity.clone(),late.clone()]} else {vec![late.clone(),infinity.clone()]};
                let result=mul(&args);
                assert!(!result.is_canonical());
                let Kind::Mul(coefficient,factors)=result.kind() else {panic!("late zoo must remain a factor")};
                assert!(coefficient.is_one());
                assert_eq!(factors.len(),2);
                assert!(factors.contains(&(if negative {Expr::neg_oo()} else {Expr::oo()},Expr::one())));
                assert!(factors.contains(&(Expr::zoo(),Expr::one())));
            }
        }
        let x=p("x");
        assert_eq!(mul(&[Expr::zoo(),late.clone(),x.clone()]),p("zoo*x"));
        assert_eq!(mul(&[power(Expr::zero(),-1),Expr::oo(),x]),p("zoo*x"));
        assert_eq!(mul(&[Expr::zoo(),late,Expr::i()]),Expr::zoo());
        // Do not reproduce SymPy's order-sensitive zoo-overwrites-zero bug.
        for reversed in [false,true] {
            let mut args=vec![Expr::zero(),power(Expr::oo(),1),Expr::zoo()];
            if reversed {args.reverse();}
            assert_eq!(mul(&args),Expr::nan());
        }
    }
    #[test]
    fn shared_active_product_dag_is_compressed_and_memoized() {
        let x=p("x");let mut value=x.clone();
        for _ in 0..80 {value=raw_mul(vec![value.clone(),value]);}
        let mut cache=ProductMemo::default();
        let result=mul_with_cache(&[Expr::int(2),value.clone()],0,&mut cache);
        let exponent=Expr::num(Number::from_bigint(BigInt::one()<<80));
        assert_eq!(result,Expr::from_kind(Kind::Mul(Number::Int(2),vec![(x,exponent)])));
        assert_eq!(cache.parts.len(),80);
        assert_eq!(cache.parts.values().map(|(_,parts)|parts.factors.len()).sum::<usize>(),80);
        assert!(!value.is_canonical());
    }

    #[test]
    fn expression_pivot_simplification_exposes_zero_complex_pivots() {
        let a=p("2*I*x+1/3");
        let opposite=neg(&a);
        let candidate=sub(&a,&div(&pow(&opposite,&Expr::int(2)),&a));
        assert_eq!(crate::matrix::linalg::pivot_simp(&candidate),Expr::zero());
    }
}

#[cfg(test)]
mod power_of_power_tests {
    use super::*;
    fn p(source: &str) -> Expr { crate::parser::parse_str(source).unwrap() }
    fn raw_pow(base: &str, exponent: &str) -> Expr {
        Expr::from_unevaluated_kind(Kind::Pow(p(base), p(exponent)))
    }
    #[test]
    fn exponents_below_one_in_magnitude_combine_like_sympy() {
        // Pow._eval_power: a real inner exponent a with |a| < 1 has no
        // branch correction, so (b**a)**e = b**(a*e).
        assert_eq!(pow(&raw_pow("t", "1/2"), &p("1/2")), p("t**(1/4)"));
        assert_eq!(pow(&raw_pow("t", "1/2"), &p("y")), p("t**(y/2)"));
        assert_eq!(pow(&raw_pow("t", "2/3"), &p("1/2")), p("t**(1/3)"));
        assert_eq!(pow(&raw_pow("t", "-1/2"), &p("1/3")), p("t**(-1/6)"));
        assert_eq!(pow(&raw_pow("x+1", "1/2"), &p("1/3")), p("(x+1)**(1/6)"));
        // |a| >= 1 keeps the nested power for a base of unknown sign.
        for (base, exponent) in [("t", "2"), ("t", "3/2"), ("t", "-1"), ("t", "y")] {
            let inner = raw_pow(base, exponent);
            assert!(matches!(pow(&inner, &p("1/2")).kind(), Kind::Pow(b, _) if matches!(b.kind(), Kind::Pow(..))));
        }
    }
}
