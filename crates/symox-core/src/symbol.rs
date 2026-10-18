//! Symbols and their assumptions (`Symbol('x', positive=True)`).

use std::fmt;
use std::hash::{BuildHasher, Hash, Hasher};
use rustc_hash::FxHashSet;
use std::ops::Deref;
use std::sync::{Arc, Mutex, OnceLock};
use num_bigint::BigInt;

/// Assumption flags attached to a symbol. Only *known true* facts are
/// stored; implied facts are added by [`Assumptions::closure`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Assumptions(pub u32);

macro_rules! flags {
    ($($name:ident = $bit:expr, $py:literal;)*) => {
        impl Assumptions {
            $(pub const $name: Assumptions = Assumptions(1 << $bit);)*

            /// SymPy keyword names, in canonical order.
            pub const NAMES: &'static [(&'static str, Assumptions)] = &[$(($py, Assumptions::$name)),*];
        }
    };
}

flags! {
    REAL = 0, "real";
    POSITIVE = 1, "positive";
    NEGATIVE = 2, "negative";
    NONNEGATIVE = 3, "nonnegative";
    NONPOSITIVE = 4, "nonpositive";
    NONZERO = 5, "nonzero";
    ZERO = 6, "zero";
    INTEGER = 7, "integer";
    RATIONAL = 8, "rational";
    EVEN = 9, "even";
    ODD = 10, "odd";
    PRIME = 11, "prime";
    FINITE = 12, "finite";
    COMPLEX = 13, "complex";
    IMAGINARY = 14, "imaginary";
    IRRATIONAL = 15, "irrational";
    // negated facts
    NOT_REAL = 16, "!real";
    NOT_INTEGER = 17, "!integer";
    NOT_POSITIVE = 18, "!positive";
    NOT_NEGATIVE = 19, "!negative";
    NOT_ZERO = 20, "!zero";
    NOT_NONZERO = 21, "!nonzero";
    NOT_RATIONAL = 22, "!rational";
    NOT_NONNEGATIVE = 23, "!nonnegative";
    NOT_NONPOSITIVE = 24, "!nonpositive";
}

impl Assumptions {
    pub const NONE: Assumptions = Assumptions(0);

    #[inline]
    pub fn has(self, f: Assumptions) -> bool {
        self.0 & f.0 == f.0 && f.0 != 0
    }

    #[inline]
    pub fn with(self, f: Assumptions) -> Assumptions {
        Assumptions(self.0 | f.0)
    }

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Add every implied fact. Returns an error for contradictions.
    pub fn closure(self) -> Result<Assumptions, String> {
        let mut a = self;
        loop {
            let before = a;
            let imp = |a: &mut Assumptions, cond: Assumptions, then: Assumptions| {
                if a.has(cond) {
                    *a = a.with(then);
                }
            };
            imp(&mut a, Self::PRIME, Self::POSITIVE.with(Self::INTEGER));
            imp(&mut a, Self::EVEN, Self::INTEGER);
            imp(&mut a, Self::ODD, Self::INTEGER.with(Self::NONZERO));
            imp(&mut a, Self::INTEGER, Self::RATIONAL);
            imp(&mut a, Self::RATIONAL, Self::REAL);
            imp(&mut a, Self::NONZERO, Self::REAL.with(Self::NOT_ZERO));
            imp(&mut a, Self::NOT_REAL, Self::NOT_NONZERO.with(Self::NOT_RATIONAL).with(Self::NOT_INTEGER).with(Self::NOT_POSITIVE).with(Self::NOT_NEGATIVE).with(Self::NOT_NONNEGATIVE).with(Self::NOT_NONPOSITIVE));
            imp(&mut a, Self::NOT_NONNEGATIVE, Self::NOT_POSITIVE.with(Self::NOT_ZERO));
            imp(&mut a, Self::NOT_NONPOSITIVE, Self::NOT_NEGATIVE.with(Self::NOT_ZERO));
            imp(&mut a, Self::NOT_RATIONAL, Self::NOT_INTEGER);
            imp(&mut a, Self::NOT_INTEGER, Self::NOT_ZERO);
            imp(&mut a, Self::IRRATIONAL, Self::REAL.with(Self::NONZERO));
            imp(&mut a, Self::POSITIVE, Self::REAL.with(Self::NONZERO).with(Self::NONNEGATIVE).with(Self::NOT_NEGATIVE).with(Self::NOT_NONPOSITIVE));
            imp(&mut a, Self::NEGATIVE, Self::REAL.with(Self::NONZERO).with(Self::NONPOSITIVE).with(Self::NOT_POSITIVE).with(Self::NOT_NONNEGATIVE));
            imp(&mut a, Self::ZERO, Self::REAL.with(Self::NONNEGATIVE).with(Self::NONPOSITIVE).with(Self::EVEN).with(Self::NOT_POSITIVE).with(Self::NOT_NEGATIVE).with(Self::NOT_NONZERO));
            imp(&mut a, Self::NONNEGATIVE, Self::REAL.with(Self::NOT_NEGATIVE));
            imp(&mut a, Self::NONPOSITIVE, Self::REAL.with(Self::NOT_POSITIVE));
            imp(&mut a, Self::REAL, Self::COMPLEX.with(Self::FINITE));
            imp(&mut a, Self::COMPLEX, Self::FINITE);
            imp(&mut a, Self::IMAGINARY, Self::COMPLEX.with(Self::NOT_REAL).with(Self::NOT_ZERO).with(Self::NOT_NONZERO));
            if a.has(Self::REAL) && a.has(Self::NOT_ZERO) { a = a.with(Self::NONZERO); }
            if a.has(Self::REAL) && a.has(Self::NOT_NONZERO) { a = a.with(Self::ZERO); }
            if a.has(Self::NOT_ZERO) && a.has(Self::NOT_NONZERO) { a = a.with(Self::NOT_REAL); }
            if a.has(Self::REAL) && a.has(Self::NOT_RATIONAL) { a = a.with(Self::IRRATIONAL); }
            if a.has(Self::REAL) && a.has(Self::NOT_NONNEGATIVE) { a = a.with(Self::NEGATIVE); }
            if a.has(Self::REAL) && a.has(Self::NOT_NONPOSITIVE) { a = a.with(Self::POSITIVE); }
            if a.has(Self::REAL) && a.has(Self::NOT_NEGATIVE) { a = a.with(Self::NONNEGATIVE); }
            if a.has(Self::REAL) && a.has(Self::NOT_POSITIVE) { a = a.with(Self::NONPOSITIVE); }
            if a.has(Self::NOT_NONNEGATIVE) && a.has(Self::NOT_NEGATIVE) { a = a.with(Self::NOT_REAL); }
            if a.has(Self::NOT_NONPOSITIVE) && a.has(Self::NOT_POSITIVE) { a = a.with(Self::NOT_REAL); }
            if a.has(Self::NONNEGATIVE) && a.has(Self::NONZERO) {
                a = a.with(Self::POSITIVE);
            }
            if a.has(Self::NONPOSITIVE) && a.has(Self::NONZERO) {
                a = a.with(Self::NEGATIVE);
            }
            if a.has(Self::NONNEGATIVE) && a.has(Self::NONPOSITIVE) {
                a = a.with(Self::ZERO);
            }
            if a == before {
                break;
            }
        }
        let bad = |x: Assumptions, y: Assumptions| a.has(x) && a.has(y);
        if bad(Self::POSITIVE, Self::NOT_POSITIVE)
            || bad(Self::NEGATIVE, Self::NOT_NEGATIVE)
            || bad(Self::REAL, Self::NOT_REAL)
            || bad(Self::INTEGER, Self::NOT_INTEGER)
            || bad(Self::EVEN, Self::ODD)
            || bad(Self::ZERO, Self::NONZERO)
            || bad(Self::ZERO, Self::NOT_ZERO)
            || bad(Self::NONZERO, Self::NOT_NONZERO)
            || bad(Self::RATIONAL, Self::NOT_RATIONAL)
            || bad(Self::RATIONAL, Self::IRRATIONAL)
            || bad(Self::NONNEGATIVE, Self::NOT_NONNEGATIVE)
            || bad(Self::NONPOSITIVE, Self::NOT_NONPOSITIVE)
        {
            return Err("inconsistent assumptions".into());
        }
        Ok(a)
    }

    /// Build from SymPy keyword arguments (`positive=True`, `real=False`).
    /// Returns the given (non-closed) flags after checking consistency.
    pub fn from_kwargs<'a>(kw: impl IntoIterator<Item = (&'a str, bool)>) -> Result<Assumptions, String> {
        let mut a = Assumptions::NONE;
        for (k, v) in kw {
            let f = match (k, v) {
                ("commutative", true) => continue,
                ("real", false) | ("extended_real", false) => Self::NOT_REAL,
                ("integer", false) => Self::NOT_INTEGER,
                ("positive", false) => Self::NOT_POSITIVE,
                ("negative", false) => Self::NOT_NEGATIVE,
                ("zero", false) => Self::NOT_ZERO,
                ("nonzero", false) => Self::NOT_NONZERO,
                ("rational", false) => Self::NOT_RATIONAL,
                ("nonnegative", false) => Self::NOT_NONNEGATIVE,
                ("nonpositive", false) => Self::NOT_NONPOSITIVE,
                ("finite", false) | ("complex", false) | ("even", false)
                | ("odd", false) | ("prime", false) | ("irrational", false) | ("imaginary", false) => {
                    continue
                }
                ("extended_real", true) => Self::REAL,
                ("extended_positive", true) => Self::POSITIVE,
                ("extended_negative", true) => Self::NEGATIVE,
                ("extended_nonnegative", true) => Self::NONNEGATIVE,
                ("extended_nonpositive", true) => Self::NONPOSITIVE,
                (k, true) => match Self::NAMES.iter().find(|(n, _)| *n == k) {
                    Some((_, f)) => *f,
                    None => return Err(format!("unknown assumption {k}")),
                },
                (k, false) => return Err(format!("unknown assumption {k}")),
            };
            a = a.with(f);
        }
        a.closure()?;
        Ok(a)
    }

    /// SymPy keyword form (`["positive=True", "real=False"]`), sorted by
    /// name like SymPy's `srepr`.
    pub fn keywords(self) -> Vec<String> {
        let mut v: Vec<String> = Vec::new();
        for (n, f) in Self::NAMES {
            if self.has(*f) {
                match n.strip_prefix('!') {
                    Some(m) => v.push(format!("{m}=False")),
                    None => v.push(format!("{n}=True")),
                }
            }
        }
        v.sort();
        v
    }
}

/// Fresh indices keep a process namespace independent of foreign generators.
/// Explicit restoration reserves an index but never pulls this sequence into a
/// foreign namespace. Only restored/imported indices occupy the reservation set.
struct DummyState {
    pid: u32,
    next: BigInt,
    anonymous: BigInt,
    reserved: FxHashSet<BigInt>,
}
impl DummyState {
    fn seed(pid: u32) -> BigInt {
        let nanos=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default().as_nanos();
        let random=std::collections::hash_map::RandomState::new().build_hasher().finish();
        (BigInt::from(1u8)<<192usize) + (BigInt::from(nanos)<<96usize)
            + (BigInt::from(pid)<<64usize) + BigInt::from(random)
    }
    fn new() -> Self {
        let pid=std::process::id();
        Self{pid,next:Self::seed(pid),anonymous:BigInt::from(0),reserved:FxHashSet::default()}
    }
    fn fresh(&mut self) -> BigInt {
        let pid=std::process::id();
        if pid!=self.pid {
            // A fork copies global state: the child must start its own namespace.
            self.pid=pid;self.next=Self::seed(pid);self.anonymous=BigInt::from(0);
        }
        loop {
            let result=self.next.clone();self.next+=1;
            if !self.reserved.contains(&result) { return result; }
        }
    }
}
fn with_dummy_state<T>(operation:impl FnOnce(&mut DummyState)->T) -> T {
    static STATE:OnceLock<Mutex<DummyState>>=OnceLock::new();
    operation(&mut STATE.get_or_init(||Mutex::new(DummyState::new())).lock().unwrap())
}
fn next_dummy_index(reserve: Option<&BigInt>) -> BigInt {
    with_dummy_state(|state| match reserve {
        Some(index)=>{state.reserved.insert(index.clone());index.clone()},
        None=>state.fresh(),
    })
}

/// A symbol: a name plus assumptions. Dereferences to its name.
///
/// Identity uses the name, closed assumptions, and optional Dummy index.
/// The assumptions originally provided remain available for round trips.
#[derive(Clone, Debug)]
pub struct Symbol {
    name: Arc<str>,
    given: Assumptions,
    known: Assumptions,
    given_facts: crate::assumption_facts::FactSet,
    facts: crate::assumption_facts::FactSet,
    // Preserve historical hashes where the legacy projection is lossless.
    // New distinctions contribute their complete canonical facts.
    hash_extra: Option<crate::assumption_facts::FactSet>,
    /// `Some(exclude)` for pattern wildcards (`Wild('a', exclude=[x])`).
    wild: Option<Arc<Vec<crate::expr::Expr>>>,
    /// Explicit identity of a Dummy; ordinary symbols never infer this from a name.
    dummy_index: Option<Arc<BigInt>>,
}

impl Symbol {
    pub fn new(name: &str) -> Symbol {
        Self::with_assumptions(name,Assumptions::NONE)
    }

    /// `given` are the user's assumptions; their closure is computed here.
    pub fn with_assumptions(name: &str, given: Assumptions) -> Symbol {
        Self::with_facts(name,crate::assumption_facts::FactSet::from_legacy(given)).expect("consistent Symbol facts")
    }

    pub fn with_facts(name:&str,given_facts:crate::assumption_facts::FactSet)->Result<Self,String>{
        let facts=given_facts.closure()?;
        let known=facts.legacy();
        let hash_extra=(crate::assumption_facts::FactSet::from_legacy(known).closure()?!=facts).then_some(facts);
        Ok(Self{name:Arc::from(name),given:given_facts.legacy(),known,given_facts,facts,hash_extra,wild:None,dummy_index:None})
    }
    pub fn facts(&self)->crate::assumption_facts::FactSet{self.facts}
    pub fn given_facts(&self)->crate::assumption_facts::FactSet{self.given_facts}
    pub(crate) fn cmp_basic(&self,other:&Self)->std::cmp::Ordering{
        (1+self.facts.len()+usize::from(self.is_dummy())).cmp(&(1+other.facts.len()+usize::from(other.is_dummy())))
            .then_with(||self.name.cmp(&other.name))
            .then_with(||self.facts.iter().cmp(other.facts.iter()))
            .then_with(||self.dummy_index.cmp(&other.dummy_index))
            .then_with(||self.wild.cmp(&other.wild))
    }
    pub fn dummy_with_facts(name:&str,given:crate::assumption_facts::FactSet,index:Option<BigInt>)->Result<Self,String>{
        let mut symbol=Self::with_facts(name,given)?;
        symbol.dummy_index=Some(Arc::new(next_dummy_index(index.as_ref())));Ok(symbol)
    }
    pub fn fresh_dummy_like(&self,name:&str)->Self{
        Self::dummy_with_facts(name,self.given_facts,None).expect("existing Symbol facts")
    }
    pub fn anonymous_dummy_with_facts(given:crate::assumption_facts::FactSet)->Result<Self,String>{
        let mut symbol=Self::with_facts("",given)?;
        let(index,name)=with_dummy_state(|state|{let index=state.fresh();let name=format!("Dummy_{}",state.anonymous);state.anonymous+=1;(index,name)});
        symbol.name=Arc::from(name);symbol.dummy_index=Some(Arc::new(index));Ok(symbol)
    }

    /// A fresh Dummy, whose identity is independent of its displayed name.
    pub fn dummy(name: &str, given: Assumptions) -> Symbol {
        let index = next_dummy_index(None);
        Self::dummy_registered(name, given, index)
    }

    /// Restore a Dummy identity, reserving its index against fresh construction.
    pub fn dummy_with_index(name: &str, given: Assumptions, index: BigInt) -> Symbol {
        next_dummy_index(Some(&index));
        Self::dummy_registered(name, given, index)
    }

    pub fn anonymous_dummy(given: Assumptions) -> Symbol {
        let (index,name)=with_dummy_state(|state| {
            let index=state.fresh();
            let name=format!("Dummy_{}",state.anonymous);state.anonymous+=1;
            (index,name)
        });
        Self::dummy_registered(&name, given, index)
    }

    fn dummy_registered(name: &str, given: Assumptions, index: BigInt) -> Symbol {
        let mut symbol = Self::with_assumptions(name, given);
        symbol.dummy_index = Some(Arc::new(index));
        symbol
    }

    #[inline]
    pub fn dummy_index(&self) -> Option<&BigInt> {
        self.dummy_index.as_deref()
    }

    #[inline]
    pub fn is_dummy(&self) -> bool {
        self.dummy_index.is_some()
    }

    /// A pattern wildcard that must not match expressions containing any
    /// of `exclude`.
    pub fn wild(name: &str, exclude: Vec<crate::expr::Expr>) -> Symbol {
        let mut symbol=Self::new(name);symbol.wild=Some(Arc::new(exclude));symbol
    }

    pub fn wild_with_facts(name: &str, exclude: Vec<crate::expr::Expr>, given: crate::assumption_facts::FactSet) -> Result<Symbol, String> {
        let mut symbol = Self::with_facts(name, given)?;
        symbol.wild = Some(Arc::new(exclude));
        Ok(symbol)
    }

    #[inline]
    pub fn is_wild(&self) -> bool {
        self.wild.is_some()
    }

    /// Exclusion list of a wildcard.
    pub fn exclude(&self) -> &[crate::expr::Expr] {
        match &self.wild {
            Some(v) => v,
            None => &[],
        }
    }

    #[inline]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// All known facts (closure of the given assumptions).
    #[inline]
    pub fn assumptions(&self) -> Assumptions {
        self.known
    }

    /// The assumptions as given by the user.
    #[inline]
    pub fn given(&self) -> Assumptions {
        self.given
    }
}

impl PartialEq for Symbol {
    fn eq(&self, o: &Symbol) -> bool {
        self.name == o.name
            && self.facts == o.facts
            && self.dummy_index == o.dummy_index
            && match (&self.wild, &o.wild) {
                (None, None) => true,
                (Some(a), Some(b)) => a == b,
                _ => false,
            }
    }
}
impl Eq for Symbol {}

impl PartialOrd for Symbol {
    fn partial_cmp(&self, o: &Symbol) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Symbol {
    fn cmp(&self, o: &Symbol) -> std::cmp::Ordering {
        self.name
            .cmp(&o.name)
            .then(self.known.cmp(&o.known))
            .then(self.facts.cmp(&o.facts))
            .then(self.dummy_index.cmp(&o.dummy_index))
            .then_with(|| self.wild.cmp(&o.wild))
    }
}

impl Hash for Symbol {
    fn hash<H: Hasher>(&self, h: &mut H) {
        self.name.hash(h);
        if !self.known.is_empty(){self.known.hash(h);}
        if let Some(extra)=self.hash_extra{0xfac7u16.hash(h);extra.hash(h);}
        if let Some(index) = &self.dummy_index {
            0xd00du16.hash(h);
            index.hash(h);
        }
        if let Some(w) = &self.wild {
            0xa11du16.hash(h);
            w.len().hash(h);
        }
    }
}

impl Deref for Symbol {
    type Target = str;
    fn deref(&self) -> &str {
        &self.name
    }
}

impl AsRef<str> for Symbol {
    fn as_ref(&self) -> &str {
        &self.name
    }
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.name)
    }
}

impl PartialEq<str> for Symbol {
    fn eq(&self, o: &str) -> bool {
        &*self.name == o
    }
}

impl PartialEq<&str> for Symbol {
    fn eq(&self, o: &&str) -> bool {
        &*self.name == *o
    }
}

#[cfg(test)]
mod tests {
    use super::Assumptions as A;
    #[test]
    fn real_nonzero_and_complex_nonzero_are_distinct() {
        let nonzero=A::from_kwargs([("nonzero",true)]).unwrap().closure().unwrap();
        assert!(nonzero.has(A::REAL));assert!(nonzero.has(A::NOT_ZERO));
        let not_zero=A::from_kwargs([("zero",false)]).unwrap().closure().unwrap();
        assert!(!not_zero.has(A::REAL));assert!(!not_zero.has(A::NONZERO));
        let not_nonzero=A::from_kwargs([("nonzero",false)]).unwrap().closure().unwrap();
        assert!(!not_nonzero.has(A::ZERO));assert!(!not_nonzero.has(A::REAL));
        let imaginary=A::from_kwargs([("imaginary",true)]).unwrap().closure().unwrap();
        assert!(imaginary.has(A::NOT_ZERO));assert!(imaginary.has(A::NOT_NONZERO));assert!(!imaginary.has(A::NONZERO));
        assert!(A::from_kwargs([("imaginary",true),("nonzero",true)]).is_err());
        assert!(A::from_kwargs([("zero",true),("zero",false)]).is_err());
        assert!(A::from_kwargs([("rational",true),("rational",false)]).is_err());
    }
}

#[cfg(test)]
mod dummy_tests {
    use super::*;
    use crate::{expr::Expr, ops, parser, printer, traverse};
    use std::collections::{HashSet, BTreeSet};

    #[test]
    fn dummy_identity_is_explicit_and_independent_of_display_name() {
        let a=Expr::dummy("x",Assumptions::POSITIVE);
        let b=Expr::dummy("x",Assumptions::POSITIVE);
        let index=a.as_sym().unwrap().dummy_index().unwrap().clone();
        let restored=Expr::dummy_with_index("x",Assumptions::POSITIVE,index.clone());
        assert_eq!(a,restored);assert_ne!(a,b);
        assert_ne!(a,Expr::symbol_with("x",Assumptions::POSITIVE));
        assert_ne!(a,Expr::dummy_with_index("y",Assumptions::POSITIVE,index.clone()));
        assert_ne!(a,Expr::dummy_with_index("x",Assumptions::NONE,index));
        assert_eq!(HashSet::from([a.clone(),b.clone(),restored.clone()]).len(),2);
        assert_eq!(BTreeSet::from([a.clone(),b.clone(),restored]).len(),2);
        assert_eq!(a.class_name(),"Dummy");
        assert_eq!(a.to_string(),"_x");
        assert_eq!(Expr::symbol("_x").class_name(),"Symbol");
    }

    #[test]
    fn dummy_arithmetic_has_substitution_and_derivatives_do_not_capture_names() {
        let a=Expr::dummy("x",Assumptions::NONE);let b=Expr::dummy("x",Assumptions::NONE);
        let x=Expr::symbol("x");let underscore=Expr::symbol("_x");
        let e=ops::add(&[a.clone(),b.clone(),x.clone(),underscore.clone()]);
        assert_eq!(traverse::free_symbols(&e).len(),4);
        assert_eq!(crate::diff::diff(&e,&a),Expr::one());
        let replaced=crate::subs::subs1(&e,&a,&Expr::int(2));
        assert!(!traverse::has(&replaced,&a));assert!(traverse::has(&replaced,&b));
        assert!(traverse::has(&replaced,&x));assert!(traverse::has(&replaced,&underscore));
        let f=crate::functions::undefined("f",vec![x.powi(2)]);
        let first=crate::diff::diff(&f,&x);let second=crate::diff::diff(&f,&x);
        assert_eq!(first,second); // held Subs alpha identity ignores fresh bound indices
        assert_eq!(traverse::free_symbols(&first),vec![x]);
    }

    #[test]
    fn dummy_srepr_roundtrips_arbitrary_indices_and_assumptions() {
        for index in [BigInt::from(-7),BigInt::from(1)<<180] {
            let e=Expr::dummy_with_index("a'\\b",Assumptions::REAL,index);
            let encoded=printer::srepr(&e);
            assert!(encoded.starts_with("Dummy("));assert_eq!(parser::parse_str(&encoded).unwrap(),e);
            let pair=ops::add(&[e.clone(),Expr::dummy("a'\\b",Assumptions::REAL)]);
            assert_eq!(parser::parse_str(&printer::srepr(&pair)).unwrap(),pair);
        }
    }

    #[test]
    fn explicit_restoration_reserves_indices_for_future_fresh_atoms() {
        let high:BigInt=BigInt::from(1)<<200;
        let restored=Symbol::dummy_with_index("x",Assumptions::NONE,high.clone());
        let fresh=Symbol::dummy("x",Assumptions::NONE);
        assert_ne!(fresh.dummy_index().unwrap(),&high);assert_ne!(restored,fresh);
    }

    #[test]
    fn fresh_namespace_skips_explicit_indices_and_deduplicates_reservations() {
        with_dummy_state(|state|{
            let index=state.next.clone();
            let before=state.reserved.len();
            state.reserved.insert(index.clone());state.reserved.insert(index.clone());
            assert_eq!(state.reserved.len(),before+1);
            assert_ne!(state.fresh(),index);
        });
    }

    #[test]
    fn dummy_cxx_names_contain_the_exact_identity() {
        let d=Expr::dummy_with_index("x",Assumptions::NONE,BigInt::from(57));
        assert_eq!(crate::printing::code::cxxcode(&d).unwrap(),"x_57");
        assert_eq!(crate::printing::code::cxxcode(&Expr::symbol("_x")).unwrap(),"_x");
    }
}
