//! Symbols and their assumptions (`Symbol('x', positive=True)`).

use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::sync::Arc;

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
            imp(&mut a, Self::IRRATIONAL, Self::REAL.with(Self::NONZERO));
            imp(&mut a, Self::POSITIVE, Self::REAL.with(Self::NONZERO).with(Self::NONNEGATIVE).with(Self::NOT_NEGATIVE));
            imp(&mut a, Self::NEGATIVE, Self::REAL.with(Self::NONZERO).with(Self::NONPOSITIVE).with(Self::NOT_POSITIVE));
            imp(&mut a, Self::ZERO, Self::REAL.with(Self::NONNEGATIVE).with(Self::NONPOSITIVE).with(Self::EVEN).with(Self::NOT_POSITIVE).with(Self::NOT_NEGATIVE));
            imp(&mut a, Self::NONNEGATIVE, Self::REAL.with(Self::NOT_NEGATIVE));
            imp(&mut a, Self::NONPOSITIVE, Self::REAL.with(Self::NOT_POSITIVE));
            imp(&mut a, Self::REAL, Self::COMPLEX.with(Self::FINITE));
            imp(&mut a, Self::IMAGINARY, Self::COMPLEX.with(Self::NOT_REAL).with(Self::NONZERO));
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
            || bad(Self::RATIONAL, Self::IRRATIONAL)
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
                ("zero", false) => Self::NONZERO,
                ("nonzero", false) => Self::ZERO,
                ("nonnegative", false) => Self::NEGATIVE,
                ("nonpositive", false) => Self::POSITIVE,
                ("finite", false) | ("complex", false) | ("rational", false) | ("even", false)
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

/// A symbol: a name plus assumptions. Dereferences to its name.
///
/// Identity (equality, hashing) uses the name and the assumptions *as
/// given* (like SymPy's `assumptions0`); queries use their closure.
#[derive(Clone, Debug)]
pub struct Symbol {
    name: Arc<str>,
    given: Assumptions,
    known: Assumptions,
}

impl Symbol {
    pub fn new(name: &str) -> Symbol {
        Symbol { name: Arc::from(name), given: Assumptions::NONE, known: Assumptions::NONE }
    }

    /// `given` are the user's assumptions; their closure is computed here.
    pub fn with_assumptions(name: &str, given: Assumptions) -> Symbol {
        let known = given.closure().unwrap_or(given);
        Symbol { name: Arc::from(name), given, known }
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
        self.name == o.name && self.known == o.known
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
        self.name.cmp(&o.name).then(self.known.cmp(&o.known))
    }
}

impl Hash for Symbol {
    fn hash<H: Hasher>(&self, h: &mut H) {
        self.name.hash(h);
        if !self.known.is_empty() {
            self.known.hash(h);
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
