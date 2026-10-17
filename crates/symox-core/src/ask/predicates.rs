//! Predicate keys of the new assumptions system (SymPy's `Q`).

use std::collections::HashSet;
use std::sync::Mutex;

use super::facts_data::PREDICATE_TABLE;
use crate::expr::RelOp;

/// A predicate of the `Q` namespace, or a user `Predicate(name)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Pred {
    Algebraic,
    Antihermitian,
    Commutative,
    Complex,
    ComplexElements,
    Composite,
    Diagonal,
    Eq,
    Even,
    ExtendedNegative,
    ExtendedNonnegative,
    ExtendedNonpositive,
    ExtendedNonzero,
    ExtendedPositive,
    ExtendedReal,
    Finite,
    Fullrank,
    Ge,
    Gt,
    Hermitian,
    Imaginary,
    Infinite,
    Integer,
    IntegerElements,
    Invertible,
    Irrational,
    IsTrue,
    Le,
    LowerTriangular,
    Lt,
    Ne,
    Negative,
    NegativeInfinite,
    Noninteger,
    Nonnegative,
    Nonpositive,
    Nonzero,
    Normal,
    Odd,
    Orthogonal,
    Positive,
    PositiveDefinite,
    PositiveInfinite,
    Prime,
    Rational,
    Real,
    RealElements,
    Singular,
    Square,
    Symmetric,
    Transcendental,
    Triangular,
    UnitTriangular,
    Unitary,
    UpperTriangular,
    Zero,
    /// `Predicate("name")` (SymPy's `UndefinedPredicate`), interned.
    Undefined(&'static str),
}

static INTERNED: Mutex<Option<HashSet<&'static str>>> = Mutex::new(None);

impl Pred {
    /// Every defined predicate, in `Q` attribute order.
    pub fn all() -> impl Iterator<Item = Pred> {
        PREDICATE_TABLE.iter().map(|(p, _, _)| *p)
    }
    fn row(self) -> Option<&'static (Pred, &'static str, &'static str)> {
        PREDICATE_TABLE.iter().find(|(p, _, _)| *p == self)
    }
    /// SymPy name (`positive`, `eq`, ...) or the undefined predicate name.
    pub fn name(self) -> &'static str {
        match self {
            Pred::Undefined(n) => n,
            _ => self.row().map(|r| r.1).unwrap_or("?"),
        }
    }
    /// SymPy class name (`PositivePredicate`, `UndefinedPredicate`, ...).
    pub fn class_name(self) -> &'static str {
        match self {
            Pred::Undefined(_) => "UndefinedPredicate",
            _ => self.row().map(|r| r.2).unwrap_or("Predicate"),
        }
    }
    /// A defined `Q` predicate by name.
    pub fn from_name(name: &str) -> Option<Pred> {
        PREDICATE_TABLE
            .iter()
            .find(|(_, n, _)| *n == name)
            .map(|r| r.0)
    }
    /// `Predicate(name)`: never equal to a `Q` predicate.
    pub fn undefined(name: &str) -> Pred {
        let mut guard = INTERNED.lock().unwrap_or_else(|e| e.into_inner());
        let set = guard.get_or_insert_with(HashSet::new);
        if let Some(s) = set.get(name) {
            return Pred::Undefined(s);
        }
        let leaked: &'static str = Box::leak(name.to_string().into_boxed_str());
        set.insert(leaked);
        Pred::Undefined(leaked)
    }
    pub fn is_undefined(self) -> bool {
        matches!(self, Pred::Undefined(_))
    }
    pub fn is_binrel(self) -> bool {
        self.relop().is_some()
    }
    pub fn relop(self) -> Option<RelOp> {
        Some(match self {
            Pred::Eq => RelOp::Eq,
            Pred::Ne => RelOp::Ne,
            Pred::Gt => RelOp::Gt,
            Pred::Ge => RelOp::Ge,
            Pred::Lt => RelOp::Lt,
            Pred::Le => RelOp::Le,
            _ => return None,
        })
    }
    pub fn from_relop(op: RelOp) -> Pred {
        match op {
            RelOp::Eq => Pred::Eq,
            RelOp::Ne => Pred::Ne,
            RelOp::Gt => Pred::Gt,
            RelOp::Ge => Pred::Ge,
            RelOp::Lt => Pred::Lt,
            RelOp::Le => Pred::Le,
        }
    }
    /// Predicate of the reversed relation (`Q.gt` -> `Q.lt`).
    pub fn reversed(self) -> Pred {
        match self {
            Pred::Gt => Pred::Lt,
            Pred::Lt => Pred::Gt,
            Pred::Ge => Pred::Le,
            Pred::Le => Pred::Ge,
            p => p,
        }
    }
    /// Predicate of the negated relation (`Q.gt` -> `Q.le`).
    pub fn negated(self) -> Option<Pred> {
        Some(match self {
            Pred::Eq => Pred::Ne,
            Pred::Ne => Pred::Eq,
            Pred::Gt => Pred::Le,
            Pred::Le => Pred::Gt,
            Pred::Ge => Pred::Lt,
            Pred::Lt => Pred::Ge,
            _ => return None,
        })
    }
    /// Matrix predicates (their facts are instantiated for matrix arguments).
    pub fn is_matrix(self) -> bool {
        matches!(
            self,
            Pred::Square
                | Pred::Symmetric
                | Pred::Invertible
                | Pred::Orthogonal
                | Pred::Unitary
                | Pred::PositiveDefinite
                | Pred::UpperTriangular
                | Pred::LowerTriangular
                | Pred::Diagonal
                | Pred::Fullrank
                | Pred::Singular
                | Pred::Normal
                | Pred::Triangular
                | Pred::UnitTriangular
                | Pred::IntegerElements
                | Pred::RealElements
                | Pred::ComplexElements
        )
    }
    /// Name of the matching old-style assumption (`expr.is_<name>`).
    pub fn old_assumption(self) -> Option<&'static str> {
        Some(match self {
            Pred::Algebraic => "algebraic",
            Pred::Antihermitian => "antihermitian",
            Pred::Commutative => "commutative",
            Pred::Complex => "complex",
            Pred::Composite => "composite",
            Pred::Even => "even",
            Pred::ExtendedNegative => "extended_negative",
            Pred::ExtendedNonnegative => "extended_nonnegative",
            Pred::ExtendedNonpositive => "extended_nonpositive",
            Pred::ExtendedNonzero => "extended_nonzero",
            Pred::ExtendedPositive => "extended_positive",
            Pred::ExtendedReal => "extended_real",
            Pred::Finite => "finite",
            Pred::Hermitian => "hermitian",
            Pred::Imaginary => "imaginary",
            Pred::Infinite => "infinite",
            Pred::Integer => "integer",
            Pred::Irrational => "irrational",
            Pred::Negative => "negative",
            Pred::Noninteger => "noninteger",
            Pred::Nonnegative => "nonnegative",
            Pred::Nonpositive => "nonpositive",
            Pred::Nonzero => "nonzero",
            Pred::Odd => "odd",
            Pred::Positive => "positive",
            Pred::Prime => "prime",
            Pred::Rational => "rational",
            Pred::Real => "real",
            Pred::Transcendental => "transcendental",
            Pred::Zero => "zero",
            _ => return None,
        })
    }
}
