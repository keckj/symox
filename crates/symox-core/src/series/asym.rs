//! Auxiliary functions with an asymptotic power series in `1/z` at an
//! infinite argument, the counterparts of SymPy's `_erfs` and `_eis`
//! helpers. The Gruntz preparation writes the special functions whose
//! argument tends to `±oo` as an elementary part (exponentials,
//! oscillating `sin`/`cos`, constants) times these functions, so that the
//! mrv set sees the exponentials and the series engine expands the rest
//! (`_eval_rewrite_as_tractable` and `_eval_aseries` of erf, erfc, erfi,
//! Ei, li, expint, uppergamma, Si, Ci, Shi, Chi, fresnels and fresnelc).
//!
//! The functions are internal undefined heads: they never leave the limit
//! engine, and the series engine expands them only when the argument
//! tends to an infinite value.

use num_bigint::BigInt;

use crate::expr::{Expr, Func, Kind, UndefinedFunction};
use crate::number::Number;
use crate::ops::{add, mul2, pow};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Asym {
    /// `exp(z**2)*erfc(z)`: `Σ (-1)**j (2j-1)!!/2**j / sqrt(pi) * z**-(2j+1)`
    Erfs,
    /// `exp(-z**2)*erfi(z)` (real `z`): `Σ (2j-1)!!/2**j / sqrt(pi) * z**-(2j+1)`
    Erfis,
    /// `exp(-z)*Ei(z)`: `Σ j! z**-(j+1)`
    Eis,
    /// `exp(z)*expint(nu, z)`: `Σ (-1)**j rf(nu, j) z**-(j+1)` (two arguments)
    Expints,
    /// `f` of `Si(z) = pi/2 - f cos(z) - g sin(z)`: `Σ (-1)**j (2j)! z**-(2j+1)`
    SiF,
    /// `g`: `Σ (-1)**j (2j+1)! z**-(2j+2)`
    SiG,
    /// `f` of `S(z) = 1/2 - f cos(pi z**2/2) - g sin(pi z**2/2)`:
    /// `Σ (-1)**m (4m-1)!!/pi**(2m+1) z**-(4m+1)`
    FresF,
    /// `g`: `Σ (-1)**m (4m+1)!!/pi**(2m+2) z**-(4m+3)`
    FresG,
}

const NAMES: [(Asym, &str); 8] = [
    (Asym::Erfs, "_symox_erfs"),
    (Asym::Erfis, "_symox_erfis"),
    (Asym::Eis, "_symox_eis"),
    (Asym::Expints, "_symox_expints"),
    (Asym::SiF, "_symox_sif"),
    (Asym::SiG, "_symox_sig"),
    (Asym::FresF, "_symox_fresf"),
    (Asym::FresG, "_symox_fresg"),
];

impl Asym {
    fn name(self) -> &'static str {
        NAMES.iter().find(|(k, _)| *k == self).unwrap().1
    }

    pub fn of(f: &Func) -> Option<Asym> {
        match f {
            Func::Undefined(h) if h.internal => NAMES.iter().find(|(_, n)| *n == h.name.as_ref()).map(|(k, _)| *k),
            _ => None,
        }
    }

    /// The node `self(args)`.
    pub fn apply(self, args: Vec<Expr>) -> Expr {
        let mut head = UndefinedFunction::new(self.name(), crate::symbol::Assumptions::NONE, Some(vec![args.len()]));
        head.internal = true;
        Expr::from_kind(Kind::Func(Func::Undefined(head), args))
    }

    /// Coefficient of `z**-k` (`k >= 1`) in the asymptotic series; `args`
    /// are the arguments of the node (the order `nu` first for `Expints`).
    pub fn coeff(self, args: &[Expr], k: i64) -> Expr {
        let zero = Expr::zero();
        let int = |b: BigInt| Expr::num(Number::from_bigint(b));
        let sgn = |j: i64| if j % 2 == 0 { 1 } else { -1 };
        let rsqrtpi = pow(&Expr::pi(), &Expr::rational(-1, 2));
        match self {
            Asym::Erfs | Asym::Erfis => {
                if k % 2 == 0 {
                    return zero;
                }
                let j = (k - 1) / 2;
                let num = int(odd_factorial(j));
                let s = if self == Asym::Erfs { sgn(j) } else { 1 };
                let c = mul2(&Expr::int(s), &crate::ops::div(&num, &pow(&Expr::int(2), &Expr::int(j))));
                mul2(&c, &rsqrtpi)
            }
            Asym::Eis => int(factorial(k - 1)),
            Asym::Expints => {
                let nu = &args[0];
                let j = k - 1;
                let mut acc = Expr::int(sgn(j));
                for i in 0..j {
                    acc = mul2(&acc, &add(&[nu.clone(), Expr::int(i)]));
                }
                crate::expand::expand(&acc)
            }
            Asym::SiF => {
                if k % 2 == 0 {
                    return zero;
                }
                let j = (k - 1) / 2;
                mul2(&Expr::int(sgn(j)), &int(factorial(2 * j)))
            }
            Asym::SiG => {
                if k % 2 == 1 {
                    return zero;
                }
                let j = (k - 2) / 2;
                mul2(&Expr::int(sgn(j)), &int(factorial(2 * j + 1)))
            }
            Asym::FresF | Asym::FresG => {
                let (r, d) = if self == Asym::FresF { (1, 0) } else { (3, 1) };
                if k % 4 != r {
                    return zero;
                }
                let m = (k - r) / 4;
                // (4m-1)!! = odd_factorial(2m) and (4m+1)!! = odd_factorial(2m+1)
                let num = int(odd_factorial(2 * m + d));
                let p = pow(&Expr::pi(), &Expr::int(-(2 * m + 1 + d)));
                mul2(&Expr::int(sgn(m)), &mul2(&num, &p))
            }
        }
    }
}

fn factorial(n: i64) -> BigInt {
    (1..=n.max(0)).fold(BigInt::from(1), |a, i| a * i)
}

/// `(2j-1)!! = 1*3*...*(2j-1)` (`1` for `j = 0`).
fn odd_factorial(j: i64) -> BigInt {
    (1..=j.max(0)).fold(BigInt::from(1), |a, i| a * (2 * i - 1))
}

/// Does `e` contain one of the auxiliary functions?
pub fn contains_asym(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(f, _) if Asym::of(f).is_some() => true,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) | Kind::Const(_) => false,
        _ => e.args().iter().any(contains_asym),
    }
}
