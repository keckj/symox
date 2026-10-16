//! Polynomials: sparse multivariate polynomials over `ZZ` and `QQ`, dense
//! univariate fast paths, gcd, square-free decomposition, factorisation,
//! rational function normalisation and root finding.
//!
//! Layers:
//!
//! * [`int`]: coefficient rings (`Int` with an `i64` fast path, `Rat`);
//! * [`mpoly`]: sparse multivariate polynomials, [`dense`] and [`gf`]:
//!   dense univariate polynomials over `ZZ`/`QQ` and `GF(p)`;
//! * algorithms: [`gcd`], [`sqf`], [`factor`] (Zassenhaus with Hensel
//!   lifting, multivariate lifting), [`resultant`];
//! * [`convert`]/[`poly`]: conversion between expressions and polynomials
//!   (automatic generator detection) and the high level `Poly` object;
//! * [`rational`], [`roots`]: `cancel`, `apart`, `collect`, `roots`,
//!   `nroots` working on expressions.

pub mod dense;
pub mod gf;
pub mod int;
pub mod mpoly;

pub use int::{Coeff, Field, Int, Rat};
pub use mpoly::MPoly;
pub mod factor_multi;
pub mod factor_uni;
pub mod gcd;
pub mod sqf;
pub mod convert;
pub mod polyobj;
pub mod api;
pub mod resultant;
pub mod roots;
pub mod modgcd;
pub mod ratfunc;
pub mod groebner;
pub mod systems;
pub mod parametric;
pub mod radical_cert;
pub mod fieldalg;
pub mod fraction;
pub mod algebraic;

pub mod f5b;
