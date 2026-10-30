//! The compiled program with SciPy's float64 gamma (`crate::cephes_gamma`)
//! in place of `gamma_real`, for the array calls where SymPy's `lambdify`
//! calls `scipy.special.gamma`.

use super::{Compiled, Ins};
use crate::special::numeric::gamma::gamma_real;

impl Compiled {
    /// The program with `scipy.special.gamma` for gamma, or None when it
    /// does not call gamma.
    pub fn with_scipy_gamma(&self) -> Option<Compiled> {
        let real = gamma_real as fn(f64) -> f64 as usize;
        let mut found = false;
        let code = self
            .code
            .iter()
            .map(|i| match i {
                Ins::Call1(f, a) if *f as usize == real => {
                    found = true;
                    Ins::Call1(crate::cephes_gamma::gamma, *a)
                }
                other => other.clone(),
            })
            .collect();
        found.then(|| Compiled { code, outputs: self.outputs.clone(), n_args: self.n_args })
    }
}

#[cfg(test)]
mod tests {
    use crate::lambdify::Compiled;
    use crate::prelude::*;

    #[test]
    fn gamma_is_replaced() {
        let x = Expr::symbol("x");
        let e = add(&[func(Func::Gamma, vec![x.clone()]), x.clone()]);
        let c = Compiled::new(&[x.clone()], &[e]).unwrap();
        let s = c.with_scipy_gamma().unwrap();
        let mut out = [0.0];
        s.call(&[0.5], &mut out);
        assert_eq!(out[0], crate::cephes_gamma::gamma(0.5) + 0.5);
        let c = Compiled::new(&[x.clone()], &[sin(&x)]).unwrap();
        assert!(c.with_scipy_gamma().is_none());
    }
}
