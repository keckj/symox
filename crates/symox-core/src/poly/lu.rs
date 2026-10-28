//! Exact field in-place LU (QQ, GF(p), QQ_I, QQ<a>), with the public
//! routine's publication order.
use super::Field;
#[cfg(test)]
use super::Rat;

pub enum Event<'a, C = super::Rat> {
    Swap(usize, usize),
    Set(usize, usize, &'a C),
}

/// Return None only for a ragged shape, before publishing any mutation.
/// Zero columns are skipped at their original diagonal index: this is LU,
/// not an echelon reduction or a rank-revealing pivoting substitution.
pub fn in_place<C: Field, E>(
    a: &mut [Vec<C>],
    check: &mut impl FnMut() -> Result<(), E>,
    publish: &mut impl FnMut(Event<'_, C>) -> Result<(), E>,
) -> Result<Option<Vec<(usize, usize)>>, E> {
    check()?;
    let m = a.len();
    let Some(first) = a.first() else { return Ok(Some(Vec::new())) };
    let n = first.len();
    if a.iter().any(|row| row.len() != n) { return Ok(None); }
    let mut swaps = Vec::new();
    for i in 0..m.min(n) {
        check()?;
        if a[i][i].is_zero() {
            let replacement = ((i + 1)..m).find(|&ip| !a[ip][i].is_zero());
            let Some(ip) = replacement else { continue };
            swaps.push((i, ip));
            publish(Event::Swap(i, ip))?;
            a.swap(i, ip);
        }
        // The pivot is nonzero whenever a row below is eliminated; its exact
        // inverse is a pure function of it, so l_ji = a[j][i] / a[i][i].
        let inverse = if i + 1 < m && !a[i][i].is_zero() { Some(a[i][i].inv()) } else { None };
        for j in (i + 1)..m {
            check()?;
            let multiplier = match &inverse { Some(inverse) => a[j][i].mul(inverse), None => a[j][i].mul(&a[i][i].inv()) };
            publish(Event::Set(j, i, &multiplier))?;
            a[j][i] = multiplier.clone();
            for k in (i + 1)..n {
                check()?;
                let value = a[j][k].sub(&multiplier.mul(&a[i][k]));
                publish(Event::Set(j, k, &value))?;
                a[j][k] = value;
            }
        }
    }
    check()?;
    Ok(Some(swaps))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Int;
    fn q(n: i64, d: i64) -> Rat { Rat::new(Int::from(n), Int::from(d)) }

    #[test]
    fn rational_multipliers_and_row_swaps_follow_diagonal_pivots() {
        let mut a = vec![vec![q(0,1),q(2,3),q(1,4)],vec![q(2,5),q(3,7),q(8,9)],vec![q(1,3),q(0,1),q(4,11)]];
        let mut events = Vec::new();
        let swaps = in_place(&mut a, &mut || Ok::<_,()>(()), &mut |event| {
            events.push(match event { Event::Swap(i,j) => (i,j,None),Event::Set(i,j,v) => (i,j,Some(v.clone())) });
            Ok(())
        }).unwrap().unwrap();
        assert_eq!(swaps, vec![(0,1)]);
        assert_eq!(a[0],vec![q(2,5),q(3,7),q(8,9)]);
        assert_eq!(a[1],vec![q(0,1),q(2,3),q(1,4)]);
        assert_eq!(a[2],vec![q(5,6),q(-15,28),q(-8089,33264)]);
        assert_eq!(events[0],(0,1,None));
        assert_eq!(events[1],(1,0,Some(q(0,1))));
        assert_eq!(events.len(),9);
    }

    #[test]
    fn finite_field_lu_matches_sympy_gf7_factors() {
        use crate::poly::f5b::{Fp, PrimeScope};
        let _scope = PrimeScope::new(7).unwrap();
        let f = |n: i64| Fp::from_rat(&q(n, 1), 7).unwrap();
        let mut a = vec![vec![f(0),f(2),f(3)],vec![f(4),f(5),f(6)],vec![f(1),f(1),f(1)]];
        let swaps = in_place(&mut a, &mut || Ok::<_,()>(()), &mut |_: Event<'_, Fp>| Ok(())).unwrap().unwrap();
        // SymPy 1.14: L=[[1,0,0],[0,1,0],[2,-1,1]], U=[[-3,-2,-1],[0,2,3],[0,0,-1]] over GF(7).
        assert_eq!(swaps, vec![(0,1)]);
        assert_eq!(a, vec![vec![f(4),f(5),f(6)],vec![f(0),f(2),f(3)],vec![f(2),f(6),f(6)]]);
    }

    #[test]
    fn zero_diagonal_columns_and_empty_rectangles_keep_their_positions() {
        let mut a = vec![vec![q(0,1),q(3,1),q(4,1)],vec![q(0,1),q(2,1),q(7,1)],vec![q(0,1),q(1,1),q(5,1)]];
        assert_eq!(in_place(&mut a,&mut || Ok::<_,()>(()),&mut |_| Ok(())).unwrap(),Some(vec![]));
        assert_eq!(a[2],vec![q(0,1),q(1,2),q(3,2)]);
        for mut a in [Vec::<Vec<Rat>>::new(),vec![vec![]],vec![vec![],vec![]]] {
            assert_eq!(in_place(&mut a,&mut || Ok::<_,()>(()),&mut |_| panic!("no empty publication")).unwrap(),Some(vec![]));
        }
    }

    #[test]
    fn ragged_refusal_and_interruption_preserve_completed_publication_prefix() {
        let mut ragged = vec![vec![q(1,1),q(2,1)],vec![q(3,1)]];
        let before = ragged.clone();
        assert_eq!(in_place(&mut ragged,&mut || Ok::<_,()>(()),&mut |_| panic!("ragged publication")).unwrap(),None);
        assert_eq!(ragged,before);
        let mut a = vec![vec![q(2,1),q(2,1)],vec![q(6,1),q(4,1)]];
        let mut count=0;
        let result=in_place(&mut a,&mut || Ok(()),&mut |_| {count+=1;if count==2 {Err("interrupt")}else {Ok(())}});
        assert_eq!(result,Err("interrupt"));
        assert_eq!(a,vec![vec![q(2,1),q(2,1)],vec![q(3,1),q(4,1)]]);
        assert_eq!(count,2);
    }
}
