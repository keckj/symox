//! Exact modulo-D Hermite recurrence with ordered publication and modulus provenance.
use super::Int;

pub enum Event<'a> {
    Work(usize, usize, &'a Int),
    Output(usize, usize, &'a Int),
    WorkModulus(usize, usize),
    OutputModulus(usize, usize),
    ModulusUpdate(&'a Int),
}

fn gcdex(a: &Int, b: &Int) -> (Int, Int, Int) {
    let (g, mut u, mut v) = a.gcdex(b);
    if !a.is_zero() && b.div_rem_trunc(a).1.is_zero() {
        u = Int::from(if a.is_negative() { -1 } else { 1 });
        v = Int::zero();
    }
    (u, v, g)
}

fn symmetric(value: Int, modulus: &Int) -> Int {
    let value = value.div_mod_floor(modulus).1;
    // symmetric_residue keeps +R/2 on the positive side for even R.
    if value.add(&value) > *modulus { value.sub(modulus) } else { value }
}

/// Original modulo-D recurrence, with arbitrary precision scalar arithmetic.
/// The caller keeps the actual domain/error/shape/conversion prefixes and W
/// construction. Invalid structural input declines before any publication.
/// A valid determinant multiple remains the public algorithm's assumption;
/// this routine does not invent an empty answer for a deficient-rank input.
pub fn modular<E>(a: &mut [Vec<Int>], n: usize, modulus: &Int,
    check: &mut impl FnMut() -> Result<(), E>,
    publish: &mut impl FnMut(Event<'_>) -> Result<(), E>,
) -> Result<Option<Vec<Vec<Int>>>, E> {
    check()?;
    let m = a.len();
    if n < m || modulus <= &Int::zero() || a.iter().any(|row| row.len() != n) {
        return Ok(None);
    }
    let mut w = vec![vec![Int::zero(); m]; m];
    let mut r = modulus.clone();
    let mut k = n;
    for i in (0..m).rev() {
        check()?;
        k -= 1;
        for j in (0..k).rev() {
            check()?;
            if a[i][j].is_zero() { continue; }
            let (u, v, d) = gcdex(&a[i][k], &a[i][j]);
            let left = a[i][k].div_exact_unchecked(&d);
            let right = a[i][j].div_exact_unchecked(&d);
            for (ii, row) in a.iter_mut().enumerate() {
                check()?;
                let old = row[k].clone();
                let value = symmetric(u.mul(&old).add(&v.mul(&row[j])), &r);
                publish(Event::Work(ii, k, &value))?;
                row[k] = value;
                let value = symmetric(right.neg().mul(&old).add(&left.mul(&row[j])), &r);
                publish(Event::Work(ii, j, &value))?;
                row[j] = value;
            }
        }
        let mut b = a[i][k].clone();
        if b.is_zero() {
            b = r.clone();
            publish(Event::WorkModulus(i, k))?;
            a[i][k] = b.clone();
        }
        let (u, _, d) = gcdex(&b, &r);
        for ii in 0..m {
            check()?;
            let value = u.mul(&a[ii][k]).div_mod_floor(&r).1;
            publish(Event::Output(ii, i, &value))?;
            w[ii][i] = value;
        }
        if w[i][i].is_zero() {
            publish(Event::OutputModulus(i, i))?;
            w[i][i] = r.clone();
        }
        for j in (i+1)..m {
            check()?;
            let q = w[i][j].div_mod_floor(&w[i][i]).0;
            for (ii, row) in w.iter_mut().enumerate() {
                check()?;
                let old = row[j].clone();
                let value = old.sub(&q.mul(&row[i]));
                publish(Event::Output(ii, j, &value))?;
                row[j] = value;
                // Original add_columns writes the retained i cell too.
                publish(Event::Output(ii, i, &row[i]))?;
            }
        }
        let next_modulus = r.div_exact_unchecked(&d);
        publish(Event::ModulusUpdate(&next_modulus))?;
        r = next_modulus;
    }
    check()?;
    Ok(Some(w))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn matrix(a: &[&[i64]]) -> Vec<Vec<Int>> {
        a.iter().map(|r| r.iter().copied().map(Int::from).collect()).collect()
    }
    fn run(a: &mut [Vec<Int>], n: usize, d: Int) -> Vec<Vec<Int>> {
        modular(a,n,&d,&mut || Ok::<_,()>(()),&mut |_| Ok(())).unwrap().unwrap()
    }
    #[test]
    fn textbook_and_arbitrary_determinant_multiples() {
        let a=matrix(&[&[12,6,4],&[3,9,6],&[2,16,14]]);
        let expected=matrix(&[&[10,0,2],&[0,15,3],&[0,0,2]]);
        for d in [300,600,900,1200] { assert_eq!(run(&mut a.clone(),3,Int::from(d)),expected); }
        let huge=Int::from(1_000_000_000).pow(20).mul(&Int::from(300));
        assert_eq!(run(&mut a.clone(),3,huge),expected);
    }
    #[test]
    fn full_rank_rectangles_match_unimodular_ordinary_reduction() {
        // First m columns have triangular determinant, the remaining columns
        // enlarge the same integer lattice. No heuristic bound is needed.
        for m in 1..=6 {
            for extra in 0..=4 {
                for seed in 0..=7 {
                    let n=m+extra;
                    let mut a=vec![vec![Int::zero();n];m];
                    let mut d=Int::one();
                    for i in 0..m {
                        let pivot=Int::from((i+2+seed) as i64);
                        a[i][i]=pivot.clone();d=d.mul(&pivot);
                        for j in (i+1)..n {a[i][j]=Int::from(((seed*7+i*3+j*5)%17) as i64-8);}
                    }
                    // A negative pivot and full row permutation are deliberate.
                    a[0][0]=a[0][0].neg();
                    a.rotate_left(seed % m);
                    let mut ordinary=a.clone();
                    let k=super::super::hnf::ordinary(&mut ordinary,n,&mut ||Ok::<_,()>(()),&mut |_|Ok(())).unwrap().unwrap();
                    assert_eq!(k,extra);
                    let expected=ordinary.iter().map(|row|row[k..].to_vec()).collect::<Vec<_>>();
                    for multiple in [1,2,7] {assert_eq!(run(&mut a.clone(),n,d.mul(&Int::from(multiple))),expected);}
                }
            }
        }
    }
    #[test]
    fn symmetric_ties_empty_shapes_and_before_publication_refusals() {
        assert_eq!(symmetric(Int::from(2),&Int::from(4)),Int::from(2));
        assert_eq!(symmetric(Int::from(3),&Int::from(4)),Int::from(-1));
        assert_eq!(run(&mut [],17,Int::from(3)),Vec::<Vec<Int>>::new());
        for (mut a,n,d) in [(matrix(&[&[1],&[2]]),1,3),(matrix(&[&[1],&[2,3]]),2,3),(matrix(&[&[1]]),1,0)] {
            let old=a.clone();
            assert!(modular(&mut a,n,&Int::from(d),&mut ||Ok::<_,()>(()),&mut |_|panic!("refusal wrote")).unwrap().is_none());
            assert_eq!(a,old);
        }
    }
    #[test]
    fn interrupted_publication_preserves_completed_assignments() {
        let mut a=matrix(&[&[2,3]]);let mut writes=0;
        let result=modular(&mut a,2,&Int::from(6),&mut ||Ok(()),&mut |event|{
            if let Event::Work(..)=event {writes+=1;if writes==2{return Err("interrupted");}}
            Ok(())
        });
        assert!(matches!(result,Err("interrupted")));
        assert_eq!(writes,2);
        assert_eq!(a,matrix(&[&[2,1]]));
    }
}
