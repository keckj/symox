//! Exact row LLL with the reference nearest-integer and swap conventions.
//! Input is immutable; interruption/rank failure publish no partial result.
use super::{Int, Rat};

#[derive(Debug, PartialEq, Eq)]
pub enum Error<E> {
    Delta,
    Shape,
    Rank,
    Interrupted(E),
}

#[derive(Debug, PartialEq, Eq)]
pub struct Reduced {
    pub basis: Vec<Vec<Int>>,
    pub transform: Option<Vec<Vec<Int>>>,
    /// Source row when only permutation occurred; changed rows have no source.
    /// Basis and transformation undergo the same row operations.
    pub row_sources: Vec<Option<usize>>,
}

fn checkpoint<E>(check: &mut impl FnMut() -> Result<(), E>) -> Result<(), Error<E>> {
    check().map_err(Error::Interrupted)
}

fn dot<E>(
    a: &[Rat],
    b: &[Rat],
    check: &mut impl FnMut() -> Result<(), E>,
) -> Result<Rat, Error<E>> {
    let mut result = Rat::zero();
    for (i, (a, b)) in a.iter().zip(b).enumerate() {
        if i % 64 == 0 {
            checkpoint(check)?;
        }
        result = result.add(&a.mul(b));
    }
    Ok(result)
}

fn reduce<E>(
    basis: &mut [Vec<Int>],
    mu: &mut [Vec<Rat>],
    transform: &mut Option<Vec<Vec<Int>>>,
    k: usize,
    j: usize,
    half: &Rat,
    check: &mut impl FnMut() -> Result<(), E>,
) -> Result<(), Error<E>> {
    // Exactly floor(mu + 1/2), including negative half-integer ties.
    let r = mu[k][j].add(half).floor();
    for col in 0..basis[k].len() {
        if col % 64 == 0 {
            checkpoint(check)?;
        }
        basis[k][col] = basis[k][col].sub(&r.mul(&basis[j][col]));
    }
    for col in 0..j {
        checkpoint(check)?;
        mu[k][col] = mu[k][col].sub(&mu[j][col].mul_int(&r));
    }
    mu[k][j] = mu[k][j].sub(&Rat::from(r.clone()));
    if let Some(t) = transform {
        for col in 0..t.len() {
            if col % 64 == 0 {
                checkpoint(check)?;
            }
            t[k][col] = t[k][col].sub(&r.mul(&t[j][col]));
        }
    }
    Ok(())
}

/// Row reduction; columns also represents zero-row rectangular inputs.
/// The public adapter preserves its domain/protocol/error ordering.
pub fn checked<E>(
    input: &[Vec<Int>],
    columns: usize,
    delta: &Rat,
    with_transform: bool,
    check: &mut impl FnMut() -> Result<(), E>,
) -> Result<Reduced, Error<E>> {
    checkpoint(check)?;
    if delta <= &Rat::new(Int::one(), Int::from(4)) || delta >= &Rat::one() {
        return Err(Error::Delta);
    }
    let m = input.len();
    if m > columns || input.iter().any(|row| row.len() != columns) {
        return Err(Error::Shape);
    }
    let half = Rat::new(Int::one(), Int::from(2));
    let mut basis = input.to_vec();
    let mut row_sources: Vec<_> = (0..m).map(Some).collect();
    let mut stars: Vec<Vec<Rat>> = Vec::with_capacity(m);
    let mut mu = vec![vec![Rat::zero(); m]; m];
    let mut norms = vec![Rat::zero(); m];
    let mut transform = with_transform.then(|| {
        (0..m)
            .map(|i| (0..m).map(|j| Int::from(i64::from(i == j))).collect())
            .collect()
    });
    for i in 0..m {
        checkpoint(check)?;
        let mut original = Vec::with_capacity(columns);
        for (col, value) in basis[i].iter().enumerate() {
            if col % 64 == 0 {
                checkpoint(check)?;
            }
            original.push(Rat::from(value.clone()));
        }
        let mut star = original.clone();
        for j in 0..i {
            checkpoint(check)?;
            let product = dot(&original, &stars[j], check)?;
            if norms[j].is_zero() {
                return Err(Error::Rank);
            }
            mu[i][j] = product.div(&norms[j]);
            for (col, (value, previous)) in star.iter_mut().zip(&stars[j]).enumerate() {
                if col % 64 == 0 {
                    checkpoint(check)?;
                }
                *value = value.sub(&mu[i][j].mul(previous));
            }
        }
        norms[i] = dot(&star, &star, check)?;
        stars.push(star);
    }
    let mut k = 1;
    while k < m {
        checkpoint(check)?;
        if mu[k][k - 1].abs() > half {
            reduce(&mut basis, &mut mu, &mut transform, k, k - 1, &half, check)?;
            row_sources[k] = None;
        }
        let boundary = delta
            .sub(&mu[k][k - 1].mul(&mu[k][k - 1]))
            .mul(&norms[k - 1]);
        if norms[k] >= boundary {
            for j in (0..k - 1).rev() {
                checkpoint(check)?;
                if mu[k][j].abs() > half {
                    reduce(&mut basis, &mut mu, &mut transform, k, j, &half, check)?;
                    row_sources[k] = None;
                }
            }
            k += 1;
        } else {
            let nu = mu[k][k - 1].clone();
            let alpha = norms[k].add(&nu.mul(&nu).mul(&norms[k - 1]));
            if alpha.is_zero() {
                return Err(Error::Rank);
            }
            let beta = norms[k - 1].div(&alpha);
            mu[k][k - 1] = nu.mul(&beta);
            norms[k] = norms[k].mul(&beta);
            norms[k - 1] = alpha;
            basis.swap(k, k - 1);
            row_sources.swap(k, k - 1);
            for j in 0..k - 1 {
                let old = mu[k][j].clone();
                mu[k][j] = mu[k - 1][j].clone();
                mu[k - 1][j] = old;
            }
            for i in k + 1..m {
                checkpoint(check)?;
                let xi = mu[i][k].clone();
                mu[i][k] = mu[i][k - 1].sub(&nu.mul(&xi));
                mu[i][k - 1] = mu[k][k - 1].mul(&mu[i][k]).add(&xi);
            }
            if let Some(t) = &mut transform {
                t.swap(k, k - 1);
            }
            k = k.saturating_sub(1).max(1);
        }
    }
    checkpoint(check)?;
    Ok(Reduced {
        basis,
        transform,
        row_sources,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ints(rows: &[&[i64]]) -> Vec<Vec<Int>> {
        rows.iter()
            .map(|r| r.iter().copied().map(Int::from).collect())
            .collect()
    }
    #[test]
    fn exact_basis_and_unimodular_transform() {
        let input = ints(&[&[1, 1, 1], &[-1, 0, 2], &[3, 5, 6]]);
        let output = checked(
            &input,
            3,
            &Rat::new(Int::from(3), Int::from(4)),
            true,
            &mut || Ok::<_, ()>(()),
        )
        .unwrap();
        assert_eq!(output.basis, ints(&[&[0, 1, 0], &[1, 0, 1], &[-1, 0, 2]]));
        let t = output.transform.unwrap();
        for (i, source) in output.row_sources.iter().enumerate() {
            if let Some(j) = source {
                assert_eq!(output.basis[i], input[*j]);
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let value = (0..3).fold(Int::zero(), |v, k| v.add(&t[i][k].mul(&input[k][j])));
                assert_eq!(value, output.basis[i][j]);
            }
        }
        assert_eq!(input, ints(&[&[1, 1, 1], &[-1, 0, 2], &[3, 5, 6]]));
    }
    #[test]
    fn validation_rank_empty_single_zero_and_interruptions() {
        let delta = Rat::new(Int::from(3), Int::from(4));
        let mut check = || Ok::<_, ()>(());
        assert_eq!(
            checked(&ints(&[&[1], &[2]]), 1, &delta, false, &mut check),
            Err(Error::Shape)
        );
        assert_eq!(
            checked(&[], 0, &Rat::one(), false, &mut check),
            Err(Error::Delta)
        );
        assert_eq!(
            checked(&ints(&[&[1, 2], &[2, 4]]), 2, &delta, false, &mut check),
            Err(Error::Rank)
        );
        assert_eq!(
            checked(&[], 99, &delta, true, &mut check).unwrap(),
            Reduced {
                basis: vec![],
                transform: Some(vec![]),
                row_sources: vec![]
            }
        );
        assert_eq!(
            checked(&ints(&[&[0, 0]]), 2, &delta, false, &mut check)
                .unwrap()
                .basis,
            ints(&[&[0, 0]])
        );
        let input = ints(&[&[1, 1, 1], &[-1, 0, 2], &[3, 5, 6]]);
        let original = input.clone();
        let mut calls = 0;
        assert_eq!(
            checked(&input, 3, &delta, true, &mut || {
                calls += 1;
                if calls == 16 {
                    Err("stop")
                } else {
                    Ok(())
                }
            }),
            Err(Error::Interrupted("stop"))
        );
        assert_eq!(input, original);
    }
    #[test]
    fn untouched_rows_and_swaps_keep_source_identity_metadata() {
        let delta = Rat::new(Int::from(3), Int::from(4));
        for (rows, sources) in [
            (ints(&[&[1, 0], &[0, 2]]), vec![Some(0), Some(1)]),
            (ints(&[&[2, 0], &[0, 1]]), vec![Some(1), Some(0)]),
            (ints(&[&[1, 0], &[2, 1]]), vec![Some(0), None]),
        ] {
            let out = checked(&rows, 2, &delta, true, &mut || Ok::<_, ()>(())).unwrap();
            assert_eq!(out.row_sources, sources);
        }
    }
    #[test]
    fn nearest_integer_negative_ties_are_floor_after_adding_half() {
        let half = Rat::new(Int::one(), Int::from(2));
        for (n, expected) in [(-5, -2), (-3, -1), (-1, 0), (1, 1), (3, 2)] {
            assert_eq!(
                Rat::new(Int::from(n), Int::from(2)).add(&half).floor(),
                Int::from(expected)
            );
        }
    }
}
