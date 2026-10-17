//! Bounded structural transforms. Permutations store only their finite support,
//! including when a sparse caller has dimensions larger than machine integers.
use super::functions::budget::Budget;
use super::{MResult, Matrix, MatrixError};
use crate::expr::Kind;
use crate::{Expr, Number};
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use std::collections::{BTreeMap, BTreeSet};

pub const MAX_ITEMS: usize = 100_000;
fn limit() -> MatrixError {
    MatrixError::NotImplemented("Matrix structural transformation resource limit exceeded".into())
}
pub fn check_count(count: &BigInt) -> MResult<()> {
    if count < &BigInt::zero() || count > &BigInt::from(MAX_ITEMS) {
        Err(limit())
    } else {
        Ok(())
    }
}
fn coordinate(x: &BigInt) -> MResult<()> {
    if x.bits() > 1_000_000 {
        Err(limit())
    } else {
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct IndexPermutation {
    source: BTreeMap<BigInt, BigInt>,
}
impl IndexPermutation {
    pub fn new(
        dimension: &BigInt,
        specification: &[Vec<BigInt>],
        cycles: bool,
        backward: bool,
    ) -> MResult<Self> {
        coordinate(dimension)?;
        if dimension < &BigInt::zero() {
            return Err(MatrixError::Value("Negative permutation dimension".into()));
        }
        let count = specification
            .iter()
            .try_fold(0usize, |n, x| n.checked_add(x.len()).ok_or_else(limit))?;
        if count > MAX_ITEMS || specification.len() > MAX_ITEMS {
            return Err(limit());
        }
        let mut source = BTreeMap::new();
        for row in specification {
            for x in row {
                coordinate(x)?;
                if x < &BigInt::zero() || x > dimension {
                    return Err(MatrixError::Index("swap indices out of range".into()));
                }
            }
        }
        if cycles {
            let mut apply = |cycle: &Vec<BigInt>| -> MResult<()> {
                let unique: BTreeSet<_> = cycle.iter().collect();
                if unique.len() != cycle.len() {
                    return Err(MatrixError::Value(
                        "All elements must be unique in a cycle.".into(),
                    ));
                }
                if cycle.len() > 1 {
                    let values: Vec<_> = cycle
                        .iter()
                        .map(|i| source.get(i).cloned().unwrap_or_else(|| i.clone()))
                        .collect();
                    for (position, i) in cycle.iter().enumerate() {
                        let value = values[(position + 1) % values.len()].clone();
                        if i == &value {
                            source.remove(i);
                        } else {
                            source.insert(i.clone(), value);
                        }
                    }
                }
                Ok(())
            };
            if backward {
                for cycle in specification.iter().rev() {
                    apply(cycle)?;
                }
            } else {
                for cycle in specification {
                    apply(cycle)?;
                }
            }
        } else {
            if specification.len() != 1 {
                return Err(MatrixError::Value("Expected one permutation array".into()));
            }
            let values = &specification[0];
            let unique: BTreeSet<_> = values.iter().collect();
            if unique.len() != values.len()
                || values
                    .iter()
                    .any(|x| x.to_usize().is_none_or(|i| i >= values.len()))
            {
                return Err(MatrixError::Value(
                    "Integers 0 through n-1 must appear exactly once in permutation array".into(),
                ));
            }
            for (i, value) in values.iter().enumerate() {
                let i = BigInt::from(i);
                if &i != value {
                    source.insert(i, value.clone());
                }
            }
        }
        if source.iter().any(|(i, j)| i >= dimension || j >= dimension) {
            return Err(MatrixError::Index("swap indices out of range".into()));
        }
        Ok(Self { source })
    }
    /// Old index -> new index, suitable for dictionary-backed sparse entries.
    pub fn destinations(&self) -> BTreeMap<BigInt, BigInt> {
        self.source
            .iter()
            .map(|(new, old)| (old.clone(), new.clone()))
            .collect()
    }
    pub fn apply(&self, m: &Matrix, rows: bool) -> MResult<Matrix> {
        check_count(&BigInt::from(m.data().len()))?;
        if m.data().is_empty() {
            return Ok(m.clone());
        }
        let mapping: BTreeMap<usize, usize> = self
            .source
            .iter()
            .map(|(a, b)| {
                Ok((
                    a.to_usize().ok_or_else(limit)?,
                    b.to_usize().ok_or_else(limit)?,
                ))
            })
            .collect::<MResult<_>>()?;
        let dimension = if rows { m.rows() } else { m.cols() };
        if mapping
            .iter()
            .any(|(a, b)| *a >= dimension || *b >= dimension)
        {
            return Err(MatrixError::Index("swap indices out of range".into()));
        }
        let data = (0..m.data().len())
            .map(|position| {
                let (i, j) = (position / m.cols(), position % m.cols());
                m.get(
                    if rows {
                        *mapping.get(&i).unwrap_or(&i)
                    } else {
                        i
                    },
                    if rows {
                        j
                    } else {
                        *mapping.get(&j).unwrap_or(&j)
                    },
                )
                .clone()
            })
            .collect();
        Matrix::new(m.rows(), m.cols(), data)
    }
}

/// Transform compact vectors in one bounded native arithmetic operation.
/// Coefficients multiply on the LEFT, including noncommutative coefficients.
pub fn elementary_values(
    add: bool,
    coefficient: &Expr,
    targets: &[Expr],
    sources: &[Expr],
) -> MResult<Vec<Expr>> {
    if targets.len() != sources.len() {
        return Err(MatrixError::Shape(
            "Elementary operation vector lengths differ".into(),
        ));
    }
    check_count(&BigInt::from(targets.len()))?;
    let mut budget = Budget::new();
    budget.expression(coefficient)?;
    budget.allocate(targets.len())?;
    let mut values = Vec::with_capacity(targets.len());
    for (target, source) in targets.iter().zip(sources) {
        budget.work(4)?;
        budget.expression(source)?;
        if add {
            budget.expression(target)?;
        }
        let raw_product = Expr::from_unevaluated_kind(Kind::Mul(
            Number::int(1),
            vec![
                (coefficient.clone(), Expr::one()),
                (source.clone(), Expr::one()),
            ],
        ));
        let raw = if add {
            Expr::from_unevaluated_kind(Kind::Add(
                Number::int(0),
                vec![
                    (target.clone(), Number::int(1)),
                    (raw_product, Number::int(1)),
                ],
            ))
        } else {
            raw_product
        };
        budget.scalar_preflight(&raw, None)?;
        let product = super::nonfinite_product(coefficient, source);
        let value = if add {
            crate::ops::add(&[target.clone(), product])
        } else {
            product
        };
        budget.expression(&value)?;
        values.push(value);
    }
    Ok(values)
}

impl Matrix {
    pub fn elementary_transform(
        &self,
        rows: bool,
        operation: u8,
        target: usize,
        source: usize,
        coefficient: &Expr,
    ) -> MResult<Matrix> {
        let dimension = if rows { self.rows() } else { self.cols() };
        let width = if rows { self.cols() } else { self.rows() };
        if target >= dimension || (operation != 0 && source >= dimension) {
            return Err(MatrixError::Value(
                "Elementary operation index out of range".into(),
            ));
        }
        if operation > 2 {
            return Err(MatrixError::Value("Unknown elementary operation".into()));
        }
        if operation != 0 && target == source {
            return Err(MatrixError::Value(
                "Elementary operation indices must be distinct".into(),
            ));
        }
        check_count(&BigInt::from(self.data().len()))?;
        if operation == 1 {
            let mut result = self.clone();
            if rows {
                result.row_swap(target, source);
            } else {
                result.col_swap(target, source);
            }
            return Ok(result);
        }
        let target_values: Vec<_> = (0..width)
            .map(|p| {
                if rows {
                    self.get(target, p)
                } else {
                    self.get(p, target)
                }
                .clone()
            })
            .collect();
        let source_values = if operation == 0 {
            target_values.clone()
        } else {
            (0..width)
                .map(|p| {
                    if rows {
                        self.get(source, p)
                    } else {
                        self.get(p, source)
                    }
                    .clone()
                })
                .collect()
        };
        let values =
            elementary_values(operation == 2, coefficient, &target_values, &source_values)?;
        let mut result = self.clone();
        for (p, value) in values.into_iter().enumerate() {
            if rows {
                result.set(target, p, value);
            } else {
                result.set(p, target, value);
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ints(v: &[i64]) -> Vec<BigInt> {
        v.iter().map(|n| BigInt::from(*n)).collect()
    }
    #[test]
    fn cycle_directions_and_sparse_support() {
        let p = vec![ints(&[0, 1]), ints(&[0, 2])];
        let a = Matrix::new(3, 1, vec![Expr::int(1), Expr::int(2), Expr::int(3)]).unwrap();
        let f = IndexPermutation::new(&3.into(), &p, true, false)
            .unwrap()
            .apply(&a, true)
            .unwrap();
        assert_eq!(f.data(), &[Expr::int(3), Expr::int(1), Expr::int(2)]);
        let b = IndexPermutation::new(&3.into(), &p, true, true)
            .unwrap()
            .apply(&f, true)
            .unwrap();
        assert_eq!(b, a);
        let huge = BigInt::from(10u8).pow(100);
        assert_eq!(
            IndexPermutation::new(&huge, &[ints(&[0, 1])], true, false)
                .unwrap()
                .destinations()
                .len(),
            2
        );
    }
    #[test]
    fn elementary_inverse_and_nonfinite_zero() {
        let a = Matrix::new(
            2,
            2,
            vec![Expr::int(1), Expr::int(2), Expr::int(3), Expr::int(4)],
        )
        .unwrap();
        let b = a
            .elementary_transform(true, 2, 0, 1, &Expr::int(2))
            .unwrap();
        assert_eq!(
            b.elementary_transform(true, 2, 0, 1, &Expr::int(-2))
                .unwrap(),
            a
        );
        assert_eq!(
            elementary_values(
                false,
                &Expr::from_kind(Kind::Const(crate::expr::Constant::Infinity)),
                &[Expr::zero()],
                &[Expr::zero()]
            )
            .unwrap(),
            vec![Expr::nan()]
        );
    }
}
