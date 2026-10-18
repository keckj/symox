//! Finite permutation matrices and held row/column permutations.
use super::{finite_dim, MatrixExpr, MatrixKind};
use crate::expr::{Expr, Func};
use crate::functions::func;
use crate::matrix::{MatrixError, MResult};
use std::sync::OnceLock;

pub(super) fn axis_expr(axis: u8) -> Expr {
    // scalar_map caches by DAG identity. The typed axis must expose the same
    // scalar owner during collection and reconstruction, like stored fields.
    static AXES: OnceLock<[Expr;2]> = OnceLock::new();
    AXES.get_or_init(||[Expr::zero(),Expr::one()])[axis as usize].clone()
}

pub(super) fn validate(images: &[usize]) -> MResult<()> {
    let mut seen = vec![false; images.len()];
    for &image in images {
        if image >= images.len() || seen[image] {
            return Err(MatrixError::Value("Permutation images must be a bijection".into()));
        }
        seen[image] = true;
    }
    Ok(())
}

pub(super) fn inverse(images: &[usize]) -> Vec<usize> {
    let mut output = vec![0; images.len()];
    for (i, &image) in images.iter().enumerate() { output[image] = i; }
    output
}

pub(super) fn compose(first: &[usize], second: &[usize]) -> Vec<usize> {
    first.iter().map(|&i| second[i]).collect()
}

pub(super) fn signature(images: &[usize]) -> i64 {
    let mut seen = vec![false; images.len()];
    let mut sign = 1;
    for start in 0..images.len() {
        if seen[start] { continue; }
        let mut index = start;
        let mut length = 0;
        while !seen[index] {
            seen[index] = true;
            length += 1;
            index = images[index];
        }
        if length % 2 == 0 { sign = -sign; }
    }
    sign
}

pub(crate) fn label(images: &[usize]) -> String {
    let mut seen = vec![false; images.len()];
    let mut output = String::new();
    if let Some(last) = images.len().checked_sub(1) {
        if images[last] == last { output.push_str(&format!("({last})")); }
    }
    for start in 0..images.len() {
        if seen[start] { continue; }
        let mut cycle = Vec::new();
        let mut index = start;
        while !seen[index] {
            seen[index] = true;
            cycle.push(index.to_string());
            index = images[index];
        }
        if cycle.len() > 1 { output.push_str(&format!("({})", cycle.join(" "))); }
    }
    if output.is_empty() { output.push_str("()"); }
    output
}

pub(super) fn power(images: &[usize], exponent: i64) -> Vec<usize> {
    let mut base = if exponent < 0 { inverse(images) } else { images.to_vec() };
    let mut result = (0..images.len()).collect::<Vec<_>>();
    let mut count = exponent.unsigned_abs();
    while count != 0 {
        if count & 1 != 0 { result = compose(&result, &base); }
        count >>= 1;
        if count != 0 { base = compose(&base, &base); }
    }
    result
}

pub(super) fn apply(images: &[usize], index: &Expr) -> MResult<Expr> {
    if images.is_empty() {
        return Err(MatrixError::Index("Cannot index an empty permutation".into()));
    }
    let size = Expr::int(images.len() as i64);
    MatrixExpr::validate_index(index, &size)?;
    let index = normalize_index(index, &size);
    if let Some(i) = index.as_i64() { return Ok(Expr::int(images[i as usize] as i64)); }
    if images.iter().enumerate().all(|(i, &j)| i == j) { return Ok(index); }
    let cases = images.iter().enumerate().map(|(i, &image)| {
        let condition = if i+1 == images.len() { Expr::boolean(true) }
                        else { crate::relational::eq(&index, &Expr::int(i as i64)) };
        (Expr::int(image as i64), condition)
    }).collect();
    Ok(crate::relational::piecewise(cases))
}

fn normalize_index(index: &Expr, size: &Expr) -> Expr {
    if crate::assumptions::is_nonnegative(index)==Some(true) { return index.clone(); }
    let shifted=crate::ops::add2(index,size);
    if crate::assumptions::is_negative(index)==Some(true) { return shifted; }
    crate::relational::piecewise(vec![
        (shifted,crate::relational::lt(index,&Expr::zero())),
        (index.clone(),Expr::boolean(true)),
    ])
}

impl MatrixExpr {
    pub fn permutation(images: Vec<usize>, evaluate: bool) -> MResult<Self> {
        validate(&images)?;
        let size = Expr::int(images.len() as i64);
        if evaluate && images.iter().enumerate().all(|(i, &j)| i == j) {
            return Self::identity(size);
        }
        Ok(Self::new(MatrixKind::Permutation(images), size.clone(), size))
    }

    pub fn permute(&self, mut images: Vec<usize>, axis: u8, evaluate: bool) -> MResult<Self> {
        validate(&images)?;
        if axis > 1 { return Err(MatrixError::Value("Permutation axis must be zero or one".into())); }
        let size = finite_dim(if axis == 0 { self.rows() } else { self.cols() })?;
        if images.len() > size {
            if images[size..].iter().enumerate().any(|(i, &j)| j != i+size) {
                return Err(MatrixError::Shape("Permutation cannot shrink across a moved index".into()));
            }
            images.truncate(size);
        } else { images.extend(images.len()..size); }
        if evaluate {
            if images.iter().enumerate().all(|(i, &j)| i == j) ||
                matches!(self.kind(), MatrixKind::Zero | MatrixKind::One) {
                return Ok(self.clone());
            }
            match self.kind() {
                MatrixKind::Identity => return Self::permutation(
                    if axis == 0 { images } else { inverse(&images) }, true),
                MatrixKind::Permute(parent, inner, inner_axis) if axis == *inner_axis =>
                    return parent.permute(compose(&images, inner), axis, true),
                _ => {},
            }
        }
        Ok(Self::new(MatrixKind::Permute(self.clone(), images, axis),
                     self.rows().clone(), self.cols().clone()))
    }

    pub(super) fn permutation_entry(images: &[usize], i: &Expr, j: &Expr) -> MResult<Expr> {
        let image = apply(images, i)?;
        let j = normalize_index(j, &Expr::int(images.len() as i64));
        Ok(func(Func::Special(crate::special::SpecialFn::KroneckerDelta), vec![image, j]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::Matrix;
    use crate::ops;

    fn permutations(values: &mut [usize], start: usize, output: &mut Vec<Vec<usize>>) {
        if start == values.len() { output.push(values.to_vec()); return; }
        for i in start..values.len() {
            values.swap(i,start);
            permutations(values,start+1,output);
            values.swap(i,start);
        }
    }

    #[test]
    fn finite_permutations_preserve_composition_powers_and_determinants() {
        for size in 0..=5 {
            let mut images = (0..size).collect::<Vec<_>>();
            let mut cases = Vec::new();
            permutations(&mut images,0,&mut cases);
            for images in cases {
                let p = MatrixExpr::permutation(images.clone(),false).unwrap();
                let matrix = p.as_explicit().unwrap();
                assert_eq!(p.determinant().unwrap(),matrix.det().unwrap());
                assert_eq!(p.trace().unwrap(),matrix.trace().unwrap());
                assert_eq!(p.transpose(true).unwrap().as_explicit().unwrap(),matrix.transpose());
                assert_eq!(p.adjoint(true).unwrap().as_explicit().unwrap(),matrix.transpose());
                assert_eq!(p.conjugate(true).unwrap(),p);
                for exponent in [-3,-1,0,1,2,5] {
                    assert_eq!(p.pow(Expr::int(exponent),true).unwrap().as_explicit().unwrap(),
                               matrix.powi(exponent).unwrap());
                }
                let q_images = (0..size).rev().collect::<Vec<_>>();
                let q = MatrixExpr::permutation(q_images,false).unwrap();
                let product = MatrixExpr::mul(Expr::one(),&[p.clone(),q.clone()],true).unwrap();
                assert_eq!(product.as_explicit().unwrap(),matrix.matmul(&q.as_explicit().unwrap()).unwrap());
                assert_eq!(p.permute(images.clone(),0,true).unwrap().as_explicit().unwrap(),
                           matrix.matmul(&matrix).unwrap());
            }
        }
    }

    #[test]
    fn row_column_permutations_resize_and_compose_without_densifying() {
        let matrix = Matrix::from_fn(4,3,|i,j|Expr::int((i*10+j) as i64));
        let expression = MatrixExpr::explicit(matrix.clone());
        let rows = vec![2,0,1]; // Extend by the fixed point 3.
        let row_result = expression.permute(rows.clone(),0,false).unwrap();
        let expected = Matrix::from_fn(4,3,|i,j|matrix.get([2,0,1,3][i],j).clone());
        assert_eq!(row_result.as_explicit().unwrap(),expected);
        let columns = vec![1,0,2,3]; // Drop the fixed point 3.
        let result = row_result.permute(columns,1,false).unwrap();
        assert_eq!(result.as_explicit().unwrap(),Matrix::from_fn(4,3,|i,j|expected.get(i,[1,0,2][j]).clone()));
        assert_eq!(result.transpose(true).unwrap().as_explicit().unwrap(),result.as_explicit().unwrap().transpose());
        let nested = row_result.permute(vec![1,2,0,3],0,true).unwrap();
        assert_eq!(nested,expression);
        assert!(expression.permute(vec![3,1,2,0],1,false).is_err());
        assert!(expression.permute(vec![0,1],2,false).is_err());
        assert!(MatrixExpr::permutation(vec![1,1],false).is_err());
        let large = MatrixExpr::symbol("A",Expr::int(10000),Expr::int(10000)).unwrap();
        let permuted = large.permute(vec![1,0],0,false).unwrap();
        assert_eq!(permuted.matrix_symbols(),vec![large]);
        assert!(permuted.scalar_free_symbols().is_empty());
    }

    #[test]
    fn symbolic_permutation_entries_and_scalar_derivatives_are_exact() {
        let p = MatrixExpr::permutation(vec![1,2,0],false).unwrap();
        let i = Expr::symbol("i");
        let j = Expr::symbol("j");
        let entry = p.element(&i,&j).unwrap();
        for row in -3..3 {
            for col in -3..3 {
                let result = crate::subs::subs(&entry,&[(i.clone(),Expr::int(row)),(j.clone(),Expr::int(col))]);
                assert_eq!(result,p.element(&Expr::int(row),&Expr::int(col)).unwrap());
            }
        }
        assert_eq!(p.element(&Expr::int(-1),&Expr::zero()).unwrap(),Expr::one());
        assert!(p.element(&Expr::int(-4),&Expr::zero()).is_err());
        let x = Expr::symbol("x");
        let a = MatrixExpr::explicit(Matrix::from_fn(3,2,|row,col|
            ops::mul2(&Expr::int((1+row+col) as i64),&ops::pow(&x,&Expr::int(2)))));
        let permuted = a.permute(vec![2,0,1],0,false).unwrap();
        let derivative = permuted.diff(&x).unwrap().as_explicit().unwrap();
        assert_eq!(derivative,permuted.as_explicit().unwrap().applyfunc(|e|crate::diff::diff(e,&x)));
        let mapped = permuted.scalar_map(&mut |e|crate::subs::subs(e,&[(x.clone(),Expr::int(3))])).unwrap();
        assert_eq!(mapped.as_explicit().unwrap(),permuted.as_explicit().unwrap().applyfunc(|e|crate::subs::subs(e,&[(x.clone(),Expr::int(3))])));
        let square = MatrixExpr::symbol("A",Expr::int(3),Expr::int(3)).unwrap();
        let row = square.permute(vec![1,2,0],0,false).unwrap();
        let column = square.permute(vec![1,2,0],1,false).unwrap();
        assert_eq!(row.scalar_map(&mut |e|crate::subs::subs(e,&[(Expr::zero(),Expr::one())])).unwrap(),column);
    }

    #[test]
    fn native_and_sympy_permutation_text_roundtrips() {
        for (text,images) in [
            ("Permutation()",vec![]),
            ("Permutation(2)",vec![0,1,2]),
            ("Permutation(0, 1, 2)",vec![1,2,0]),
            ("Permutation(3)(1, 2)",vec![0,2,1,3]),
            ("Permutation([2,0,1])",vec![2,0,1]),
        ] {
            let matrix = MatrixExpr::permutation(images,false).unwrap();
            for source in [format!("PermutationMatrix({text})"),matrix.to_repr()] {
                let parsed = crate::parser::parse_str(&source).unwrap();
                assert!(matches!(parsed.kind(),crate::expr::Kind::Func(Func::MatrixProxy(parent),_) if **parent==matrix));
            }
            let parent = MatrixExpr::symbol("A",matrix.rows().clone(),Expr::int(2)).unwrap();
            let permuted = parent.permute(if let MatrixKind::Permutation(images)=matrix.kind(){images.clone()}else{unreachable!()},0,false).unwrap();
            assert!(crate::parser::parse_str(&permuted.to_repr()).is_ok());
        }
        for invalid in ["PermutationMatrix(Permutation([1,1]))",
                        "PermutationMatrix(Permutation(100000000000))",
                        "PermutationMatrix(Permutation(0,0))"] {
            assert!(crate::parser::parse_str(invalid).is_err());
        }
    }
}
