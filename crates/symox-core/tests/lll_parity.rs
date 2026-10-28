//! Exact reference fixtures plus independently recomputed lattice certificates.
//! The fixture provenance records reference floating-floor failures separately.
use num_bigint::BigInt;
use symox::poly::{
    lll::{checked, Error},
    Int, Rat,
};
fn integers(text: &str) -> Vec<Int> {
    if text.is_empty() {
        return vec![];
    }
    text.split(',')
        .map(|v| Int::from(v.parse::<BigInt>().unwrap()))
        .collect()
}
#[test]
fn all_recorded_exact_lll_bases_and_transformations() {
    let mut cases = 0;
    for line in include_str!("fixtures/lll_sympy_1_14.tsv").lines() {
        let f: Vec<_> = line.split('\t').collect();
        let m: usize = f[0].parse().unwrap();
        let n: usize = f[1].parse().unwrap();
        let d = Rat::new(
            Int::from(f[2].parse::<BigInt>().unwrap()),
            Int::from(f[3].parse::<BigInt>().unwrap()),
        );
        let flat = integers(f[5]);
        let rows: Vec<Vec<Int>> = (0..m).map(|i| flat[i * n..(i + 1) * n].to_vec()).collect();
        let output = checked(&rows, n, &d, true, &mut || Ok::<_, ()>(()));
        match f[4] {
            "rank" => assert_eq!(output, Err(Error::Rank), "case {cases}"),
            "shape" => assert_eq!(output, Err(Error::Shape), "case {cases}"),
            "delta" => assert_eq!(output, Err(Error::Delta), "case {cases}"),
            "ok" => {
                let output = output.unwrap();
                assert_eq!(
                    output.basis.iter().flatten().cloned().collect::<Vec<_>>(),
                    integers(f[6]),
                    "basis case {cases}"
                );
                assert_eq!(
                    output
                        .transform
                        .as_ref()
                        .unwrap()
                        .iter()
                        .flatten()
                        .cloned()
                        .collect::<Vec<_>>(),
                    integers(f[7]),
                    "transform case {cases}"
                );
                let no_t = checked(&rows, n, &d, false, &mut || Ok::<_, ()>(())).unwrap();
                assert_eq!(no_t.basis, output.basis);
                assert!(no_t.transform.is_none());
                for (i, row) in output.basis.iter().enumerate() {
                    for (j, value) in row.iter().enumerate() {
                        let actual = (0..m).fold(Int::zero(), |sum, k| {
                            sum.add(&output.transform.as_ref().unwrap()[i][k].mul(&rows[k][j]))
                        });
                        assert_eq!(&actual, value, "T*A=B case {cases}");
                    }
                }
            }
            _ => panic!("unknown fixture"),
        }
        cases += 1;
    }
    assert_eq!(cases, 651);
}
