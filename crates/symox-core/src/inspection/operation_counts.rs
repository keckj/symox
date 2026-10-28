//! SymPy operation labels over native scalar DAGs, with exact occurrence weights.
//! Inspection must preserve held arithmetic: synthesised public arguments use
//! held constructors, and fraction normalization runs only after resource checks.
use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::{One, Zero};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::{
    assumptions,
    expr::{Constant, Expr, Func, Kind, RelOp},
    logic::LogicOp,
    number::Number,
    ops,
    sets::SetOp,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CountError {
    NumericPower,
    MultiplicativeFactors,
    DistinctNodes,
    Cycle,
}
impl CountError {
    pub fn message(&self) -> &'static str {
        match self {
            Self::NumericPower => "Fraction normalization exceeds the numeric power resource limit",
            Self::MultiplicativeFactors => {
                "Fraction normalization exceeds 10000 multiplicative factors"
            }
            Self::DistinctNodes => "Operation counting exceeds 100000 distinct nodes",
            Self::Cycle => "Cyclic operation-count input",
        }
    }
}

#[derive(Default, Debug, PartialEq, Eq)]
pub struct OperationCounts {
    pub labels: BTreeMap<String, BigInt>,
}
impl OperationCounts {
    pub fn total(&self) -> BigInt {
        self.labels.values().sum()
    }
    pub fn visual(&self) -> Expr {
        ops::add(
            &self
                .labels
                .iter()
                .filter(|(_, n)| !n.is_zero())
                .map(|(name, n)| {
                    ops::mul2(
                        &Expr::num(Number::from_bigint(n.clone())),
                        &Expr::symbol(&name.to_uppercase()),
                    )
                })
                .collect::<Vec<_>>(),
        )
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Mode {
    Top,
    Expr,
    Basic,
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Key {
    value: Expr,
    mode: Mode,
}
#[derive(Default)]
struct Instruction {
    labels: Vec<(String, BigInt)>,
    children: Vec<Key>,
}
impl Instruction {
    fn label(&mut self, name: &str, amount: usize) {
        if amount != 0 {
            self.labels.push((name.to_owned(), BigInt::from(amount)));
        }
    }
    fn children(&mut self, children: impl IntoIterator<Item = Expr>, mode: Mode) {
        self.children
            .extend(children.into_iter().map(|value| Key { value, mode }));
    }
}

fn tuple(e: &Expr) -> bool {
    matches!(e.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple"))
}
fn basic_function(e: &Expr) -> bool {
    matches!(
        e.kind(),
        Kind::Func(Func::Predicate(_) | Func::AppliedPredicate, _)
    )
}
fn nonexpr_set(e: &Expr) -> bool {
    matches!(e.kind(),Kind::Func(Func::Set(op),_) if !matches!(op,SetOp::Lambda|SetOp::IdentityFunction|SetOp::Contains|SetOp::FiniteSet))
}
fn function(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Undefined(head), _) => !presentation_head(head),
        Kind::Func(
            Func::GenericExpr
            | Func::SymbolicDerivative
            | Func::Subs
            | Func::Limit
            | Func::Integral
            | Func::Sum
            | Func::Product
            | Func::Order
            | Func::Predicate(_)
            | Func::AppliedPredicate,
            _,
        ) => false,
        Kind::Func(Func::Set(op), _) => matches!(op, SetOp::Lambda),
        Kind::Func(_, _) | Kind::Piecewise(_) => true,
        _ => false,
    }
}

/// Constructor-only views of standard Python scalar frontends. These heads
/// are structural containers, rather than applied undefined functions.
pub fn presentation(name: &str, args: Vec<Expr>) -> Option<Expr> {
    if !matches!(name, "Indexed" | "IndexedBase" | "Idx" | "ExprCondPair") {
        return None;
    }
    let mut head =
        crate::expr::UndefinedFunction::new(name, crate::symbol::Assumptions::NONE, None);
    head.internal = true;
    Some(Expr::from_kind(Kind::Func(Func::Undefined(head), args)))
}
fn presentation_head(head: &crate::expr::UndefinedFunction) -> bool {
    head.internal
        && matches!(
            head.as_ref(),
            "Tuple" | "UnevaluatedExpr" | "Indexed" | "IndexedBase" | "Idx" | "ExprCondPair"
        )
}
fn indexed_view(identity: &Expr) -> Option<Expr> {
    if !tuple(identity) {
        return None;
    }
    let metadata = identity.args();
    if metadata.len() != 3 || !tuple(&metadata[2]) {
        return None;
    }
    let mut base = vec![metadata[0].clone()];
    if !matches!(metadata[1].kind(), Kind::Bool(false)) {
        base.push(metadata[1].clone());
    }
    let mut args = vec![presentation("IndexedBase", base)?];
    for entry in metadata[2].args() {
        if !tuple(&entry) {
            return None;
        }
        let fields = entry.args();
        if fields.len() == 2 && matches!(fields[0].kind(), Kind::Bool(false)) {
            args.push(fields[1].clone());
        } else if matches!(fields.first().map(Expr::kind), Some(Kind::Bool(true)))
            && (fields.len() == 2 || fields.len() == 4)
        {
            let mut index = vec![fields[1].clone()];
            if fields.len() == 4 {
                index.push(crate::functions::undefined("Tuple", fields[2..].to_vec()));
            }
            args.push(presentation("Idx", index)?);
        } else {
            return None;
        }
    }
    presentation("Indexed", args)
}
fn children(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        // These helpers construct payloads without evaluating a power. In
        // particular add_term also preserves the public tail of a held Mul.
        Kind::Add(..) | Kind::Mul(..) => crate::traverse::ordered_args(e),
        Kind::Func(Func::Heaviside, args) if args.len() == 1 => {
            vec![args[0].clone(), Expr::rational(1, 2)]
        }
        Kind::Func(Func::Order, args) if !args.is_empty() => {
            let mut result = vec![args[0].clone()];
            result.extend(
                args[1..]
                    .chunks(2)
                    .map(|pair| crate::functions::undefined("Tuple", pair.to_vec())),
            );
            result
        }
        Kind::Piecewise(branches) => branches
            .iter()
            .map(|(value, condition)| {
                presentation("ExprCondPair", vec![value.clone(), condition.clone()]).unwrap()
            })
            .collect(),
        _ => super::scalar_children(e),
    }
}
fn negative_atom(e: &Expr) -> bool {
    e.as_number().is_some_and(Number::is_negative) || e.is_constant(Constant::NegInfinity)
}
fn negative(e: &Expr) -> bool {
    if e.is_mul() {
        children(e).first().is_some_and(negative_atom)
    } else {
        negative_atom(e)
    }
}
fn held_product(values: &[Expr]) -> Expr {
    match values.len() {
        0 => Expr::one(),
        1 => values[0].clone(),
        _ => Expr::from_unevaluated_kind(Kind::Mul(
            Number::ONE,
            values.iter().cloned().map(|v| (v, Expr::one())).collect(),
        )),
    }
}
fn positive_part(e: &Expr) -> Expr {
    if e.is_constant(Constant::NegInfinity) {
        return Expr::oo();
    }
    match e.kind() {
        Kind::Num(n) if n.is_negative() => Expr::num(n.abs()),
        Kind::Mul(_, _) => {
            let mut values = children(e);
            if values.first().is_some_and(negative_atom) {
                let abs = positive_part(&values[0]);
                if abs.is_one() {
                    values.remove(0);
                } else {
                    values[0] = abs;
                }
                held_product(&values)
            } else {
                e.clone()
            }
        }
        _ => e.clone(),
    }
}
// Mul.__neg__ extracts only a rational leading coefficient, then preserves
// the literal tail. In particular it must not normalize every numeric factor.
fn negated(e: &Expr) -> Expr {
    if e.is_constant(Constant::NegInfinity) {
        return Expr::oo();
    }
    if let Some(number) = e.as_number() {
        return Expr::num(number.neg());
    }
    if !e.is_mul() {
        return ops::neg(e);
    }
    let values = children(e);
    let mut coefficient = Expr::one();
    let mut tail = values.clone();
    if let Some(first) = values.first() {
        if first.is_rational_number() {
            coefficient = first.clone();
            tail = values[1..].to_vec();
        } else if negative_atom(first) {
            coefficient = Expr::minus_one();
            tail[0] = negated(first);
        }
    }
    if !tail
        .first()
        .is_some_and(|value| value.is_constant(Constant::ComplexInfinity))
    {
        coefficient = ops::neg(&coefficient);
    }
    if !coefficient.is_one() {
        if tail.first().is_some_and(|value| {
            value.as_number().is_some()
                || matches!(
                    value.kind(),
                    Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::NaN)
                )
        }) {
            tail[0] = if coefficient.is_minus_one() {
                negated(&tail[0])
            } else {
                ops::mul2(&coefficient, &tail[0])
            };
        } else {
            tail.insert(0, coefficient);
        }
    }
    held_product(&tail)
}

fn denominator_possible(e: &Expr) -> bool {
    if let Some(n) = e.as_number().filter(|n| n.is_rational()) {
        return !n.is_integer();
    }
    let exponent = match e.kind() {
        Kind::Pow(_, exponent) => exponent,
        Kind::Func(Func::Exp, args) if args.len() == 1 => &args[0],
        _ => return false,
    };
    assumptions::is_negative(exponent) == Some(true)
        || assumptions::is_positive(exponent) != Some(true) && exponent.is_mul()
}

// Only the extra numeric-base normalization below can materialize arithmetic
// hidden inside a held power. Bound that optimization, keeping the original
// symbolic power when its integer storage or depth would be excessive.
fn bounded_numeric_base(root: &Expr) -> bool {
    const BITS: u64 = 1 << 20;
    let mut memo: FxHashMap<Expr, Option<(u64, u16)>> = FxHashMap::default();
    let mut pending = vec![(root.clone(), false)];
    while let Some((value, visited)) = pending.pop() {
        if memo.contains_key(&value) {
            continue;
        }
        if memo.len() + pending.len() > 10000 {
            return false;
        }
        let args = match value.kind() {
            Kind::Mul(..) => children(&value),
            Kind::Pow(base, exponent) if exponent.is_integer_number() => vec![base.clone()],
            _ => Vec::new(),
        };
        if !args.is_empty() && !visited {
            pending.push((value, true));
            pending.extend(args.into_iter().map(|arg| (arg, false)));
            continue;
        }
        let estimate = match value.kind() {
            Kind::Num(number) => Some((
                match number {
                    Number::Int(n) => u64::from(64 - n.unsigned_abs().leading_zeros()).max(1),
                    Number::BigInt(n) => n.bits().max(1),
                    Number::Rational(q) => q.numer().bits().saturating_add(q.denom().bits()).max(1),
                    Number::Float(_) => 53,
                    Number::MpFloat(_, precision) => u64::from(*precision),
                },
                1,
            )),
            Kind::Const(Constant::I) => Some((1, 1)),
            Kind::Mul(..) => args.iter().try_fold((0u64, 1u16), |(bits, depth), arg| {
                let (child_bits, child_depth) = memo[arg]?;
                Some((
                    bits.saturating_add(child_bits),
                    depth.max(child_depth.saturating_add(1)),
                ))
            }),
            Kind::Pow(_, exponent) if exponent.is_integer_number() => {
                match (args.first().and_then(|arg| memo[arg]), exponent.as_i64()) {
                    (Some((bits, depth)), Some(exponent)) => Some((
                        bits.saturating_mul(exponent.unsigned_abs()).max(1),
                        depth.saturating_add(1),
                    )),
                    _ => None,
                }
            }
            _ => None,
        }
        .filter(|(bits, depth)| *bits <= BITS && *depth <= 128);
        memo.insert(value, estimate);
    }
    memo[root].is_some()
}

fn normalized_product(values: &[Expr]) -> Expr {
    // AssocOp removes its identity before deciding whether a single held
    // argument can be returned unchanged.
    let values = values
        .iter()
        .filter(|value| !value.is_one())
        .cloned()
        .collect::<Vec<_>>();
    if values.len() < 2 {
        return held_product(&values);
    }
    let mut factors = Vec::new();
    let mut pending = values.into_iter().rev().collect::<Vec<_>>();
    while let Some(value) = pending.pop() {
        if value.is_mul() {
            pending.extend(children(&value).into_iter().rev());
        } else if let Kind::Pow(base, exponent) = value.kind() {
            if base.is_mul() && exponent.is_integer_number() && bounded_numeric_base(&value) {
                factors.push(ops::pow(&normalized_product(&children(base)), exponent));
            } else {
                factors.push(value);
            }
        } else {
            factors.push(value);
        }
    }
    ops::mul(&factors)
}

// fraction(exact=False) rebuilds both parts with an evaluated Mul. Keeping a
// raw nested numeric product here changes DIV/NEG counts after normalization.
fn fraction(value: &Expr) -> (Expr, Expr) {
    let terms = if value.is_mul() {
        children(value)
    } else {
        vec![value.clone()]
    };
    let mut numerator = Vec::new();
    let mut denominator = Vec::new();
    for term in terms {
        let powlike = term.is_pow() || matches!(term.kind(), Kind::Func(Func::Exp, _));
        if powlike && crate::decompose::commutativity(&term) == Some(true) {
            let (base, exponent) = crate::decompose::as_base_exp(&term);
            if assumptions::is_negative(&exponent) == Some(true) {
                let positive = ops::neg(&exponent);
                if exponent.is_minus_one() {
                    denominator.push(base);
                } else {
                    let base = if base.is_mul()
                        && positive.is_integer_number()
                        && bounded_numeric_base(&Expr::from_unevaluated_kind(Kind::Pow(
                            base.clone(),
                            positive.clone(),
                        ))) {
                        normalized_product(&children(&base))
                    } else {
                        base
                    };
                    denominator.push(ops::pow(&base, &positive));
                }
            } else if assumptions::is_positive(&exponent) == Some(true) {
                numerator.push(term);
            } else if exponent.is_mul() {
                let (n, d) = crate::simptools::as_numer_denom(&term);
                if !n.is_one() {
                    numerator.push(n);
                }
                denominator.push(d);
            } else {
                numerator.push(term);
            }
        } else if let Some(number) = term
            .as_number()
            .filter(|n| n.is_rational() && !n.is_integer())
        {
            let (n, d) = number.numer_denom().unwrap();
            if n != BigInt::one() {
                numerator.push(Expr::num(Number::from_bigint(n)));
            }
            denominator.push(Expr::num(Number::from_bigint(d)));
        } else {
            numerator.push(term);
        }
    }
    (
        normalized_product(&numerator),
        normalized_product(&denominator),
    )
}

#[derive(Default)]
struct Builder {
    product_sizes: FxHashMap<Expr, usize>,
}
impl Builder {
    fn multiplication_size(&mut self, e: &Expr) -> Result<usize, CountError> {
        let mut pending = vec![(e.clone(), false)];
        while let Some((node, visited)) = pending.pop() {
            if self.product_sizes.contains_key(&node) {
                continue;
            }
            if let Kind::Pow(base, exponent) = node.kind() {
                if crate::functions::is_numeric_value(base)
                    && exponent.is_integer_number()
                    && exponent.as_i64().is_none_or(|n| n.unsigned_abs() > 100000)
                {
                    return Err(CountError::NumericPower);
                }
            }
            if !node.is_mul() {
                self.product_sizes.insert(node, 1);
            } else if visited {
                let size = children(&node).iter().fold(0usize, |sum, arg| {
                    sum.saturating_add(self.product_sizes[arg]).min(10001)
                });
                self.product_sizes.insert(node, size);
            } else {
                pending.push((node.clone(), true));
                pending.extend(
                    children(&node)
                        .into_iter()
                        .filter(|arg| !self.product_sizes.contains_key(arg))
                        .map(|arg| (arg, false)),
                );
            }
        }
        Ok(self.product_sizes[e])
    }

    // None asks the adapter to retain the external Basic/container protocol.
    fn instruction(&mut self, key: &Key) -> Result<Option<Instruction>, CountError> {
        let mut instruction = Instruction::default();
        let mut value = if let Some(identity) = key
            .value
            .as_sym()
            .and_then(|symbol| symbol.scalar_identity())
        {
            let Some(value) = indexed_view(identity) else {
                return Ok(None);
            };
            value
        } else {
            key.value.clone()
        };
        let mut mode = key.mode;
        if mode == Mode::Top {
            match value.kind() {
                Kind::Rel(relation, a, b) => {
                    let label = match relation {
                        RelOp::Eq => "EQ",
                        RelOp::Ne => "NE",
                        RelOp::Lt => "LT",
                        RelOp::Le => "LE",
                        RelOp::Gt => "GT",
                        RelOp::Ge => "GE",
                    };
                    instruction.label(label, 1);
                    instruction.children([a.clone(), b.clone()], Mode::Top);
                    return Ok(Some(instruction));
                }
                Kind::Func(Func::Logic(op), args) => {
                    instruction.label(&op.name().to_uppercase(), 1);
                    instruction.children(args.clone(), Mode::Top);
                    return Ok(Some(instruction));
                }
                Kind::Func(Func::Set(SetOp::FiniteSet), args) => {
                    instruction.children(args.clone(), Mode::Top);
                    return Ok(Some(instruction));
                }
                _ if tuple(&value) => {
                    instruction.children(children(&value), Mode::Top);
                    return Ok(Some(instruction));
                }
                _ if nonexpr_set(&value) => return Ok(None),
                _ if basic_function(&value)
                    || matches!(value.kind(), Kind::Func(Func::Set(SetOp::Contains), _)) =>
                {
                    mode = Mode::Basic
                }
                _ => mode = Mode::Expr,
            }
        }
        if mode == Mode::Basic {
            let args = children(&value);
            if !args.is_empty() {
                let count = if matches!(value.kind(), Kind::Func(Func::Logic(_), _)) {
                    args.len().saturating_sub(1)
                } else {
                    1
                };
                instruction.label(&value.class_name().to_uppercase(), count);
                instruction.children(args, Mode::Basic);
            }
            return Ok(Some(instruction));
        }
        if let Some(number) = value.as_number().filter(|n| n.is_rational()) {
            instruction.label("NEG", usize::from(number.is_negative()));
            instruction.label("DIV", usize::from(!number.is_integer()));
            return Ok(Some(instruction));
        }
        if value.is_mul() {
            if negative(&value) {
                instruction.label("NEG", 1);
                // SymPy drops an explicit leading -1 without evaluating the
                // remaining held product; other coefficients use Mul.__neg__ coefficient extraction.
                value = if children(&value).first().is_some_and(Expr::is_minus_one) {
                    positive_part(&value)
                } else {
                    negated(&value)
                };
            }
            let potential = if value.is_mul() {
                children(&value).iter().any(denominator_possible)
            } else {
                denominator_possible(&value)
            };
            let normalize = if potential {
                if self.multiplication_size(&value)? > 10000 {
                    return Err(CountError::MultiplicativeFactors);
                }
                true
            } else {
                // SymPy always rebuilds fraction parts, even with no explicit
                // denominator. Keep existing support for enormous held powers
                // and shared products by inspecting their original shape when
                // that reconstruction would exceed the normalization budget.
                matches!(self.multiplication_size(&value),Ok(size) if size<=10000)
            };
            let (numerator, denominator) = if normalize {
                fraction(&value)
            } else {
                (value.clone(), Expr::one())
            };
            if numerator.is_integer_number() {
                instruction.label("DIV", 1);
                instruction.label("NEG", usize::from(negative(&numerator)));
                instruction.children([denominator], Mode::Expr);
                return Ok(Some(instruction));
            }
            if !denominator.is_one() {
                instruction.label("DIV", 1);
                instruction.children([numerator], Mode::Expr);
                if !denominator.is_integer_number() {
                    instruction.children([denominator], Mode::Expr);
                }
                return Ok(Some(instruction));
            }
        } else if value.is_add() {
            let args = children(&value);
            let negatives = args.iter().filter(|arg| negative(arg)).count();
            if negatives == args.len() {
                instruction.label("NEG", 1);
                instruction.label("SUB", args.len().saturating_sub(1));
            } else {
                instruction.label("SUB", negatives);
                instruction.label("ADD", args.len() - negatives - 1);
            }
            instruction.children(
                args.into_iter()
                    .map(|arg| if negative(&arg) { negated(&arg) } else { arg }),
                Mode::Expr,
            );
            return Ok(Some(instruction));
        }
        if let Kind::Pow(base, exponent) = value.kind() {
            if exponent.is_minus_one() {
                instruction.label("DIV", 1);
                instruction.children([base.clone()], Mode::Expr);
                return Ok(Some(instruction));
            }
            if base.is_constant(Constant::E) {
                instruction.label("EXP", 1);
                instruction.children([exponent.clone()], Mode::Expr);
                return Ok(Some(instruction));
            }
        }
        if value.is_constant(Constant::E) {
            instruction.label("EXP", 1);
            return Ok(Some(instruction));
        }
        let args = children(&value);
        if value.is_mul()
            || matches!(
                value.kind(),
                Kind::Func(
                    Func::Logic(LogicOp::And | LogicOp::Or) | Func::Max | Func::Min,
                    _
                )
            )
        {
            instruction.label(
                &value.class_name().to_uppercase(),
                args.len().saturating_sub(1),
            );
        } else if !args.is_empty()
            && (value.is_pow()
                || function(&value)
                || matches!(
                    value.kind(),
                    Kind::Derivative(..)
                        | Kind::Func(Func::Integral | Func::Sum | Func::SymbolicDerivative, _)
                        | Kind::Piecewise(_)
                ))
        {
            let label = match value.kind() {
                Kind::Func(Func::Undefined(head), _) if !head.internal => {
                    format!("FUNC_{}", head.as_ref())
                }
                _ => value.class_name().to_uppercase(),
            };
            instruction.label(&label, 1);
        }
        if !value.is_symbol() {
            instruction.children(args, Mode::Expr);
        }
        Ok(Some(instruction))
    }
}

/// None retains host/matrix and external container protocols at the API layer.
/// A successful result counts tree occurrences without unfolding a shared DAG.
pub fn count(e: &Expr) -> Result<Option<OperationCounts>, CountError> {
    if e.has_host() || e.has_matrix() {
        return Ok(None);
    }
    let root = Key {
        value: e.clone(),
        mode: Mode::Top,
    };
    let mut builder = Builder::default();
    let mut instructions = FxHashMap::default();
    let mut active = FxHashSet::default();
    let mut order = Vec::new();
    let mut pending = vec![(root.clone(), false)];
    while let Some((key, visited)) = pending.pop() {
        if visited {
            active.remove(&key);
            order.push(key);
            continue;
        }
        if active.contains(&key) {
            return Err(CountError::Cycle);
        }
        if key.value.has_host() || key.value.has_matrix() {
            return Ok(None);
        }
        if instructions.contains_key(&key) {
            continue;
        }
        let Some(instruction) = builder.instruction(&key)? else {
            return Ok(None);
        };
        let children = instruction.children.clone();
        instructions.insert(key.clone(), instruction);
        if instructions.len() > 100000 {
            return Err(CountError::DistinctNodes);
        }
        active.insert(key.clone());
        pending.push((key, true));
        pending.extend(children.into_iter().rev().map(|child| (child, false)));
    }
    let mut weights = FxHashMap::default();
    weights.insert(root, BigInt::one());
    let mut result = OperationCounts::default();
    for key in order.into_iter().rev() {
        let weight = weights.remove(&key).unwrap_or_else(BigInt::zero);
        let instruction = &instructions[&key];
        for (label, number) in &instruction.labels {
            *result.labels.entry(label.clone()).or_default() += &weight * number;
        }
        for child in &instruction.children {
            *weights.entry(child.clone()).or_default() += &weight;
        }
    }
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_graph_weights_are_unbounded() {
        let mut e = Expr::from_unevaluated_kind(Kind::Mul(
            Number::Int(2),
            vec![(Expr::symbol("x"), Expr::one())],
        ));
        for _ in 0..80 {
            e = Expr::from_unevaluated_kind(Kind::Mul(
                Number::ONE,
                vec![(e.clone(), Expr::one()), (e, Expr::one())],
            ));
        }
        let counts = count(&e).unwrap().unwrap();
        assert_eq!(counts.total(), (BigInt::one() << 81) - 1);
        assert_eq!(counts.labels["MUL"], counts.total());
    }
    #[test]
    fn numeric_reconstruction_keeps_large_towers_symbolic() {
        let mut value = Expr::int(2);
        for _ in 0..200 {
            let base = Expr::from_unevaluated_kind(Kind::Mul(
                Number::ONE,
                vec![(Expr::one(), Expr::one()), (value, Expr::one())],
            ));
            value = Expr::from_unevaluated_kind(Kind::Pow(base, Expr::int(2)));
        }
        assert!(!bounded_numeric_base(&value));
    }

    #[test]
    fn undefined_heads_do_not_collide_with_builtins() {
        for name in ["Order", "Tuple", "FiniteSet", "Max", "And"] {
            let head =
                crate::expr::UndefinedFunction::new(name, crate::symbol::Assumptions::NONE, None);
            let e = Expr::from_kind(Kind::Func(Func::Undefined(head), vec![Expr::symbol("x")]));
            let counts = count(&e).unwrap().unwrap();
            assert_eq!(
                counts.labels.get(&format!("FUNC_{}", name)),
                Some(&BigInt::one())
            );
            assert_eq!(counts.total(), BigInt::one());
        }
    }
}
