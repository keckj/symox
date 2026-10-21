//! Targeted simplification and expansion routines: `trigsimp`,
//! `expand_trig`, `powsimp`, `expand_log`, `logcombine`, `expand_mul`,
//! `expand_power_base`, `expand_power_exp`, `nsimplify`.

use rustc_hash::FxHashMap;

use crate::assumptions::is_positive;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::*;
use crate::number::Number;
use crate::ops::{add, add2, div, mul, mul2, mul_factor, neg, pow, sqrt, sub};
use crate::parser::parse_str;
use crate::subs::xreplace;
use crate::traverse::{map_children, node_count};

// ----------------------------------------------------------------------
// expand_trig
// ----------------------------------------------------------------------

/// Chebyshev polynomial T_n(t) (first kind).
fn cheb_t(n: u32, t: &Expr) -> Expr {
    let (mut a, mut b) = (Expr::one(), t.clone());
    if n == 0 {
        return a;
    }
    for _ in 1..n {
        let c = expand(&sub(&mul(&[Expr::int(2), t.clone(), b.clone()]), &a));
        a = b;
        b = c;
    }
    b
}

/// Chebyshev polynomial U_n(t) (second kind).
fn cheb_u(n: u32, t: &Expr) -> Expr {
    let (mut a, mut b) = (Expr::one(), mul2(&Expr::int(2), t));
    if n == 0 {
        return a;
    }
    for _ in 1..n {
        let c = expand(&sub(&mul(&[Expr::int(2), t.clone(), b.clone()]), &a));
        a = b;
        b = c;
    }
    b
}

/// Expand trigonometric functions of sums and integer multiples.
pub fn expand_trig(e: &Expr) -> Expr { expand_trig_deep(e,true) }

pub fn expand_trig_deep(e: &Expr,deep:bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut expand_trig) };
    let Kind::Func(f, args) = e2.kind() else { return e2 };
    use crate::special::SpecialFn;
    // the reciprocal functions expand as the reciprocal of their partner
    // (the reference `_calculate_reciprocal`)
    let partner = match f {
        Func::Sec => Some(cos(&args[0])),
        Func::Csc => Some(sin(&args[0])),
        Func::Special(SpecialFn::Sech) => Some(cosh(&args[0])),
        Func::Special(SpecialFn::Csch) => Some(sinh(&args[0])),
        _ => None,
    };
    if let Some(partner) = partner {
        return pow(&expand_trig_deep(&partner, false), &Expr::minus_one());
    }
    match f {
        Func::Tan | Func::Cot | Func::Tanh | Func::Coth => {
            return expand_trig_ratio(&e2, f, &args[0]);
        }
        Func::Sinh | Func::Cosh => return expand_trig_hyperbolic(&e2, f, &args[0]),
        Func::Sin | Func::Cos => {}
        _ => return e2,
    }
    let u = &args[0];
    match u.kind() {
        Kind::Add(..) => {
            let terms = crate::sortkey::ordered_terms(u);
            let (a, rest) = (terms[0].clone(), add(&terms[1..]));
            let (sa, ca, sb, cb) = (
                expand_trig(&sin(&a)),
                expand_trig(&cos(&a)),
                expand_trig(&sin(&rest)),
                expand_trig(&cos(&rest)),
            );
            match f {
                Func::Sin => add2(&mul2(&sa, &cb), &mul2(&ca, &sb)),
                _ => sub(&mul2(&ca, &cb), &mul2(&sa, &sb)),
            }
        }
        Kind::Mul(c, facs) if c.is_integer() && c.is_positive() && !c.is_one() => {
            let n = c.to_i64().unwrap_or(1);
            if n > 64 {
                return e2;
            }
            let n = n as u32;
            let x = crate::ops::add_term(&mul_from(facs), &Number::ONE);
            match f {
                Func::Sin => {
                    let s = sin(&x);
                    if n % 2 == 1 {
                        let sign = if ((n - 1) / 2) % 2 == 0 { 1 } else { -1 };
                        expand(&mul2(&Expr::int(sign), &cheb_t(n, &s)))
                    } else {
                        let sign = if (n / 2 - 1) % 2 == 0 { 1 } else { -1 };
                        expand(&mul(&[Expr::int(sign), cos(&x), cheb_u(n - 1, &s)]))
                    }
                }
                _ => expand(&cheb_t(n, &cos(&x))),
            }
        }
        _ => e2,
    }
}

/// `expr.expand(trig=True)` with the default hints, as the hyperbolic
/// expansions of the reference call it on their sum formulas.
fn expand_trig_full(e: &Expr) -> Expr {
    let options = crate::expand::ExpandOptions { trig: true, ..Default::default() };
    crate::expand::expand_with_options(e, &options, &mut Default::default())
        .unwrap_or_else(|_| e.clone())
}

/// The reference `sinh._eval_expand_trig` and `cosh._eval_expand_trig`: a
/// sum `x + y` (its first term and the rest) or `n*x` as `x + (n - 1)*x`
/// through the sum formula, expanded again with `trig=True`.
fn expand_trig_hyperbolic(e: &Expr, f: &Func, u: &Expr) -> Expr {
    let (a, b) = match u.kind() {
        Kind::Add(..) => {
            let terms = crate::sortkey::ordered_terms(u);
            (terms[0].clone(), add(&terms[1..]))
        }
        Kind::Mul(c, facs) if c.is_integer() && !c.is_one() => {
            let x = crate::ops::add_term(&mul_from(facs), &Number::ONE);
            let n = Expr::num(c.clone());
            (x.clone(), mul2(&sub(&n, &Expr::one()), &x))
        }
        _ => return e.clone(),
    };
    let r = match f {
        Func::Sinh => add2(&mul2(&sinh(&a), &cosh(&b)), &mul2(&sinh(&b), &cosh(&a))),
        _ => add2(&mul2(&cosh(&a), &cosh(&b)), &mul2(&sinh(&a), &sinh(&b))),
    };
    expand_trig_full(&r)
}

/// Elementary symmetric polynomial `e_k` of `values`, as the reference
/// `symmetric_poly` builds it: the sum of the products of the `k`-subsets.
fn elementary_symmetric(k: usize, values: &[Expr]) -> Expr {
    if k == 0 {
        return Expr::one();
    }
    let n = values.len();
    if k > n {
        return Expr::zero();
    }
    let mut terms = Vec::new();
    let mut index: Vec<usize> = (0..k).collect();
    loop {
        terms.push(mul(&index.iter().map(|&i| values[i].clone()).collect::<Vec<_>>()));
        let mut i = k;
        while i > 0 && index[i - 1] == n - k + i - 1 {
            i -= 1;
        }
        if i == 0 {
            break;
        }
        index[i - 1] += 1;
        for j in i..k {
            index[j] = index[j - 1] + 1;
        }
    }
    add(&terms)
}

fn binomial_expr(n: u64, k: u64) -> Expr {
    let mut value = num_bigint::BigInt::from(1u32);
    for i in 0..k {
        value = value * (n - i) / (i + 1);
    }
    Expr::num(Number::from_bigint(value))
}

/// The reference `_eval_expand_trig` of `tan`, `cot`, `tanh` and `coth`: a
/// sum through the symmetric functions of the expanded terms, an integer
/// multiple through the binomial coefficients.
fn expand_trig_ratio(e: &Expr, f: &Func, u: &Expr) -> Expr {
    let own = |a: &Expr| match f {
        Func::Tan => tan(a),
        Func::Cot => cot(a),
        Func::Tanh => tanh(a),
        _ => coth(a),
    };
    let sign = |i: usize| if (i % 4) / 2 == 0 { Expr::one() } else { Expr::minus_one() };
    match u.kind() {
        Kind::Add(..) => {
            let values: Vec<Expr> =
                u.args().iter().map(|a| expand_trig_deep(&own(a), false)).collect();
            let n = values.len();
            let (mut num, mut den) = (Vec::new(), Vec::new());
            for i in 0..=n {
                let e_i = elementary_symmetric(i, &values);
                match f {
                    // p[1 - i % 2] += (-1)**((i % 4)//2)*e_i
                    Func::Tan => {
                        let term = mul2(&sign(i), &e_i);
                        if i % 2 == 1 { num.push(term) } else { den.push(term) }
                    }
                    // p[i % 2] += e_i, the numerator at odd i
                    Func::Tanh => {
                        if i % 2 == 1 { num.push(e_i) } else { den.push(e_i) }
                    }
                    // p[(n - i) % 2] += (-1)**(((n - i) % 4)//2)*e_i
                    Func::Cot => {
                        let term = mul2(&sign(n - i), &e_i);
                        if (n - i) % 2 == 0 { num.push(term) } else { den.push(term) }
                    }
                    _ => {
                        if (n - i) % 2 == 0 { num.push(e_i) } else { den.push(e_i) }
                    }
                }
            }
            div(&add(&num), &add(&den))
        }
        Kind::Mul(c, facs) if c.is_integer() && c.is_positive() && !c.is_one() => {
            let Some(n) = c.to_i64().filter(|n| *n <= 64) else { return e.clone() };
            let n = n as u64;
            let x = crate::ops::add_term(&mul_from(facs), &Number::ONE);
            let t = own(&x);
            let (mut num, mut den) = (Vec::new(), Vec::new());
            for k in 0..=n {
                let term = mul2(&binomial_expr(n, k), &pow(&t, &Expr::int(k as i64)));
                match f {
                    // im((1 + I*t)**n)/re((1 + I*t)**n)
                    Func::Tan => {
                        let term = mul2(&sign(k as usize), &term);
                        if k % 2 == 1 { num.push(term) } else { den.push(term) }
                    }
                    // re((t + I)**n)/im((t + I)**n)
                    Func::Cot => {
                        let term = mul2(&sign((n - k) as usize), &term);
                        if (n - k) % 2 == 0 { num.push(term) } else { den.push(term) }
                    }
                    Func::Tanh => {
                        if k % 2 == 1 { num.push(term) } else { den.push(term) }
                    }
                    _ => {
                        if (n - k) % 2 == 0 { num.push(term) } else { den.push(term) }
                    }
                }
            }
            div(&add(&num), &add(&den))
        }
        _ => e.clone(),
    }
}

fn mul_from(facs: &[(Expr, Expr)]) -> Expr {
    mul(&facs.iter().map(|(b, x)| mul_factor(b, x)).collect::<Vec<_>>())
}

// ----------------------------------------------------------------------
// trigsimp
// ----------------------------------------------------------------------

struct Rule {
    pat: Expr,
    rep: Expr,
    /// Trigonometric heads of `pat` (see `head_bit`).
    heads: u8,
}

fn wilds() -> FxHashMap<&'static str, Expr> {
    let mut m = FxHashMap::default();
    for n in ["a", "b", "c", "d", "n"] {
        m.insert(n, Expr::wild(n, vec![]));
    }
    m
}

fn rule(p: &str, r: &str, w: &FxHashMap<&'static str, Expr>) -> Rule {
    let sub_w = |s: &str| -> Expr {
        let e = parse_str(s).expect("rule parses");
        let map: FxHashMap<Expr, Expr> = w.iter().map(|(k, v)| (Expr::symbol(k), v.clone())).collect();
        xreplace(&e, &map)
    };
    let pat = sub_w(p);
    Rule { heads: heads(&pat), pat, rep: sub_w(r) }
}

fn trig_rules() -> &'static Vec<Rule> {
    use std::sync::OnceLock;
    static RULES: OnceLock<Vec<Rule>> = OnceLock::new();
    RULES.get_or_init(|| {
        let w = wilds();
        vec![
            // Pythagorean identities
            rule("a*sin(b)**2 + a*cos(b)**2 + c", "a + c", &w),
            rule("a*cosh(b)**2 - a*sinh(b)**2 + c", "a + c", &w),
            rule("a*tan(b)**2 + a + c", "a/cos(b)**2 + c", &w),
            // double angle
            rule("2*a*sin(b)*cos(b) + c", "a*sin(2*b) + c", &w),
            rule("a*cos(b)**2 - a*sin(b)**2 + c", "a*cos(2*b) + c", &w),
            rule("2*a*cos(b)**2 - a + c", "a*cos(2*b) + c", &w),
            rule("a - 2*a*sin(b)**2 + c", "a*cos(2*b) + c", &w),
            // complements
            rule("a - a*sin(b)**2 + c", "a*cos(b)**2 + c", &w),
            rule("a - a*cos(b)**2 + c", "a*sin(b)**2 + c", &w),
            rule("a*cosh(b)**2 - a + c", "a*sinh(b)**2 + c", &w),
            // sum formulas
            rule("a*sin(b)*cos(d) + a*sin(d)*cos(b) + c", "a*sin(b + d) + c", &w),
            rule("a*cos(b)*cos(d) - a*sin(b)*sin(d) + c", "a*cos(b + d) + c", &w),
            rule("a*cos(b)*cos(d) + a*sin(b)*sin(d) + c", "a*cos(b - d) + c", &w),
            // quotients
            rule("a*sin(b)**n*cos(b)**(-n)", "a*tan(b)**n", &w),
            rule("a*sinh(b)**n*cosh(b)**(-n)", "a*tanh(b)**n", &w),
        ]
    })
}

/// Bit of a trigonometric head in the `heads` masks below.
fn head_bit(e: &Expr) -> u8 {
    match e.kind() {
        Kind::Func(Func::Sin, _) => 1,
        Kind::Func(Func::Cos, _) => 2,
        Kind::Func(Func::Tan, _) => 4,
        Kind::Func(Func::Sinh, _) => 8,
        Kind::Func(Func::Cosh, _) => 16,
        _ => 0,
    }
}

/// Trigonometric heads occurring anywhere in `e`.
fn heads(e: &Expr) -> u8 {
    let mut h = 0;
    let mut pending = vec![e.clone()];
    let mut seen = FxHashMap::default();
    while let Some(node) = pending.pop() {
        if seen.insert(node.dag_id(), node.clone()).is_some() { continue; }
        h |= head_bit(&node);
        pending.extend(node.args());
    }
    h
}

/// Memo tables shared by every pass of one trigsimp call. Each rewrite is a
/// pure function of the rule and the expression, so results computed in one
/// pass (or for a sub-expression in a nested call) stay valid for the
/// whole search: a node that reappears unchanged is not matched again.
/// Keys compare structurally; a hit is used only when the stored key has
/// the same presentation (held arithmetic and bound variables compare equal
/// to differently presented nodes).
#[derive(Default)]
pub(crate) struct TrigMemo {
    /// (rule index, expression) -> (rewritten expression, its heads).
    rewrite: FxHashMap<(usize, Expr), (Expr, u8)>,
    /// apply_rules fixpoints.
    rules: FxHashMap<Expr, Expr>,
    /// trigsimp results.
    trigsimp: FxHashMap<Expr, Expr>,
    /// Node id -> (node, tree size, height), see `TrigMemo::node_count`.
    sizes: NodeCache<(usize, usize)>,
    /// Node id -> (node, trigonometric heads), see `TrigMemo::heads`.
    head_sets: NodeCache<u8>,
}

/// Values attached to live nodes by identity. The weak reference keeps the
/// node's allocation, so its address cannot be reused by another node while
/// the entry exists, but not its children: the many rewrites tried and
/// rejected on a large sum are freed, and their entries are pruned when the
/// table has doubled. Atoms and nodes without another owner are not stored:
/// `args()` builds the terms and powers of sums and products afresh, and
/// those copies are never looked up again.
struct NodeCache<T> {
    entries: FxHashMap<usize, (crate::expr::WeakExpr, T)>,
    prune_at: usize,
}

impl<T> Default for NodeCache<T> {
    fn default() -> Self {
        NodeCache { entries: FxHashMap::default(), prune_at: 1 << 16 }
    }
}

impl<T: Copy> NodeCache<T> {
    fn get(&self, e: &Expr) -> Option<T> {
        self.entries.get(&e.dag_id()).filter(|(node, _)| node.is(e)).map(|(_, value)| *value)
    }

    fn insert(&mut self, e: &Expr, value: T) {
        if self.entries.len() >= self.prune_at {
            self.entries.retain(|_, (node, _)| node.upgrade().is_some());
            self.prune_at = (2 * self.entries.len()).max(1 << 16);
        }
        self.entries.insert(e.dag_id(), (e.downgrade(), value));
    }
}

/// Height bound of `traverse::node_count`: a node at this depth counts as
/// `usize::MAX`.
const NODE_COUNT_DEPTH: usize = 512;

impl TrigMemo {
    /// `node_count(e)`, reusing the sizes of the sub-expressions seen
    /// earlier in the search: a rewritten sum shares almost all of its terms
    /// with the previous one, so only its new nodes are visited. Each node
    /// keeps its height, and the shared sizes are used only when no node
    /// reaches node_count's depth bound; otherwise node_count itself runs,
    /// so the result is always node_count's.
    fn node_count(&mut self, e: &Expr) -> usize {
        fn count(e: &Expr, depth: usize, sizes: &mut NodeCache<(usize, usize)>) -> Option<(usize, usize)> {
            if e.is_atom() {
                return (depth < NODE_COUNT_DEPTH).then_some((1, 0));
            }
            if let Some((size, height)) = sizes.get(e) {
                return (depth + height < NODE_COUNT_DEPTH).then_some((size, height));
            }
            if depth >= NODE_COUNT_DEPTH {
                return None;
            }
            let mut height = 0;
            let mut child = |c: &Expr, size: usize, sizes: &mut NodeCache<(usize, usize)>| {
                let (s, h) = count(c, depth + 1, sizes)?;
                height = height.max(h + 1);
                Some(size.saturating_add(s))
            };
            let size = if let Kind::Derivative(body, variables) = e.kind() {
                let mut size = child(body, 1, sizes)?;
                for (variable, _) in variables {
                    size = child(variable, size, sizes)?.saturating_add(2);
                }
                size
            } else {
                let mut size = 1usize;
                for c in e.args().iter() {
                    size = child(c, size, sizes)?;
                }
                size
            };
            if e.strong_count() > 1 {
                sizes.insert(e, (size, height));
            }
            Some((size, height))
        }
        match count(e, 0, &mut self.sizes) {
            Some((size, _)) => size,
            None => node_count(e),
        }
    }

    /// `heads(e)`, reusing the heads of the sub-expressions seen earlier in
    /// the search (the heads of a node are its own and those of its
    /// children).
    fn heads(&mut self, e: &Expr) -> u8 {
        fn collect(e: &Expr, depth: usize, head_sets: &mut NodeCache<u8>) -> u8 {
            if e.is_atom() {
                return 0;
            }
            if let Some(h) = head_sets.get(e) {
                return h;
            }
            // Deep chains go through the iterative walk instead of the stack.
            let h = if depth >= NODE_COUNT_DEPTH {
                heads(e)
            } else {
                e.args().iter().fold(head_bit(e), |h, c| h | collect(c, depth + 1, head_sets))
            };
            if e.strong_count() > 1 {
                head_sets.insert(e, h);
            }
            h
        }
        collect(e, 0, &mut self.head_sets)
    }
}

/// `replace_with` for one trig rule, skipping the sub-expressions that lack
/// one of the rule's function heads: every function of a pattern has to
/// match a function of the same head inside the expression, so those nodes
/// cannot match and the result is unchanged. Returns the heads of the result.
fn rewrite_cached(e: &Expr, index: usize, r: &Rule, memo: &mut TrigMemo) -> (Expr, u8) {
    let key = (index, e.clone());
    let reusable = match memo.rewrite.get_key_value(&key) {
        Some(((_, stored), (result, heads))) if stored.same_presentation(e) => return (result.clone(), *heads),
        Some(_) => false,
        None => true,
    };
    let result = rewrite_uncached(e, index, r, memo);
    if reusable {
        memo.rewrite.insert(key, result.clone());
    }
    result
}

fn rewrite_uncached(e: &Expr, index: usize, r: &Rule, memo: &mut TrigMemo) -> (Expr, u8) {
    if e.is_atom() {
        return (e.clone(), 0);
    }
    let mut h = head_bit(e);
    let rebuilt = map_children(e, &mut |c| {
        let (nc, hc) = rewrite_cached(c, index, r, memo);
        h |= hc;
        nc
    });
    if h & r.heads == r.heads {
        if let Some(bind) = crate::pattern::pattern_match(&rebuilt, &r.pat)
            .expect("built-in trigonometric patterns have no host predicates")
        {
            let v = xreplace(&r.rep, &bind);
            let (new_size, old_size) = (memo.node_count(&v), memo.node_count(&rebuilt));
            if new_size < old_size || (new_size == old_size && v != rebuilt) {
                let v = share_terms(&v, &rebuilt);
                let hv = memo.heads(&v);
                return (v, hv);
            }
        }
    }
    (rebuilt, h)
}

/// The sum `v` with each term that also occurs in the sum `old` taken from
/// `old`. The rewritten sum gets the rest of its terms back through the
/// pattern bindings, which rebuild them: without this every later rule
/// would compare them in full against the memo keys of the original terms.
/// Only terms with the same presentation are exchanged, so the result is
/// the same expression.
fn share_terms(v: &Expr, old: &Expr) -> Expr {
    let (Kind::Add(c, terms), Kind::Add(_, old_terms)) = (v.kind(), old.kind()) else {
        return v.clone();
    };
    if !v.is_canonical() {
        return v.clone();
    }
    let originals: FxHashMap<&Expr, &Expr> = old_terms.iter().map(|(t, _)| (t, t)).collect();
    let mut changed = false;
    let shared: Vec<(Expr, Number)> = terms
        .iter()
        .map(|(t, k)| match originals.get(t) {
            Some(&o) if !o.ptr_eq(t) && o.same_presentation(t) => {
                changed = true;
                (o.clone(), k.clone())
            }
            _ => (t.clone(), k.clone()),
        })
        .collect();
    if !changed {
        return v.clone();
    }
    Expr::from_kind(Kind::Add(c.clone(), shared))
}

fn apply_rules(e: &Expr, memo: &mut TrigMemo) -> Expr {
    let reusable = match memo.rules.get_key_value(e) {
        Some((stored, result)) if stored.same_presentation(e) => return result.clone(),
        Some(_) => false,
        None => true,
    };
    let mut cur = e.clone();
    let mut cur_heads = memo.heads(&cur);
    for _ in 0..16 {
        let mut changed = false;
        for (index, r) in trig_rules().iter().enumerate() {
            if cur_heads & r.heads != r.heads {
                continue;
            }
            let (next, _) = rewrite_cached(&cur, index, r, memo);
            if next != cur {
                cur = next;
                cur_heads = memo.heads(&cur);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    if reusable {
        memo.rules.insert(e.clone(), cur.clone());
    }
    cur
}

/// sin(u)**2 -> 1 - cos(u)**2 (or the reverse) everywhere.
fn rewrite_squares(e: &Expr, to_cos: bool) -> Expr {
    rewrite_squares_cached(e, to_cos, &mut FxHashMap::default())
}

fn rewrite_squares_cached(e: &Expr, to_cos: bool, memo: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
    if let Some((_, result)) = memo.get(&e.dag_id()) { return result.clone(); }
    let result = rewrite_squares_uncached(e, to_cos, memo);
    memo.insert(e.dag_id(), (e.clone(), result.clone()));
    result
}

fn rewrite_squares_uncached(e: &Expr, to_cos: bool, memo: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
    let e2 = if e.is_atom() { e.clone() } else { map_children(e, &mut |c| rewrite_squares_cached(c, to_cos, memo)) };
    if let Kind::Pow(b, x) = e2.kind() {
        if let (Kind::Func(f, a), Some(k)) = (b.kind(), x.as_i64()) {
            if k >= 2 && k % 2 == 0 {
                let target = if to_cos { Func::Sin } else { Func::Cos };
                if *f == target {
                    let other = if to_cos { cos(&a[0]) } else { sin(&a[0]) };
                    return pow(&sub(&Expr::one(), &other.powi(2)), &Expr::int(k / 2));
                }
            }
        }
    }
    e2
}

pub(crate) fn has_trig(e: &Expr) -> bool {
    let mut pending = vec![e.clone()];
    let mut seen = FxHashMap::default();
    while let Some(node) = pending.pop() {
        if seen.insert(node.dag_id(), node.clone()).is_some() { continue; }
        if matches!(node.kind(), Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc | Func::Sinh | Func::Cosh | Func::Tanh, _)) {
            return true;
        }
        pending.extend(node.args());
    }
    false
}

/// Simplify trigonometric expressions using identities.
pub fn trigsimp(e: &Expr) -> Expr {
    trigsimp_cached(e, &mut TrigMemo::default())
}

/// `trigsimp` sharing its memo tables with other calls of one search.
pub(crate) fn trigsimp_cached(e: &Expr, memo: &mut TrigMemo) -> Expr {
    let reusable = match memo.trigsimp.get_key_value(e) {
        Some((stored, result)) if stored.same_presentation(e) => return result.clone(),
        Some(_) => false,
        None => true,
    };
    let result = trigsimp_uncached(e, memo);
    if reusable {
        memo.trigsimp.insert(e.clone(), result.clone());
    }
    result
}

fn trigsimp_uncached(e: &Expr, memo: &mut TrigMemo) -> Expr {
    if !has_trig(e) {
        return e.clone();
    }
    // simplify inside first
    let e1 = if e.is_atom() { e.clone() } else { map_children(e, &mut |c| trigsimp_cached(c, memo)) };
    let mut best = apply_rules(&e1, memo);
    let mut cands = vec![];
    let ex = expand(&e1);
    if ex != e1 {
        cands.push(apply_rules(&ex, memo));
    }
    for to_cos in [true, false] {
        let r = expand(&rewrite_squares(&e1, to_cos));
        cands.push(apply_rules(&r, memo));
    }
    for c in cands {
        if node_count(&c) < node_count(&best) {
            best = c;
        }
    }
    best
}

// ----------------------------------------------------------------------
// powers and logarithms
// ----------------------------------------------------------------------

/// Combine powers with equal bases (and exponentials), SymPy's `powsimp`.
pub fn powsimp(e: &Expr) -> Expr {
    powsimp_deep(e, true)
}

pub fn powsimp_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| powsimp_deep(c, deep)) };
    let Kind::Mul(c, facs) = e2.kind() else { return e2 };
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let mut exps: Vec<Expr> = Vec::new();
    for (b, x) in facs {
        if let Kind::Func(Func::Exp, a) = b.kind() {
            exps.push(mul2(&a[0], x));
            continue;
        }
        if b.is_constant(Constant::E) {
            exps.push(x.clone());
            continue;
        }
        match groups.iter_mut().find(|(g, _)| g == b) {
            Some((_, v)) => v.push(x.clone()),
            None => groups.push((b.clone(), vec![x.clone()])),
        }
    }
    let mut out = vec![Expr::num(c.clone())];
    for (b, xs) in groups {
        out.push(Expr::from_kind(Kind::Pow(b.clone(), add(&xs))).normalize_pow());
    }
    if !exps.is_empty() {
        out.push(exp(&add(&exps)));
    }
    mul(&out)
}

impl Expr {
    /// Re-canonicalise a raw `Pow` node.
    fn normalize_pow(self) -> Expr {
        match self.kind() {
            Kind::Pow(b, x) => pow(b, x),
            _ => self,
        }
    }
}

/// `log(x*y) -> log(x) + log(y)`, `log(x**n) -> n*log(x)` for positive
/// arguments (or everything with `force`).
pub fn expand_log(e: &Expr, force: bool) -> Expr {
    expand_log_deep(e, force, true)
}

pub fn expand_log_deep(e: &Expr, force: bool, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_log_deep(c, force, deep)) };
    let Kind::Func(Func::Log, a) = e2.kind() else { return e2 };
    if a.len() != 1 {
        return e2;
    }
    expand_log_arg(&a[0], force).unwrap_or(e2)
}

/// SymPy's `log._eval_expand_log` (without the `factor` hint); `None` keeps
/// the logarithm unchanged.
fn expand_log_arg(arg: &Expr, force: bool) -> Option<Expr> {
    // the logarithm of a factor, expanded again when it stays a logarithm
    let log_of = |x: &Expr| -> Expr {
        let a = log(x);
        match a.kind() {
            Kind::Func(Func::Log, b) if b.len() == 1 => expand_log_arg(&b[0], force).unwrap_or(a),
            // log(1/4) evaluates to -log(4), whose argument is expanded too
            Kind::Mul(..) => map_children(&a, &mut |c| match c.kind() {
                Kind::Func(Func::Log, b) if b.len() == 1 => expand_log_arg(&b[0], force).unwrap_or_else(|| c.clone()),
                _ => c.clone(),
            }),
            _ => a,
        }
    };
    let polar = |x: &Expr| crate::special::polar::is_polar(x);
    match arg.kind() {
        Kind::Num(n) if n.is_integer() => {
            // remove perfect powers
            let (p, _) = n.numer_denom()?;
            let (base, k) = crate::ntheory::factor::perfect_power(&p, None, true)?;
            Some(mul2(&Expr::int(k as i64), &log(&Expr::num(Number::from_bigint(base)))))
        }
        Kind::Num(n) if n.is_rational() => {
            let (p, q) = n.numer_denom()?;
            // SymPy's expand repeats the log hint to a fixed point.
            Some(sub(&log_of(&Expr::num(Number::from_bigint(p))), &log_of(&Expr::num(Number::from_bigint(q)))))
        }
        Kind::Mul(..) => {
            let mut parts = Vec::new();
            let mut nonpos = Vec::new();
            for x in arg.args() {
                if force || is_positive(&x) == Some(true) || polar(&x) {
                    parts.push(log_of(&x));
                } else if crate::assumptions::is_negative(&x) == Some(true) {
                    parts.push(log(&neg(&x)));
                    nonpos.push(Expr::minus_one());
                } else {
                    nonpos.push(x);
                }
            }
            parts.push(log(&mul(&nonpos)));
            Some(add(&parts))
        }
        Kind::Pow(b, x) => expand_log_power(b, x, force, &log_of),
        Kind::Func(Func::Exp, x) if x.len() == 1 => expand_log_power(&Expr::constant(Constant::E), &x[0], force, &log_of),
        _ => None,
    }
}

fn expand_log_power(b: &Expr, x: &Expr, force: bool, log_of: &dyn Fn(&Expr) -> Expr) -> Option<Expr> {
    let real = crate::assumptions::fact(x, "extended_real") == Some(true);
    let near = || {
        is_positive(&add2(x, &Expr::one())) == Some(true)
            && crate::assumptions::fact(&sub(x, &Expr::one()), "nonpositive") == Some(true)
    };
    if force || (real && (is_positive(b) == Some(true) || near())) || crate::special::polar::is_polar(b) {
        return Some(mul2(&crate::special::polar::unpolarify(x, false), &log_of(b)));
    }
    None
}

/// SymPy's `logcombine`: `a*log(x) + b*log(y) -> log(x**a * y**b)` for
/// positive arguments and real coefficients (or everything with `force`).
/// Ported from sympy/simplify/simplify.py (see LICENSE-SYMPY).
pub fn logcombine(e: &Expr, force: bool) -> Expr {
    crate::simptools::bottom_up(e, &mut |rv| logcombine_node(rv, force))
}

fn single_log(e: &Expr) -> Option<&Expr> {
    match e.kind() {
        Kind::Func(Func::Log, a) if a.len() == 1 => Some(&a[0]),
        _ => None,
    }
}

fn unevaluated_log(a: Expr) -> Expr {
    Expr::from_unevaluated_kind(Kind::Func(Func::Log, vec![a]))
}

fn logcombine_node(rv: &Expr, force: bool) -> Expr {
    use crate::assumptions::fact;
    use crate::simptools::{add_args, args, ordered};
    if !(rv.is_add() || rv.is_mul()) {
        return rv.clone();
    }
    let gooda = |a: &Expr| {
        !a.is_minus_one() && {
            let real = fact(a, "extended_real");
            real == Some(true) || (force && real != Some(false))
        }
    };
    let goodlog = |l: &Expr| {
        let a = single_log(l).unwrap();
        is_positive(a) == Some(true) || (force && fact(a, "nonpositive") != Some(false))
    };
    // log1: insertion-ordered map from a tuple of "other" factors to
    // (coefficients, log) pairs
    let mut log1: Vec<(Vec<Expr>, Vec<(Vec<Expr>, Expr)>)> = Vec::new();
    fn push(log1: &mut Vec<(Vec<Expr>, Vec<(Vec<Expr>, Expr)>)>, key: Vec<Expr>, item: (Vec<Expr>, Expr)) {
        match log1.iter_mut().find(|(k, _)| *k == key) {
            Some((_, v)) => v.push(item),
            None => log1.push((key, vec![item])),
        }
    }
    let mut other: Vec<Expr> = Vec::new();
    let mut logs: Vec<(Vec<Expr>, Vec<Expr>, Vec<Expr>)> = Vec::new();
    for a in add_args(rv) {
        if single_log(&a).is_some() && goodlog(&a) {
            push(&mut log1, vec![], (vec![], a));
        } else if !a.is_mul() {
            other.push(a);
        } else {
            let (mut ot, mut co, mut lo) = (Vec::new(), Vec::new(), Vec::new());
            for ai in args(&a) {
                if matches!(ai.kind(), Kind::Num(n) if n.is_rational() && n.is_negative()) {
                    ot.push(Expr::minus_one());
                    co.push(neg(&ai));
                } else if single_log(&ai).is_some() && goodlog(&ai) {
                    lo.push(ai);
                } else if gooda(&ai) {
                    co.push(ai);
                } else {
                    ot.push(ai);
                }
            }
            if lo.len() > 1 {
                logs.push((ot, co, lo));
            } else if let Some(l) = lo.pop() {
                push(&mut log1, ot, (co, l));
            } else {
                other.push(a);
            }
        }
    }
    // a single log in other joins the good logs
    if other.len() == 1 && single_log(&other[0]).is_some() {
        let l = other.pop().unwrap();
        push(&mut log1, vec![], (vec![], l));
    }
    // SymPy compares the (coefficients, log) tuple with [], which never
    // holds: only an absence of logs leaves the expression alone.
    if logs.is_empty() && log1.is_empty() {
        return rv.clone();
    }
    for (o, e, l) in logs {
        let mut l = ordered(l);
        let first = l.remove(0);
        let mut combined = log(&pow(single_log(&first).unwrap(), &mul(&e)));
        for li in l {
            combined = log(&pow(single_log(&li).unwrap(), &combined));
        }
        let c = mul(&o);
        if single_log(&combined).is_some() {
            push(&mut log1, vec![c], (vec![], combined));
        } else {
            other.push(mul2(&c, &combined));
        }
    }
    // logs with the same coefficient multiply
    let mut merged: Vec<(Expr, Expr)> = Vec::new();
    for (k, items) in log1 {
        let product = mul(&items.iter().map(|(c, l)| pow(single_log(l).unwrap(), &mul(c))).collect::<Vec<_>>());
        let value = unevaluated_log(logcombine(&product, force));
        let key = mul(&k);
        match merged.iter_mut().find(|(k2, _)| *k2 == key) {
            Some(slot) => slot.1 = value,
            None => merged.push((key, value)),
        }
    }
    // logs with opposite coefficients divide
    let keys = ordered(merged.iter().map(|(k, _)| k.clone()).collect());
    let take = |merged: &mut Vec<(Expr, Expr)>, k: &Expr| -> Option<Expr> {
        let i = merged.iter().position(|(k2, _)| k2 == k)?;
        Some(merged.remove(i).1)
    };
    for k in keys {
        if !merged.iter().any(|(k2, _)| *k2 == k) {
            continue;
        }
        let nk = neg(&k);
        if merged.iter().any(|(k2, _)| *k2 == nk) {
            let (mut num, mut den) = (k.clone(), nk);
            if crate::simptools::methods::count_ops(&num) > crate::simptools::methods::count_ops(&den) {
                std::mem::swap(&mut num, &mut den);
            }
            let ln = take(&mut merged, &num).unwrap();
            let ld = take(&mut merged, &den).unwrap();
            let ratio = div(single_log(&ln).unwrap(), single_log(&ld).unwrap());
            other.push(mul2(&num, &unevaluated_log(ratio)));
        } else {
            let l = take(&mut merged, &k).unwrap();
            other.push(mul2(&k, &l));
        }
    }
    add(&other)
}

/// `(x*y)**a -> x**a * y**a` (for positive factors, or with `force`).
pub fn expand_power_base(e: &Expr, force: bool) -> Expr {
    expand_power_base_deep(e, force, true)
}

pub fn expand_power_base_deep(e: &Expr, force: bool, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_power_base_deep(c, force, deep)) };
    if let Kind::Pow(b, x) = e2.kind() {
        if let Kind::Mul(..) = b.kind() {
            let fs = b.args();
            if force || fs.iter().all(|f| is_positive(f) == Some(true)) {
                return mul(&fs.iter().map(|f| pow(f, x)).collect::<Vec<_>>());
            }
        }
    }
    e2
}

/// `x**(a+b) -> x**a * x**b` when the base is known to be non-zero.
pub fn expand_power_exp(e: &Expr) -> Expr {
    expand_power_exp_deep(e, true)
}

pub fn expand_power_exp_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else if let Kind::Add(c, terms) = e.kind() {
        // Scaling an expanded child with mul2 would merge its factors again.
        let transformed: Vec<_> = terms.iter().map(|(t,k)| (expand_power_exp_deep(t, true), k.clone())).collect();
        if transformed.iter().zip(terms).all(|((a,_),(b,_))| a == b) { e.clone() }
        else { Expr::from_unevaluated_kind(Kind::Add(c.clone(), transformed)) }
    } else { map_children(e, &mut |c| expand_power_exp_deep(c, deep)) };
    match e2.kind() {
        Kind::Pow(b, x) if x.is_add() && crate::assumptions::is_zero(b) == Some(false) => {
            expanded_product(&x.args().iter().map(|t| pow(b, t)).collect::<Vec<_>>())
        }
        // ExpBase._eval_expand_power_exp splits only a commutative exponent.
        Kind::Func(Func::Exp, a) if a[0].is_add() && crate::decompose::commutes(&a[0]) => expanded_product(&a[0].args().iter().map(exp).collect::<Vec<_>>()),
        Kind::Mul(..) if deep => {
            let args = e2.args();
            let changed: Vec<_> = args.iter().map(|f| expand_power_exp_deep(f, false)).collect();
            if args.iter().zip(&changed).any(|(a,b)| a != b) { expanded_product(&changed) } else { e2 }
        }
        _ => e2,
    }
}

/// The ordinary product constructor merges equal bases, undoing this
/// expansion. Preserve its explicit factors just as for evaluate=False.
pub(crate) fn expanded_product(args: &[Expr]) -> Expr {
    let mut coefficient = Number::ONE;
    let mut factors = Vec::new();
    for arg in args {
        match arg.kind() {
            Kind::Num(n) => coefficient = coefficient.mul(n),
            Kind::Mul(c, fs) => { coefficient = coefficient.mul(c); factors.extend(fs.iter().cloned()); },
            Kind::Pow(base, exponent) => factors.push((base.clone(), exponent.clone())),
            _ => factors.push((arg.clone(), Expr::one())),
        }
    }
    if factors.is_empty() { return Expr::num(coefficient); }
    Expr::from_unevaluated_kind(Kind::Mul(coefficient, factors))
}

/// Distribute products over sums (but do not expand powers).
pub fn expand_mul(e: &Expr) -> Expr {
    expand_mul_deep(e, true)
}

pub fn expand_mul_deep(e: &Expr, deep: bool) -> Expr {
    let e2 = if !deep || e.is_atom() { e.clone() } else { map_children(e, &mut |c| expand_mul_deep(c, deep)) };
    if let Kind::Mul(..) = e2.kind() {
        let fs = e2.args();
        if fs.iter().any(|f| f.is_add()) {
            let mut acc = vec![Expr::one()];
            for f in fs {
                let terms = if f.is_add() { f.args() } else { vec![f] };
                let mut next = Vec::with_capacity(acc.len() * terms.len());
                for a in &acc {
                    for t in &terms {
                        next.push(mul2(a, t));
                    }
                }
                acc = next;
            }
            return add(&acc);
        }
    }
    e2
}

// ----------------------------------------------------------------------
// nsimplify
// ----------------------------------------------------------------------

/// Best rational approximation of `x` with denominator <= `max_den`.
pub(crate) fn rational_approx(x: f64, max_den: i64, tol: f64) -> Option<(i64, i64)> {
    let (mut h0, mut h1, mut k0, mut k1) = (0i64, 1i64, 1i64, 0i64);
    let mut v = x;
    for _ in 0..64 {
        let a = v.floor();
        if a.abs() > 1e15 {
            return None;
        }
        let a = a as i64;
        let h2 = a.checked_mul(h1)?.checked_add(h0)?;
        let k2 = a.checked_mul(k1)?.checked_add(k0)?;
        if k2 > max_den {
            return None;
        }
        if ((h2 as f64) / (k2 as f64) - x).abs() <= tol * x.abs().max(1.0) {
            return Some((h2, k2));
        }
        (h0, h1, k0, k1) = (h1, h2, k1, k2);
        let frac = v - a as f64;
        if frac.abs() < 1e-300 {
            return None;
        }
        v = 1.0 / frac;
    }
    None
}

/// Find a simple exact expression close to a float (rationals, square
/// roots of rationals, rational multiples of the given constants).
pub fn nsimplify_float(x: f64, constants: &[Expr]) -> Option<Expr> {
    let tol = 1e-12;
    if x == 0.0 {
        return Some(Expr::zero());
    }
    // simplest candidates first: small rationals, multiples of the given
    // constants, square roots, then rationals with larger denominators
    if let Some((p, q)) = rational_approx(x, 10_000, tol) {
        return Some(Expr::rational(p, q));
    }
    for c in constants {
        let cv = crate::evalf::eval_f64(c)?;
        if let Some((p, q)) = rational_approx(x / cv, 10_000, tol) {
            return Some(mul2(&Expr::rational(p, q), c));
        }
    }
    let sq = x * x;
    if let Some((p, q)) = rational_approx(sq, 10_000, 1e-11) {
        let r = sqrt(&Expr::rational(p, q));
        return Some(if x < 0.0 { neg(&r) } else { r });
    }
    if let Some((p, q)) = rational_approx(x, 1_000_000_000, tol) {
        return Some(Expr::rational(p, q));
    }
    None
}

/// Replace floats in `e` by simple exact numbers.
pub fn nsimplify(e: &Expr, constants: &[Expr]) -> Expr {
    // An exact algebraic expression must not be replaced by a nearby
    // rational obtained by rounding it to a machine float first.
    if e.float_info().max_precision.is_none() { return e.clone(); }
    match e.kind() {
        Kind::Num(Number::Float(f)) => nsimplify_float(*f, constants).unwrap_or_else(|| e.clone()),
        _ if e.is_atom() => e.clone(),
        _ => {
            if crate::functions::is_numeric_value(e) {
                if let Some(v) = crate::evalf::eval_f64(e) {
                    if let Some(r) = nsimplify_float(v, constants) {
                        return r;
                    }
                }
            }
            match e.kind() {
                // numeric coefficients are not children: go through args()
                Kind::Add(..) => add(&e.args().iter().map(|c| nsimplify(c, constants)).collect::<Vec<_>>()),
                Kind::Mul(..) => mul(&e.args().iter().map(|c| nsimplify(c, constants)).collect::<Vec<_>>()),
                _ => map_children(e, &mut |c| nsimplify(c, constants)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(f: fn(&Expr) -> Expr, s: &str) -> String {
        f(&parse_str(s).unwrap()).to_string()
    }

    #[test]
    fn trig() {
        assert_eq!(t(expand_trig, "sin(2*x)"), "2*sin(x)*cos(x)");
        assert_eq!(t(expand_trig, "cos(2*x)"), "2*cos(x)**2 - 1");
        assert_eq!(t(expand_trig, "sin(x+y)"), "sin(x)*cos(y) + sin(y)*cos(x)");
        assert_eq!(t(expand_trig, "cos(x+y)"), "-sin(x)*sin(y) + cos(x)*cos(y)");
        assert_eq!(t(expand_trig, "sin(3*x)"), "-4*sin(x)**3 + 3*sin(x)");
        assert_eq!(t(expand_trig, "cos(3*x)"), "4*cos(x)**3 - 3*cos(x)");
        assert_eq!(t(trigsimp, "sin(x)**2+cos(x)**2"), "1");
        assert_eq!(t(trigsimp, "2*sin(x)**2+2*cos(x)**2"), "2");
        assert_eq!(t(trigsimp, "sin(x)/cos(x)"), "tan(x)");
        assert_eq!(t(trigsimp, "1-sin(x)**2"), "cos(x)**2");
        assert_eq!(t(trigsimp, "2*sin(x)*cos(x)"), "sin(2*x)");
        assert_eq!(t(trigsimp, "cos(x)**2-sin(x)**2"), "cos(2*x)");
        assert_eq!(t(trigsimp, "sin(x)**4 - cos(x)**4"), "-cos(2*x)");
        assert_eq!(t(trigsimp, "1+tan(x)**2"), "cos(x)**(-2)");
        assert_eq!(t(trigsimp, "sin(x)*cos(y)+cos(x)*sin(y)"), "sin(x + y)");
        assert_eq!(t(trigsimp, "cosh(x)**2-sinh(x)**2"), "1");
        assert_eq!(t(trigsimp, "x*sin(y)**2 + x*cos(y)**2"), "x");
    }

    #[test]
    fn trig_passes_visit_shared_subexpressions_once_per_pass() {
        let x = Expr::symbol("x");
        let mut expression = sin(&x);
        for _ in 0..64 {
            expression = undefined("shared", vec![expression.clone(), expression]);
        }
        assert_eq!(heads(&expression), 1);
        assert!(has_trig(&expression));
        assert_eq!(trigsimp(&expression), expression);
    }

    #[test]
    fn shared_rewrites_keep_every_occurrence_and_exact_identity() {
        let x = Expr::symbol("x");
        let mut expression = add2(&sin(&x).powi(2), &cos(&x).powi(2));
        let mut expected = Expr::one();
        for _ in 0..64 {
            expression = undefined("shared", vec![expression.clone(), expression]);
            expected = undefined("shared", vec![expected.clone(), expected]);
        }
        assert_eq!(trigsimp(&expression), expected);
        assert_eq!(heads(&expected), 0);
        assert!(!has_trig(&expected));
    }

    #[test]
    fn powers_logs() {
        assert_eq!(t(powsimp, "x**a*x**b"), "x**(a + b)");
        assert_eq!(t(powsimp, "exp(x)*exp(y)"), "exp(x + y)");
        assert_eq!(t(powsimp, "2**x*2**y"), "2**(x + y)");
        assert_eq!(t(powsimp, "x**2*x**a"), "x**(a + 2)");
        assert_eq!(t(expand_mul, "x*(y+z)"), "x*y + x*z");
        assert_eq!(nsimplify(&Expr::float(0.5), &[]).to_string(), "1/2");
        assert_eq!(nsimplify(&Expr::float(1.0 / 3.0), &[]).to_string(), "1/3");
        assert_eq!(nsimplify(&Expr::float(2f64.sqrt()), &[]).to_string(), "sqrt(2)");
        assert_eq!(nsimplify(&Expr::float(std::f64::consts::PI), &[Expr::pi()]).to_string(), "pi");
    }
}

