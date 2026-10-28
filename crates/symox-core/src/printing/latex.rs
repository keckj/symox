//! LaTeX printer reproducing SymPy's `latex()` output.

use rustc_hash::FxHashMap;

use super::{
    float_str, is_known_negative, minus_sign, ordered_factors, precedence_traditional, PREC_ATOM,
    PREC_MUL, PREC_POW,
};
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::Number;
use crate::ops::{mul_factor, neg, pow};

/// Output environment, like SymPy's `mode` setting.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum LatexMode {
    #[default]
    Plain,
    Inline,
    Equation,
    EquationStar,
}

/// How inverse trigonometric functions are written.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum InvTrigStyle {
    /// `\operatorname{asin}`
    #[default]
    Abbreviated,
    /// `\arcsin`
    Full,
    /// `\sin^{-1}`
    Power,
}

/// Settings of the LaTeX printer (a subset of SymPy's).
#[derive(Clone, Debug)]
pub struct LatexSettings {
    pub mode: LatexMode,
    /// `None` (juxtaposition), `"dot"`, `"times"`, `"ldot"` or a literal.
    pub mul_symbol: Option<String>,
    /// Print `a / b` instead of `\frac{a}{b}` for short fractions
    /// (default: only in inline mode).
    pub fold_short_frac: Option<bool>,
    /// SymPy's `long_frac_ratio`: the allowed ratio of the numerator width
    /// to the denominator width before a long fraction is broken off
    /// (`None`, the default, never breaks).
    pub long_frac_ratio: Option<f64>,
    pub fold_frac_powers: bool,
    pub fold_func_brackets: bool,
    pub inv_trig_style: InvTrigStyle,
    pub ln_notation: bool,
    pub root_notation: bool,
    pub itex: bool,
    /// `"i"`, `"j"`, `"ri"`, `"rj"`, `"ti"`, `"tj"` or a literal.
    pub imaginary_unit: String,
    /// `"d"`, `"rd"`, `"td"` or a literal.
    pub diff_operator: String,
    /// Custom LaTeX for symbols, keyed by symbol name.
    pub symbol_names: FxHashMap<String, String>,
    pub parenthesize_super: bool,
    /// SymPy's `mat_symbol_style`: `"plain"` or `"bold"`.
    pub mat_symbol_bold: bool,
    /// SymPy's `adjoint_style`: `"dagger"`, `"star"` or `"hermitian"`.
    pub adjoint_style: String,
    /// SymPy's `mat_delim`: `"["`, `"("` or empty (no delimiters).
    pub mat_delim: String,
    /// SymPy's `mat_str` (`None` chooses `matrix`, `smallmatrix` or `array`).
    pub mat_str: Option<String>,
}

impl Default for LatexSettings {
    fn default() -> Self {
        LatexSettings {
            mode: LatexMode::Plain,
            mul_symbol: None,
            fold_short_frac: None,
            long_frac_ratio: None,
            fold_frac_powers: false,
            fold_func_brackets: false,
            inv_trig_style: InvTrigStyle::Abbreviated,
            ln_notation: false,
            root_notation: true,
            itex: false,
            imaginary_unit: "i".into(),
            diff_operator: "d".into(),
            symbol_names: FxHashMap::default(),
            parenthesize_super: true,
            mat_symbol_bold: false,
            adjoint_style: "dagger".into(),
            mat_delim: "[".into(),
            mat_str: None,
        }
    }
}

/// SymPy's `latex(expr)` with default settings.
pub fn latex(e: &Expr) -> String {
    latex_with(e, &LatexSettings::default())
}

/// `latex(expr, **settings)`.
pub fn latex_with(e: &Expr, settings: &LatexSettings) -> String {
    let _session = super::order::session();
    let p = LatexPrinter::new(settings);
    let tex = p.print(e);
    match settings.mode {
        LatexMode::Plain => tex,
        LatexMode::Inline => format!("${tex}$"),
        _ if settings.itex => format!("$${tex}$$"),
        LatexMode::Equation => format!("\\begin{{equation}}{tex}\\end{{equation}}"),
        LatexMode::EquationStar => format!("\\begin{{equation*}}{tex}\\end{{equation*}}"),
    }
}

// ----------------------------------------------------------------------
// symbol names
// ----------------------------------------------------------------------

const GREEKS: [&str; 24] = [
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
    "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi",
    "psi", "omega",
];

fn tex_greek(s: &str) -> Option<&'static str> {
    Some(match s {
        "Alpha" => r"\mathrm{A}",
        "Beta" => r"\mathrm{B}",
        "Gamma" => r"\Gamma",
        "Delta" => r"\Delta",
        "Epsilon" => r"\mathrm{E}",
        "Zeta" => r"\mathrm{Z}",
        "Eta" => r"\mathrm{H}",
        "Theta" => r"\Theta",
        "Iota" => r"\mathrm{I}",
        "Kappa" => r"\mathrm{K}",
        "Lambda" => r"\Lambda",
        "Mu" => r"\mathrm{M}",
        "Nu" => r"\mathrm{N}",
        "Xi" => r"\Xi",
        "omicron" => "o",
        "Omicron" => r"\mathrm{O}",
        "Pi" => r"\Pi",
        "Rho" => r"\mathrm{P}",
        "Sigma" => r"\Sigma",
        "Tau" => r"\mathrm{T}",
        "Upsilon" => r"\Upsilon",
        "Phi" => r"\Phi",
        "Chi" => r"\mathrm{X}",
        "Psi" => r"\Psi",
        "Omega" => r"\Omega",
        "lamda" => r"\lambda",
        "Lamda" => r"\Lambda",
        "khi" => r"\chi",
        "Khi" => r"\mathrm{X}",
        "varepsilon" => r"\varepsilon",
        "varkappa" => r"\varkappa",
        "varphi" => r"\varphi",
        "varpi" => r"\varpi",
        "varrho" => r"\varrho",
        "varsigma" => r"\varsigma",
        "vartheta" => r"\vartheta",
        _ => return None,
    })
}

const OTHER_SYMBOLS: [&str; 10] = [
    "aleph", "beth", "daleth", "gimel", "ell", "eth", "hbar", "hslash", "mho", "wp",
];

/// Modifiers sorted by decreasing length.
const MODIFIERS: &[(&str, &str, &str)] = &[
    ("mathring", r"\mathring{", "}"),
    ("ddddot", r"\ddddot{", "}"),
    ("dddot", r"\dddot{", "}"),
    ("check", r"\check{", "}"),
    ("breve", r"\breve{", "}"),
    ("acute", r"\acute{", "}"),
    ("grave", r"\grave{", "}"),
    ("tilde", r"\tilde{", "}"),
    ("prime", "{", "}'"),
    ("ddot", r"\ddot{", "}"),
    ("bold", r"\boldsymbol{", "}"),
    ("frak", r"\mathfrak{", "}"),
    ("norm", r"\left\|{", r"}\right\|"),
    ("dot", r"\dot{", "}"),
    ("hat", r"\hat{", "}"),
    ("bar", r"\bar{", "}"),
    ("vec", r"\vec{", "}"),
    ("prm", "{", "}'"),
    ("cal", r"\mathcal{", "}"),
    ("scr", r"\mathscr{", "}"),
    ("avg", r"\left\langle{", r"}\right\rangle"),
    ("abs", r"\left|{", r"}\right|"),
    ("mag", r"\left|{", r"}\right|"),
    ("bm", r"\boldsymbol{", "}"),
];

/// `sympy.printing.latex.translate`.
pub fn translate(s: &str) -> String {
    if let Some(t) = tex_greek(s) {
        return t.to_string();
    }
    let lower = s.to_lowercase();
    if GREEKS.contains(&lower.as_str()) {
        return format!("\\{lower}");
    }
    if OTHER_SYMBOLS.contains(&s) {
        return format!("\\{s}");
    }
    for (key, pre, post) in MODIFIERS.iter() {
        if lower.ends_with(key) && s.chars().count() > key.len() {
            // keys are ASCII, so the suffix has the same byte length
            let head = &s[..s.len() - key.len()];
            return format!("{pre}{}{post}", translate(head));
        }
    }
    s.to_string()
}

/// `sympy.printing.conventions.split_super_sub`.
pub(crate) fn split_super_sub(text: &str) -> (String, Vec<String>, Vec<String>) {
    if text.is_empty() {
        return (String::new(), vec![], vec![]);
    }
    let b = text.as_bytes();
    let mut pos = 0usize;
    let mut name: Option<String> = None;
    let mut supers = Vec::new();
    let mut subs = Vec::new();
    let find = |c: u8, start: usize| -> usize {
        if start >= b.len() {
            return b.len();
        }
        b[start..].iter().position(|&x| x == c).map(|p| p + start).unwrap_or(b.len())
    };
    while pos < b.len() {
        let mut start = pos + 1;
        if text[pos..].starts_with("__") {
            start += 1;
        }
        let pos_next = find(b'^', start).min(find(b'_', start));
        let part = &text[pos..pos_next];
        pos = pos_next;
        if name.is_none() {
            name = Some(part.to_string());
        } else if let Some(p) = part.strip_prefix('^') {
            supers.push(p.to_string());
        } else if let Some(p) = part.strip_prefix("__") {
            supers.push(p.to_string());
        } else if let Some(p) = part.strip_prefix('_') {
            subs.push(p.to_string());
        }
    }
    let mut name = name.unwrap_or_default();
    // a name ending with digits: treat the digits as a subscript
    let digits_start = name
        .char_indices()
        .rev()
        .take_while(|(_, c)| c.is_ascii_digit())
        .last()
        .map(|(i, _)| i);
    if let Some(i) = digits_start {
        let head = &name[..i];
        if !head.is_empty() && head.chars().all(|c| c.is_alphabetic()) {
            subs.insert(0, name[i..].to_string());
            name = head.to_string();
        }
    }
    (name, supers, subs)
}

/// `LatexPrinter._deal_with_super_sub`.
pub(crate) fn deal_with_super_sub(s: &str) -> String {
    if s.contains('{') {
        return s.to_string();
    }
    let (name, supers, subs) = split_super_sub(s);
    let mut out = translate(&name);
    if !supers.is_empty() {
        let v: Vec<String> = supers.iter().map(|x| translate(x)).collect();
        out.push_str(&format!("^{{{}}}", v.join(" ")));
    }
    if !subs.is_empty() {
        let v: Vec<String> = subs.iter().map(|x| translate(x)).collect();
        out.push_str(&format!("_{{{}}}", v.join(" ")));
    }
    out
}

const ACCEPTED_LATEX_FUNCTIONS: [&str; 19] = [
    "arcsin", "arccos", "arctan", "sin", "cos", "tan", "sinh", "cosh", "tanh", "sqrt", "ln", "log",
    "sec", "csc", "cot", "coth", "re", "im", "frac",
];

fn is_accepted(f: &str) -> bool {
    ACCEPTED_LATEX_FUNCTIONS.contains(&f) || f == "root" || f == "arg"
}

const INV_TRIG: [&str; 12] = [
    "asin", "acos", "atan", "acsc", "asec", "acot", "asinh", "acosh", "atanh", "acsch", "asech",
    "acoth",
];

// ----------------------------------------------------------------------
// printer
// ----------------------------------------------------------------------

struct LatexPrinter<'a> {
    s: &'a LatexSettings,
    sep: String,
    numsep: String,
    iu: String,
    diff_op: String,
    fold_short_frac: bool,
}

fn mul_symbol_table(k: &str) -> Option<&'static str> {
    Some(match k {
        "ldot" => r" \,.\, ",
        "dot" => r" \cdot ",
        "times" => r" \times ",
        _ => return None,
    })
}

impl<'a> LatexPrinter<'a> {
    fn new(s: &'a LatexSettings) -> Self {
        let sep = match &s.mul_symbol {
            None => " ".to_string(),
            Some(m) => mul_symbol_table(m).map(|x| x.to_string()).unwrap_or_else(|| m.clone()),
        };
        let numsep = match &s.mul_symbol {
            None => r" \cdot ".to_string(),
            Some(m) => match mul_symbol_table(m) {
                Some(x) => x.to_string(),
                None => {
                    if ["", " ", "\\", "\\,", "\\:", "\\;", "\\quad"].contains(&m.trim()) {
                        r" \cdot ".to_string()
                    } else {
                        m.clone()
                    }
                }
            },
        };
        let iu = match s.imaginary_unit.as_str() {
            "i" => "i".to_string(),
            "ri" => r"\mathrm{i}".to_string(),
            "ti" => r"\text{i}".to_string(),
            "j" => "j".to_string(),
            "rj" => r"\mathrm{j}".to_string(),
            "tj" => r"\text{j}".to_string(),
            o => o.to_string(),
        };
        let diff_op = match s.diff_operator.as_str() {
            "d" => "d".to_string(),
            "rd" => r"\mathrm{d}".to_string(),
            "td" => r"\text{d}".to_string(),
            o => o.to_string(),
        };
        let fold_short_frac = s.fold_short_frac.unwrap_or(s.mode == LatexMode::Inline);
        LatexPrinter { s, sep, numsep, iu, diff_op, fold_short_frac }
    }

    fn parens(s: &str) -> String {
        format!("\\left({s}\\right)")
    }

    fn parenthesize(&self, e: &Expr, level: i32, strict: bool) -> String {
        let p = precedence_traditional(e);
        if p < level || (!strict && p <= level) {
            Self::parens(&self.print(e))
        } else {
            self.print(e)
        }
    }

    fn parenthesize_super(&self, s: String) -> String {
        if s.contains('^') {
            if self.s.parenthesize_super {
                Self::parens(&s)
            } else {
                format!("{{{s}}}")
            }
        } else {
            s
        }
    }

    fn print(&self, e: &Expr) -> String {
        if e.has_host() {
            if let Some(text) = crate::host::print(e, crate::host::HostFormat::Latex) { return text; }
        }
        match e.kind() {
            Kind::Num(n) => self.print_number(n),
            Kind::Sym(symbol) => {
                let name = crate::host::atom_name(e, symbol);
                match self.s.symbol_names.get(&*name) {
                    Some(t) => t.clone(),
                    None => deal_with_super_sub(&name),
                }
            },
            Kind::Const(c) => match c {
                Constant::Pi => r"\pi".into(),
                Constant::E => "e".into(),
                Constant::I => self.iu.clone(),
                Constant::EulerGamma => r"\gamma".into(),
                Constant::Catalan => "G".into(),
                Constant::GoldenRatio => r"\phi".into(),
                Constant::TribonacciConstant => r"\text{TribonacciConstant}".into(),
                Constant::Infinity => r"\infty".into(),
                Constant::NegInfinity => r"-\infty".into(),
                Constant::ComplexInfinity => r"\tilde{\infty}".into(),
                Constant::NaN => r"\text{NaN}".into(),
            },
            Kind::Bool(b) => format!("\\text{{{}}}", if *b { "True" } else { "False" }),
            Kind::Add(..) => self.print_add(e),
            Kind::Mul(..) => self.print_mul(e),
            Kind::Pow(b, x) => self.print_pow(e, b, x),
            Kind::Func(Func::SymbolicDerivative,_) => {let(body,vars)=crate::diff::symbolic::parts(e).unwrap();self.print_symbolic_derivative(body,&vars)},
            Kind::Func(..) => self.print_func(e, None),
            Kind::Derivative(f, vars) => self.print_derivative(f, vars),
            Kind::Rel(op, a, b) => {
                let o = match op {
                    RelOp::Eq => "=",
                    RelOp::Ne => r"\neq",
                    RelOp::Lt => {
                        if self.s.itex {
                            r"\lt"
                        } else {
                            "<"
                        }
                    }
                    RelOp::Gt => {
                        if self.s.itex {
                            r"\gt"
                        } else {
                            ">"
                        }
                    }
                    RelOp::Le => r"\leq",
                    RelOp::Ge => r"\geq",
                };
                format!("{} {} {}", self.print(a), o, self.print(b))
            }
            Kind::Piecewise(p) => {
                let mut pairs: Vec<String> = Vec::with_capacity(p.len());
                for (i, (x, c)) in p.iter().enumerate() {
                    if i + 1 == p.len() && matches!(c.kind(), Kind::Bool(true)) {
                        pairs.push(format!("{} & \\text{{otherwise}}", self.print(x)));
                    } else {
                        pairs.push(format!("{} & \\text{{for}}\\: {}", self.print(x), self.print(c)));
                    }
                }
                format!("\\begin{{cases}} {} \\end{{cases}}", pairs.join(" \\\\"))
            }
        }
    }

    fn print_number(&self, n: &Number) -> String {
        match n {
            Number::Int(_) | Number::BigInt(_) => n.to_string(),
            Number::Rational(r) => {
                let (mut p, q) = (r.numer().clone(), r.denom().clone());
                let mut sign = "";
                if p < num_bigint::BigInt::from(0) {
                    sign = "- ";
                    p = -p;
                }
                if self.fold_short_frac {
                    format!("{sign}{p} / {q}")
                } else {
                    format!("{sign}\\frac{{{p}}}{{{q}}}")
                }
            }
            Number::Float(f) => {
                let s = float_str(*f);
                if let Some((m, x)) = s.split_once('e') {
                    let x = x.strip_prefix('+').unwrap_or(x);
                    format!("{m}{}10^{{{x}}}", self.numsep)
                } else if s == "inf" {
                    r"\infty".into()
                } else if s == "-inf" {
                    r"- \infty".into()
                } else {
                    s
                }
            }
            Number::MpFloat(..) => {
                let s = n.to_string();
                match s.split_once('e') {
                    Some((mantissa, exponent)) => format!("{mantissa}{}10^{{{}}}", self.numsep, exponent.trim_start_matches('+')),
                    None => s,
                }
            }
        }
    }

    fn print_add(&self, e: &Expr) -> String {
        let terms = crate::printer::ordered_terms(e);
        let mut tex = String::new();
        for (i, t) in terms.iter().enumerate() {
            let term = if i == 0 {
                t.clone()
            } else if minus_sign(t) {
                tex.push_str(" - ");
                neg(t)
            } else {
                tex.push_str(" + ");
                t.clone()
            };
            let tt = self.print(&term);
            if matches!(term.kind(), Kind::Rel(..) | Kind::Add(..)) || has_mod(&term) {
                tex.push_str(&Self::parens(&tt));
            } else {
                tex.push_str(&tt);
            }
        }
        tex
    }

    fn needs_mul_brackets(&self, e: &Expr, first: bool) -> bool {
        let base = match e.kind() {
            Kind::Mul(..) => !first && minus_sign(e),
            Kind::Rel(..) | Kind::Piecewise(_) => true,
            _ => precedence_traditional(e) < PREC_MUL,
        };
        base || has_mod(e)
    }

    fn convert_args(&self, args: &[Expr]) -> String {
        let mut tex = String::new();
        let mut last = String::new();
        for (i, term) in args.iter().enumerate() {
            let mut tt = self.print(term);
            if self.needs_mul_brackets(term, i == 0) {
                tt = Self::parens(&tt);
            }
            if ends_with_number(&last) && starts_like_number(&tt) {
                tex.push_str(&self.numsep);
            } else if !tex.is_empty() {
                tex.push_str(&self.sep);
            }
            tex.push_str(&tt);
            last = tt;
        }
        tex
    }

    /// `convert()` of SymPy's `_print_Mul` on a product of `factors`.
    fn convert_list(&self, factors: &[Expr]) -> String {
        match factors.len() {
            0 => "1".into(),
            1 => self.print(&factors[0]),
            _ => {
                let mut v = factors.to_vec();
                super::order::sort_commutative_factors(&mut v);
                // numbers first
                v.sort_by_key(|f| !f.is_number());
                self.convert_args(&v)
            }
        }
    }

    /// `fraction(expr, exact=True)` as factor lists.
    fn fraction(&self, e: &Expr) -> (Vec<Expr>, Vec<Expr>) {
        let mut numer = Vec::new();
        let mut denom = Vec::new();
        let push = |f: Expr, numer: &mut Vec<Expr>, denom: &mut Vec<Expr>| match f.kind() {
            Kind::Pow(..) | Kind::Func(Func::Exp, _) if !crate::decompose::commutes(&f) => numer.push(f),
            Kind::Pow(b, x) => {
                if is_known_negative(x) {
                    if x.is_minus_one() {
                        denom.push(b.clone());
                    } else {
                        denom.push(pow(b, &neg(x)));
                    }
                } else {
                    numer.push(f.clone());
                }
            }
            Kind::Func(Func::Exp, a) if is_known_negative(&a[0]) => {
                denom.push(crate::functions::exp(&neg(&a[0])));
            }
            Kind::Num(Number::Rational(r)) => {
                let p = r.numer().clone();
                if p != num_bigint::BigInt::from(1) {
                    numer.push(Expr::num(Number::from_bigint(p)));
                }
                denom.push(Expr::num(Number::from_bigint(r.denom().clone())));
            }
            _ => numer.push(f),
        };
        match e.kind() {
            Kind::Mul(c, facs) => {
                if !c.is_one() {
                    push(Expr::num(c.clone()), &mut numer, &mut denom);
                }
                for (b, x) in facs {
                    push(mul_factor(b, x), &mut numer, &mut denom);
                }
            }
            _ => push(e.clone(), &mut numer, &mut denom),
        }
        (numer, denom)
    }

    fn print_mul(&self, e: &Expr) -> String {
        let mut tex = String::new();
        let mut include_parens = false;
        let mut expr = e.clone();
        if minus_sign(e) {
            expr = neg(e);
            tex.push_str("- ");
            if expr.is_add() {
                tex.push('(');
                include_parens = true;
            }
        }
        let (numer, denom) = self.fraction(&expr);
        if denom.is_empty() {
            if expr.is_mul() {
                tex.push_str(&self.convert_args(&ordered_factors(&expr)));
            } else {
                tex.push_str(&self.print(&expr));
            }
        } else {
            let snumer = self.convert_list(&numer);
            let sdenom = self.convert_list(&denom);
            let ldenom = sdenom.split_whitespace().count();
            if self.fold_short_frac && ldenom <= 2 && !sdenom.contains('^') {
                let numer_needs = numer.len() == 1 && self.needs_mul_brackets(&numer[0], false);
                if numer_needs {
                    tex.push_str(&format!("\\left({snumer}\\right) / {sdenom}"));
                } else {
                    tex.push_str(&format!("{snumer} / {sdenom}"));
                }
            } else if let Some(ratio) =
                self.s.long_frac_ratio.filter(|r| snumer.split_whitespace().count() as f64 > r * ldenom as f64)
            {
                tex.push_str(&self.long_fraction(&numer, &snumer, &sdenom, ratio * ldenom as f64));
            } else {
                tex.push_str(&format!("\\frac{{{snumer}}}{{{sdenom}}}"));
            }
        }
        if include_parens {
            tex.push(')');
        }
        tex
    }

    /// `convert()` of SymPy's `_print_Mul` on one expression.
    fn convert_expr(&self, e: &Expr) -> String {
        if e.is_mul() {
            self.convert_args(&ordered_factors(e))
        } else {
            self.print(e)
        }
    }

    /// SymPy's `_needs_mul_brackets(expr, last=True)` for a numerator.
    fn needs_mul_brackets_last(&self, e: &Expr) -> bool {
        self.needs_mul_brackets(e, false)
    }

    /// The long-fraction branch of SymPy's `_print_Mul` (`long_frac_ratio`).
    fn long_fraction(&self, numer: &[Expr], snumer: &str, sdenom: &str, width: f64) -> String {
        let sep = &self.sep;
        let numer_expr = crate::ops::mul(numer);
        if self.needs_mul_brackets_last(&numer_expr) {
            return format!("\\frac{{1}}{{{sdenom}}}{sep}\\left({snumer}\\right)");
        }
        if numer_expr.is_mul() {
            // Split a long numerator: the factors that keep it short enough
            // stay above the bar.
            let mut a = Expr::one();
            let mut b = Expr::one();
            for x in ordered_factors(&numer_expr) {
                let ax = crate::ops::mul2(&a, &x);
                let noncommuting = !crate::decompose::commutes(&b) && !crate::decompose::commutes(&x);
                if self.needs_mul_brackets(&x, false)
                    || self.convert_expr(&ax).split_whitespace().count() as f64 > width
                    || noncommuting
                {
                    b = crate::ops::mul2(&b, &x);
                } else {
                    a = ax;
                }
            }
            let (sa, sb) = (self.convert_expr(&a), self.convert_expr(&b));
            if self.needs_mul_brackets_last(&b) {
                return format!("\\frac{{{sa}}}{{{sdenom}}}{sep}\\left({sb}\\right)");
            }
            return format!("\\frac{{{sa}}}{{{sdenom}}}{sep}{sb}");
        }
        format!("\\frac{{1}}{{{sdenom}}}{sep}{snumer}")
    }

    fn print_pow(&self, e: &Expr, b: &Expr, x: &Expr) -> String {
        if let Some((p, q)) = small_rational(x) {
            if p.abs() == 1 && q != 1 && self.s.root_notation {
                let base = self.print(b);
                let tex = if q == 2 {
                    format!("\\sqrt{{{base}}}")
                } else if self.s.itex {
                    format!("\\root{{{q}}}{{{base}}}")
                } else {
                    format!("\\sqrt[{q}]{{{base}}}")
                };
                return if p < 0 { format!("\\frac{{1}}{{{tex}}}") } else { tex };
            } else if self.s.fold_frac_powers && q != 1 {
                let mut base = self.parenthesize(b, PREC_POW, false);
                if b.is_symbol() {
                    base = self.parenthesize_super(base);
                }
                if matches!(b.kind(), Kind::Func(..)) {
                    return self.print_func(b, Some(format!("{p}/{q}")));
                }
                return format!("{base}^{{{p}/{q}}}");
            } else if p < 0 && crate::decompose::commutes(b) {
                return self.print_mul(e);
            }
        } else if x.is_rational_number() && is_known_negative(x) && crate::decompose::commutes(b) {
            // big rational exponents
            return self.print_mul(e);
        }
        if matches!(b.kind(), Kind::Func(..)) {
            // A host `_latex(printer, exp)` such as `adjoint._latex`:
            // the printed base in parentheses raised to the exponent.
            if b.has_host() && matches!(b.kind(), Kind::Func(Func::Undefined(_), _)) {
                if let Some(text) = crate::host::print(b, crate::host::HostFormat::Latex) {
                    return format!("\\left({text}\\right)^{{{}}}", self.print(x));
                }
            }
            return self.print_func(b, Some(self.print(x)));
        }
        let exp = self.print(x);
        let mut base = self.parenthesize(b, PREC_POW, false);
        if b.is_symbol() {
            base = self.parenthesize_super(base);
        } else if matches!(b.kind(), Kind::Num(Number::Float(_))) {
            base = format!("{{{base}}}");
        }
        format!("{base}^{{{exp}}}")
    }

    fn do_exponent(tex: String, exp: Option<String>) -> String {
        match exp {
            Some(x) => format!("\\left({tex}\\right)^{{{x}}}"),
            None => tex,
        }
    }

    fn sup(tex: String, exp: Option<String>) -> String {
        match exp {
            Some(x) => format!("{tex}^{{{x}}}"),
            None => tex,
        }
    }

    fn hprint_function(&self, func: &str) -> String {
        let func = deal_with_super_sub(func);
        let sup = func.find('^');
        let sub = func.find('_');
        if is_accepted(&func) {
            format!("\\{func}")
        } else if func.chars().count() == 1 || func.starts_with('\\') || sub == Some(1) || sup == Some(1) {
            func
        } else {
            let cut = match (sup, sub) {
                (Some(a), Some(b)) if a > 0 && b > 0 => Some(a.min(b)),
                (Some(a), _) if a > 0 => Some(a),
                (_, Some(b)) if b > 0 => Some(b),
                _ => None,
            };
            match cut {
                Some(i) => format!("\\operatorname{{{}}}{}", &func[..i], &func[i..]),
                None => format!("\\operatorname{{{func}}}"),
            }
        }
    }

    /// The reference latex set printers.
    fn print_set(&self, op: crate::sets::SetOp, args: &[Expr]) -> Option<String> {
        use crate::sets::SetOp as O;
        let is_op = |e: &Expr, ops: &[O]| matches!(e.kind(), Kind::Func(Func::Set(o), _) if ops.contains(o));
        let traditional = [O::Union, O::Intersection, O::Complement, O::ProductSet];
        let joined = |delim: &str| {
            args.iter()
                .map(|a| if is_op(a, &traditional) { Self::parens(&self.print(a)) } else { self.print(a) })
                .collect::<Vec<_>>()
                .join(delim)
        };
        Some(match op {
            O::FiniteSet => {
                let mut v = args.to_vec();
                v.sort_by_cached_key(crate::sortkey::sort_key);
                format!("\\left\\{{{}\\right\\}}", v.iter().map(|a| self.print(a)).collect::<Vec<_>>().join(", "))
            }
            O::Interval => {
                let e = Expr::from_kind(Kind::Func(Func::Set(op), args.to_vec()));
                let (a, b, lo, ro) = crate::sets::interval_parts(&e)?;
                if a == b {
                    return Some(format!("\\left\\{{{}\\right\\}}", self.print(&a)));
                }
                format!("\\left{}{}, {}\\right{}", if lo { "(" } else { "[" }, self.print(&a), self.print(&b), if ro { ")" } else { "]" })
            }
            O::Union => joined(" \\cup "),
            O::Intersection => joined(" \\cap "),
            O::Complement => joined(" \\setminus "),
            O::ProductSet if args.len() > 1 && args.iter().any(|a| *a != args[0]) => joined(" \\times "),
            O::Integers => "\\mathbb{Z}".into(),
            O::Naturals => "\\mathbb{N}".into(),
            O::Naturals0 => "\\mathbb{N}_0".into(),
            O::Reals => "\\mathbb{R}".into(),
            O::Complexes => "\\mathbb{C}".into(),
            O::Rationals => "\\mathbb{Q}".into(),
            O::EmptySet => "\\emptyset".into(),
            O::UniversalSet => "\\mathbb{U}".into(),
            O::ImageSet => {
                let Kind::Func(Func::Set(O::Lambda), lam) = args.first()?.kind() else { return None };
                if lam.len() != 2 {
                    return None;
                }
                let sig = crate::special::tuple_items(&lam[0]).map(|v| v.to_vec()).unwrap_or_else(|| vec![lam[0].clone()]);
                let sets = super::pretty::image_base_sets(&sig, &args[1..]);
                let xinys = sig.iter().zip(&sets)
                    .map(|(x, y)| format!("{} \\in {}", self.print(x), self.print(y)))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("\\left\\{{{}\\; \\middle|\\; {}\\right\\}}", self.print(&lam[1]), xinys)
            }
            O::ConditionSet if args.len() == 3 => {
                let syms = crate::special::tuple_items(&args[0]).map(|v| v.to_vec()).unwrap_or_else(|| vec![args[0].clone()]);
                let vars = syms.iter().map(|a| self.print(a)).collect::<Vec<_>>().join(", ");
                if is_op(&args[2], &[O::UniversalSet]) {
                    return Some(format!("\\left\\{{{}\\; \\middle|\\; {} \\right\\}}", vars, self.print(&args[1])));
                }
                format!("\\left\\{{{}\\; \\middle|\\; {} \\in {} \\wedge {} \\right\\}}", vars, vars, self.print(&args[2]), self.print(&args[1]))
            }
            O::Range => {
                let items = super::pretty::range_display(args, true)?;
                let body = items.iter()
                    .map(|i| match i { Some(e) => self.print(e), None => "\\ldots".to_string() })
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("\\left\\{{{}\\right\\}}", body)
            }
            _ => return None,
        })
    }

    fn print_func(&self, e: &Expr, exp: Option<String>) -> String {
        let Kind::Func(f, args) = e.kind() else {
            return self.print(e);
        };
        // `_print_Lambda`: `\left( signature \mapsto body \right)`.
        if let (Func::Set(crate::sets::SetOp::Lambda), [signature, body]) = (f, args.as_slice()) {
            return Self::sup(format!("\\left( {} \\mapsto {} \\right)", self.lambda_signature(signature), self.print(body)), exp);
        }
        if let (Func::Set(op), None) = (f, &exp) {
            if let Some(tex) = self.print_set(*op, args) {
                return tex;
            }
        }
        if let Some((_, polynomial, fun)) = crate::printer::root_sum_parts(e) {
            // The reference's _print_RootSum: the polynomial expression, then the function
            // unless it is the identity.
            let mut parts = vec![self.print(&polynomial)];
            let identity = crate::sets::lambda_parts(&fun)
                .is_some_and(|(variables, body)| variables.len() == 1 && variables[0] == body);
            if !identity {
                parts.push(self.print(&fun));
            }
            let body = format!("\\operatorname{{RootSum}} {{\\left({}\\right)}}", parts.join(", "));
            return match exp { Some(exp) => format!("{body}^{{{exp}}}"), None => body };
        }
        if let Func::Undefined(head) = f {
            if crate::poly::rootof::is_rootof_head(head, args.len()) {
                // The reference's _print_ComplexRootOf.
                let body = format!(
                    "\\operatorname{{CRootOf}} {{\\left({}, {}\\right)}}",
                    self.print(&args[0]),
                    self.print(&args[1])
                );
                return match exp { Some(exp) => format!("{body}^{{{exp}}}"), None => body };
            }
        }
        if let (Some(_), Func::Undefined(head)) = (crate::host::expression_head(e), f) {
            // SymPy's _print_Basic for an Expr subclass that is not a Function.
            let name = deal_with_super_sub(&head.name);
            let body = if args.is_empty() {
                format!("\\text{{{name}}}")
            } else {
                let args = args.iter().map(|arg| self.print(arg)).collect::<Vec<_>>().join(", ");
                format!("\\operatorname{{{name}}}\\left({args}\\right)")
            };
            return match exp { Some(exp) => format!("{body}^{{{exp}}}"), None => body };
        }
        if let Func::Predicate(p) = f {
            return format!("\\operatorname{{Q}}_{{\\text{{{}}}}}", p.name().replace('_', "\\_"));
        }
        if matches!(f, Func::AppliedPredicate) {
            let tail = args[1..].iter().map(|arg|self.print(arg)).collect::<Vec<_>>().join(", ");
            return format!("{}({tail})", self.print(&args[0]));
        }
        if matches!(f,Func::Sum|Func::Product) {
            // SymPy's _print_Sum/_print_Product: one limit as sub/superscripts,
            // several as a substack of inequalities, innermost first.
            let limits=crate::binders::limits(args);
            let op=if *f==Func::Sum {"\\sum"} else {"\\prod"};
            let mut value=if limits.len()==1 {
                format!("{op}_{{{}={}}}^{{{}}} ",self.print(&limits[0][0]),self.print(&limits[0][1]),self.print(&limits[0][2]))
            } else {
                let rows:Vec<String>=limits.iter().map(|l|format!("{} \\leq {} \\leq {}",self.print(&l[1]),self.print(&l[0]),self.print(&l[2]))).collect();
                format!("{op}_{{\\substack{{{}}}}} ",rows.join("\\\\"))
            };
            if matches!(args[0].kind(),Kind::Add(..)) {value.push_str(&format!("\\left({}\\right)",self.print(&args[0])));}
            else {value.push_str(&self.print(&args[0]));}
            return Self::do_exponent(value,exp);
        }
        if matches!(f,Func::Integral) {
            // The reference _print_Integral: `\iint` for up to four
            // indefinite variables, else one `\int\limits` per limit; the
            // integrand is parenthesized strictly below Mul, or always when
            // its coefficient is negative.
            let limits=crate::binders::limits(args);
            let mut value=String::new();
            if limits.len()<=4 && limits.iter().all(|l|l.len()==1) {
                value.push_str(&format!("\\i{}nt","i".repeat(limits.len().saturating_sub(1))));
            } else {
                for limit in limits.iter().rev() {
                    value.push_str("\\int");
                    if limit.len()>1 && self.s.mode!=LatexMode::Inline && !self.s.itex {value.push_str("\\limits");}
                    match limit.len(){
                        2=>value.push_str(&format!("^{{{}}}",self.print(&limit[1]))),
                        3=>value.push_str(&format!("_{{{}}}^{{{}}}",self.print(&limit[1]),self.print(&limit[2]))),_=>{},
                    }
                }
            }
            value.push(' ');
            let negative=matches!(args[0].kind(),Kind::Mul(c,_) if c.is_negative())
                || matches!(args[0].as_number(),Some(n) if n.is_negative());
            if negative {value.push_str(&Self::parens(&self.print(&args[0])));}
            else {value.push_str(&self.parenthesize(&args[0],PREC_MUL,true));}
            for limit in limits {value.push_str(&format!("\\, d{}",self.print(&limit[0])));}
            return Self::do_exponent(value,exp);
        }
        if matches!(f,Func::Order) {
            // The reference _print_Order.
            let mut s=self.print(&args[0]);
            let vars:Vec<&Expr>=args[1..].iter().step_by(2).collect();
            let points:Vec<&Expr>=args[1..].iter().skip(1).step_by(2).collect();
            let tuple=|items:&[&Expr]|format!("\\left( {}\\right)",items.iter().map(|e|self.print(e)).collect::<Vec<_>>().join(", \\  "));
            if points.iter().any(|p|!p.is_zero()) || vars.len()>1 {
                s.push_str("; ");
                if vars.len()>1 {s.push_str(&tuple(&vars));} else if let Some(v)=vars.first() {s.push_str(&self.print(v));}
                s.push_str("\\rightarrow ");
                if points.len()>1 {s.push_str(&tuple(&points));} else if let Some(p)=points.first() {s.push_str(&self.print(p));}
            }
            return Self::sup(format!("O\\left({s}\\right)"),exp);
        }
        if let Func::Special(sf @ (crate::special::SpecialFn::Hyper | crate::special::SpecialFn::MeijerG)) = f {
            if args.len()==3 {
                // The reference _print_hyper and _print_meijerg.
                let vec=|items:&[Expr]|items.iter().map(|e|self.print(e)).collect::<Vec<_>>().join(", ");
                let z=self.print(&args[2]);
                let tex=if *sf==crate::special::SpecialFn::Hyper {
                    let (ap,bq)=(args[0].args(),args[1].args());
                    format!("{{{{}}_{{{}}}F_{{{}}}\\left(\\begin{{matrix}} {} \\\\ {} \\end{{matrix}}\\middle| {{{}}} \\right)}}",ap.len(),bq.len(),vec(&ap),vec(&bq),z)
                } else {
                    let (a,b)=(args[0].args(),args[1].args());
                    if a.len()==2 && b.len()==2 {
                        let (an,ao,bm,bo)=(a[0].args(),a[1].args(),b[0].args(),b[1].args());
                        format!("{{G_{{{}, {}}}^{{{}, {}}}\\left(\\begin{{matrix}} {} & {} \\\\{} & {} \\end{{matrix}} \\middle| {{{}}} \\right)}}",
                            an.len()+ao.len(),bm.len()+bo.len(),bm.len(),an.len(),vec(&an),vec(&ao),vec(&bm),vec(&bo),z)
                    } else { String::new() }
                };
                if !tex.is_empty() {
                    return match exp { Some(x)=>format!("{{{tex}}}^{{{x}}}"), None=>tex };
                }
            }
        }
        if matches!(f,Func::Limit){
            let direction=args[3].as_sym().map(|s|s.name()).unwrap_or("+");
            let suffix=if !matches!(args[2].kind(),Kind::Const(crate::expr::Constant::Infinity|crate::expr::Constant::NegInfinity)) && direction!="+-" {format!("^{direction}")}else{String::new()};
            let head=format!("\\lim_{{{} \\to {}{suffix}}}",self.print(&args[1]),self.print(&args[2]));
            // The reference wraps an AssocOp (Add/Mul) argument in \left(...\right) without a space.
            let tex=if matches!(args[0].kind(),Kind::Add(..)|Kind::Mul(..)){format!("{head}\\left({}\\right)",self.print(&args[0]))}else{format!("{head} {}",self.print(&args[0]))};
            return Self::do_exponent(tex,exp);
        }
        if matches!(f,Func::Subs){
            let variables=args[1].args();let points=args[2].args();
            let entries=variables.iter().zip(&points).map(|(v,p)|format!("{}={}",self.print(v),self.print(p))).collect::<Vec<_>>();
            let lower=format!("\\substack{{ {} }}",entries.join("\\\\ "));
            let value=format!("\\left. {} \\right|_{{{lower}}}",self.print(&args[0]));
            return exp.map_or_else(||value.clone(),|power|format!("\\left({value}\\right)^{{{power}}}"));
        }
        if let Func::Special(sf) = f {
            if sf.is_transform() && args.len() >= 3 {
                let name = sf.name();
                let inverse = name.starts_with("Inverse");
                let letter = match name.trim_start_matches("Inverse") {
                    "LaplaceTransform" => Some("L"),
                    "FourierTransform" => Some("F"),
                    "SineTransform" => Some("SIN"),
                    "CosineTransform" => Some("COS"),
                    "MellinTransform" => Some("M"),
                    _ => None,
                };
                if let Some(letter) = letter {
                    let tex = format!(
                        "\\mathcal{{{letter}}}{}_{{{}}}\\left[{}\\right]\\left({}\\right)",
                        if inverse { "^{-1}" } else { "" },
                        self.print(&args[1]),
                        self.print(&args[0]),
                        self.print(&args[2])
                    );
                    return Self::do_exponent(tex, exp);
                }
            }
        }
        let a0 = || self.print(&args[0]);
        match f {
            // SymPy's _print_DotProduct: precedence Mul - 1, not strict.
            Func::MatrixDotProduct(a,b)=>return Self::sup(format!("{} \\cdot {}",self.mat_parenthesize(a,PREC_MUL-1,false),self.mat_parenthesize(b,PREC_MUL-1,false)),exp),
            Func::MatrixProxy(parent) => return Self::sup(self.print_matrix(parent), exp),
            Func::MatrixTrace(parent) => return Self::sup(format!("\\operatorname{{tr}}\\left({} \\right)", self.print_matrix(parent)), exp),
            Func::MatrixDeterminant(parent) => {
                let inner = match parent.kind() {
                    crate::matrixexpr::MatrixKind::Block(grid) => { let (r, c) = grid.shape(); self.matrix_contents(r, c, (0..r * c).map(|k| self.print_matrix(grid.get(k / c, k % c))).collect()) }
                    _ => self.print_matrix(parent),
                };
                return Self::sup(format!("\\left|{{{inner}}}\\right|"), exp);
            }
            Func::MatrixElement(parent) => return Self::sup(format!("{{{}}}_{{{},{}}}", self.mat_parenthesize(parent, PREC_ATOM, true), self.print(&args[0]), self.print(&args[1])), exp),
            Func::Exp => return Self::do_exponent(format!("e^{{{}}}", a0()), exp),
            Func::Log => {
                let t = if self.s.ln_notation {
                    format!("\\ln{{\\left({} \\right)}}", a0())
                } else {
                    format!("\\log{{\\left({} \\right)}}", a0())
                };
                return Self::sup(t, exp);
            }
            Func::Abs => return Self::sup(format!("\\left|{{{}}}\\right|", a0()), exp),
            Func::Re | Func::Im => {
                let name = if *f == Func::Re { "re" } else { "im" };
                let t = format!("\\operatorname{{{name}}}{{{}}}", self.parenthesize(&args[0], PREC_ATOM, false));
                return Self::do_exponent(t, exp);
            }
            Func::Conjugate => return Self::sup(format!("\\overline{{{}}}", a0()), exp),
            Func::Floor => return Self::sup(format!("\\left\\lfloor{{{}}}\\right\\rfloor", a0()), exp),
            Func::Ceiling => return Self::sup(format!("\\left\\lceil{{{}}}\\right\\rceil", a0()), exp),
            Func::Max | Func::Min => {
                let mut v = args.clone();
                v.sort_by_cached_key(crate::sortkey::sort_key);
                let texargs: Vec<String> = v.iter().map(|a| self.print(a)).collect();
                let name = if *f == Func::Max { "max" } else { "min" };
                return Self::sup(format!("\\{name}\\left({}\\right)", texargs.join(", ")), exp);
            }
            Func::Gamma => {
                let t = match &exp {
                    Some(x) => format!("\\Gamma^{{{x}}}\\left({}\\right)", a0()),
                    None => format!("\\Gamma\\left({}\\right)", a0()),
                };
                return t;
            }
            // SymPy's _print_beta: `\operatorname{B}^{exp}\left(x, y\right)`,
            // the one-argument form printed as beta(x, x).
            Func::Special(crate::special::SpecialFn::Beta) if !args.is_empty() && args.len() <= 2 => {
                let x = a0();
                let y = if args.len() == 1 { x.clone() } else { self.print(&args[1]) };
                return match &exp {
                    Some(e) => format!("\\operatorname{{B}}^{{{e}}}\\left({x}, {y}\\right)"),
                    None => format!("\\operatorname{{B}}\\left({x}, {y}\\right)"),
                };
            }
            Func::Special(sp) => {
                use crate::special::SpecialFn as S;
                match sp {
                    S::Mod => {
                        let t = format!(
                            "{} \\bmod {}",
                            self.parenthesize(&args[0], PREC_MUL, true),
                            self.parenthesize(&args[1], PREC_MUL, true)
                        );
                        return Self::do_exponent(t, exp);
                    }
                    S::KroneckerDelta | S::LeviCivita => {
                        let shown = if *sp == S::KroneckerDelta { &args[..2] } else { &args[..] };
                        let idx: Vec<String> = shown.iter().map(|a| self.print(a)).collect();
                        let sym = if *sp == S::KroneckerDelta { "\\delta" } else { "\\varepsilon" };
                        let sep = if shown.iter().all(|a| a.is_atom()) { " " } else { ", " };
                        return Self::do_exponent(format!("{sym}_{{{}}}", idx.join(sep)), exp);
                    }
                    S::Frac => return Self::sup(format!("\\operatorname{{frac}}{{\\left({}\\right)}}", a0()), exp),
                    _ => {}
                }
            }
            Func::Heaviside | Func::DiracDelta => {
                let name = if *f == Func::Heaviside { "\\theta" } else { "\\delta" };
                // SymPy's Heaviside prints its pargs: H0 unless it is 1/2.
                let t = if *f == Func::Heaviside && args.len() > 1 && args[1] != Expr::rational(1, 2) {
                    format!("{name}\\left({}, {}\\right)", a0(), self.print(&args[1]))
                } else {
                    format!("{name}\\left({}\\right)", a0())
                };
                return match exp {
                    Some(x) => format!("\\left({t}\\right)^{{{x}}}"),
                    None => t,
                };
            }
            _ => {}
        }
        // generic function
        let mut func = f.name().to_string();
        let texargs: Vec<String> = args.iter().map(|a| self.print(a)).collect();
        let mut can_fold = self.s.fold_func_brackets && args.len() == 1 && !self.needs_function_brackets(&args[0]);
        let mut inv_power = false;
        if !matches!(f, Func::Undefined(_)) && INV_TRIG.contains(&func.as_str()) {
            match self.s.inv_trig_style {
                InvTrigStyle::Abbreviated => {}
                InvTrigStyle::Full => {
                    func = format!("{}{}", if func.ends_with('h') { "ar" } else { "arc" }, &func[1..]);
                }
                InvTrigStyle::Power => {
                    func = func[1..].to_string();
                    inv_power = true;
                    if exp.is_some() {
                        can_fold = false;
                    }
                }
            }
        }
        let mut name = if inv_power {
            if is_accepted(&func) {
                format!("\\{func}^{{-1}}")
            } else {
                format!("\\operatorname{{{func}}}^{{-1}}")
            }
        } else if let Some(x) = &exp {
            let ft = self.parenthesize_super(self.hprint_function(&func));
            format!("{ft}^{{{x}}}")
        } else {
            self.hprint_function(&func)
        };
        let joined = texargs.join(",");
        if can_fold {
            if is_accepted(&func) {
                name.push_str(&format!(" {{{joined}}}"));
            } else {
                name.push_str(&joined);
            }
        } else {
            name.push_str(&format!("{{\\left({joined} \\right)}}"));
        }
        if inv_power {
            if let Some(x) = exp {
                name.push_str(&format!("^{{{x}}}"));
            }
        }
        name
    }

    fn needs_brackets(e: &Expr) -> bool {
        match e.kind() {
            Kind::Num(n) => !(n.is_integer() && !n.is_negative()) && !n.is_float(),
            Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
            _ => true,
        }
    }

    fn needs_function_brackets(&self, e: &Expr) -> bool {
        if !Self::needs_brackets(e) {
            return false;
        }
        match e.kind() {
            Kind::Mul(_, facs) => facs.iter().any(|(b, x)| x.is_one() && matches!(b.kind(), Kind::Func(..))),
            Kind::Pow(b, _) => Self::needs_brackets(b),
            Kind::Add(..) | Kind::Func(..) => true,
            _ => false,
        }
    }

    fn print_symbolic_derivative(&self,f:&Expr,vars:&[(Expr,Expr)])->String {
        let d=if crate::traverse::free_symbols(f).len()>1{r"\partial".to_string()}else{self.diff_op.clone()};
        let mut denominator=String::new();
        for(v,n)in vars.iter().rev(){
            if crate::diff::symbolic::one(n){denominator.push_str(&format!("{d} {}",self.print(v)));}
            else{denominator.push_str(&format!("{d} {}^{{{}}}",self.parenthesize_super(self.print(v)),self.print(n)));}
        }
        let total=match crate::diff::symbolic::total(vars){Ok(total)=>total,Err(_)=>return crate::printer::to_str(&crate::diff::symbolic::raw(f,vars))};
        let numerator=if crate::diff::symbolic::one(&total){d}else{format!("{d}^{{{}}}",self.print(&total))};
        let body=if minus_sign(f){Self::parens(&self.print(f))}else{self.parenthesize(f,PREC_MUL,true)};
        format!("\\frac{{{numerator}}}{{{denominator}}} {body}")
    }
    fn print_derivative(&self, f: &Expr, vars: &[(Expr, u32)]) -> String {
        let partial = crate::traverse::free_symbols(f).len() > 1;
        let d = if partial { r"\partial".to_string() } else { self.diff_op.clone() };
        let mut tex = String::new();
        let mut dim = 0u32;
        for (x, n) in vars.iter().rev() {
            dim += n;
            if *n == 1 {
                tex.push_str(&format!("{d} {}", self.print(x)));
            } else {
                tex.push_str(&format!("{d} {}^{{{n}}}", self.parenthesize_super(self.print(x))));
            }
        }
        let head = if dim == 1 {
            format!("\\frac{{{d}}}{{{tex}}}")
        } else {
            format!("\\frac{{{d}^{{{dim}}}}}{{{tex}}}")
        };
        let body = if minus_sign(f) {
            Self::parens(&self.print(f))
        } else {
            self.parenthesize(f, PREC_MUL, true)
        };
        format!("{head} {body}")
    }
}

/// `(p, q)` of a machine size rational number.
// ----------------------------------------------------------------------
// matrix expressions (SymPy's LatexPrinter matrix handlers)
// ----------------------------------------------------------------------

use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand, SliceSpec};

/// `could_extract_minus_sign` of a MatMul: a negative Number coefficient.
fn matmul_negative(c: &Expr) -> bool {
    c.as_number().is_some_and(|n| n.is_negative())
}

/// SymPy's `precedence_traditional` of a matrix expression.
pub(crate) fn matrix_precedence(m: &MatrixExpr) -> i32 {
    matrix_precedence_with(m, true)
}

/// SymPy's `precedence` (`traditional == false`, where Transpose and
/// Adjoint are atoms) or `precedence_traditional` of a matrix expression.
pub(crate) fn matrix_precedence_with(m: &MatrixExpr, traditional: bool) -> i32 {
    match m.kind() {
        MatrixKind::Transpose(_) | MatrixKind::Adjoint(_) | MatrixKind::Conjugate(_) if !traditional => PREC_ATOM,
        MatrixKind::Add(_) => super::PREC_ADD,
        MatrixKind::Mul(c, _) => if matmul_negative(c) { super::PREC_ADD } else { PREC_MUL },
        MatrixKind::Pow(..) | MatrixKind::Inverse(_) | MatrixKind::HadamardPower(..)
        | MatrixKind::Transpose(_) | MatrixKind::Adjoint(_) | MatrixKind::Conjugate(_) => PREC_POW,
        MatrixKind::HadamardProduct(_) | MatrixKind::Kronecker(_) => PREC_MUL,
        _ => PREC_ATOM,
    }
}

/// SymPy's `latexslice`/`ppslice`: the printed parts of a slice, with a
/// unit step dropped and a zero start or full stop left empty.
pub(crate) fn slice_parts(spec: &SliceSpec, dim: &Expr) -> Vec<Option<Expr>> {
    let start = if spec.start.is_zero() { None } else { Some(spec.start.clone()) };
    let full = spec.stop_is_omitted || &spec.stop == dim
        || (spec.stop.is_number() && dim.is_number()
            && crate::ops::sub(&spec.stop, dim).as_number().is_some_and(|d| !d.is_negative()));
    let stop = if full { None } else { Some(spec.stop.clone()) };
    let mut parts = vec![start, stop];
    if !spec.step.is_one() { parts.push(Some(spec.step.clone())); }
    parts
}


/// `precedence_traditional`, whose table holds the class names Add, Mul
/// and Pow at fixed levels (so `-x` is a Mul, not a negative term).
fn class_precedence_traditional(e: &Expr) -> i32 {
    match e.kind() {
        Kind::Add(..) => super::PREC_ADD,
        Kind::Mul(..) => PREC_MUL,
        Kind::Pow(..) => PREC_POW,
        _ => precedence_traditional(e),
    }
}

impl LatexPrinter<'_> {
    fn mat_parenthesize(&self, m: &MatrixExpr, level: i32, strict: bool) -> String {
        let p = matrix_precedence(m);
        if p < level || (!strict && p <= level) { Self::parens(&self.print_matrix(m)) } else { self.print_matrix(m) }
    }

    fn symbol_style(&self, name: &str) -> String {
        if !self.s.mat_symbol_bold || name.contains('{') { return deal_with_super_sub(name); }
        let (head, supers, subs) = split_super_sub(name);
        let mut out = format!("\\mathbf{{{}}}", translate(&head));
        if !supers.is_empty() { out.push_str(&format!("^{{{}}}", supers.iter().map(|x| translate(x)).collect::<Vec<_>>().join(" "))); }
        if !subs.is_empty() { out.push_str(&format!("_{{{}}}", subs.iter().map(|x| translate(x)).collect::<Vec<_>>().join(" "))); }
        out
    }

    /// `_print_matrix_contents`.
    fn matrix_contents(&self, rows: usize, cols: usize, cells: Vec<String>) -> String {
        let lines: Vec<String> = (0..rows).map(|i| cells[i * cols..(i + 1) * cols].join(" & ")).collect();
        let mat_str = match &self.s.mat_str {
            Some(m) => m.clone(),
            None if self.s.mode == LatexMode::Inline => "smallmatrix".into(),
            None if cols <= 10 => "matrix".into(),
            None => "array".into(),
        };
        let spec = if mat_str == "array" { format!("{{{}}}", "c".repeat(cols)) } else { String::new() };
        format!("\\begin{{{mat_str}}}{spec}{}\\end{{{mat_str}}}", lines.join("\\\\"))
    }

    /// `_print_MatrixBase`.
    fn matrix_base(&self, rows: usize, cols: usize, cells: Vec<String>) -> String {
        let body = self.matrix_contents(rows, cols, cells);
        match self.s.mat_delim.as_str() {
            "" => body,
            "(" => format!("\\left({body}\\right)"),
            _ => format!("\\left[{body}\\right]"),
        }
    }

    fn zero_matrix(&self) -> String { if self.s.mat_symbol_bold { "\\mathbf{0}".into() } else { "0".into() } }

    fn operand(&self, o: &MatrixOperand) -> (String, i32) {
        match o {
            MatrixOperand::Scalar(e) => (self.print(e), class_precedence_traditional(e)),
            MatrixOperand::Matrix(m) => (self.print_matrix(m), matrix_precedence(m)),
        }
    }

    /// `_print_MatPow`, `_print_Transpose` and `_print_Adjoint` share this
    /// rule: a compound base gets parentheses, a symbol only when it
    /// already carries a superscript.
    fn matrix_super(&self, base: &MatrixExpr, sup: &str, block_plain: bool) -> String {
        let plain = matches!(base.kind(), MatrixKind::Symbol(_)) || (block_plain && matches!(base.kind(), MatrixKind::Block(_) | MatrixKind::BlockDiag(_)));
        if !plain { return format!("\\left({}\\right)^{{{sup}}}", self.print_matrix(base)); }
        let s = if block_plain { self.mat_parenthesize(base, PREC_POW, true) } else { self.print_matrix(base) };
        if s.contains('^') { format!("\\left({s}\\right)^{{{sup}}}") } else { format!("{s}^{{{sup}}}") }
    }

    fn adjoint_symbol(&self) -> &'static str {
        match self.s.adjoint_style.as_str() { "star" => "\\ast", "hermitian" => "\\mathsf{H}", _ => "\\dagger" }
    }

    /// `_print_MatMul` with the coefficient and matrix factors.
    fn print_matmul(&self, c: &Expr, factors: &[MatrixExpr]) -> String {
        let negative = matmul_negative(c);
        let level = if negative { super::PREC_ADD } else { PREC_MUL };
        let scalar = |e: &Expr| if matches!(e.kind(), Kind::Mul(..)) { self.print(e) } else { self.parenthesize(e, level, false) };
        let mut parts: Vec<String> = Vec::with_capacity(factors.len() + 1);
        let coefficient = if negative { neg(c) } else { c.clone() };
        if !coefficient.is_one() || (!negative && !c.is_one()) { parts.push(scalar(&coefficient)); }
        parts.extend(factors.iter().map(|m| self.mat_parenthesize(m, level, false)));
        let body = parts.join(" ");
        if negative { format!("- {body}") } else { body }
    }

    /// SymPy's LaTeX form of a matrix expression.
    pub(crate) fn print_matrix(&self, m: &MatrixExpr) -> String {
        match m.kind() {
            MatrixKind::Symbol(name) => self.symbol_style(name),
            MatrixKind::Identity | MatrixKind::GenericIdentity => if self.s.mat_symbol_bold { "\\mathbf{I}".into() } else { "\\mathbb{I}".into() },
            MatrixKind::Zero | MatrixKind::GenericZero => self.zero_matrix(),
            MatrixKind::Dft(inverse) => format!("\\text{{{}}}_{{{}}}", if *inverse { "IDFT" } else { "DFT" }, self.print(m.rows())),
            MatrixKind::One => if self.s.mat_symbol_bold { "\\mathbf{1}".into() } else { "1".into() },
            MatrixKind::Explicit(x) => self.matrix_base(x.rows(), x.cols(), x.data().iter().map(|e| self.print(e)).collect()),
            MatrixKind::Block(grid) => { let (r, c) = grid.shape(); self.matrix_base(r, c, (0..r * c).map(|k| self.print_matrix(grid.get(k / c, k % c))).collect()) }
            MatrixKind::BlockDiag(blocks) => {
                let k = blocks.len();
                self.matrix_base(k, k, (0..k * k).map(|p| if p / k == p % k { self.print_matrix(&blocks[p / k]) } else { self.zero_matrix() }).collect())
            }
            // `_print_Add` through the ordered terms.
            MatrixKind::Add(terms) => {
                let mut tex = String::new();
                for (i, term) in super::order::matrix_ordered(terms).iter().enumerate() {
                    let text = match term.kind() {
                        MatrixKind::Mul(c, factors) if i > 0 && matmul_negative(c) => { tex.push_str(" - "); self.print_matmul(&neg(c), factors) }
                        _ => { if i > 0 { tex.push_str(" + "); } self.print_matrix(term) }
                    };
                    if matches!(term.kind(), MatrixKind::Add(_)) { tex.push_str(&Self::parens(&text)); } else { tex.push_str(&text); }
                }
                tex
            }
            MatrixKind::Mul(c, factors) => self.print_matmul(c, factors),
            MatrixKind::Pow(base, exp) => self.matrix_super(base, &self.print(exp), false),
            MatrixKind::Inverse(base) => self.matrix_super(base, "-1", false),
            MatrixKind::Transpose(base) => self.matrix_super(base, "T", true),
            MatrixKind::Adjoint(base) => self.matrix_super(base, self.adjoint_symbol(), true),
            // SymPy holds the conjugate of a matrix expression as Adjoint(Transpose(base)).
            MatrixKind::Conjugate(base) => format!("\\left({}\\right)^{{{}}}", self.matrix_super(base, "T", true), self.adjoint_symbol()),
            MatrixKind::HadamardProduct(factors) => factors.iter().map(|f| self.mat_parenthesize(f, PREC_POW, true)).collect::<Vec<_>>().join(" \\circ "),
            MatrixKind::Kronecker(factors) => factors.iter().map(|f| self.mat_parenthesize(f, PREC_POW, true)).collect::<Vec<_>>().join(" \\otimes "),
            // `_print_HadamardPower` through `_helper_print_standard_power`.
            MatrixKind::HadamardPower(base, exp) => {
                let (exp_tex, exp_prec) = self.operand(exp);
                let (base_tex, base_prec) = self.operand(base);
                let base_tex = if base_prec <= PREC_POW { Self::parens(&base_tex) } else { base_tex };
                if exp_prec < PREC_MUL { format!("{base_tex}^{{\\circ \\left({{{exp_tex}}}\\right)}}") } else { format!("{base_tex}^{{\\circ {{{exp_tex}}}}}") }
            }
            MatrixKind::Slice(parent, rows, cols) => {
                let part = |spec: &SliceSpec, dim: &Expr| slice_parts(spec, dim).iter().map(|x| x.as_ref().map(|e| self.print(e)).unwrap_or_default()).collect::<Vec<_>>().join(":");
                format!("{}\\left[{}, {}\\right]", self.mat_parenthesize(parent, PREC_ATOM, true), part(rows, parent.rows()), part(cols, parent.cols()))
            }
            MatrixKind::ApplyFunction(function, parent) => format!("{{{}}}_{{\\circ}}\\left({{{}}}\\right)", self.print(function), self.print_matrix(parent)),
            _ => m.to_str(),
        }
    }
}


impl LatexPrinter<'_> {
    /// A Lambda signature: a lone symbol, or `_print_tuple` of the
    /// (possibly nested) variables.
    fn lambda_signature(&self, signature: &Expr) -> String {
        match crate::special::tuple_items(signature) {
            Some([single]) if crate::special::tuple_items(single).is_none() => self.print(single),
            Some(items) => self.latex_tuple(items),
            None => self.print(signature),
        }
    }

    fn latex_tuple(&self, items: &[Expr]) -> String {
        let print = |e: &Expr| match crate::special::tuple_items(e) { Some(inner) => self.latex_tuple(inner), None => self.print(e) };
        let sep = ",";
        if items.len() == 1 { return format!("\\left( {}{sep}\\right)", print(&items[0])); }
        format!("\\left( {}\\right)", items.iter().map(print).collect::<Vec<_>>().join(&format!("{sep} \\  ")))
    }
}

fn has_mod(e: &Expr) -> bool {
    crate::special::has_head(e, crate::special::SpecialFn::Mod)
}

fn small_rational(x: &Expr) -> Option<(i64, i64)> {
    match x.kind() {
        Kind::Num(Number::Int(i)) => Some((*i, 1)),
        Kind::Num(Number::Rational(r)) => {
            use num_traits::ToPrimitive;
            Some((r.numer().to_i64()?, r.denom().to_i64()?))
        }
        _ => None,
    }
}

/// `re.search(r'[0-9][} ]*$', s)`
fn ends_with_number(s: &str) -> bool {
    let t = s.trim_end_matches(['}', ' ']);
    t.ends_with(|c: char| c.is_ascii_digit())
}

/// `re.match(r'(\d|\\frac{\d+}{\d+})', tex)`
fn starts_like_number(tex: &str) -> bool {
    if tex.starts_with(|c: char| c.is_ascii_digit()) {
        return true;
    }
    let Some(rest) = tex.strip_prefix("\\frac{") else { return false };
    let d1 = rest.chars().take_while(|c| c.is_ascii_digit()).count();
    if d1 == 0 {
        return false;
    }
    let Some(rest) = rest[d1..].strip_prefix("}{") else { return false };
    let d2 = rest.chars().take_while(|c| c.is_ascii_digit()).count();
    d2 > 0 && rest[d2..].starts_with('}')
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn l(s: &str) -> String {
        latex(&parse_str(s).unwrap())
    }

    #[test]
    fn basics() {
        assert_eq!(l("x/2"), r"\frac{x}{2}");
        assert_eq!(l("-x/2"), r"- \frac{x}{2}");
        assert_eq!(l("sqrt(x)"), r"\sqrt{x}");
        assert_eq!(l("x**(1/3)"), r"\sqrt[3]{x}");
        assert_eq!(l("sin(x)"), r"\sin{\left(x \right)}");
        assert_eq!(l("sin(x)**2"), r"\sin^{2}{\left(x \right)}");
        assert_eq!(l("alpha + beta"), r"\alpha + \beta");
        assert_eq!(l("x1"), r"x_{1}");
        assert_eq!(l("2*3**x"), r"2 \cdot 3^{x}");
        assert_eq!(l("x**2 - 2*x + 1"), r"x^{2} - 2 x + 1");
        assert_eq!(l("x*y/(2*z)"), r"\frac{x y}{2 z}");
        assert_eq!(l("(x+1)/(x-1)"), r"\frac{x + 1}{x - 1}");
        assert_eq!(l("1/sqrt(x)"), r"\frac{1}{\sqrt{x}}");
        assert_eq!(l("exp(x)**y"), r"\left(e^{x}\right)^{y}");
        assert_eq!(l("atan2(y, x)"), r"\operatorname{atan}_{2}{\left(y,x \right)}");
    }

    #[test]
    fn names() {
        assert_eq!(deal_with_super_sub("alphahat"), r"\hat{\alpha}");
        assert_eq!(deal_with_super_sub("theta_i_j"), r"\theta_{i j}");
        assert_eq!(deal_with_super_sub("x^2"), r"x^{2}");
        assert_eq!(deal_with_super_sub("xdot"), r"\dot{x}");
        assert_eq!(deal_with_super_sub("Gamma"), r"\Gamma");
    }
}
