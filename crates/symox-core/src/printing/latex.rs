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
}

impl Default for LatexSettings {
    fn default() -> Self {
        LatexSettings {
            mode: LatexMode::Plain,
            mul_symbol: None,
            fold_short_frac: None,
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
pub(crate) fn translate(s: &str) -> String {
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
        match e.kind() {
            Kind::Num(n) => self.print_number(n),
            Kind::Sym(name) => match self.s.symbol_names.get(&**name) {
                Some(t) => t.clone(),
                None => deal_with_super_sub(name),
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
                crate::printer::sort_factors(&mut v);
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
            } else {
                tex.push_str(&format!("\\frac{{{snumer}}}{{{sdenom}}}"));
            }
        }
        if include_parens {
            tex.push(')');
        }
        tex
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
            } else if p < 0 {
                return self.print_mul(e);
            }
        } else if x.is_rational_number() && is_known_negative(x) {
            // big rational exponents
            return self.print_mul(e);
        }
        if matches!(b.kind(), Kind::Func(..)) {
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

    fn print_func(&self, e: &Expr, exp: Option<String>) -> String {
        let Kind::Func(f, args) = e.kind() else {
            return self.print(e);
        };
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
            let limits=crate::binders::limits(args);
            let mut value=String::new();
            for limit in limits.iter().rev() {
                value.push_str("\\int");match limit.len(){
                    2=>value.push_str(&format!("^{{{}}}",self.print(&limit[1]))),
                    3=>value.push_str(&format!("_{{{}}}^{{{}}}",self.print(&limit[1]),self.print(&limit[2]))),_=>{},
                }value.push(' ');
            }
            value.push_str(&self.parenthesize(&args[0],PREC_MUL,false));
            for limit in limits {value.push_str(&format!("\\, d{}",self.print(&limit[0])));}
            return Self::do_exponent(value,exp);
        }
        if matches!(f,Func::Limit){
            let direction=args[3].as_sym().map(|s|s.name()).unwrap_or("+");
            let suffix=if !matches!(args[2].kind(),Kind::Const(crate::expr::Constant::Infinity|crate::expr::Constant::NegInfinity)) && direction!="+-" {format!("^{{{direction}}}")}else{String::new()};
            let tex=format!("\\lim_{{{} \\to {}{suffix}}} {}",self.print(&args[1]),self.print(&args[2]),self.parenthesize(&args[0],PREC_MUL,false));
            return Self::do_exponent(tex,exp);
        }
        if matches!(f,Func::Subs){
            let variables=args[1].args();let points=args[2].args();
            let entries=variables.iter().zip(&points).map(|(v,p)|format!("{}={}",self.print(v),self.print(p))).collect::<Vec<_>>();
            let lower=if entries.len()==1{entries[0].clone()}else{format!("\\substack{{{}}}",entries.join(" \\\\ "))};
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
            Func::MatrixDotProduct(a,b)=>{let parent=|m:&crate::matrixexpr::MatrixExpr|if matches!(m.kind(),crate::matrixexpr::MatrixKind::Add(_)){format!("\\left({}\\right)",m)}else{m.to_str()};return Self::sup(format!("{} \\cdot {}",parent(a),parent(b)),exp);},
            Func::MatrixProxy(parent) => return Self::sup(parent.to_str(), exp),
            Func::MatrixTrace(parent) => return Self::sup(format!("\\operatorname{{tr}}\\left({}\\right)", parent.to_str()), exp),
            Func::MatrixDeterminant(parent) => return Self::sup(format!("\\det\\left({}\\right)", parent.to_str()), exp),
            Func::MatrixElement(parent) => return Self::sup(format!("{}_{{{}, {}}}", parent.to_str(), self.print(&args[0]), self.print(&args[1])), exp),
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
                v.sort_by_cached_key(|a| (!a.is_number(), crate::printer::to_str(a)));
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
                let t = format!("{name}\\left({}\\right)", a0());
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
