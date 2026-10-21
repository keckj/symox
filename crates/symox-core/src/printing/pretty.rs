//! 2D unicode pretty printer reproducing SymPy's `pretty()`/`pprint()`.
//!
//! Pictures are rectangular blocks of text with a baseline and a binding
//! power, combined with the same box algebra as SymPy's
//! `stringPict`/`prettyForm` (horizontal juxtaposition aligned on
//! baselines, centred stacking, growing parentheses, ...).

use std::rc::Rc;

use super::{float_str, ordered_factors};
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::Number;
use crate::ops::neg;

/// Settings of the pretty printer.
#[derive(Clone, Debug)]
pub struct PrettySettings {
    /// Unicode box drawing / greek letters (default) or plain ASCII.
    pub use_unicode: bool,
    /// Break lines longer than `num_columns` like SymPy (default `true`).
    pub wrap_line: bool,
    /// Terminal width used for line wrapping (SymPy falls back to 80).
    pub num_columns: usize,
    pub root_notation: bool,
    pub use_unicode_sqrt_char: bool,
}

impl Default for PrettySettings {
    fn default() -> Self {
        PrettySettings {
            use_unicode: true,
            wrap_line: true,
            num_columns: 80,
            root_notation: true,
            use_unicode_sqrt_char: true,
        }
    }
}

/// SymPy's `pretty(expr)` (unicode).
pub fn pretty(e: &Expr) -> String {
    pretty_with(e, &PrettySettings::default())
}

/// `pretty(expr, **settings)`.
pub fn pretty_with(e: &Expr, s: &PrettySettings) -> String {
    let _session = super::order::session();
    let _ascii = AsciiGuard::set(!s.use_unicode);
    let p = PrettyPrinter { s };
    let pict = p.print_top(e);
    pict.render(s)
}

/// Render explicit matrix cells using the scalar printer's baseline-aware boxes.
pub fn pretty_matrix_with(matrix: &crate::matrix::Matrix, s: &PrettySettings) -> String {
    let _session = super::order::session();
    let _ascii = AsciiGuard::set(!s.use_unicode);
    let printer = PrettyPrinter { s };
    if matrix.rows() == 0 || matrix.cols() == 0 {
        return "[]".to_owned();
    }
    let cells: Vec<Vec<Pict>> = (0..matrix.rows()).map(|i|
        (0..matrix.cols()).map(|j| printer.print(matrix.get(i, j))).collect()
    ).collect();
    let widths: Vec<usize> = (0..matrix.cols()).map(|j|
        cells.iter().map(|row| row[j].width()).max().unwrap_or(0)
    ).collect();
    let mut lines = Vec::new();
    for row in cells {
        let mut parts = Vec::new();
        for (j, mut cell) in row.into_iter().enumerate() {
            if j != 0 { parts.push(Pict::new("  ")); }
            cell.lines = cell.lines.iter().map(|line| center(line, widths[j])).collect();
            parts.push(cell);
        }
        if !lines.is_empty() { lines.push(String::new()); }
        lines.extend(Pict::next(&parts).lines);
    }
    let width = lines.iter().map(|line| line_width(line)).max().unwrap_or(0);
    let baseline = (lines.len() / 2) as isize;
    let body = Pict {
        lines: lines.iter().map(|line| center(line, width)).collect(),
        baseline, binding: Binding::Atom, func: None,
    };
    body.parens("[", "]").render(s)
}

/// `pprint(expr)`: print the pretty form to stdout.
pub fn pprint(e: &Expr) {
    println!("{}", pretty(e));
}

// ----------------------------------------------------------------------
// pictures
// ----------------------------------------------------------------------

thread_local! {
    static ASCII: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Sets the ASCII mode of the current thread for one printing call.
struct AsciiGuard {
    old: bool,
}

impl AsciiGuard {
    fn set(v: bool) -> AsciiGuard {
        AsciiGuard { old: ASCII.with(|c| c.replace(v)) }
    }
}

impl Drop for AsciiGuard {
    fn drop(&mut self) {
        ASCII.with(|c| c.set(self.old));
    }
}

fn ascii() -> bool {
    ASCII.with(|c| c.get())
}

/// `xsym`: unicode symbol or its ASCII fallback.
fn xsym(u: &'static str, a: &'static str) -> &'static str {
    if ascii() {
        a
    } else {
        u
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Binding {
    Atom = 0,
    Func = 1,
    Div = 2,
    Pow = 3,
    Mul = 4,
    Add = 5,
    Neg = 6,
    Open = 7,
}

#[derive(Clone, Debug)]
struct Pict {
    lines: Vec<String>,
    baseline: isize,
    binding: Binding,
    /// `(prettyFunc, prettyArgs)` of a function application.
    func: Option<Rc<(Pict, Pict)>>,
}

fn is_combining(c: char) -> bool {
    matches!(c as u32, 0x0300..=0x036F | 0x1AB0..=0x1AFF | 0x1DC0..=0x1DFF | 0x20D0..=0x20FF | 0xFE20..=0xFE2F)
}

fn line_width(s: &str) -> usize {
    s.chars().filter(|c| !is_combining(*c)).count()
}

fn center(s: &str, width: usize) -> String {
    let w = line_width(s);
    let d = width.saturating_sub(w);
    let l = d / 2;
    let mut out = String::with_capacity(s.len() + d);
    out.extend(std::iter::repeat(' ').take(l));
    out.push_str(s);
    out.extend(std::iter::repeat(' ').take(d - l));
    out
}

enum Item {
    P(Pict),
    Line,
}

impl Pict {
    fn new(s: &str) -> Pict {
        Pict::with_baseline(s, 0)
    }

    fn with_baseline(s: &str, baseline: isize) -> Pict {
        let raw: Vec<&str> = s.split('\n').collect();
        let width = raw.iter().map(|l| line_width(l)).max().unwrap_or(0);
        let lines = raw.iter().map(|l| center(l, width)).collect();
        Pict { lines, baseline, binding: Binding::Atom, func: None }
    }

    fn bind(mut self, b: Binding) -> Pict {
        self.binding = b;
        self
    }

    fn height(&self) -> usize {
        self.lines.len()
    }

    fn width(&self) -> usize {
        self.lines.first().map(|l| line_width(l)).unwrap_or(0)
    }

    fn text(&self) -> String {
        self.lines.join("\n")
    }

    /// `stringPict.next`
    fn next(objs: &[Pict]) -> Pict {
        let base = objs.iter().map(|o| o.baseline).max().unwrap_or(0);
        let below = objs.iter().map(|o| o.height() as isize - o.baseline).max().unwrap_or(0);
        let h = (base + below).max(0) as usize;
        let mut lines = vec![String::new(); h];
        for o in objs {
            let empty = " ".repeat(o.width());
            let pad = (base - o.baseline) as usize;
            let total = h - o.height();
            let mut col: Vec<&str> = Vec::with_capacity(h);
            col.extend(std::iter::repeat(empty.as_str()).take(pad));
            col.extend(o.lines.iter().map(|s| s.as_str()));
            col.extend(std::iter::repeat(empty.as_str()).take(total.saturating_sub(pad)));
            for (l, c) in lines.iter_mut().zip(col) {
                l.push_str(c);
            }
        }
        Pict { lines, baseline: base, binding: Binding::Atom, func: None }
    }

    fn right(&self, others: &[Pict]) -> Pict {
        let mut v = Vec::with_capacity(others.len() + 1);
        v.push(self.clone());
        v.extend(others.iter().cloned());
        Pict::next(&v)
    }

    fn left(&self, others: &[Pict]) -> Pict {
        let mut v: Vec<Pict> = others.to_vec();
        v.push(self.clone());
        Pict::next(&v)
    }

    /// `stringPict.stack`: baseline is that of the second element.
    fn stack(items: Vec<Item>) -> Pict {
        let width = items
            .iter()
            .filter_map(|i| match i {
                Item::P(p) => Some(p.width()),
                Item::Line => None,
            })
            .max()
            .unwrap_or(0);
        let line = Pict::new(&xsym("─", "-").repeat(width));
        let objs: Vec<Pict> = items
            .into_iter()
            .map(|i| match i {
                Item::P(p) => p,
                Item::Line => line.clone(),
            })
            .collect();
        let mut lines = Vec::new();
        for o in &objs {
            for l in &o.lines {
                lines.push(center(l, width));
            }
        }
        let baseline = objs[0].height() as isize + objs.get(1).map(|o| o.baseline).unwrap_or(0);
        Pict { lines, baseline, binding: Binding::Atom, func: None }
    }

    fn below(&self, items: Vec<Item>) -> Pict {
        let mut v = vec![Item::P(self.clone())];
        v.extend(items);
        let mut p = Pict::stack(v);
        p.baseline = self.baseline;
        p
    }

    fn above(&self, items: Vec<Item>) -> Pict {
        let mut v = items;
        v.push(Item::P(self.clone()));
        let mut p = Pict::stack(v);
        p.baseline = p.height() as isize - self.height() as isize + self.baseline;
        p
    }

    fn parens(&self, left: &str, right: &str) -> Pict {
        let h = self.height();
        let b = self.baseline;
        let mut res = self.clone();
        if !left.is_empty() {
            let l = Pict::with_baseline(&vobj(left, h), b);
            res = l.right(&[res]);
        }
        if !right.is_empty() {
            let r = Pict::with_baseline(&vobj(right, h), b);
            res = res.right(&[r]);
        }
        res.binding = Binding::Atom;
        res.func = None;
        res
    }

    fn paren(&self) -> Pict {
        self.parens("(", ")")
    }

    // --------------------------------------------------------------
    // prettyForm algebra
    // --------------------------------------------------------------

    fn add(forms: Vec<Pict>) -> Pict {
        let mut result = Vec::with_capacity(forms.len() * 2);
        for (i, mut a) in forms.into_iter().enumerate() {
            if a.binding > Binding::Neg {
                a = a.paren();
            }
            if i > 0 && a.binding != Binding::Neg {
                result.push(Pict::new(" + "));
            }
            result.push(a);
        }
        Pict::next(&result).bind(Binding::Add)
    }

    fn mul(forms: Vec<Pict>) -> Pict {
        if forms.len() == 1 {
            return forms.into_iter().next().unwrap();
        }
        let mut result: Vec<Pict> = Vec::with_capacity(forms.len() * 2);
        for (i, mut a) in forms.into_iter().enumerate() {
            if i > 0 {
                result.push(Pict::new(xsym("⋅", "*")));
            }
            if a.binding > Binding::Mul && a.binding != Binding::Neg {
                a = a.paren();
            }
            result.push(a);
        }
        if result.len() > 1 && result[0].text() == "-1" && result[1].text() == xsym("⋅", "*") {
            result.remove(0);
            result[0] = Pict::new("-");
        }
        let first_neg = result[0].lines.first().map(|l| l.starts_with('-')).unwrap_or(false);
        let bin = if first_neg {
            if result[0].text() == "-" && result.len() > 1 {
                let r = &result[1];
                let bl = r.baseline.max(0) as usize;
                if r.lines.get(bl).map(|l| l.starts_with('-')).unwrap_or(false) {
                    result[0] = Pict::new("- ");
                }
            }
            Binding::Neg
        } else {
            Binding::Mul
        };
        Pict::next(&result).bind(bin)
    }

    fn div(num: Pict, den: Pict) -> Pict {
        let mut num = num;
        let mut den = den;
        if num.binding == Binding::Div {
            num = num.paren();
        }
        if den.binding == Binding::Div {
            den = den.paren();
        }
        if num.binding == Binding::Neg {
            num = num.right(&[Pict::new(" ")]);
        }
        Pict::stack(vec![Item::P(num), Item::Line, Item::P(den)]).bind(Binding::Div)
    }

    fn pow(a: Pict, b: Pict) -> Pict {
        let mut a = a;
        let mut b = b;
        if b.binding == Binding::Pow {
            b = b.paren();
        }
        let mut inline = false;
        if a.binding > Binding::Func {
            a = a.paren();
        } else if a.binding == Binding::Func {
            if b.height() > 1 {
                a = a.paren();
            } else {
                inline = true;
            }
        }
        if inline {
            let parts = a.func.clone().unwrap();
            let (pf, pa) = (&parts.0, &parts.1);
            b.baseline = pf.baseline + b.height() as isize;
            let f = pf.right(&[b]);
            return f.right(&[pa.clone()]);
        }
        let top = b.left(&[Pict::new(&" ".repeat(a.width()))]);
        let bot = a.right(&[Pict::new(&" ".repeat(b.width()))]);
        bot.above(vec![Item::P(top)]).bind(Binding::Pow)
    }

    fn render(&self, s: &PrettySettings) -> String {
        if !s.wrap_line || self.width() <= s.num_columns {
            return self.text();
        }
        let ncols = if s.num_columns == 0 { 80 } else { s.num_columns };
        let (ostart, oend) = if ascii() { ("> ", " >") } else { ("↪ ", " ↪") };
        let chunks = |line: &str| -> Vec<String> {
            let mut out = Vec::new();
            let mut prefix = String::new();
            let mut width = line_width(&format!("{prefix}{oend}"));
            let chars: Vec<char> = line.chars().collect();
            let mut start = 0usize;
            for (i, c) in chars.iter().enumerate() {
                let wx = if is_combining(*c) { 0 } else { 1 };
                if width + wx > ncols {
                    let piece: String = chars[start..i].iter().collect();
                    out.push(format!("{prefix}{piece}{oend}"));
                    prefix = ostart.to_string();
                    width = line_width(&format!("{prefix}{oend}"));
                    start = i;
                }
                width += wx;
            }
            let piece: String = chars[start..].iter().collect();
            out.push(format!("{prefix}{piece}"));
            out
        };
        let per_line: Vec<Vec<String>> = self.lines.iter().map(|l| chunks(l)).collect();
        let n = per_line.iter().map(|c| c.len()).min().unwrap_or(0);
        let mut screens = Vec::with_capacity(n);
        for k in 0..n {
            let screen: Vec<&str> = per_line.iter().map(|c| c[k].as_str()).collect();
            screens.push(screen.join("\n"));
        }
        screens.join("\n\n")
    }
}

/// Vertical delimiter of height `h` (`xobj`).
fn vobj(sym: &str, h: usize) -> String {
    // (ext, top, bot, mid, single)
    let (ext, top, bot, mid, c1): (&str, &str, &str, Option<&str>, &str) = if ascii() {
        match sym {
            "(" => ("|", "/", "\\", None, "("),
            ")" => ("|", "\\", "/", None, ")"),
            "[" => ("[", "[", "[", None, "["),
            "]" => ("]", "]", "]", None, "]"),
            "{" => ("|", "/", "\\", Some("<"), "{"),
            "}" => ("|", "\\", "/", Some(">"), "}"),
            "|" => ("|", "|", "|", None, "|"),
            other => (other, other, other, None, other),
        }
    } else {
        match sym {
        "(" => ("⎜", "⎛", "⎝", None, "("),
        ")" => ("⎟", "⎞", "⎠", None, ")"),
        "[" => ("⎢", "⎡", "⎣", None, "["),
        "]" => ("⎥", "⎤", "⎦", None, "]"),
        "{" => ("⎪", "⎧", "⎩", Some("⎨"), "{"),
        "}" => ("⎪", "⎫", "⎭", Some("⎬"), "}"),
        "|" => ("│", "│", "│", None, "│"),
        "lfloor" => ("⎢", "⎢", "⎣", None, "⌊"),
        "rfloor" => ("⎥", "⎥", "⎦", None, "⌋"),
        "lceil" => ("⎢", "⎡", "⎢", None, "⌈"),
        "rceil" => ("⎥", "⎤", "⎥", None, "⌉"),
        other => (other, other, other, None, other),
        }
    };
    let mut length = h.max(1);
    let midc = match mid {
        Some(m) => {
            if length % 2 == 0 {
                length += 1;
            }
            m
        }
        None => ext,
    };
    if length == 1 {
        return c1.to_string();
    }
    let next = (length - 2) / 2;
    let nmid = (length - 2) - next * 2;
    let mut res: Vec<&str> = vec![top];
    res.extend(std::iter::repeat(ext).take(next));
    res.extend(std::iter::repeat(midc).take(nmid));
    res.extend(std::iter::repeat(ext).take(next));
    res.push(bot);
    res.join("\n")
}

// ----------------------------------------------------------------------
// symbol names
// ----------------------------------------------------------------------

fn greek_unicode(s: &str) -> Option<&'static str> {
    Some(match s {
        "alpha" => "α",
        "beta" => "β",
        "gamma" => "γ",
        "delta" => "δ",
        "epsilon" => "ε",
        "zeta" => "ζ",
        "eta" => "η",
        "theta" => "θ",
        "iota" => "ι",
        "kappa" => "κ",
        "lambda" | "lamda" => "λ",
        "mu" => "μ",
        "nu" => "ν",
        "xi" => "ξ",
        "omicron" => "ο",
        "pi" => "π",
        "rho" => "ρ",
        "sigma" => "σ",
        "tau" => "τ",
        "upsilon" => "υ",
        "phi" => "φ",
        "chi" => "χ",
        "psi" => "ψ",
        "omega" => "ω",
        "varsigma" => "ς",
        "Alpha" => "Α",
        "Beta" => "Β",
        "Gamma" => "Γ",
        "Delta" => "Δ",
        "Epsilon" => "Ε",
        "Zeta" => "Ζ",
        "Eta" => "Η",
        "Theta" => "Θ",
        "Iota" => "Ι",
        "Kappa" => "Κ",
        "Lambda" | "Lamda" => "Λ",
        "Mu" => "Μ",
        "Nu" => "Ν",
        "Xi" => "Ξ",
        "Omicron" => "Ο",
        "Pi" => "Π",
        "Rho" => "Ρ",
        "Sigma" => "Σ",
        "Tau" => "Τ",
        "Upsilon" => "Υ",
        "Phi" => "Φ",
        "Chi" => "Χ",
        "Psi" => "Ψ",
        "Omega" => "Ω",
        _ => return None,
    })
}

fn sub_char(s: &str) -> Option<&'static str> {
    Some(match s {
        "(" => "₍",
        ")" => "₎",
        "+" => "₊",
        "-" => "₋",
        "=" => "₌",
        "0" => "₀",
        "1" => "₁",
        "2" => "₂",
        "3" => "₃",
        "4" => "₄",
        "5" => "₅",
        "6" => "₆",
        "7" => "₇",
        "8" => "₈",
        "9" => "₉",
        "a" => "ₐ",
        "e" => "ₑ",
        "h" => "ₕ",
        "i" => "ᵢ",
        "k" => "ₖ",
        "l" => "ₗ",
        "m" => "ₘ",
        "n" => "ₙ",
        "o" => "ₒ",
        "p" => "ₚ",
        "r" => "ᵣ",
        "s" => "ₛ",
        "t" => "ₜ",
        "u" => "ᵤ",
        "v" => "ᵥ",
        "x" => "ₓ",
        "beta" => "ᵦ",
        "gamma" => "ᵧ",
        "rho" => "ᵨ",
        "phi" => "ᵩ",
        "chi" => "ᵪ",
        _ => return None,
    })
}

fn sup_char(s: &str) -> Option<&'static str> {
    Some(match s {
        "(" => "⁽",
        ")" => "⁾",
        "+" => "⁺",
        "-" => "⁻",
        "=" => "⁼",
        "0" => "⁰",
        "1" => "¹",
        "2" => "²",
        "3" => "³",
        "4" => "⁴",
        "5" => "⁵",
        "6" => "⁶",
        "7" => "⁷",
        "8" => "⁸",
        "9" => "⁹",
        "i" => "ⁱ",
        "n" => "ⁿ",
        _ => return None,
    })
}

fn center_accent(s: &str, accent: char) -> String {
    let chars: Vec<char> = s.chars().collect();
    let mid = (chars.len() / 2 + 1).min(chars.len());
    let mut out: String = chars[..mid].iter().collect();
    out.push(accent);
    out.extend(chars[mid..].iter());
    out
}

const PRETTY_MODIFIERS: [&str; 19] = [
    "mathring", "ddddot", "dddot", "check", "breve", "acute", "grave", "tilde", "prime", "ddot",
    "norm", "dot", "hat", "bar", "vec", "prm", "avg", "abs", "mag",
];

fn apply_modifier(key: &str, s: String) -> String {
    match key {
        "mathring" => center_accent(&s, '\u{030A}'),
        "ddddot" => center_accent(&s, '\u{20DC}'),
        "dddot" => center_accent(&s, '\u{20DB}'),
        "ddot" => center_accent(&s, '\u{0308}'),
        "dot" => center_accent(&s, '\u{0307}'),
        "check" => center_accent(&s, '\u{030C}'),
        "breve" => center_accent(&s, '\u{0306}'),
        "acute" => center_accent(&s, '\u{0301}'),
        "grave" => center_accent(&s, '\u{0300}'),
        "tilde" => center_accent(&s, '\u{0303}'),
        "hat" => center_accent(&s, '\u{0302}'),
        "bar" => center_accent(&s, '\u{0305}'),
        "vec" => center_accent(&s, '\u{20D7}'),
        "prime" | "prm" => format!("{s}′"),
        "norm" => format!("‖{s}‖"),
        "avg" => format!("⟨{s}⟩"),
        "abs" | "mag" => format!("|{s}|"),
        _ => s,
    }
}

fn translate(s: &str) -> String {
    if let Some(g) = greek_unicode(s) {
        return g.to_string();
    }
    let lower = s.to_lowercase();
    for key in PRETTY_MODIFIERS {
        if lower.ends_with(key) && s.chars().count() > key.len() {
            return apply_modifier(key, translate(&s[..s.len() - key.len()]));
        }
    }
    s.to_string()
}

/// `pretty_symbol(name)`.
fn pretty_symbol(name: &str) -> String {
    if ascii() {
        return name.to_string();
    }
    let (n, sups, subs) = super::latex::split_super_sub(name);
    let n = translate(&n);
    let pretty_list = |l: &[String], f: fn(&str) -> Option<&'static str>| -> Option<Vec<String>> {
        let mut out = Vec::new();
        for s in l {
            if let Some(p) = f(s) {
                out.push(p.to_string());
                continue;
            }
            let mut acc = String::new();
            for c in s.chars() {
                acc.push_str(f(&c.to_string())?);
            }
            out.push(acc);
        }
        Some(out)
    };
    let psups = pretty_list(&sups, sup_char);
    let psubs = if psups.is_some() { pretty_list(&subs, sub_char) } else { None };
    match (psups, psubs) {
        (Some(a), Some(b)) => format!("{n}{}{}", a.join(" "), b.join(" ")),
        _ => {
            let mut out = n;
            if !subs.is_empty() {
                let v: Vec<String> = subs.iter().map(|s| translate(s)).collect();
                out.push('_');
                out.push_str(&v.join("_"));
            }
            if !sups.is_empty() {
                let v: Vec<String> = sups.iter().map(|s| translate(s)).collect();
                out.push_str("__");
                out.push_str(&v.join("__"));
            }
            out
        }
    }
}

// ----------------------------------------------------------------------
// the printer
// ----------------------------------------------------------------------

struct PrettyPrinter<'a> {
    s: &'a PrettySettings,
}

impl PrettyPrinter<'_> {
    fn print_top(&self, e: &Expr) -> Pict {
        if let Kind::Num(Number::Float(f)) = e.kind() {
            return Pict::new(&crate::number::format_float_full(*f));
        }
        self.print(e)
    }

    fn print(&self, e: &Expr) -> Pict {
        match e.kind() {
            Kind::Num(n) => self.print_number(n),
            Kind::Sym(s) => Pict::new(&pretty_symbol(&crate::host::atom_name(e, s))),
            Kind::Const(c) => Pict::new(match c {
                Constant::Pi => xsym("π", "pi"),
                Constant::E => xsym("ℯ", "e"),
                Constant::I => xsym("ⅈ", "I"),
                Constant::EulerGamma => xsym("γ", "EulerGamma"),
                Constant::Catalan => "G",
                Constant::GoldenRatio => xsym("φ", "GoldenRatio"),
                Constant::TribonacciConstant => "TribonacciConstant",
                Constant::Infinity => xsym("∞", "oo"),
                Constant::NegInfinity => xsym("-∞", "-oo"),
                Constant::ComplexInfinity => "zoo",
                Constant::NaN => "nan",
            }),
            Kind::Bool(b) => Pict::new(if *b { "True" } else { "False" }),
            Kind::Add(..) => self.print_add(e),
            Kind::Mul(..) => self.print_mul(e),
            Kind::Pow(b, x) => self.print_pow(b, x),
            Kind::Func(f, args) => self.print_func(f, args),
            Kind::Derivative(f, vars) => self.print_derivative(f, vars),
            Kind::Rel(op, a, b) => {
                let o = match op {
                    RelOp::Eq => "=",
                    RelOp::Ne => xsym("≠", "!="),
                    RelOp::Lt => "<",
                    RelOp::Le => xsym("≤", "<="),
                    RelOp::Gt => ">",
                    RelOp::Ge => xsym("≥", ">="),
                };
                Pict::next(&[self.print(a), Pict::new(&format!(" {o} ")), self.print(b)]).bind(Binding::Open)
            }
            Kind::Piecewise(p) => self.print_piecewise(p),
        }
    }

    fn print_number(&self, n: &Number) -> Pict {
        match n {
            Number::Int(_) | Number::BigInt(_) => {
                let p = Pict::new(&n.to_string());
                if n.is_negative() {
                    p.bind(Binding::Neg)
                } else {
                    p
                }
            }
            Number::Rational(r) => {
                use num_traits::Signed;
                let (p, q) = (r.numer(), r.denom());
                let ten = num_bigint::BigInt::from(10);
                if p.abs() >= ten && q.abs() >= ten {
                    let mut num = Pict::new(&p.to_string());
                    if p.is_negative() {
                        num = num.bind(Binding::Neg);
                    }
                    Pict::div(num, Pict::new(&q.to_string()))
                } else {
                    Pict::new(&format!("{p}/{q}"))
                }
            }
            Number::Float(f) => Pict::new(&float_str(*f)),
            Number::MpFloat(..) => Pict::new(&n.to_string()),
        }
    }

    fn pretty_negative(p: Pict, index: usize) -> Pict {
        let neg = if index == 0 {
            if p.height() > 1 {
                "- "
            } else {
                "-"
            }
        } else {
            " - "
        };
        let p = if p.binding > Binding::Neg || p.binding == Binding::Add { p.paren() } else { p };
        Pict::next(&[Pict::new(neg), p]).bind(Binding::Neg)
    }

    fn print_add(&self, e: &Expr) -> Pict {
        let terms = crate::printer::ordered_terms(e);
        let mut forms: Vec<Option<Pict>> = Vec::with_capacity(terms.len());
        let mut deferred = Vec::new();
        for (i, t) in terms.iter().enumerate() {
            match t.kind() {
                Kind::Mul(c, _) if c.is_negative() => {
                    forms.push(Some(Self::pretty_negative(self.print(&neg(t)), i)));
                }
                Kind::Num(Number::Rational(_)) => {
                    forms.push(None);
                    deferred.push(i);
                }
                Kind::Num(n) if n.is_negative() => {
                    forms.push(Some(Self::pretty_negative(self.print_number(&n.neg()), i)));
                }
                Kind::Const(Constant::NegInfinity) => {
                    forms.push(Some(Self::pretty_negative(Pict::new(xsym("∞", "oo")), i)));
                }
                Kind::Rel(..) => forms.push(Some(self.print(t).paren())),
                _ => forms.push(Some(self.print(t))),
            }
        }
        if !deferred.is_empty() {
            let large = forms.iter().flatten().any(|p| p.height() > 1);
            for i in deferred {
                let Kind::Num(n) = terms[i].kind() else { unreachable!() };
                let (n, negative) = if n.is_negative() { (n.neg(), true) } else { (n.clone(), false) };
                let mut p = if large {
                    let (pp, qq) = n.numer_denom().unwrap();
                    Pict::div(Pict::new(&pp.to_string()), Pict::new(&qq.to_string()))
                } else {
                    self.print_number(&n)
                };
                if negative {
                    p = Self::pretty_negative(p, i);
                }
                forms[i] = Some(p);
            }
        }
        Pict::add(forms.into_iter().flatten().collect())
    }

    fn print_mul(&self, e: &Expr) -> Pict {
        let mut a: Vec<Expr> = Vec::new();
        let mut b: Vec<Den> = Vec::new();
        for item in ordered_factors(e) {
            match item.kind() {
                Kind::Pow(base, x) if matches!(x.as_number(), Some(n) if n.is_rational() && n.is_negative()) => {
                    if x.is_minus_one() {
                        b.push(Den::E(base.clone()));
                    } else {
                        b.push(Den::Pow(base.clone(), neg(x)));
                    }
                }
                Kind::Num(n) if n.is_rational() => {
                    let (p, q) = n.numer_denom().unwrap();
                    if p != num_bigint::BigInt::from(1) {
                        a.push(Expr::num(Number::from_bigint(p)));
                    }
                    if q != num_bigint::BigInt::from(1) {
                        b.push(Den::E(Expr::num(Number::from_bigint(q))));
                    }
                }
                _ => a.push(item),
            }
        }
        let mut pa: Vec<Pict> = a.iter().map(|x| self.print(x)).collect();
        let pb: Vec<Pict> = b
            .iter()
            .map(|d| match d {
                Den::E(x) => self.print(x),
                Den::Pow(base, x) => self.print_pow(base, x),
            })
            .collect();
        if pb.is_empty() {
            if pa.is_empty() {
                return Pict::new("1");
            }
            Pict::mul(pa)
        } else {
            if pa.is_empty() {
                pa.push(Pict::new("1"));
            }
            Pict::div(Pict::mul(pa), Pict::mul(pb))
        }
    }

    fn nth_root(&self, base: &Expr, root: &Expr) -> Pict {
        let bp = self.print(base);
        if self.s.use_unicode_sqrt_char
            && !ascii()
            && root.as_i64() == Some(2)
            && bp.height() == 1
            && (bp.width() == 1 || matches!(base.as_number(), Some(n) if n.is_integer() && !n.is_negative()))
        {
            return bp.left(&[Pict::new("√")]);
        }
        let rp = self.print(root);
        if rp.height() != 1 {
            return Pict::pow(bp, self.print(&crate::ops::inv(root)));
        }
        let exp = if root.as_i64() == Some(2) {
            String::new()
        } else {
            let mut s = rp.text();
            while line_width(&s) < 2 {
                s.push(' ');
            }
            s
        };
        let mut rootsign = xsym("╲╱", "\\/").to_string();
        if line_width(&exp) > 2 {
            rootsign = format!("{}{rootsign}", " ".repeat(line_width(&exp) - 2));
        }
        let mut rs = Pict::with_baseline(&format!("{exp}\n{rootsign}"), 0);
        let linelength = bp.height() as isize - 1;
        let diag_text: Vec<String> = (0..linelength.max(0))
            .map(|i| format!("{}{}{}", " ".repeat((linelength - i - 1) as usize), xsym("╱", "/"), " ".repeat(i as usize)))
            .collect();
        let diag = Pict::with_baseline(&diag_text.join("\n"), linelength - 1);
        rs = rs.right(&[diag]);
        rs.baseline = 1.max(bp.baseline);
        let bar = Pict::new(&"_".repeat(2 + bp.width()));
        let s = bp.above(vec![Item::P(bar)]);
        s.left(&[rs])
    }

    fn print_pow(&self, b: &Expr, x: &Expr) -> Pict {
        // as_base_exp: (1/q)**e -> q**(-e)
        if let Kind::Num(Number::Rational(r)) = b.kind() {
            if r.numer() == &num_bigint::BigInt::from(1) {
                let q = Expr::num(Number::from_bigint(r.denom().clone()));
                return self.print_pow(&q, &neg(x));
            }
        }
        if x.is_minus_one() {
            return Pict::div(Pict::new("1"), self.print(b));
        }
        if self.s.root_notation {
            // fraction(e) == (1, d) with d an atom
            match x.kind() {
                Kind::Num(Number::Rational(r)) if r.numer() == &num_bigint::BigInt::from(1) => {
                    return self.nth_root(b, &Expr::num(Number::from_bigint(r.denom().clone())));
                }
                Kind::Pow(d, m) if m.is_minus_one() && d.is_symbol() => return self.nth_root(b, d),
                _ => {}
            }
        }
        if let Some(n) = x.as_number() {
            if n.is_rational() && n.is_negative() {
                let pos = Expr::from_kind(Kind::Pow(b.clone(), neg(x)));
                return Pict::div(Pict::new("1"), self.print(&pos));
            }
        }
        let base = if matches!(b.kind(), Kind::Rel(..)) { self.print(b).paren() } else { self.print(b) };
        Pict::pow(base, self.print(x))
    }

    fn seq(&self, args: &[Expr]) -> Pict {
        let mut forms = Vec::with_capacity(args.len() * 2);
        for (i, a) in args.iter().enumerate() {
            if i > 0 {
                forms.push(Pict::new(", "));
            }
            forms.push(self.print(a));
        }
        if forms.is_empty() {
            return Pict::new("");
        }
        Pict::next(&forms)
    }

    fn function(&self, name: &str, args: &[Expr]) -> Pict {
        let pf = Pict::new(&pretty_symbol(name));
        let pa = self.seq(args).paren();
        let mut p = Pict::next(&[pf.clone(), pa.clone()]).bind(Binding::Func);
        p.func = Some(Rc::new((pf, pa)));
        p
    }

    /// SymPy's `_print_Product`: one box sign per limit, innermost first.
    fn print_product(&self, args: &[Expr]) -> Pict {
        let mut pretty_func = self.print(&args[0]);
        let (horizontal, corner, vertical) = if ascii() { ("_", "_", "|") } else { ("─", "┬", "│") };
        let func_height = pretty_func.height();
        let (mut first, mut max_upper, mut sign_height) = (true, 0usize, 0usize);
        for limit in crate::binders::limits(args) {
            let lower = Pict::next(&[self.print(&limit[0]), Pict::new(" = "), self.print(&limit[1])]);
            let upper = self.print(&limit[2]);
            let width = (func_height + 2) * 5 / 3 - 2;
            let mut lines = vec![Item::P(Pict::new(&format!("{horizontal}{corner}{}{corner}{horizontal}", horizontal.repeat(width - 2))))];
            for _ in 0..=func_height {
                lines.push(Item::P(Pict::new(&format!(" {vertical}{}{vertical} ", " ".repeat(width - 2)))));
            }
            let mut sign = Pict::stack(lines);
            max_upper = max_upper.max(upper.height());
            if first { sign_height = sign.height(); }
            sign = sign.above(vec![Item::P(upper)]);
            sign = sign.below(vec![Item::P(lower)]);
            if first { pretty_func.baseline = 0; first = false; }
            let padding = Pict::stack((0..sign.height() - 1).map(|_| Item::P(Pict::new(" "))).collect());
            sign = sign.right(&[padding]);
            pretty_func = sign.right(&[pretty_func]);
        }
        pretty_func.baseline = (max_upper + sign_height / 2) as isize;
        pretty_func.bind(Binding::Mul)
    }

    fn print_func(&self, f: &Func, args: &[Expr]) -> Pict {
        if *f == Func::Product { return self.print_product(args); }
        if let Func::Undefined(head) = f {
            if crate::poly::rootof::is_rootof_head(head, args.len()) {
                // The reference's _print_ComplexRootOf: an atom.
                return Pict::next(&[Pict::new("CRootOf"), self.seq(args).paren()]);
            }
        }
        if matches!(f,Func::Predicate(_) | Func::AppliedPredicate | Func::Subs | Func::Limit | Func::Integral | Func::Sum){return Pict::new(&crate::printer::to_str(&Expr::from_kind(Kind::Func(f.clone(),args.to_vec()))));}
        if let Func::MatrixElement(parent) = f {
            return Pict::new(&format!("{}[{}, {}]",parent.to_str(),args[0],args[1]));
        }
        if let Func::MatrixDotProduct(a,b)=f{return Pict::new(&format!("DotProduct({}, {})",a,b));}
        if let Func::MatrixProxy(parent) = f { return Pict::new(&parent.to_str()); }
        if let Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent) = f {
            return Pict::new(&format!("{}({})",f.name(),parent.to_str()));
        }
        match f {
            Func::Exp => Pict::pow(Pict::new(xsym("ℯ", "e")), self.print(&args[0])),
            Func::Abs => self.print(&args[0]).parens("|", "|"),
            Func::Floor if !ascii() => self.print(&args[0]).parens("lfloor", "rfloor"),
            Func::Ceiling if !ascii() => self.print(&args[0]).parens("lceil", "rceil"),
            Func::Conjugate => {
                let p = self.print(&args[0]);
                let bar = Pict::new(&"_".repeat(p.width()));
                p.above(vec![Item::P(bar)])
            }
            Func::Atan2 => self.seq(args).paren().left(&[Pict::new("atan2")]),
            Func::Gamma => {
                let pf = Pict::new(xsym("Γ", "Gamma"));
                let pa = self.seq(args).paren();
                let mut p = Pict::next(&[pf.clone(), pa.clone()]).bind(Binding::Func);
                p.func = Some(Rc::new((pf, pa)));
                p
            }
            Func::Heaviside => self.print(&args[0]).paren().left(&[Pict::new(xsym("θ", "Heaviside"))]),
            Func::Special(crate::special::SpecialFn::Mod) => {
                let mut p = self.print(&args[0]);
                if p.binding > Binding::Mul {
                    p = p.paren();
                }
                p.right(&[Pict::new(" mod "), self.print(&args[1])]).bind(Binding::Open)
            }
            Func::Special(crate::special::SpecialFn::KroneckerDelta) => {
                let b = self.print(&args[0]).right(&[Pict::new(","), self.print(&args[1])]);
                let a = Pict::new(xsym("δ", "d"));
                let top = b.left(&[Pict::new(&" ".repeat(a.width()))]);
                let bot = a.right(&[Pict::new(&" ".repeat(b.width()))]);
                bot.below(vec![Item::P(top)]).bind(Binding::Pow)
            }
            Func::DiracDelta if !ascii() => self.print(&args[0]).paren().left(&[Pict::new("δ")]),
            // Max/Min are not Functions in SymPy: printed with str()
            Func::Max | Func::Min => {
                let mut v = args.to_vec();
                v.sort_by_cached_key(crate::sortkey::sort_key);
                Pict::new(&crate::printer::to_str(&Expr::from_kind(Kind::Func(f.clone(), v))))
            }
            _ => self.function(f.name(), args),
        }
    }

    fn print_derivative(&self, f: &Expr, vars: &[(Expr, u32)]) -> Pict {
        let partial = crate::traverse::free_symbols(f).len() > 1;
        let d = if partial && !ascii() { "∂" } else { "d" };
        let mut x: Option<Pict> = None;
        let mut total = 0u32;
        for (v, n) in vars.iter().rev() {
            let s = self.print(v);
            let mut ds = s.left(&[Pict::new(d)]);
            total += n;
            if *n > 1 {
                ds = Pict::pow(ds, Pict::new(&n.to_string()));
            }
            x = Some(match x {
                None => ds,
                Some(p) => p.right(&[Pict::new(" "), ds]),
            });
        }
        let fp = self.print(f).paren().bind(Binding::Func);
        let mut p = Pict::new(d);
        if total > 1 {
            p = Pict::pow(p, Pict::new(&total.to_string()));
        }
        let mut p = p.below(vec![Item::Line, Item::P(x.unwrap())]);
        p.baseline += 1;
        Pict::next(&[p, fp]).bind(Binding::Mul)
    }

    fn print_piecewise(&self, pieces: &[(Expr, Expr)]) -> Pict {
        let n = pieces.len();
        let mut cells: Vec<[Pict; 2]> = Vec::with_capacity(n);
        for (e, c) in pieces {
            let pe = self.print(e);
            let pc = if matches!(c.kind(), Kind::Bool(true)) {
                Pict::new("otherwise")
            } else {
                Pict::new("for ").right(&[self.print(c)])
            };
            cells.push([pe, pc]);
        }
        let maxw = [
            cells.iter().map(|c| c[0].width()).max().unwrap_or(0),
            cells.iter().map(|c| c[1].width()).max().unwrap_or(0),
        ];
        let mut d: Option<Pict> = None;
        for row in cells {
            let mut drow: Option<Pict> = None;
            for (j, p) in row.into_iter().enumerate() {
                let wd = maxw[j] - p.width();
                let wl = wd / 2;
                let wr = wd - wl;
                let p = p.right(&[Pict::new(&" ".repeat(wr))]).left(&[Pict::new(&" ".repeat(wl))]);
                drow = Some(match drow {
                    None => p,
                    Some(r) => r.right(&[Pict::new("  "), p]),
                });
            }
            let drow = drow.unwrap();
            d = Some(match d {
                None => drow,
                Some(dd) => dd.below(vec![Item::P(Pict::new(" "))]).below(vec![Item::P(drow)]),
            });
        }
        let mut d = d.unwrap().parens("{", "");
        d.baseline = d.height() as isize / 2;
        d.binding = Binding::Open;
        d
    }
}

enum Den {
    E(Expr),
    Pow(Expr, Expr),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> String {
        pretty(&parse_str(s).unwrap())
    }

    #[test]
    fn basics() {
        assert_eq!(p("x**2"), " 2\nx ");
        assert_eq!(p("x**2 + 1/x"), " 2   1\nx  + ─\n     x");
        assert_eq!(p("(x+1)/(x-1)"), "x + 1\n─────\nx - 1");
        assert_eq!(p("sqrt(x)"), "√x");
        assert_eq!(p("sqrt(x+1)"), "  _______\n╲╱ x + 1 ");
        assert_eq!(p("x**(1/3)"), "3 ___\n╲╱ x ");
        assert_eq!(p("sin(x)**2"), "   2   \nsin (x)");
        assert_eq!(p("sin(x/y)"), "   ⎛x⎞\nsin⎜─⎟\n   ⎝y⎠");
        assert_eq!(p("-x/2"), "-x \n───\n 2 ");
        assert_eq!(p("x1"), "x₁");
        assert_eq!(p("Abs(x)"), "│x│");
        assert_eq!(p("Piecewise((x, x > 0), (0, True))"), "⎧x  for x > 0\n⎨            \n⎩0  otherwise");
        let a = |s: &str| {
            let st = PrettySettings { use_unicode: false, ..Default::default() };
            pretty_with(&parse_str(s).unwrap(), &st)
        };
        assert_eq!(a("sqrt(x)"), "  ___\n\\/ x ");
        assert_eq!(a("x/y + pi"), "x     \n- + pi\ny     ");
    }
}
