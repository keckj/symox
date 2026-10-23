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

    /// `stringPict(s, baseline)`: Python `splitlines` drops one trailing
    /// empty line.
    fn from_splitlines(s: &str, baseline: isize) -> Pict {
        if s.is_empty() {
            return Pict { lines: Vec::new(), baseline, binding: Binding::Atom, func: None };
        }
        Pict::with_baseline(s.strip_suffix('\n').unwrap_or(s), baseline)
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

/// The elements the reference `_print_Range` shows, `None` for the dots;
/// `None` when the range is not numeric. `latex` selects the latex rule
/// (`size < 4` instead of `len <= 4`).
pub(crate) fn range_display(args: &[Expr], latex: bool) -> Option<Vec<Option<Expr>>> {
    if args.len() != 3 {
        return None;
    }
    let (start, stop, step) = (&args[0], &args[1], &args[2]);
    let inf = |e: &Expr| *e == Expr::oo() || *e == Expr::neg_oo();
    let st = step.as_i64()?;
    if inf(start) && inf(stop) {
        let v = if st > 0 { [-1, 0, 1] } else { [1, 0, -1] };
        return Some(vec![None, Some(Expr::int(v[0])), Some(Expr::int(v[1])), Some(Expr::int(v[2])), None]);
    }
    if inf(start) {
        let last = stop.as_i64()? - st;
        return Some(vec![None, Some(Expr::int(last - st)), Some(Expr::int(last))]);
    }
    let a = start.as_i64()?;
    if inf(stop) {
        return Some(vec![Some(Expr::int(a)), Some(Expr::int(a + st)), None]);
    }
    let b = stop.as_i64()?;
    if st == 0 {
        return None;
    }
    let (d, m) = if st > 0 { (b - a, st) } else { (a - b, -st) };
    let size = if d <= 0 { 0 } else { (d + m - 1) / m };
    let elem = |i: i64| Some(Expr::int(a + i * st));
    let few = if latex { size < 4 } else { size <= 4 };
    if few {
        return Some((0..size).map(elem).collect());
    }
    Some(vec![elem(0), elem(1), None, elem(size - 1)])
}

/// The reference `ImageSet.base_sets`: symox keeps several base sets as
/// one `ProductSet`, unpacked here for a signature of the same length.
pub(crate) fn image_base_sets(sig: &[Expr], sets: &[Expr]) -> Vec<Expr> {
    if sig.len() > 1 && sets.len() == 1 {
        if let Kind::Func(Func::Set(crate::sets::SetOp::ProductSet), parts) = sets[0].kind() {
            if parts.len() == sig.len() {
                return parts.clone();
            }
        }
    }
    sets.to_vec()
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
            "int" => (" | ", "  /", "/  ", None, " \n|\n "),
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
        "int" => ("⎮", "⌠", "⌡", None, "∫"),
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

    /// The reference `_print_Order`.
    fn print_order(&self, args: &[Expr]) -> Pict {
        let mut p = self.print(&args[0]);
        let vars: Vec<Expr> = args[1..].iter().step_by(2).cloned().collect();
        let points: Vec<Expr> = args[1..].iter().skip(1).step_by(2).cloned().collect();
        if points.iter().any(|q| !q.is_zero()) || vars.len() > 1 {
            let mut parts = vec![p, Pict::new("; ")];
            parts.push(if vars.len() > 1 { self.seq(&vars).paren() } else { self.print(&vars[0]) });
            parts.push(Pict::new(xsym(" → ", " -> ")));
            parts.push(if points.len() > 1 { self.seq(&points).paren() } else { self.print(&points[0]) });
            p = Pict::next(&parts);
        }
        p.paren().left(&[Pict::new("O")])
    }

    /// `_hprint_vec`: the items joined by `, `, or a blank.
    fn hprint_vec(&self, items: &[Expr]) -> Pict {
        if items.is_empty() {
            return Pict::new(" ");
        }
        self.seq(items)
    }

    /// `_hprint_vseparator`: `p1 | p2` with a bar of their joint height.
    fn hprint_vseparator(p1: Pict, p2: Pict) -> Pict {
        let tmp = Pict::next(&[p1.clone(), p2.clone()]);
        let sep = Pict::with_baseline(&vobj("|", tmp.height()), tmp.baseline);
        Pict::next(&[p1, sep, p2])
    }

    /// The `annotated` letter drawing `(sz, t, b, add, img)`.
    fn annotated(letter: char) -> (isize, isize, isize, isize, &'static str) {
        match (letter, ascii()) {
            ('F', false) => (2, 0, 2, 0, "┌─\n├─\n╵"),
            ('G', false) => (3, 0, 3, 1, "╭─╮\n│╶┐\n╰─╯"),
            ('F', true) => (3, 0, 3, 0, " _\n|_\n|\n"),
            _ => (3, 0, 3, 1, " __\n/__\n\\_|"),
        }
    }

    /// The letter of `_print_hyper`/`_print_meijerg` sized for `d`, and
    /// `above` and `add`.
    fn hyper_letter(letter: char, d: &Pict) -> (Pict, isize, isize) {
        let above = d.height() as isize / 2 - 1;
        let below = d.height() as isize - above - 1;
        let (sz, t, b, add, img) = Self::annotated(letter);
        let text = format!("{}{}{}", "\n".repeat((above - t).max(0) as usize), img, "\n".repeat((below - b).max(0) as usize));
        (Pict::from_splitlines(&text, above + sz), above, add)
    }

    /// The reference `_print_hyper`.
    fn print_hyper(&self, ap: &[Expr], bq: &[Expr], z: &Expr) -> Pict {
        let mut p = self.print(z);
        p.baseline = p.height() as isize / 2;
        let mut d = self.hprint_vec(ap).below(vec![Item::P(Pict::new(" "))]);
        d = d.below(vec![Item::P(self.hprint_vec(bq))]);
        d.baseline = d.height() as isize / 2;
        let p = p.left(&[Pict::new(" ")]);
        let d = d.right(&[Pict::new(" ")]);
        let d = Self::hprint_vseparator(d, p).paren();
        let (f, above, _) = Self::hyper_letter('F', &d);
        // `_print_hyper` replaces the annotated offset by `(sz + 1)//2`.
        let add = (Self::annotated('F').0 + 1) / 2;
        let mut f = f.left(&[Pict::new(&ap.len().to_string())]);
        f = f.right(&[Pict::new(&bq.len().to_string())]);
        f.baseline = above + add;
        f.right(&[Pict::new(" "), d])
    }

    /// The reference `_print_meijerg`.
    fn print_meijerg(&self, an: &[Expr], aother: &[Expr], bm: &[Expr], bother: &[Expr], z: &Expr) -> Pict {
        let mut p = self.print(z);
        p.baseline = p.height() as isize / 2;
        let mut vp = [[self.hprint_vec(an), self.hprint_vec(aother)], [self.hprint_vec(bm), self.hprint_vec(bother)]];
        for i in 0..2 {
            let maxw = vp[0][i].width().max(vp[1][i].width());
            for row in vp.iter_mut() {
                let s = &row[i];
                let left = (maxw - s.width()) / 2;
                let right = maxw - left - s.width();
                row[i] = s.left(&[Pict::new(&" ".repeat(left))]).right(&[Pict::new(&" ".repeat(right))]);
            }
        }
        let d1 = vp[0][0].right(&[Pict::new("  "), vp[0][1].clone()]).below(vec![Item::P(Pict::new(" "))]);
        let d2 = vp[1][0].right(&[Pict::new("  "), vp[1][1].clone()]);
        let mut d = d1.below(vec![Item::P(d2)]);
        d.baseline = d.height() as isize / 2;
        let p = p.left(&[Pict::new(" ")]);
        let d = d.right(&[Pict::new(" ")]);
        let d = Self::hprint_vseparator(d, p).paren();
        let (mut f, above, add) = Self::hyper_letter('G', &d);
        let num = |n: usize| Pict::new(&n.to_string());
        let adjust = |p1: Pict, p2: Pict| -> (Pict, Pict) {
            let (w1, w2) = (p1.width(), p2.width());
            if w1 > w2 {
                let p2 = p2.left(&[Pict::new(&" ".repeat(w1 - w2))]);
                (p1, p2)
            } else if w2 > w1 {
                (p1.left(&[Pict::new(&" ".repeat(w2 - w1))]), p2)
            } else {
                (p1, p2)
            }
        };
        let (pp, pm) = adjust(num(an.len() + aother.len()), num(bm.len()));
        let (pq, pn) = adjust(num(bm.len() + bother.len()), num(an.len()));
        let mut pu = pm.right(&[Pict::new(", "), pn]);
        let pl = pp.right(&[Pict::new(", "), pq]);
        let ht = f.baseline - above - 2;
        if ht > 0 {
            pu = pu.below(vec![Item::P(Pict::from_splitlines(&"\n".repeat(ht as usize), 0))]);
        }
        let pr = pu.below(vec![Item::P(pl)]);
        f.baseline = above;
        let mut f = f.right(&[pr]);
        f.baseline = above + add;
        f.right(&[Pict::new(" "), d])
    }

    /// `parens(left, right, ifascii_nougly=True)`: single line brackets in
    /// ASCII mode.
    fn parens_nougly(p: Pict, left: &str, right: &str) -> Pict {
        if !ascii() {
            return p.parens(left, right);
        }
        let mut parts = Vec::with_capacity(3);
        if !left.is_empty() {
            parts.push(Pict::new(left));
        }
        parts.push(p);
        if !right.is_empty() {
            parts.push(Pict::new(right));
        }
        let mut res = Pict::next(&parts);
        res.binding = Binding::Atom;
        res
    }

    /// The reference `_print_seq` over already printed forms.
    fn seq_forms(forms: Vec<Pict>, left: &str, right: &str, delimiter: &str) -> Pict {
        let mut parts = Vec::with_capacity(forms.len() * 2);
        for (i, f) in forms.into_iter().enumerate() {
            if i > 0 {
                parts.push(Pict::new(delimiter));
            }
            parts.push(f);
        }
        let p = if parts.is_empty() { Pict::new("") } else { Pict::next(&parts) };
        Self::parens_nougly(p, left, right)
    }

    /// `_hprint_vseparator(p1, p2, left, right, delimiter=' ',
    /// ifascii_nougly=True)`.
    fn vseparator_braced(p1: Pict, p2: Pict) -> Pict {
        if ascii() {
            return Self::seq_forms(vec![p1, Pict::new("|"), p2], "{", "}", " ");
        }
        let tmp = Self::seq_forms(vec![p1.clone(), p2.clone()], "{", "}", " ");
        let sep = Pict::with_baseline(&vobj("|", tmp.height()), tmp.baseline);
        Self::seq_forms(vec![p1, sep, p2], "{", "}", " ")
    }

    /// The reference set printers (`_print_FiniteSet`, `_print_Interval`,
    /// `_print_Union`, `_print_ImageSet`, `_print_ConditionSet`, ...).
    fn print_set(&self, op: crate::sets::SetOp, args: &[Expr]) -> Option<Pict> {
        use crate::sets::SetOp as O;
        let is_op = |e: &Expr, ops: &[O]| matches!(e.kind(), Kind::Func(Func::Set(o), _) if ops.contains(o));
        let joined = |delim: &str, paren: &[O]| {
            let forms = args.iter().map(|a| {
                let p = self.print(a);
                if is_op(a, paren) { p.paren() } else { p }
            }).collect();
            Self::seq_forms(forms, "", "", delim)
        };
        Some(match op {
            O::FiniteSet => {
                let mut v = args.to_vec();
                v.sort_by_cached_key(crate::sortkey::sort_key);
                Self::seq_forms(v.iter().map(|a| self.print(a)).collect(), "{", "}", ", ")
            }
            O::Interval => {
                let e = Expr::from_kind(Kind::Func(Func::Set(op), args.to_vec()));
                let (a, b, lo, ro) = crate::sets::interval_parts(&e)?;
                if a == b {
                    return Some(Self::seq_forms(vec![self.print(&a)], "{", "}", ", "));
                }
                Self::seq_forms(vec![self.print(&a), self.print(&b)], if lo { "(" } else { "[" }, if ro { ")" } else { "]" }, ", ")
            }
            O::Union => joined(xsym(" ∪ ", " U "), &[O::ProductSet, O::Intersection, O::Complement]),
            O::Intersection => joined(xsym(" ∩ ", " n "), &[O::ProductSet, O::Union, O::Complement]),
            O::Complement => joined(" \\ ", &[O::ProductSet, O::Intersection, O::Union]),
            O::ProductSet if args.len() > 1 && args.iter().any(|a| *a != args[0]) => {
                joined(xsym(" × ", " x "), &[O::Union, O::Intersection, O::ProductSet])
            }
            O::Integers => Pict::new(xsym("ℤ", "Integers")),
            O::Naturals => Pict::new(xsym("ℕ", "Naturals")),
            O::Naturals0 => Pict::new(xsym("ℕ₀", "Naturals0")),
            O::Complexes => Pict::new(xsym("ℂ", "Complexes")),
            O::Rationals => Pict::new(xsym("ℚ", "Rationals")),
            O::EmptySet => Pict::new(xsym("∅", "EmptySet")),
            O::UniversalSet => Pict::new(xsym("𝕌", "UniversalSet")),
            O::Reals if ascii() => Self::seq_forms(vec![Pict::new("-oo"), Pict::new("oo")], "(", ")", ", "),
            O::Reals => Pict::new("ℝ"),
            O::ImageSet => {
                let Kind::Func(Func::Set(O::Lambda), lam) = args.first()?.kind() else { return None };
                if lam.len() != 2 {
                    return None;
                }
                let sig = crate::special::tuple_items(&lam[0]).map(|v| v.to_vec()).unwrap_or_else(|| vec![lam[0].clone()]);
                let sets = image_base_sets(&sig, &args[1..]);
                let sets = &sets[..];
                if sig.len() != sets.len() {
                    return None;
                }
                let inn = xsym("∊", "in");
                let expr = self.print(&lam[1]);
                let right = if sig.len() == 1 {
                    Self::seq_forms(vec![self.print(&sig[0]), Pict::new(inn), self.print(&sets[0])], "", "", " ")
                } else {
                    let mut forms = Vec::new();
                    for (i, (v, st)) in sig.iter().zip(sets).enumerate() {
                        if i > 0 {
                            forms.push(Pict::new(", "));
                        }
                        forms.extend([self.print(v), Pict::new(" "), Pict::new(inn), Pict::new(" "), self.print(st)]);
                    }
                    Self::seq_forms(forms, "", "", "")
                };
                Self::vseparator_braced(expr, right)
            }
            O::ConditionSet if args.len() == 3 => {
                let syms = crate::special::tuple_items(&args[0]).map(|v| v.to_vec()).unwrap_or_else(|| vec![args[0].clone()]);
                let variables = Self::seq_forms(syms.iter().map(|a| self.print(a)).collect(), "", "", ", ");
                let mut cond = self.print(&args[1]);
                if !ascii() {
                    cond = cond.paren();
                }
                if is_op(&args[2], &[O::UniversalSet]) {
                    return Some(Self::vseparator_braced(variables, cond));
                }
                let c = Self::seq_forms(
                    vec![variables.clone(), Pict::new(xsym("∊", "in")), self.print(&args[2]), Pict::new(xsym("∧", "and")), cond],
                    "", "", " ");
                Self::vseparator_braced(variables, c)
            }
            O::Range => {
                let items = range_display(args, false)?;
                let dots = xsym("…", "...");
                let forms = items.iter().map(|i| match i { Some(e) => self.print(e), None => Pict::new(dots) }).collect();
                Self::seq_forms(forms, "{", "}", ", ")
            }
            _ => return None,
        })
    }

    /// The reference `_print_Integral`.
    fn print_integral(&self, args: &[Expr]) -> Pict {
        let limits = crate::binders::limits(args);
        let mut arg = self.print(&args[0]);
        if args[0].is_add() {
            arg = arg.paren();
        }
        for lim in &limits {
            let mut pa = self.print(&lim[0]);
            if pa.width() > 1 {
                pa = pa.paren();
            }
            arg = arg.right(&[Pict::new(" d"), pa]);
        }
        let asc = ascii();
        let mut s: Option<Pict> = None;
        for lim in &limits {
            let h = arg.height() as isize;
            let big = h + 2 + if asc { 2 } else { 0 };
            let mut pform = Pict::new(&vobj("int", big as usize));
            pform.baseline = arg.baseline + (big - h) / 2;
            if lim.len() > 1 {
                let (mut a, mut b) = if lim.len() == 2 {
                    (Pict::new(""), self.print(&lim[1]))
                } else {
                    (self.print(&lim[1]), self.print(&lim[2]))
                };
                if asc {
                    let spc = 3usize.saturating_sub(b.width()).max(1);
                    b = b.left(&[Pict::new(&" ".repeat(spc))]);
                    let spc = 4usize.saturating_sub(a.width()).max(1);
                    a = a.right(&[Pict::new(&" ".repeat(spc))]);
                }
                pform = pform.above(vec![Item::P(b)]);
                pform = pform.below(vec![Item::P(a)]);
            }
            if !asc {
                pform = pform.right(&[Pict::new(" ")]);
            }
            s = Some(match s {
                None => pform,
                Some(prev) => prev.left(&[pform]),
            });
        }
        let mut p = match s {
            Some(sign) => arg.left(&[sign]),
            None => arg,
        };
        p.binding = Binding::Mul;
        p.func = None;
        p
    }

    fn print_func(&self, f: &Func, args: &[Expr]) -> Pict {
        if let Func::Set(op) = f {
            if let Some(p) = self.print_set(*op, args) {
                return p;
            }
        }
        if *f == Func::Product { return self.print_product(args); }
        if *f == Func::Integral && args.len() > 1 { return self.print_integral(args); }
        if *f == Func::Order && !args.is_empty() { return self.print_order(args); }
        if let Func::Special(crate::special::SpecialFn::Hyper) = f {
            if args.len() == 3 {
                return self.print_hyper(&args[0].args(), &args[1].args(), &args[2]);
            }
        }
        if let Func::Special(crate::special::SpecialFn::MeijerG) = f {
            if args.len() == 3 {
                let (a, b) = (args[0].args(), args[1].args());
                if a.len() == 2 && b.len() == 2 {
                    return self.print_meijerg(&a[0].args(), &a[1].args(), &b[0].args(), &b[1].args(), &args[2]);
                }
            }
        }
        if let Func::Undefined(head) = f {
            if crate::poly::rootof::is_rootof_head(head, args.len()) {
                // The reference's _print_ComplexRootOf: an atom.
                return Pict::next(&[Pict::new("CRootOf"), self.seq(args).paren()]);
            }
        }
        if matches!(f,Func::Predicate(_) | Func::AppliedPredicate | Func::Subs | Func::Limit | Func::Integral | Func::Sum){return Pict::new(&crate::printer::to_str(&Expr::from_kind(Kind::Func(f.clone(),args.to_vec()))));}
        if let Some(p) = self.print_matrix_func(f, args) { return p; }
        // `_print_Lambda`: `signature ↦ body`.
        if let (Func::Set(crate::sets::SetOp::Lambda), [signature, body]) = (f, args) {
            let signature = match crate::special::tuple_items(signature) {
                Some([single]) if crate::special::tuple_items(single).is_none() => self.print(single),
                Some(items) => self.pretty_tuple(items),
                None => self.print(signature),
            };
            return Pict::next(&[signature, Pict::new(xsym(" ↦ ", " -> ")), self.print(body)]).bind(Binding::Open);
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
            // SymPy prints a non-default H0 as an ordinary function argument.
            Func::Heaviside if args.len() > 1 && args[1] != Expr::rational(1, 2) => {
                self.seq(args).paren().left(&[Pict::new(xsym("θ", "Heaviside"))])
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

// ----------------------------------------------------------------------
// matrix expressions (SymPy's PrettyPrinter matrix handlers)
// ----------------------------------------------------------------------

use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};

impl PrettyPrinter<'_> {
    /// `_print_matrix_contents`: the grid of cells, two columns between
    /// cells and an empty line between rows, baseline at mid-height.
    fn matrix_grid(&self, rows: usize, cols: usize, cells: Vec<Pict>) -> Pict {
        if rows == 0 || cols == 0 {
            return Pict::new("");
        }
        let widths: Vec<usize> = (0..cols).map(|j| (0..rows).map(|i| cells[i * cols + j].width()).max().unwrap_or(0)).collect();
        let mut lines = Vec::new();
        let mut cells = cells.into_iter();
        for _ in 0..rows {
            let mut parts = Vec::with_capacity(2 * cols);
            for (j, width) in widths.iter().enumerate() {
                let mut cell = cells.next().unwrap();
                if j != 0 { parts.push(Pict::new("  ")); }
                cell.lines = cell.lines.iter().map(|line| center(line, *width)).collect();
                parts.push(cell);
            }
            if !lines.is_empty() { lines.push(String::new()); }
            lines.extend(Pict::next(&parts).lines);
        }
        let width = lines.iter().map(|line| line_width(line)).max().unwrap_or(0);
        let baseline = (lines.len() / 2) as isize;
        Pict { lines: lines.iter().map(|line| center(line, width)).collect(), baseline, binding: Binding::Atom, func: None }
    }

    fn block_cells(&self, m: &MatrixExpr) -> Option<(usize, usize, Vec<Pict>)> {
        match m.kind() {
            MatrixKind::Block(grid) => {
                let (r, c) = grid.shape();
                Some((r, c, (0..r * c).map(|k| self.print_matrix(grid.get(k / c, k % c))).collect()))
            }
            MatrixKind::BlockDiag(blocks) => {
                let k = blocks.len();
                Some((k, k, (0..k * k).map(|p| if p / k == p % k { self.print_matrix(&blocks[p / k]) } else { Pict::new(xsym("𝟘", "0")) }).collect()))
            }
            _ => None,
        }
    }

    fn operand(&self, o: &MatrixOperand) -> (Pict, i32) {
        match o {
            MatrixOperand::Scalar(e) => (self.print(e), super::precedence(e)),
            MatrixOperand::Matrix(m) => (self.print_matrix(m), super::latex::matrix_precedence_with(m, false)),
        }
    }

    /// `_print_MatPow`, `_print_Transpose` and `_print_Adjoint`.
    fn matrix_super(&self, base: &MatrixExpr, sup: Pict, block_plain: bool) -> Pict {
        let mut p = self.print_matrix(base);
        let plain = matches!(base.kind(), MatrixKind::Symbol(_)) || (block_plain && matches!(base.kind(), MatrixKind::Block(_) | MatrixKind::BlockDiag(_)));
        if !plain { p = p.paren(); }
        Pict::pow(p, sup)
    }

    fn adjoint_symbol() -> Pict { Pict::new(xsym("†", "+")) }

    /// `_print_seq` with a delimiter, parenthesizing the selected items.
    fn matrix_seq(&self, items: &[MatrixExpr], delimiter: &str, paren: impl Fn(&MatrixExpr) -> bool) -> Pict {
        let mut forms = Vec::with_capacity(2 * items.len());
        for (i, item) in items.iter().enumerate() {
            if i > 0 { forms.push(Pict::new(delimiter)); }
            let p = self.print_matrix(item);
            forms.push(if paren(item) { p.paren() } else { p });
        }
        Pict::next(&forms)
    }

    /// SymPy's pretty form of a matrix expression.
    fn print_matrix(&self, m: &MatrixExpr) -> Pict {
        match m.kind() {
            MatrixKind::Symbol(name) => Pict::new(&pretty_symbol(name)),
            MatrixKind::Identity => Pict::new(xsym("𝕀", "I")),
            MatrixKind::Zero => Pict::new(xsym("𝟘", "0")),
            MatrixKind::One => Pict::new(xsym("𝟙", "1")),
            MatrixKind::Explicit(x) => {
                if x.rows() == 0 || x.cols() == 0 { return Pict::new("[]"); }
                self.matrix_grid(x.rows(), x.cols(), x.data().iter().map(|e| self.print(e)).collect()).parens("[", "]")
            }
            MatrixKind::Block(_) | MatrixKind::BlockDiag(_) => {
                let (r, c, mut cells) = self.block_cells(m).unwrap();
                if r == 1 && c == 1 { return cells.pop().unwrap(); }
                self.matrix_grid(r, c, cells).parens("[", "]")
            }
            // `_print_MatAdd`: the terms in argument order; a term whose
            // coefficient is negative carries its own sign.
            MatrixKind::Add(terms) => {
                let mut forms: Vec<Pict> = Vec::with_capacity(2 * terms.len());
                for (i, term) in terms.iter().enumerate() {
                    if i > 0 {
                        let negative = matches!(term.kind(), MatrixKind::Mul(c, _) if super::minus_sign(c));
                        forms.push(Pict::new(if negative { " " } else { " + " }));
                    }
                    forms.push(self.print_matrix(term));
                }
                Pict::next(&forms)
            }
            // `_print_MatMul`: compound sums and products get parentheses.
            MatrixKind::Mul(c, factors) => {
                let count = factors.len() + usize::from(!c.is_one());
                let mut forms = Vec::with_capacity(count);
                if !c.is_one() {
                    let p = self.print(c);
                    forms.push(if count > 1 && matches!(c.kind(), Kind::Add(..)) { p.paren() } else { p });
                }
                for f in factors {
                    let p = self.print_matrix(f);
                    let compound = matches!(f.kind(), MatrixKind::Add(_) | MatrixKind::HadamardProduct(_) | MatrixKind::Kronecker(_));
                    forms.push(if compound && count > 1 { p.paren() } else { p });
                }
                Pict::mul(forms)
            }
            MatrixKind::Pow(base, exp) => self.matrix_super(base, self.print(exp), false),
            MatrixKind::Inverse(base) => self.matrix_super(base, self.print(&Expr::minus_one()), false),
            MatrixKind::Transpose(base) => self.matrix_super(base, Pict::new("T"), true),
            MatrixKind::Adjoint(base) => self.matrix_super(base, Self::adjoint_symbol(), true),
            // SymPy holds the conjugate of a matrix expression as Adjoint(Transpose(base)).
            MatrixKind::Conjugate(base) => Pict::pow(self.matrix_super(base, Pict::new("T"), true).paren(), Self::adjoint_symbol()),
            MatrixKind::HadamardProduct(factors) => self.matrix_seq(factors, xsym("∘", ".*"), |f| matches!(f.kind(), MatrixKind::Add(_) | MatrixKind::Mul(..) | MatrixKind::HadamardProduct(_))),
            MatrixKind::Kronecker(factors) => self.matrix_seq(factors, xsym(" ⨂ ", " x "), |f| matches!(f.kind(), MatrixKind::Add(_) | MatrixKind::Mul(..))),
            MatrixKind::HadamardPower(base, exp) => {
                let (base, _) = self.operand(base);
                let (mut exp, prec) = self.operand(exp);
                if prec < super::PREC_MUL { exp = exp.paren(); }
                Pict::pow(base, Pict::next(&[Pict::new(xsym("∘", ".")), exp]))
            }
            MatrixKind::Slice(parent, rows, cols) => {
                let mut pf = self.print_matrix(parent);
                if !matches!(parent.kind(), MatrixKind::Symbol(_)) { pf = pf.paren(); }
                let part = |spec, dim| {
                    let mut forms = Vec::new();
                    for (i, x) in super::latex::slice_parts(spec, dim).iter().enumerate() {
                        if i > 0 { forms.push(Pict::new(":")); }
                        forms.push(x.as_ref().map(|e| self.print(e)).unwrap_or_else(|| Pict::new("")));
                    }
                    Pict::next(&forms)
                };
                let pa = Pict::next(&[part(rows, parent.rows()), Pict::new(", "), part(cols, parent.cols())]).parens("[", "]");
                let mut p = Pict::next(&[pf.clone(), pa.clone()]).bind(Binding::Func);
                p.func = Some(Rc::new((pf, pa)));
                p
            }
            // `_helper_print_function(..., elementwise=True)`.
            MatrixKind::ApplyFunction(function, parent) => {
                // A Lambda has no `__name__`, so SymPy parenthesizes it.
                let pf = self.print(function).paren();
                let pf = Pict::next(&[pf, Pict::new(xsym("˳", "."))]);
                let pa = self.print_matrix(parent).paren();
                let mut p = Pict::next(&[pf.clone(), pa.clone()]).bind(Binding::Func);
                p.func = Some(Rc::new((pf, pa)));
                p
            }
            _ => Pict::new(&m.to_str()),
        }
    }

    fn print_matrix_func(&self, f: &Func, args: &[Expr]) -> Option<Pict> {
        Some(match f {
            Func::MatrixProxy(parent) => self.print_matrix(parent),
            Func::MatrixDotProduct(a, b) => Pict::mul(vec![self.print_matrix(a), self.print_matrix(b)]),
            Func::MatrixTrace(parent) => {
                let mut d = self.print_matrix(parent).paren();
                d.baseline = (d.height() / 2) as isize;
                d.left(&[Pict::new("tr")])
            }
            Func::MatrixDeterminant(parent) => {
                if let Some((r, c, cells)) = self.block_cells(parent) {
                    return Some(self.matrix_grid(r, c, cells).parens("|", "|"));
                }
                let mut d = self.print_matrix(parent);
                d.baseline = (d.height() / 2) as isize;
                d.parens("|", "|")
            }
            Func::MatrixElement(parent) => {
                if let (MatrixKind::Symbol(name), Some(i), Some(j)) = (parent.kind(), args[0].as_number(), args[1].as_number()) {
                    if i.is_integer() && j.is_integer() {
                        return Some(Pict::new(&pretty_symbol(&format!("{name}_{i}{j}"))));
                    }
                }
                let pf = self.print_matrix(parent).paren();
                let pa = self.seq(args).parens("[", "]");
                let mut p = Pict::next(&[pf.clone(), pa.clone()]).bind(Binding::Func);
                p.func = Some(Rc::new((pf, pa)));
                p
            }
            _ => return None,
        })
    }
}


impl PrettyPrinter<'_> {
    /// `_print_tuple`: `(a, b)`, or `(a,)` for a single item.
    fn pretty_tuple(&self, items: &[Expr]) -> Pict {
        let print = |e: &Expr| match crate::special::tuple_items(e) { Some(inner) => self.pretty_tuple(inner), None => self.print(e) };
        let mut forms = Vec::with_capacity(2 * items.len());
        for (i, item) in items.iter().enumerate() {
            if i > 0 { forms.push(Pict::new(", ")); }
            forms.push(print(item));
        }
        if items.len() == 1 { forms.push(Pict::new(",")); }
        Pict::next(&forms).paren()
    }
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
