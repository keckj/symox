"""Refreshed symox vs SymPy v4, rendered only from verified frozen data.

Run refresh.py to measure and validate. See README.md for scope and provenance.
"""
import json
import math
from pathlib import Path

import numpy as np
import manimgx as m

HERE = Path(__file__).parent
ASSETS = HERE / "data"
DATA = json.loads((ASSETS / "data.json").read_text())
CASES = json.loads((ASSETS / "cases3.json").read_text())
TP = {lib: json.loads((ASSETS / f"tp_{lib}.json").read_text()) for lib in ("symox", "sympy")}
C4 = json.loads((ASSETS / "cases4.json").read_text())
SC4 = json.loads((ASSETS / "scene4.json").read_text())
EXACT = json.loads((ASSETS / "exactness.json").read_text())
DOWNSTREAM = json.loads((ASSETS / "downstream.json").read_text())
SAGE_TESTS = json.loads((ASSETS / "sage.json").read_text())
LCAPY_TESTS = json.loads((ASSETS / "lcapy-tests.json").read_text())

BG = "#1a1a19"
INK = "#f0efec"
INK2 = "#c3c2b7"
GRID = "#3a3a36"
SYMOX = "#3987e5"
SYMPY = "#d95926"
SANS = "IBM Plex Sans"
MONO = "IBM Plex Mono"

m.config.background_color = BG

INTEGRANDS = [
    r"x^5 e^x", r"\sin^6 x", r"x^3 \ln^2 x", r"x^2 \cos x \sin x", r"e^{2x} \sin 3x",
    r"\frac{x^3+2x}{(x^2+1)^2}", r"x \arctan x", r"\frac{x}{x^3-1}", r"\sin^3 x \cos^4 x",
    r"x^4 e^{-x} \cos x", r"\frac{1}{x^6-1}", r"\ln(x^2+1)", r"\sinh^4 x", r"\tan^3 x",
    r"x^7 e^{-x^2}", r"\cos^8 x", r"x \sin^2 x", r"e^x \cos^2 x", r"\frac{1}{x^3+1}",
    r"x^2 \arctan x", r"\ln^3 x", r"x^3 \sinh x", r"\arcsin x", r"\frac{1}{x^2+2x+5}",
    r"\frac{2x+3}{x^2+x+1}", r"x^2 e^{3x} \sin x", r"\sqrt{1-x^2}", r"\frac{1}{\sqrt{x^2+1}}",
    r"\sin 3x \cos 5x", r"\sin^2 x \cos^2 x", r"\tan^4 x", r"x e^x \sin x", r"e^{-x} \sinh 2x",
    r"\cosh^3 x", r"x^2 \ln(x+1)", r"\frac{\ln x}{x^2}", r"\frac{x^4}{x^2+1}", r"\frac{1}{x^4-1}",
    r"\frac{x+1}{x(x^2+1)}", r"\frac{1}{x(x+1)^3}", r"\frac{x^5}{x^3-1}", r"e^{\sqrt{x}}",
    r"x^3 \cos x^2", r"x e^{-x} \cos 2x", r"\frac{1}{x^2(x^2+1)}", r"x^9 e^x", r"x^6 \sin x",
    r"\cos^5 x \sin^5 x", r"x^5 \ln^2 x", r"x \ln(x^2+1)", r"x^2 \sinh x \cosh x",
    r"\arctan \sqrt{x}", r"\frac{1}{(x-1)^2(x+2)}", r"x^3 e^{-2x} \sin x", r"\sqrt{x^2+1}",
    r"\frac{1}{x^5+x}", r"\frac{\cos x}{1+\sin^2 x}", r"e^{3x} \cos^2 2x", r"x^4 \ln^3 x",
    r"\tan x \sec^4 x", r"x \sinh^2 x", r"\frac{x^8}{x^2+1}", r"\frac{3x^2+1}{x^3+x+1}",
]



STATS_ROWS = [
    ("normal_moment", r"E[X^4],\ X \sim \mathcal{N}(\mu, \sigma^2)"),
    ("normal_var", r"\operatorname{Var}(X^2),\ X \sim \mathcal{N}(\mu, \sigma^2)"),
    ("binom_var", r"\operatorname{Var}(B),\ B \sim \mathcal{B}(12, p)"),
    ("binom_skew", r"\operatorname{Skew}(B),\ B \sim \mathcal{B}(10, p)"),
    ("dice", r"P(D_1 + D_2 + D_3 > 12),\ \text{three dice}"),
    ("exp_P", r"P(X > 2),\ X \sim \operatorname{Exp}(\lambda)"),
    ("gamma_mgf", r"M_G(t),\ G \sim \Gamma(k, \theta)"),
    ("sum_normals", r"f_{X+Y}(z),\ X, Y \sim \mathcal{N}(0, 1)"),
    ("beta_E", r"E[B],\ B \sim \operatorname{Beta}(a, b)"),
    ("normal_entropy", r"H(X),\ X \sim \mathcal{N}(\mu, \sigma^2)"),
]

CRYPTO_ROWS = [
    ("berlekamp", "Berlekamp-Massey", "connection polynomial of a 200-bit LFSR stream"),
    ("hill", "Hill cipher", "2,100 letters, 3×3 key mod 26"),
    ("lfsr", "LFSR", "3,000 bits over GF(2)"),
    ("vigenere", "Vigenère", "24,000 letters"),
    ("rsa", "RSA", "199 round trips, 122-digit modulus"),
]

PHYSICS_ROWS = [
    ("tf", "control", "closed loop: Feedback(Series(G₁, G₂)), cancel + expand"),
    ("pauli", "Pauli algebra", "(σ·r)⁶ expanded and reduced"),
    ("lagrange2", "mechanics", "double pendulum: LagrangesMethod"),
    ("wigner9j", "Wigner 9j", "wigner_9j(4, 4, 4, 4, 4, 4, 4, 4, 4)"),
    ("cg", "Clebsch-Gordan", "⟨15 3; 12 -5 | 10 -2⟩"),
]


def fmt_s(seconds):
    if seconds < 1e-3:
        return f"{seconds*1e6:.0f} µs"
    if seconds < 1:
        return f"{seconds*1e3:.1f} ms"
    return f"{seconds:.1f} s"


def fmt_x(ratio):
    return "×" + (f"{ratio:.2f}" if ratio < 10 else f"{ratio:,.0f}")


def T(text, **kwargs):
    kwargs.setdefault("font", SANS)
    return m.Text(text, **kwargs)


def live_text(content, place, **style):
    """Text that follows content() every frame; each distinct string is typeset once."""
    cache = {}

    def build():
        text = content()
        if text not in cache:
            cache[text] = T(text, **style)
        return place(cache[text].copy())
    return m.always_redraw(build)


def timing(lib, case):
    """Seconds for a bench case, or None when it timed out or failed."""
    r = CASES[lib][case]
    if "items" in r:
        return sum(i["time"] for i in r["items"])
    return r.get("time")


def timing4(lib, case):
    """Seconds for a module case, or None when it timed out or failed."""
    return C4[case][lib].get("time")


def ratio4(case):
    return timing4("sympy", case)/timing4("symox", case)


class Showcase(m.Scene):
    def construct(self):
        self.intro()
        self.integrals()
        self.throughput()
        self.exactness()
        self.determinant()
        self.series()
        self.system()
        self.pendulum()
        self.stats()
        self.physics()
        self.crypto()
        self.downstream()
        self.integration_tests()
        self.finale()

    def clear_scene(self, run_time=0.6):
        for mob in self.mobjects:
            mob.clear_updaters()
        self.play(m.FadeOut(m.Group(*self.mobjects)), run_time=run_time)

    def heading(self, text, sub=None):
        title = T(text, font_size=44, color=INK)
        title.to_corner(m.UL, buff=0.5)
        group = [title]
        if sub:
            caption = T(sub, font_size=26, color=INK2)
            if caption.width > 12.7:
                caption.scale_to_fit_width(12.7)
            caption.next_to(title, m.DOWN, aligned_edge=m.LEFT, buff=0.22)
            group.append(caption)
        self.play(*[m.FadeIn(g) for g in group], run_time=0.5)
        return group

    # -- 1. one import apart ------------------------------------------------
    def intro(self):
        logo = T("symox", font_size=120, weight=m.BOLD, color=SYMOX)
        tag = T("SymPy, rewritten in Rust", font_size=40, color=INK2)
        tag.next_to(logo, m.DOWN, buff=0.4)
        self.play(m.Write(logo), m.FadeIn(tag, shift=0.3*m.UP), run_time=1.0)
        self.play(m.VGroup(logo, tag).animate.scale(0.5).to_edge(m.UP, buff=0.6), run_time=0.6)

        before = T("import sympy as sp", font=MONO, font_size=56, t2c={"sympy": SYMPY})
        after = T("import symox as sp", font=MONO, font_size=56, t2c={"symox": SYMOX})
        caption = T("the same code, for the examples shown", font_size=30, color=INK2)
        caption.next_to(before, m.DOWN, buff=0.7)
        self.play(m.FadeIn(before), run_time=0.4)
        self.play(m.Transform(before, after), m.FadeIn(caption, shift=0.2*m.UP), run_time=0.8)
        self.wait(0.9)
        self.clear_scene(0.5)

    def exactness(self):
        self.heading("Exact arithmetic. Requested precision.", "Three independently checked examples")
        rows = m.VGroup()
        rational = m.VGroup(
            T("Exact rationals", font_size=27, color=SYMOX),
            m.MathTex(r"\frac{1}{3}+\frac{1}{6}=\frac{1}{2}", font_size=40, color=INK),
        ).arrange(m.RIGHT, buff=0.8)
        rational.move_to(np.array([0, 1.05, 0]))
        code = T("log(1 + Rational(1, 10**40)).evalf(50)", font=MONO, font_size=25, color=INK2)
        number = T(EXACT['precision']['result'], font=MONO, font_size=26, color=INK)
        if number.width > 12.3:
            number.scale_to_fit_width(12.3)
        precision = m.VGroup(code, number,
            T("50 digits · checked against SymPy and 110-digit mpmath", font_size=21, color=SYMOX)
        ).arrange(m.DOWN, buff=0.2).move_to(np.array([0, -0.55, 0]))
        convergence = m.VGroup(
            T("Assumptions matter", font_size=27, color=SYMOX),
            m.MathTex(r"a>0\quad\Longrightarrow\quad\int_0^\infty e^{-ax}\,dx=\frac{1}{a}", font_size=34, color=INK),
        ).arrange(m.RIGHT, buff=0.65).move_to(np.array([0, -2.25, 0]))
        rows.add(rational, precision, convergence)
        self.play(m.LaggedStart(*[m.FadeIn(row, shift=0.15*m.UP) for row in rows], lag_ratio=0.4), run_time=1.5)
        self.wait(3.4)
        self.clear_scene()

    def downstream(self):
        self.heading("Inside Lcapy", "The same RC circuit, run with both symbolic backends")
        # Hand-drawn circuit symbols keep this scene independent of schematics
        # executables. The response itself comes from actual Lcapy workers.
        origin = np.array([-3.7, 0.6, 0])
        def p(a,b):
            return origin + np.array([a,b,0])
        wire = m.VGroup(
            m.Line(p(-2,-1.6),p(-2,-0.8)), m.Line(p(-2,0.2),p(-2,1)),
            m.Circle(radius=0.5).move_to(p(-2,-0.3)), m.Line(p(-2,1),p(-0.9,1)),
            m.Line(p(-2.1,-0.08),p(-1.9,-0.08)),
            m.Line(p(-2,-0.18),p(-2,0.02)),
            m.Line(p(-2.1,-0.52),p(-1.9,-0.52)),
            m.Line(p(0.9,1),p(2,1)), m.Line(p(2,1),p(2,0)),
            m.Line(p(2,-0.3),p(2,-1.6)), m.Line(p(2,-1.6),p(-2,-1.6)),
            m.Rectangle(width=1.8,height=0.5).move_to(p(0,1)),
            m.Line(p(1.5,0),p(2.5,0)), m.Line(p(1.5,-0.3),p(2.5,-0.3)),
        ).set_color(INK2)
        labels = m.VGroup(
            T("R = 1 kΩ", font_size=23, color=INK).move_to(p(0,1.65)),
            T("C = 1 µF", font_size=23, color=INK).move_to(p(3.05,-0.6)),
            T("unit step", font_size=23, color=INK).move_to(p(-2,-2.1)),
            T("Vout", font_size=24, color=SYMOX).move_to(p(2.5,1.45)),
        )
        axes = m.Axes(x_range=[0,6,1], y_range=[0,1.2,0.2], x_length=5.1, y_length=2.8,
            tips=False, axis_config={"color": GRID,"stroke_width":2}).move_to(np.array([3.6,0.05,0]))
        xs = DOWNSTREAM['milliseconds']
        curves = []
        for lib,color,width in [('sympy',SYMPY,7),('symox',SYMOX,3)]:
            curve = m.VMobject(stroke_color=color,stroke_width=width)
            curve.set_points_as_corners([axes.c2p(a,b) for a,b in zip(xs,DOWNSTREAM['records'][lib]['values'])])
            curves.append(curve)
        legend = m.VGroup(T("SymPy",font_size=24,color=SYMPY),T("symox",font_size=24,color=SYMOX))
        legend.arrange(m.RIGHT,buff=0.45).move_to(np.array([3.6,2.0,0]))
        axis_label = T("time (ms)",font_size=21,color=INK2).next_to(axes,m.DOWN,buff=0.15)
        formula = m.MathTex(r"V_{out}(t)=\big(1-e^{-t/(RC)}\big)\,\Theta(t)",font_size=34,color=INK)
        formula.move_to(np.array([0,-2.15,0]))
        sx,sp = (DOWNSTREAM['records'][lib]['time'] for lib in ('symox','sympy'))
        timing_line = T(f"symox {fmt_s(sx)}  ·  SymPy {fmt_s(sp)}  ·  ×{sp/sx:.2f}",
                        font_size=28,color=INK,t2c={'symox':SYMOX,'SymPy':SYMPY})
        timing_line.move_to(np.array([0,-2.95,0]))
        note = T(f"Lcapy {DOWNSTREAM['records']['symox']['lcapy']} · best of 3 · imports excluded · experimental routing",
                 font_size=19,color=INK2).to_edge(m.DOWN,buff=0.22)
        self.play(m.FadeIn(wire),m.FadeIn(labels),m.Create(axes),m.FadeIn(legend),m.FadeIn(axis_label),run_time=0.8)
        self.play(*[m.Create(curve) for curve in curves],m.FadeIn(formula),m.FadeIn(timing_line),m.FadeIn(note),run_time=1.5)
        self.wait(3.0)
        self.clear_scene()

    def integration_tests(self):
        self.heading('Sage & Lcapy: upstream tests', 'Unmodified assertions, the same selection on both engines')
        columns = m.VGroup(
            T('symox',font_size=30,color=SYMOX,weight=m.BOLD).move_to(np.array([1.4,1.85,0])),
            T('SymPy',font_size=30,color=SYMPY,weight=m.BOLD).move_to(np.array([4.6,1.85,0])))
        rows=m.VGroup()
        for label,scope,data,y in [('SageMath','sage/interfaces/sympy.py',SAGE_TESTS,0.7),
                                  ('Lcapy','lcapy/tests',LCAPY_TESTS,-1.05)]:
            name=T(label,font_size=34,color=INK).move_to(np.array([-6.0,y+0.15,0]),aligned_edge=m.LEFT)
            selection=T(scope,font=MONO,font_size=20,color=INK2).next_to(name,m.DOWN,aligned_edge=m.LEFT,buff=0.18)
            group=m.VGroup(name,selection)
            for lib,x in [('symox',1.4),('sympy',4.6)]:
                row=data['records'][lib]
                passed=T(f"{row['passed']}/{row['total']}",font_size=38,color=INK).move_to(np.array([x,y+0.15,0]))
                elapsed=T(fmt_s(row['wall_seconds']),font_size=25,color=INK2).next_to(passed,m.DOWN,buff=0.18)
                group.add(passed,elapsed)
            rows.add(group)
        method=T('Wall times: one complete suite run, including imports and the test runner',font_size=21,color=INK2)
        method.move_to(np.array([0,-2.65,0]))
        scope=T('Selected upstream tests · experimental process-local backend routing',font_size=20,color=INK2)
        scope.to_edge(m.DOWN,buff=0.35)
        self.play(m.FadeIn(columns),m.LaggedStart(*[m.FadeIn(row) for row in rows],lag_ratio=0.35),
                  m.FadeIn(method),m.FadeIn(scope),run_time=1.0)
        self.wait(4.5)
        self.clear_scene()

    # -- 2. integration wall ---------------------------------------------------
    def integrals(self):
        n = len(INTEGRANDS)
        self.heading(f"{n} antiderivatives", "integrate(f, x) · derivatives verified on smooth real intervals")
        times = {lib: [i["time"] for i in CASES[lib]["integrals"]["items"]] for lib in ("symox", "sympy")}
        done_at = {lib: np.cumsum(ts) for lib, ts in times.items()}
        sx_total, sp_total = done_at["symox"][-1], done_at["sympy"][-1]
        playback = max(1.0, sp_total/6.0)
        clock = m.ValueTracker(0.0)

        cols = 9
        pitch_x, pitch_y = 1.5, 0.74
        rows = math.ceil(n/cols)
        wall = m.VGroup()
        for k, tex in enumerate(INTEGRANDS):
            row, col = divmod(k, cols)
            center = np.array([(col - (cols - 1)/2)*pitch_x, -0.62 + ((rows - 1)/2 - row)*pitch_y, 0])
            tile = m.RoundedRectangle(corner_radius=0.07, width=1.38, height=0.64,
                                      stroke_color=GRID, stroke_width=2, fill_color=SYMOX, fill_opacity=0)
            tile.move_to(center)
            text = m.MathTex(tex, font_size=24, color=INK2).move_to(center + np.array([0, 0.03, 0]))
            if text.width > 1.22:
                text.scale_to_fit_width(1.22)
            if text.height > 0.46:
                text.scale_to_fit_height(0.46)
            # SymPy's mark: an orange rule along the tile's bottom edge.
            mark = m.Line(center + np.array([-0.56, -0.25, 0]), center + np.array([0.56, -0.25, 0]),
                          color=SYMPY, stroke_width=5, stroke_opacity=0)

            def light(mob, k=k):
                on = clock.get_value() >= done_at["symox"][k]
                mob.set_fill(SYMOX, opacity=0.24 if on else 0)
                mob.set_stroke(SYMOX if on else GRID)
            tile.add_updater(light)
            mark.add_updater(lambda mob, k=k: mob.set_stroke(opacity=1 if clock.get_value() >= done_at["sympy"][k] else 0))
            wall.add(m.VGroup(tile, text, mark))

        def status(lib, label):
            def content():
                count = int(np.searchsorted(done_at[lib], clock.get_value(), side="right"))
                tail = f" in {fmt_s(done_at[lib][-1])}" if count == n else ""
                return f"{label}  {count}/{n}{tail}"
            return content
        sx_status = live_text(status("symox", "symox"), lambda mob: mob.to_corner(m.UR, buff=0.5),
                              font_size=30, color=SYMOX)
        sp_status = live_text(status("sympy", "SymPy"),
                              lambda mob: mob.to_corner(m.UR, buff=0.5).shift(0.5*m.DOWN),
                              font_size=30, color=SYMPY)
        slow = 1.2/sx_total
        clock_text = live_text(
            lambda: (f"{clock.get_value()*1e3:.2f} ms · slow motion ×{slow:,.0f}"
                     if clock.get_value() <= sx_total else f"{clock.get_value():.1f} s · video ×{playback:g}"),
            lambda mob: mob.move_to(np.array([0, -3.62, 0])), font_size=24, color=INK2)
        self.play(m.FadeIn(wall, lag_ratio=0.01), m.FadeIn(sx_status), m.FadeIn(sp_status), m.FadeIn(clock_text),
                  run_time=0.8)
        # Show the short native wall in slow motion; label the reference clock's
        # playback multiplier if the refreshed run needs compression.
        self.play(clock.animate.set_value(sx_total), run_time=1.2, rate_func=m.linear)
        self.wait(0.4)
        self.play(clock.animate.set_value(sp_total + 0.02), run_time=min(sp_total, 6.0), rate_func=m.linear)
        ratio = T(fmt_x(sp_total/sx_total), font_size=96, weight=m.BOLD, color=SYMOX)
        detail = T(f"symox {fmt_s(sx_total)}  ·  SymPy {fmt_s(sp_total)}", font_size=28, color=INK)
        badge = m.VGroup(ratio, detail).arrange(m.DOWN, buff=0.25).move_to(np.array([0, -0.62, 0]))
        backdrop = m.SurroundingRectangle(badge, buff=0.35, corner_radius=0.2, color=SYMOX, stroke_width=3,
                                          fill_color=BG, fill_opacity=0.94)
        for tile, _, mark in wall:
            tile.clear_updaters()
            mark.clear_updaters()
        self.play(m.FadeOut(wall), m.FadeIn(backdrop), m.FadeIn(ratio, scale=1.5),
                  m.FadeIn(detail), run_time=0.5)
        self.wait(1.5)
        self.clear_scene()

    # -- 3. raw throughput -----------------------------------------------------
    def throughput(self):
        self.heading("Raw throughput", "expand(diff(e, x, 2)) over a stream of distinct random expressions")
        duration = 6.0
        clock = m.ValueTracker(0.0)

        def count(lib):
            samples = TP[lib]["samples"]
            return samples[min(int(clock.get_value()*60), len(samples) - 1)]

        full = TP["symox"]["count"]
        x0, width = -5.6, 11.2
        lanes = {"symox": (1.05, SYMOX, "symox"), "sympy": (-1.55, SYMPY, "SymPy")}
        mobs = []
        for lib, (y, color, label) in lanes.items():
            name = T(label, font_size=36, weight=m.BOLD, color=color)
            name.move_to(np.array([x0, y + 0.85, 0]), aligned_edge=m.LEFT)
            number = live_text(lambda lib=lib: f"{count(lib):,}",
                               lambda mob, y=y: mob.move_to(np.array([x0 + width, y + 0.85, 0]), aligned_edge=m.RIGHT),
                               font=MONO, font_size=72, color=INK)
            track = m.RoundedRectangle(corner_radius=0.08, width=width, height=0.5, stroke_width=0,
                                       fill_color=GRID, fill_opacity=0.6)
            track.move_to(np.array([x0, y, 0]), aligned_edge=m.LEFT)

            def bar(lib=lib, y=y, color=color):
                length = max(width*count(lib)/full, 0.02)
                mob = m.Rectangle(width=length, height=0.5, stroke_width=0, fill_color=color, fill_opacity=1)
                return mob.move_to(np.array([x0, y, 0]), aligned_edge=m.LEFT)
            mobs += [name, number, track, m.always_redraw(bar)]
        timer = live_text(lambda: f"{clock.get_value():.1f} s",
                          lambda mob: mob.move_to(np.array([0, -3.3, 0])), font_size=30, color=INK)
        caption = T("real time · one expression = one full computation", font_size=22, color=INK2)
        caption.move_to(np.array([0, -3.75, 0]))
        self.play(*[m.FadeIn(mob) for mob in mobs], m.FadeIn(timer), m.FadeIn(caption), run_time=0.5)
        self.play(clock.animate.set_value(duration), run_time=duration, rate_func=m.linear)
        ratio = T(fmt_x(TP["symox"]["count"]/TP["sympy"]["count"]), font_size=72, weight=m.BOLD, color=SYMOX)
        ratio.move_to(np.array([0, -0.3, 0]))
        self.play(m.FadeIn(ratio, scale=1.4), run_time=0.5)
        self.wait(1.3)
        self.clear_scene()

    # -- 4. symbolic determinant -----------------------------------------------
    def determinant(self):
        self.heading("Symbolic determinant", "expand(Matrix(n, n, a_ij).det()), one symbol per entry")
        sizes = [3, 4, 5, 6, 7]
        timeout = 120

        # Log time axis: 10 µs … 1000 s.
        lo, hi = -5, 3
        x0, x1 = -0.4, 6.4

        def xpos(seconds):
            return x0 + (x1 - x0)*(math.log10(seconds) - lo)/(hi - lo)

        axis = m.VGroup(m.Line(np.array([x0, -2.0, 0]), np.array([x1, -2.0, 0]), color=GRID, stroke_width=2))
        labels = ["10 µs", "100 µs", "1 ms", "10 ms", "100 ms", "1 s", "10 s", "100 s", "1000 s"]
        for d, label in zip(range(lo, hi + 1), labels):
            x = xpos(10.0**d)
            axis.add(m.Line(np.array([x, -2.0, 0]), np.array([x, -2.12, 0]), color=INK2, stroke_width=1.5))
            if (d - lo) % 2 == 0:
                axis.add(T(label, font_size=20, color=INK2).move_to(np.array([x, -2.4, 0])))
        limit = m.DashedLine(np.array([xpos(timeout), 1.6, 0]), np.array([xpos(timeout), -2.0, 0]),
                             color=INK2, stroke_width=1.5)
        limit_label = T("stopped at 2 min", font_size=20, color=INK2)
        limit_label.next_to(limit, m.UP, buff=0.1)

        lanes = {"symox": (0.75, SYMOX, "symox"), "sympy": (-0.75, SYMPY, "SymPy")}
        names = m.VGroup(*[
            T(label, font_size=32, weight=m.BOLD, color=color).move_to(np.array([x0, y + 0.55, 0]), aligned_edge=m.LEFT)
            for y, color, label in lanes.values()])

        def grid(n):
            side = 3.6
            cell = side/n
            group = m.VGroup()
            for i in range(n):
                for j in range(n):
                    sq = m.Square(cell*0.86, stroke_color=INK2, stroke_width=1.5,
                                  fill_color=INK, fill_opacity=0.08)
                    sq.move_to(np.array([-4.3 + (j - (n - 1)/2)*cell, -0.25 + ((n - 1)/2 - i)*cell, 0]))
                    group.add(sq)
            return group

        def caption(n):
            terms = CASES["symox"][f"det{n}"]["terms"]
            return m.VGroup(
                m.MathTex(f"{n} \\times {n}", font_size=40, color=INK),
                T(f"{terms:,} terms", font_size=26, color=INK2),
            ).arrange(m.RIGHT, buff=0.35).move_to(np.array([-4.3, -2.55, 0]))

        def lane(lib, n):
            y, color, _ = lanes[lib]
            measurement = CASES[lib].get(f"det{n}", {})
            seconds = measurement.get("time")
            capped = seconds is None
            skipped = "skipped" in measurement
            end = xpos(seconds) if not capped else xpos(timeout)
            bar = m.RoundedRectangle(corner_radius=0.06, width=max(end - x0, 0.06), height=0.42,
                                     fill_color=color, fill_opacity=0.35 if capped else 1, stroke_width=0)
            bar.move_to(np.array([x0, y, 0]), aligned_edge=m.LEFT)
            value = T("not measured" if skipped else "> 2 min" if capped else fmt_s(seconds),
                      font_size=26, color=INK)
            if skipped:
                bar.set_opacity(0)
            if capped:
                value.move_to(np.array([xpos(timeout) - 0.15, y + 0.5, 0]), aligned_edge=m.RIGHT)
            else:
                value.next_to(bar, m.RIGHT, buff=0.15)
            return m.VGroup(bar, value)

        n = sizes[0]
        cells, cap = grid(n), caption(n)
        bars = {lib: lane(lib, n) for lib in lanes}
        self.play(m.FadeIn(axis), m.FadeIn(names), m.Create(limit), m.FadeIn(limit_label),
                  m.FadeIn(cells), m.FadeIn(cap), *[m.GrowFromEdge(b, m.LEFT) for b in bars.values()],
                  run_time=0.8)
        self.wait(0.5)
        for n in sizes[1:]:
            new_bars = {lib: lane(lib, n) for lib in lanes}
            self.play(m.Transform(cells, grid(n)), m.Transform(cap, caption(n)),
                      *[m.Transform(bars[lib], new_bars[lib]) for lib in lanes], run_time=0.9)
            self.wait(0.7)
        self.wait(0.8)
        self.clear_scene()

    # -- 5. Taylor series -----------------------------------------------------
    def series(self):
        self.heading("Taylor series", "series(exp(sin(x)), x, 0, 25)")
        sx, sp = CASES["symox"]["taylor_exp"], CASES["sympy"]["taylor_exp"]
        completed = sp.get('time')
        limit = completed if completed is not None else sp['timeout']
        speed = max(1.0, limit/4.0)
        clock = m.ValueTracker(0.0)

        lines = json.load(open(ASSETS / "series_tex.json"))
        body = m.VGroup()
        for tex in lines:
            body.add(m.MathTex(tex, font_size=40, color=INK))
        body.arrange(m.DOWN, aligned_edge=m.LEFT, buff=0.28)
        if body.width > 12.6:
            body.scale_to_fit_width(12.6)
        if body.height > 3.6:
            body.scale_to_fit_height(3.6)
        body.move_to(np.array([0, -1.65, 0]))

        top = 1.35
        sx_name = T("symox", font_size=34, weight=m.BOLD, color=SYMOX)
        sx_name.move_to(np.array([-6.3, top, 0]), aligned_edge=m.LEFT)
        sx_time = T(fmt_s(sx["time"]), font_size=56, weight=m.BOLD, color=INK)
        sx_time.next_to(sx_name, m.DOWN, aligned_edge=m.LEFT, buff=0.15)
        sp_name = T("SymPy", font_size=34, weight=m.BOLD, color=SYMPY)
        sp_name.move_to(np.array([-2.2, top, 0]), aligned_edge=m.LEFT)
        sp_time = live_text(
            lambda: (f"{clock.get_value():.0f} s…" if clock.get_value() < limit
                     else fmt_s(completed) if completed is not None else "time limit"),
            lambda mob: mob.next_to(sp_name, m.DOWN, aligned_edge=m.LEFT, buff=0.15),
            font_size=56, weight=m.BOLD, color=INK)
        sp_note = T("SymPy completed" if completed is not None else f"stopped after {limit/60:g} min", font_size=22, color=INK2)
        sp_note.move_to(np.array([2.6, top, 0]), aligned_edge=m.LEFT)
        speed_label = T(f"video ×{speed:.2f}", font_size=22, color=INK2)
        speed_label.to_corner(m.DR, buff=0.4)
        measured_ratio = limit/sx['time']
        ratio = T(('> ' if completed is None else '') + fmt_x(math.floor(measured_ratio) if completed is None else measured_ratio),
                  font_size=56, weight=m.BOLD, color=SYMOX)
        ratio.next_to(sp_note, m.DOWN, aligned_edge=m.LEFT, buff=0.15)

        self.play(m.FadeIn(sx_name), m.FadeIn(sp_name), m.FadeIn(sp_time), m.FadeIn(speed_label), run_time=0.4)
        self.play(m.FadeIn(sx_time, scale=1.3), m.LaggedStart(*[m.FadeIn(l, shift=0.1*m.RIGHT) for l in body],
                                                            lag_ratio=0.25), run_time=1.0)
        self.play(clock.animate.set_value(limit), run_time=limit/speed, rate_func=m.linear)
        self.play(m.FadeIn(sp_note), m.FadeIn(ratio, scale=1.4), run_time=0.5)
        self.wait(1.6)
        self.clear_scene()

    # -- 6. polynomial system ----------------------------------------------------
    def system(self):
        self.heading("Polynomial system", "solve([x² + y² - 4, y - x⁴ + 1], [x, y])")
        sx, sp = CASES["symox"]["solve"], CASES["sympy"]["solve"]
        axes = m.Axes(x_range=[-2.5, 2.5, 1], y_range=[-2.5, 2.5, 1], x_length=4.8, y_length=4.8,
                      tips=False, axis_config={"color": GRID, "stroke_width": 2})
        axes.move_to(np.array([-3.4, -0.75, 0]))
        th = np.linspace(0, 2*np.pi, 241)
        circle = m.VMobject(stroke_color=INK2, stroke_width=3)
        circle.set_points_as_corners(axes.c2p(2*np.cos(th), 2*np.sin(th)).T)
        quartic = axes.plot(lambda x: x**4 - 1, x_range=[-1.37, 1.37], color=INK, stroke_width=3)
        self.play(m.Create(axes), m.Create(circle), m.Create(quartic), run_time=0.9)

        if 'time' not in sx or 'time' not in sp:
            statuses = m.VGroup()
            for lib,row,color,y0 in [('symox',sx,SYMOX,0.9),('SymPy',sp,SYMPY,-0.9)]:
                name = T(lib,font_size=32,weight=m.BOLD,color=color)
                if 'time' in row:
                    result = T(f"{row['count']} solutions · {fmt_s(row['time'])}",font_size=27,color=INK)
                    detail = T('completed set checked at 60 digits',font_size=19,color=INK2)
                    dots = m.VGroup(*[m.Circle(radius=0.16,stroke_color=color,stroke_width=4)
                                     .move_to(axes.c2p(a,b)) for a,b in row['real']])
                    self.add(dots)
                else:
                    result = T('worker did not complete',font_size=27,color=INK)
                    detail = T(f"{row['budget_seconds']} s process budget · operation time unavailable",font_size=18,color=INK2)
                group = m.VGroup(name,result,detail).arrange(m.DOWN,aligned_edge=m.LEFT,buff=0.18)
                if group.width > 6.0:
                    group.scale_to_fit_width(6.0)
                group.move_to(np.array([0.8,y0,0]),aligned_edge=m.LEFT)
                statuses.add(group)
            note = T('An incomplete run is not an empty solution set.',font_size=23,color=INK2)
            note.to_edge(m.DOWN,buff=0.5)
            self.play(m.FadeIn(statuses),m.FadeIn(note),run_time=0.6)
            self.wait(3.5)
            self.clear_scene()
            return

        clock = m.ValueTracker(0.0)
        speed = max(1.0, max(sx['time'],sp['time'])/4.0)

        def status(r, y, label, color):
            name = T(label, font_size=36, weight=m.BOLD, color=color)
            name.move_to(np.array([0.9, y, 0]), aligned_edge=m.LEFT)
            text = live_text(
                lambda: (f"{r['count']} solutions in {fmt_s(r['time'])}" if clock.get_value() >= r["time"]
                         else "computing…"),
                lambda mob: mob.next_to(name, m.DOWN, aligned_edge=m.LEFT, buff=0.18),
                font_size=30, color=INK)
            return name, text

        def rings(r, color, radius):
            group = m.VGroup()
            for x, y in r["real"]:
                ring = m.Circle(radius=radius, stroke_color=color, stroke_width=4)
                ring.move_to(axes.c2p(x, y))
                group.add(ring)
            return group

        timer = live_text(lambda: f"{clock.get_value():.1f} s",
                          lambda mob: mob.to_corner(m.UR, buff=0.6), font_size=34, color=INK)
        speed_label = T(f"video ×{speed:.2f}", font_size=22, color=INK2)
        speed_label.add_updater(lambda mob: mob.next_to(timer, m.DOWN, aligned_edge=m.RIGHT, buff=0.12))
        note = T("2 real, 6 complex · checked at 60 digits", font_size=22, color=INK2)
        if note.width > 5.8:
            note.scale_to_fit_width(5.8)
        note.move_to(np.array([0.9, -3.0, 0]), aligned_edge=m.LEFT)
        sx_name, sx_text = status(sx, 0.35, "symox", SYMOX)
        sp_name, sp_text = status(sp, -1.35, "SymPy", SYMPY)
        self.play(*[m.FadeIn(mob) for mob in (sx_name, sx_text, sp_name, sp_text, timer, speed_label)], run_time=0.4)
        previous = 0.0
        for index,(row,color,radius) in enumerate(sorted(
                [(sx,SYMOX,0.12),(sp,SYMPY,0.22)],key=lambda item:item[0]['time'])):
            self.play(clock.animate.set_value(row['time']),
                      run_time=max(0.05,(row['time']-previous)/speed),rate_func=m.linear)
            animations = [m.GrowFromCenter(rings(row,color,radius))]
            if index == 0:
                animations.append(m.FadeIn(note))
            self.play(*animations,run_time=0.4)
            previous = row['time']
        ratio = T(fmt_x(sp["time"]/sx["time"]), font_size=56, weight=m.BOLD, color=SYMOX)
        ratio.move_to(np.array([0.9, 1.45, 0]), aligned_edge=m.LEFT)
        self.play(m.FadeIn(ratio, scale=1.4), run_time=0.5)
        self.wait(1.4)
        self.clear_scene()

    # -- 7. N-pendulum ---------------------------------------------------------
    def pendulum(self):
        frames = np.array(DATA["frames"])  # (steps, n+1, 2) at 30 fps
        n = frames.shape[1] - 1
        ready = {lib: DATA["race"][lib][str(n)] for lib in ("symox", "sympy")}
        speed = max(1.0, max(ready.values())/5.0)
        self.heading(f"{n}-link pendulum", "Lagrangian  →  Euler-Lagrange  →  LUsolve  →  lambdify")

        clock = m.ValueTracker(0.0)
        scale = 2.6
        centers = {"symox": np.array([-3.5, 0.3, 0]), "sympy": np.array([3.5, 0.3, 0])}

        def chain(lib, color):
            pivot = centers[lib]

            def build():
                t = clock.get_value()
                live = t >= ready[lib]
                # Once its equations exist the chain swings in real time (the clock runs ×speed).
                idx = int((t - ready[lib])/speed*30) if live else 0
                pts = frames[min(idx, len(frames) - 1)]
                p = [pivot + scale*np.array([x, y, 0]) for x, y in pts]
                opacity = 1.0 if live else 0.25
                group = m.VGroup()
                if live and idx > 1:
                    tail = frames[max(0, idx - 24):idx + 1, -1]
                    trail = m.VMobject(stroke_color=color, stroke_width=3, stroke_opacity=0.45)
                    trail.set_points_as_corners([pivot + scale*np.array([x, y, 0]) for x, y in tail])
                    group.add(trail)
                for a, b in zip(p, p[1:]):
                    group.add(m.Line(a, b, color=INK, stroke_width=4, stroke_opacity=opacity))
                for k, q in enumerate(p[1:]):
                    group.add(m.Dot(q, radius=0.09 if k < n - 1 else 0.13, color=color, fill_opacity=opacity))
                group.add(m.Dot(pivot, radius=0.06, color=INK2))
                return group
            return m.always_redraw(build)

        def header(lib, label, color):
            name = T(label, font_size=36, weight=m.BOLD, color=color)
            name.move_to(centers[lib] + np.array([0, -2.85, 0]))
            status = T("deriving the equations…", font_size=26, color=INK2)
            done = T(f"ready in {fmt_s(ready[lib])}", font_size=26, color=INK)
            for mob in (status, done):
                mob.next_to(name, m.DOWN, buff=0.15)
            return name, status, done

        sx_name, sx_status, sx_done = header("symox", "symox", SYMOX)
        sp_name, sp_status, sp_done = header("sympy", "SymPy", SYMPY)
        divider = m.Line(np.array([0, 2.2, 0]), np.array([0, -2.7, 0]), color=GRID, stroke_width=2)
        timer = live_text(lambda: f"{clock.get_value():.1f} s",
                          lambda mob: mob.to_corner(m.UR, buff=0.6), font_size=34, color=INK)
        speed_label = T(f"video ×{speed:.2f}", font_size=22, color=INK2)
        speed_label.add_updater(lambda mob: mob.next_to(timer, m.DOWN, aligned_edge=m.RIGHT, buff=0.12))

        self.play(*[m.FadeIn(mob) for mob in (divider, sx_name, sp_name, sx_status, sp_status, timer, speed_label)],
                  run_time=0.5)
        self.add(chain("symox", SYMOX), chain("sympy", SYMPY))
        # symox is ready after a few ms: swap its status at once.
        self.play(clock.animate.set_value(ready["symox"]), run_time=0.05, rate_func=m.linear)
        self.remove(sx_status)
        self.add(sx_done)
        self.play(clock.animate.set_value(ready["sympy"]), run_time=(ready["sympy"] - ready["symox"])/speed,
                  rate_func=m.linear)
        self.remove(sp_status)
        self.add(sp_done)
        ratio = T(fmt_x(ready["sympy"]/ready["symox"]), font_size=64, weight=m.BOLD, color=SYMOX)
        ratio.move_to(np.array([0, -3.3, 0]))
        hold = 3.0
        self.play(clock.animate.set_value(ready["sympy"] + hold*speed), m.FadeIn(ratio, scale=1.4),
                  run_time=hold, rate_func=m.linear)
        self.clear_scene()

    # -- 9. speedups, log scale ------------------------------------------------
    def finale(self):
        title = T("Measured speedups", font_size=44, color=INK)
        title.to_edge(m.UP, buff=0.45)
        note = T("selected examples · log scale · best of up to 3 · > = timeout lower bound", font_size=20, color=INK2)
        if note.width > 12.8:
            note.scale_to_fit_width(12.8)
        note.to_edge(m.DOWN, buff=0.4)
        self.play(m.FadeIn(title), m.FadeIn(note), run_time=0.5)

        race = DATA["race"]
        series_limit = CASES["sympy"]["taylor_exp"].get("timeout", 120)
        det_sp = timing("sympy", "det6")
        det_bound = det_sp is None
        det_sp = CASES["sympy"]["det6"].get("timeout", 120) if det_bound else det_sp
        # Only timed-out operations receive a lower-bound marker.
        series_sp = timing('sympy', 'taylor_exp')
        bounded = ({"series exp(sin x), order 25"} if series_sp is None else set()) | ({"6×6 determinant"} if det_bound else set())
        items = [
            (f"{len(INTEGRANDS)} antiderivatives", timing("sympy", "integrals")/timing("symox", "integrals")),
            ("diff + expand throughput", TP["symox"]["count"]/TP["sympy"]["count"]),
            ("6×6 determinant", det_sp/timing("symox", "det6")),
            ("series exp(sin x), order 25", (series_sp if series_sp is not None else series_limit)/timing("symox", "taylor_exp")),
            ("6-link pendulum", race["sympy"]["6"]/race["symox"]["6"]),
        ]
        if timing('sympy','solve') is not None and timing('symox','solve') is not None:
            items.append(('polynomial system',timing('sympy','solve')/timing('symox','solve')))
        else:
            incomplete = T('Polynomial system: native worker incomplete (150 s process budget)',font_size=18,color=INK2)
            incomplete.move_to(np.array([0,-3.15,0]))
            self.add(incomplete)
        items += [("Binomial skewness", ratio4("binom_skew")),
                  ("Sum of two normals", ratio4("sum_normals")),
                  ("Hydrogen 5d normalization", ratio4("hydrogen")),
                  ("Berlekamp-Massey", ratio4("berlekamp"))]
        items.sort(key=lambda it: -it[1])
        decades = math.ceil(math.log10(items[0][1]))
        x0, width = -1.6, 7.2
        top = 2.05
        step = min(0.42, 5.2/len(items))

        ticks = m.VGroup(m.Line(np.array([x0, top + 0.3, 0]), np.array([x0 + width, top + 0.3, 0]),
                                color=GRID, stroke_width=1.5))
        for d in range(decades + 1):
            x = x0 + width*d/decades
            ticks.add(m.Line(np.array([x, top + 0.42, 0]), np.array([x, top + 0.3, 0]), color=INK2, stroke_width=1.5))
            mark = m.MathTex(f"10^{{{d}}}" if d > 1 else ("1" if d == 0 else "10"), font_size=24, color=INK2)
            ticks.add(mark.move_to(np.array([x, top + 0.6, 0])))
        rows = []
        for i, (name, ratio) in enumerate(items):
            y = top - i*step
            tag = T(name, font_size=22, color=INK2)
            tag.move_to(np.array([x0 - 0.2, y, 0]), aligned_edge=m.RIGHT)
            length = max(0.05, width*math.log10(ratio)/decades)
            bar = m.RoundedRectangle(corner_radius=0.06, width=length, height=step*0.62,
                                     fill_color=SYMOX, fill_opacity=1, stroke_width=0)
            bar.move_to(np.array([x0, y, 0]), aligned_edge=m.LEFT)
            value = T(("> " if name in bounded else "") + fmt_x(math.floor(ratio) if name in bounded else ratio), font_size=20, color=INK)
            value.next_to(bar, m.RIGHT, buff=0.12)
            rows.append((tag, bar, value))
        self.play(m.Create(ticks), run_time=0.4)
        self.play(m.LaggedStart(*[m.GrowFromEdge(r[1], m.LEFT) for r in rows], lag_ratio=0.1),
                  m.LaggedStart(*[m.FadeIn(m.VGroup(r[0], r[2])) for r in rows], lag_ratio=0.1),
                  run_time=2.0)
        self.wait(2.2)
        self.clear_scene()

        end = T("Try the examples.", font_size=58, weight=m.BOLD, color=INK).to_edge(m.UP,buff=0.65)
        sub = T("import symox as sp", font=MONO, font_size=35, t2c={"symox": SYMOX})
        sub.next_to(end, m.DOWN, buff=0.3)
        install = m.VGroup(*[T(line,font=MONO,font_size=24,color=INK) for line in [
            'git clone https://github.com/keckj/symox', 'cd symox', 'python -m pip install .',
        ]]).arrange(m.DOWN,aligned_edge=m.LEFT,buff=0.15).move_to(np.array([0,0.3,0]))
        requirements = T('From source · active virtualenv + Rust',font_size=18,color=INK2)
        requirements.move_to(np.array([0,-0.6,0]))
        link = T('github.com/keckj/symox',font_size=30,color=SYMOX).move_to(np.array([0,-1.15,0]))
        reproduce = T('Benchmarks & reproduction: examples/showcase_v4/README.md',font_size=22,color=INK2)
        reproduce.move_to(np.array([0,-1.85,0]))
        scope = T('Selected examples · full mathematical parity remains work in progress',font_size=21,color=INK2)
        scope.move_to(np.array([0,-2.5,0]))
        foot = T(f"symox {DATA['versions']['symox']} · SymPy {DATA['versions']['sympy']} · {DATA['cpu'].split()[2]}"
                 f" · {DATA.get('measured_at', '')[:10]}",
                 font_size=20, color=INK2)
        foot.to_edge(m.DOWN, buff=0.4)
        self.play(m.Write(end), run_time=1.0)
        self.play(*[m.FadeIn(obj) for obj in (sub,install,requirements,link,reproduce,scope,foot)], run_time=0.6)
        self.wait(5.0)

    # -- 8. sympy.stats ---------------------------------------------------------
    def stats(self):
        title = T("Symbolic probability", font_size=44, color=INK).to_corner(m.UL, buff=0.5)
        module = T("sympy.stats", font=MONO, font_size=26, color=INK2)
        module.next_to(title, m.RIGHT, buff=0.45, aligned_edge=m.DOWN)
        cols = {"query": -6.6, "result": -1.55, "symox": 4.0, "sympy": 5.4, "ratio": 6.7}
        top, step = 1.6, 0.5
        head_y = top + 0.62
        head = m.VGroup(
            T("symox", font_size=24, weight=m.BOLD, color=SYMOX)
            .move_to(np.array([cols["symox"], head_y, 0]), aligned_edge=m.RIGHT),
            T("SymPy", font_size=24, weight=m.BOLD, color=SYMPY)
            .move_to(np.array([cols["sympy"], head_y, 0]), aligned_edge=m.RIGHT),
        )
        rule = m.Line(np.array([cols["query"], top + 0.36, 0]), np.array([cols["ratio"], top + 0.36, 0]),
                      color=GRID, stroke_width=1.5)

        sx = [timing4("symox", case) for case, _ in STATS_ROWS]
        sp = [timing4("sympy", case) for case, _ in STATS_ROWS]
        limit = C4["normal_entropy"]["sympy"].get("timeout", 120)
        sp = [t if t is not None else limit for t in sp]
        capped = [timing4("sympy", case) is None for case, _ in STATS_ROWS]
        sx_done, sp_done = np.cumsum(sx), np.cumsum(sp)
        sx_total = sx_done[-1]
        # Keep the refreshed finite rows and timeout readable at explicit speeds.
        finite_end = max(d for d, c in zip(sp_done, capped) if not c)
        fast = max(1.0, finite_end/4.0)
        faster = max(1.0, (sp_done[-1]-finite_end)/2.0)
        clock = m.ValueTracker(0.0)

        def fit(mob, width, height):
            if mob.width > width:
                mob.scale_to_fit_width(width)
            if mob.height > height:
                mob.scale_to_fit_height(height)
            return mob

        queries, answers, sx_cells, sp_cells = m.VGroup(), m.VGroup(), m.VGroup(), m.VGroup()
        for i, (case, tex) in enumerate(STATS_ROWS):
            y = top - i*step
            q = fit(m.MathTex(tex, font_size=28, color=INK2), 4.85, step*0.8)
            queries.add(q.move_to(np.array([cols["query"], y, 0]), aligned_edge=m.LEFT))
            a = fit(m.MathTex("= " + SC4["stats_tex"][case], font_size=28, color=INK), 4.3, step*0.86)
            answers.add(a.move_to(np.array([cols["result"], y, 0]), aligned_edge=m.LEFT))
            sx_cells.add(T(fmt_s(sx[i]), font_size=22, color=INK)
                         .move_to(np.array([cols["symox"], y, 0]), aligned_edge=m.RIGHT))
            ratio = f"> ×{math.floor(sp[i]/sx[i]/1000):,}k" if capped[i] else fmt_x(sp[i]/sx[i])
            sp_cells.add(m.VGroup(
                T(f"> {limit/60:g} min" if capped[i] else fmt_s(sp[i]), font_size=22, color=INK)
                .move_to(np.array([cols["sympy"], y, 0]), aligned_edge=m.RIGHT),
                T(ratio, font_size=24, weight=m.BOLD, color=SYMOX)
                .move_to(np.array([cols["ratio"], y, 0]), aligned_edge=m.RIGHT)))
        for i in range(len(STATS_ROWS)):
            answers[i].add_updater(lambda mob, i=i: mob.set_opacity(1 if clock.get_value() >= sx_done[i] else 0))
            sx_cells[i].add_updater(lambda mob, i=i: mob.set_opacity(1 if clock.get_value() >= sx_done[i] else 0))
            sp_cells[i].add_updater(lambda mob, i=i: mob.set_opacity(1 if clock.get_value() >= sp_done[i] else 0))
            # The row SymPy is working on gets a marker in the margin.
        cursor = m.Triangle(fill_color=SYMPY, fill_opacity=1, stroke_width=0).scale(0.09).rotate(-math.pi/2)

        def follow(mob):
            t = clock.get_value()
            i = int(np.searchsorted(sp_done, t, side="right"))
            on = sx_total < t and i < len(STATS_ROWS)
            mob.move_to(np.array([cols["query"] - 0.22, top - min(i, len(STATS_ROWS) - 1)*step, 0]))
            mob.set_opacity(1 if on else 0)
        cursor.add_updater(follow)

        slow = 1.0/sx_total

        def clock_line():
            t = clock.get_value()
            if t <= sx_total:
                return f"{t*1e3:.2f} ms · slow motion ×{slow:,.0f}"
            return f"{t:.1f} s · video ×{fast if t <= finite_end else faster:g}"
        clock_text = live_text(clock_line, lambda mob: mob.move_to(np.array([0, -3.3, 0])), font_size=24, color=INK2)
        entropy_note = T("Normal entropy: independently checked against the analytic formula",
                         font_size=17, color=INK2).to_edge(m.DOWN, buff=0.10)

        self.play(m.FadeIn(title), m.FadeIn(module), m.FadeIn(head), m.Create(rule),
                  m.LaggedStart(*[m.FadeIn(q, shift=0.1*m.RIGHT) for q in queries], lag_ratio=0.08), run_time=1.0)
        self.add(answers, sx_cells, sp_cells, cursor, clock_text, entropy_note)
        self.play(clock.animate.set_value(sx_total), run_time=1.0, rate_func=m.linear)
        self.wait(0.5)
        self.play(clock.animate.set_value(finite_end), run_time=(finite_end - sx_total)/fast, rate_func=m.linear)
        self.play(clock.animate.set_value(sp_done[-1] + 1e-6), run_time=(sp_done[-1] - finite_end)/faster,
                  rate_func=m.linear)
        self.wait(2.0)
        self.clear_scene()

    # -- 9. sympy.physics -------------------------------------------------------
    def physics(self):
        self.heading("Physics", "sympy.physics: symbolic models")
        beam_x, hyd_x = -3.5, 3.5
        events = {"beam": (timing4("symox", "beam"), timing4("sympy", "beam")),
                  "hydrogen": (timing4("symox", "hydrogen"), timing4("sympy", "hydrogen"))}
        end = max(t for pair in events.values() for t in pair)*1.04
        span = 3.0
        slow = span/end
        clock = m.ValueTracker(0.0)

        def panel_title(x, text, code):
            head = T(text, font_size=30, weight=m.BOLD, color=INK).move_to(np.array([x, 2.05, 0]))
            sub = T(code, font=MONO, font_size=17, color=INK2).next_to(head, m.DOWN, buff=0.16)
            return m.VGroup(head, sub)

        # beam: 10 m, pinned at both ends, 8 kN at 4 m, 3 kN/m everywhere
        def bx(u):
            return beam_x - 2.5 + 0.5*u
        y_beam = 0.75
        beam = m.Rectangle(width=5.0, height=0.12, stroke_width=0, fill_color=INK2, fill_opacity=1)
        beam.move_to(np.array([beam_x, y_beam, 0]))
        supports = m.VGroup(*[
            m.Triangle(stroke_color=INK2, stroke_width=2, fill_opacity=0).scale(0.16)
            .move_to(np.array([bx(u), y_beam - 0.2, 0])) for u in (0, 10)])
        point = m.Arrow(np.array([bx(4), y_beam + 0.65, 0]), np.array([bx(4), y_beam + 0.07, 0]), buff=0,
                        color=INK, stroke_width=5, max_tip_length_to_length_ratio=0.25)
        point_label = T("8 kN", font_size=18, color=INK).next_to(point.get_start(), m.RIGHT, buff=0.1)
        udl = m.VGroup(m.Line(np.array([bx(0), y_beam + 0.45, 0]), np.array([bx(10), y_beam + 0.45, 0]),
                              color=INK2, stroke_width=1.5))
        for u in np.linspace(0, 10, 9):
            udl.add(m.Arrow(np.array([bx(u), y_beam + 0.45, 0]), np.array([bx(u), y_beam + 0.07, 0]), buff=0,
                            color=INK2, stroke_width=2, max_tip_length_to_length_ratio=0.35))
        udl_label = T("3 kN/m", font_size=16, color=INK2).next_to(udl[0], m.UP, aligned_edge=m.LEFT, buff=0.06)
        base_y, depth = -0.55, 1.25
        xs, ys = SC4["beam"]["x"], SC4["beam"]["y"]
        k = depth/max(abs(v) for v in ys)
        baseline = m.DashedLine(np.array([bx(0), base_y, 0]), np.array([bx(10), base_y, 0]), color=GRID, stroke_width=2)
        defl_label = T("deflection × EI", font_size=16, color=INK2).next_to(baseline, m.UP, aligned_edge=m.RIGHT, buff=0.08)
        beam_pts = [np.array([bx(u), base_y + k*v, 0]) for u, v in zip(xs, ys)]
        draw = 0.5/slow  # clock time for a 0.5 s draw-in

        def partial_curve(points, color, start, width=5):
            def build():
                alpha = min(max((clock.get_value() - start)/draw, 0), 1)
                count = max(2, int(alpha*len(points)))
                curve = m.VMobject(stroke_color=color, stroke_width=width, stroke_opacity=1 if alpha > 0 else 0)
                curve.set_points_as_corners(points[:count])
                return curve
            return m.always_redraw(build)
        beam_curve = partial_curve(beam_pts, SYMOX, events["beam"][0])

        # hydrogen: radial density of the 5d state, normalised by integration
        axes = m.Axes(x_range=[0, 70, 10], y_range=[0, 0.045, 0.01], x_length=5.4, y_length=2.1, tips=False,
                      axis_config={"color": GRID, "stroke_width": 2})
        axes.move_to(np.array([hyd_x, -0.6, 0]))
        rs, rho = SC4["hydrogen"]["r"], SC4["hydrogen"]["rho"]
        density = m.VMobject(stroke_color=INK, stroke_width=3)
        density.set_points_as_corners([axes.c2p(a, b) for a, b in zip(rs, rho)])
        area = m.Polygon(*[axes.c2p(a, b) for a, b in zip(rs, rho)], axes.c2p(rs[-1], 0), axes.c2p(0, 0),
                         stroke_width=0, fill_color=SYMOX, fill_opacity=0)
        area.add_updater(lambda mob: mob.set_fill(opacity=0.45*min(max((clock.get_value() - events["hydrogen"][0])/draw, 0), 1)))
        r_label = m.MathTex(r"r / a_0", font_size=22, color=INK2).next_to(axes.x_axis, m.DOWN, aligned_edge=m.RIGHT, buff=0.1)
        rho_label = m.MathTex(r"r^2 R_{52}(r)^2", font_size=24, color=INK2)
        rho_label.move_to(np.array([1.65, 1.1, 0]))
        norm = m.MathTex(r"\int_0^\infty r^2 R_{52}^2\,dr = 1", font_size=26, color=INK)
        norm.move_to(np.array([4.65, 1.1, 0]))
        norm.add_updater(lambda mob: mob.set_opacity(1 if clock.get_value() >= events["hydrogen"][0] else 0))

        def status(x, case, lib, y, color):
            name = T("symox" if lib == "symox" else "SymPy", font_size=26, weight=m.BOLD, color=color)
            name.move_to(np.array([x - 2.6, y, 0]), aligned_edge=m.LEFT)
            done = events[case][0 if lib == "symox" else 1]
            text = live_text(lambda: fmt_s(done) if clock.get_value() >= done else "computing…",
                             lambda mob: mob.move_to(np.array([x + 2.6, y, 0]), aligned_edge=m.RIGHT),
                             font_size=24, color=INK)
            return name, text

        statuses = []
        for x, case in ((beam_x, "beam"), (hyd_x, "hydrogen")):
            statuses += status(x, case, "symox", -2.3, SYMOX) + status(x, case, "sympy", -2.75, SYMPY)
        ratios = []
        for x, case in ((beam_x, "beam"), (hyd_x, "hydrogen")):
            r = T(fmt_x(ratio4(case)), font_size=44, weight=m.BOLD, color=SYMOX).move_to(np.array([x, -3.4, 0]))
            r.add_updater(lambda mob, case=case: mob.set_opacity(1 if clock.get_value() >= events[case][1] else 0))
            ratios.append(r)
        clock_text = live_text(lambda: f"{clock.get_value()*1e3:.1f} ms · slow motion ×{slow:,.0f}",
                               lambda mob: mob.to_corner(m.UR, buff=0.55), font_size=24, color=INK2)
        divider = m.Line(np.array([0, 2.3, 0]), np.array([0, -3.0, 0]), color=GRID, stroke_width=2)

        left = m.VGroup(panel_title(beam_x, "Beam deflection", "Beam(10, E, I).deflection()"),
                        beam, supports, udl, udl_label, point, point_label, baseline, defl_label)
        right = m.VGroup(panel_title(hyd_x, "Hydrogen, 5d orbital", "integrate(r**2*R_nl(5, 2, r)**2, (r, 0, oo))"),
                         axes, density, r_label, rho_label)
        self.play(m.FadeIn(left), m.FadeIn(right), m.Create(divider), run_time=0.8)
        self.add(beam_curve, area, norm, *statuses, *ratios, clock_text)
        self.play(clock.animate.set_value(end), run_time=span, rate_func=m.linear)
        self.wait(1.6)
        self.clear_scene()
        self.module_table("Physics, more", PHYSICS_ROWS)

    # -- 10. sympy.crypto ---------------------------------------------------------
    def crypto(self):
        plain = "MEETMEATTHEUSUALPLACE"*3
        cipher = C4["hill"]["symox"]["result"][:len(plain)]
        title = T("Cryptography", font_size=44, color=INK).to_corner(m.UL, buff=0.5)
        caption = T("encipher_hill(msg, K), K = [[6, 24, 1], [13, 16, 10], [20, 17, 15]] mod 26",
                    font=MONO, font_size=20, color=INK2).next_to(title, m.DOWN, aligned_edge=m.LEFT, buff=0.22)
        line = T(plain, font=MONO, font_size=34, color=INK).move_to(np.array([0, 0.3, 0]))
        out = T(cipher, font=MONO, font_size=34, color=SYMOX).move_to(line)
        if line.width > 12.8:
            line.scale_to_fit_width(12.8)
            out.scale_to_fit_width(12.8)
        self.play(m.FadeIn(title), m.FadeIn(caption), m.FadeIn(line, lag_ratio=0.05), run_time=0.8)
        # Hill works on blocks of three letters: swap block by block.
        blocks = [m.Transform(m.VGroup(*line[i:i + 3]), m.VGroup(*out[i:i + 3])) for i in range(0, len(plain), 3)]
        self.play(m.LaggedStart(*blocks, lag_ratio=0.12), run_time=1.4)
        times = T(f"2,100 letters · symox {fmt_s(timing4('symox', 'hill'))} · SymPy {fmt_s(timing4('sympy', 'hill'))}",
                  font_size=26, color=INK).move_to(np.array([0, -0.7, 0]))
        badge = T(fmt_x(ratio4("hill")), font_size=56, weight=m.BOLD, color=SYMOX).next_to(times, m.DOWN, buff=0.35)
        self.play(m.FadeIn(times), m.FadeIn(badge, scale=1.4), run_time=0.5)
        self.wait(1.2)
        self.clear_scene()
        self.module_table("Cryptography, more", [row for row in CRYPTO_ROWS if row[0] != "hill"])

    def module_table(self, heading, rows):
        cols = {"op": -6.6, "desc": -3.6, "symox": 3.0, "sympy": 4.9, "ratio": 6.7}
        top, step = 1.85, 0.6
        title = T(heading, font_size=44, color=INK).to_corner(m.UL, buff=0.5)
        head = m.VGroup(
            T("symox", font_size=26, weight=m.BOLD, color=SYMOX)
            .move_to(np.array([cols["symox"], top + 0.72, 0]), aligned_edge=m.RIGHT),
            T("SymPy", font_size=26, weight=m.BOLD, color=SYMPY)
            .move_to(np.array([cols["sympy"], top + 0.72, 0]), aligned_edge=m.RIGHT),
        )
        rule = m.Line(np.array([cols["op"], top + 0.4, 0]), np.array([cols["ratio"], top + 0.4, 0]),
                      color=GRID, stroke_width=1.5)
        cells = []
        for i, (case, label, desc) in enumerate(rows):
            y = top - i*step
            a, b = timing4("symox", case), timing4("sympy", case)
            text = T(desc, font_size=22, color=INK)
            if text.width > 5.35:
                text.scale_to_fit_width(5.35)
            cells.append(m.VGroup(
                T(label, font_size=22, color=INK2).move_to(np.array([cols["op"], y, 0]), aligned_edge=m.LEFT),
                text.move_to(np.array([cols["desc"], y, 0]), aligned_edge=m.LEFT),
                T(fmt_s(a), font_size=24, color=INK).move_to(np.array([cols["symox"], y, 0]), aligned_edge=m.RIGHT),
                T(fmt_s(b), font_size=24, color=INK).move_to(np.array([cols["sympy"], y, 0]), aligned_edge=m.RIGHT),
                T(fmt_x(b/a), font_size=26, weight=m.BOLD, color=SYMOX)
                .move_to(np.array([cols["ratio"], y, 0]), aligned_edge=m.RIGHT)))
        self.play(m.FadeIn(title), m.FadeIn(head), m.Create(rule), run_time=0.4)
        self.play(m.LaggedStart(*[m.FadeIn(c, shift=0.15*m.RIGHT) for c in cells], lag_ratio=0.3), run_time=1.8)
        self.wait(1.6)
        self.clear_scene(0.4)

