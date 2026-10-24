#!/usr/bin/env python
"""N-pendulum race: symox vs SymPy, same code, one import apart.

The equations of motion of an N-link pendulum are derived from its Lagrangian
(diff, subs, Matrix.LUsolve) and compiled with lambdify, once per library in a
fresh subprocess so that neither benefits from a warm cache. The results are
rendered as a split-screen video where each pendulum starts swinging as soon
as its library is done, above the time-to-first-frame scaling curve.

    python examples/pendulum_race.py --anim-n 6 --out pendulum_race.mp4

Benchmark results are cached in a JSON file (--results) so that the video can
be re-rendered without re-running SymPy.
"""

import argparse
import importlib
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

LIBS = ("symox", "sympy")
LABELS = {"symox": "symox", "sympy": "SymPy"}
# Dark chart surface, text tokens and categorical slots 1-2, validated
# together for the dark mode (CVD and contrast checks pass).
SURFACE = "#1a1a19"
TEXT = "#f0efec"
TEXT_2 = "#c3c2b7"
GRID = "#3a3a36"
COLORS = {"symox": "#3987e5", "sympy": "#d95926"}


def equations_of_motion(sp, n):
    """Accelerations of an N-link pendulum with unit masses and lengths.

    Returns the lambdify argument list (angles, rates, g) and the solved
    angular accelerations, using only the API shared by both libraries.
    """
    t, g = sp.symbols("t g")
    q = [sp.Function(f"q{i}")(t) for i in range(n)]
    angles = sp.symbols(f"theta0:{n}")
    rates = sp.symbols(f"omega0:{n}")
    accels = sp.symbols(f"alpha0:{n}")
    x = y = kinetic = potential = 0
    for qi in q:
        x += sp.sin(qi)
        y -= sp.cos(qi)
        kinetic += (sp.diff(x, t)**2 + sp.diff(y, t)**2)/2
        potential += g*y
    lagrangian = kinetic - potential
    eqs = [sp.diff(sp.diff(lagrangian, sp.diff(qi, t)), t) - sp.diff(lagrangian, qi) for qi in q]
    # Second derivatives first, so that q' is not replaced inside q''.
    for qi, a in zip(q, accels):
        eqs = [e.subs(sp.diff(qi, t, 2), a) for e in eqs]
    for qi, w in zip(q, rates):
        eqs = [e.subs(sp.diff(qi, t), w) for e in eqs]
    for qi, th in zip(q, angles):
        eqs = [e.subs(qi, th) for e in eqs]
    mass = sp.Matrix([[sp.diff(e, a) for a in accels] for e in eqs])
    forcing = sp.Matrix([-e.subs({a: 0 for a in accels}) for e in eqs])
    solution = mass.LUsolve(forcing)
    return list(angles) + list(rates) + [g], [solution[i] for i in range(n)]


def probe_point(n):
    """A generic state where both libraries' compiled functions are compared."""
    return [0.3 + 0.1*i for i in range(n)] + [0.2 - 0.05*i for i in range(n)] + [9.81]


def worker(lib, n):
    sp = importlib.import_module(lib)
    start = time.perf_counter()
    args, accelerations = equations_of_motion(sp, n)
    derived = time.perf_counter()
    f = sp.lambdify(args, accelerations)
    compiled = time.perf_counter()
    print(json.dumps({
        "derive": derived - start,
        "lambdify": compiled - derived,
        "probe": [float(v) for v in f(*probe_point(n))],
        "version": sp.__version__,
    }))


def run_worker(lib, n, timeout):
    try:
        done = subprocess.run(
            [sys.executable, __file__, "--worker", lib, str(n)],
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}, None
    if done.returncode:
        lines = done.stderr.strip().splitlines()
        return {"status": "error", "message": lines[-1] if lines else ""}, None
    result = json.loads(done.stdout.strip().splitlines()[-1])
    version = result.pop("version")
    result["status"] = "ok"
    return result, version


def cpu_name():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def best_of(lib, n, timeout, repeat=3, budget=5.0):
    """Fastest of up to `repeat` cold runs; long runs are measured once."""
    best, version = run_worker(lib, n, timeout)
    for _ in range(repeat - 1):
        if best["status"] != "ok" or ready_time(best) > budget:
            break
        result, _ = run_worker(lib, n, timeout)
        if result["status"] == "ok" and ready_time(result) < ready_time(best):
            best = result
    return best, version


def bench(max_n, timeout):
    results = {lib: {} for lib in LIBS}
    versions = {}
    for lib in LIBS:
        for n in range(1, max_n + 1):
            result, version = best_of(lib, n, timeout)
            results[lib][n] = result
            if version:
                versions[lib] = version
            if result["status"] != "ok":
                print(f"{LABELS[lib]:>6} N={n}: {result['status']} {result.get('message', '')}")
                break
            print(f"{LABELS[lib]:>6} N={n}: derive {result['derive']:.4f} s, "
                  f"lambdify {result['lambdify']:.4f} s")
    return {"results": results, "versions": versions, "cpu": cpu_name(), "timeout": timeout}


def load(path):
    data = json.loads(Path(path).read_text())
    data["results"] = {lib: {int(n): r for n, r in rows.items()} for lib, rows in data["results"].items()}
    return data


def ready_time(result):
    if result is None or result["status"] != "ok":
        return None
    return result["derive"] + result["lambdify"]


def report(data):
    results = data["results"]
    print(f"\n{'N':>3} {'symox (s)':>11} {'SymPy (s)':>11} {'speedup':>9} {'max |diff|':>11}")
    for n in sorted(results["symox"]):
        fast, slow = results["symox"][n], results["sympy"].get(n)
        a, b = ready_time(fast), ready_time(slow)
        if a is None:
            continue
        if b is None:
            status = slow["status"] if slow else "skipped"
            print(f"{n:>3} {a:>11.4f} {status:>11}")
            continue
        diff = max(abs(u - v) for u, v in zip(fast["probe"], slow["probe"]))
        print(f"{n:>3} {a:>11.4f} {b:>11.4f} {b/a:>8.0f}x {diff:>11.2e}")


def format_seconds(seconds):
    if seconds < 1:
        return f"{seconds*1000:.0f} ms" if seconds >= 0.01 else f"{seconds*1000:.1f} ms"
    return f"{seconds:.1f} s"


def draw_scaling(ax, data, highlight=None):
    results, timeout = data["results"], data["timeout"]
    ax.set_facecolor(SURFACE)
    top = 0
    for lib in LIBS:
        rows = sorted(results[lib].items())
        ns = [n for n, r in rows if ready_time(r) is not None]
        ts = [ready_time(results[lib][n]) for n in ns]
        top = max([top] + ts)
        ax.plot(ns, ts, "-o", color=COLORS[lib], lw=2, ms=8, mec=SURFACE, mew=2,
                label=LABELS[lib], zorder=3)
        if ns:
            ax.annotate(LABELS[lib], (ns[-1], ts[-1]), xytext=(10, 0), textcoords="offset points",
                        va="center", color=TEXT, fontsize=12)
        failed = [(n, r) for n, r in rows if r["status"] != "ok"]
        for n, r in failed:
            ax.plot([n], [timeout], marker="x", color=COLORS[lib], ms=10, mew=2, zorder=3)
            text = f"> {timeout:.0f} s" if r["status"] == "timeout" else "échec"
            ax.annotate(text, (n, timeout), xytext=(10, 0), textcoords="offset points",
                        va="center", color=TEXT_2, fontsize=11)
            top = max(top, timeout)
    if highlight is not None:
        ax.axvline(highlight, color=GRID, lw=8, alpha=0.6, zorder=1)
        a = ready_time(results["symox"].get(highlight))
        b = ready_time(results["sympy"].get(highlight))
        if a and b:
            ax.annotate(f"×{b/a:,.0f}".replace(",", " "), (highlight, math.sqrt(a*b)),
                        xytext=(-12, 0), textcoords="offset points", ha="right", va="center",
                        color=TEXT, fontsize=14, fontweight="bold")
    ax.set_yscale("log")
    ax.set_ylim(top=top*4)
    ns = sorted(results["symox"])
    ax.set_xticks(ns)
    ax.set_xlim(ns[0] - 0.4, ns[-1] + 1.2)
    ax.set_xlabel("nombre de segments N", color=TEXT_2, fontsize=11)
    ax.set_ylabel("secondes (log)", color=TEXT_2, fontsize=11)
    ax.set_title("Temps avant la première image : dérivation + lambdify", color=TEXT,
                 fontsize=13, loc="left")
    ax.grid(True, axis="y", color=GRID, lw=0.8)
    ax.tick_params(colors=TEXT_2, labelsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    legend = ax.legend(loc="upper left", frameon=False, fontsize=11)
    for text in legend.get_texts():
        text.set_color(TEXT)


def footer(fig, data):
    versions = data["versions"]
    fig.text(0.5, 0.012,
             f"symox {versions.get('symox', '?')} · SymPy {versions.get('sympy', '?')} · {data['cpu']}"
             " · même script, un processus par mesure, caches à froid",
             ha="center", color=TEXT_2, fontsize=9)


def simulate(f, n, seconds, fps, substeps=20):
    """RK4 integration; returns the joint positions sampled at fps."""
    import numpy as np
    # Unit links under g/l with l = 1/n: a chain of total length 1.
    g = 9.81*n
    state = np.array([3*math.pi/4 + 0.02*i for i in range(n)] + [0.0]*n)

    def rhs(s):
        return np.concatenate([s[n:], np.asarray(f(*s[:n], *s[n:], g), dtype=float)])

    h = 1/(fps*substeps)
    frames = []
    for _ in range(int(seconds*fps) + 1):
        theta = state[:n]
        frames.append(np.column_stack([
            np.concatenate([[0], np.cumsum(np.sin(theta))/n]),
            np.concatenate([[0], -np.cumsum(np.cos(theta))/n]),
        ]))
        for _ in range(substeps):
            k1 = rhs(state)
            k2 = rhs(state + h/2*k1)
            k3 = rhs(state + h/2*k2)
            k4 = rhs(state + h*k3)
            state = state + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return frames


def render(data, n, out, fps, hold=6.0, target=25.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
    import symox

    results, timeout = data["results"], data["timeout"]
    ready = {lib: ready_time(results[lib].get(n)) for lib in LIBS}
    if ready["symox"] is None:
        raise SystemExit(f"no symox result for N={n}")
    # The clocks run in real time unless SymPy would make the video too long;
    # the pendulums themselves always swing in real time.
    last = max(t if t is not None else timeout for t in ready.values())
    speed = max(1.0, last/target)
    length = last/speed + hold

    args, accelerations = equations_of_motion(symox, n)
    frames = simulate(symox.lambdify(args, accelerations), n, length, fps)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100, facecolor=SURFACE)
    grid = fig.add_gridspec(2, 2, height_ratios=[2.6, 1], left=0.07, right=0.97,
                            top=0.87, bottom=0.1, hspace=0.22, wspace=0.08)
    fig.text(0.07, 0.94, f"Pendule à {n} segments : même code, un import d'écart",
             color=TEXT, fontsize=20, fontweight="bold")
    fig.text(0.07, 0.9, "Lagrangien → diff → subs → Matrix.LUsolve → lambdify, puis intégration RK4",
             color=TEXT_2, fontsize=12)
    clock_note = fig.text(0.97, 0.94, f"chronos accélérés ×{speed:.1f}" if speed > 1.05 else "",
                          ha="right", color=TEXT_2, fontsize=11)

    artists = {}
    for column, lib in enumerate(LIBS):
        # The labels live in the column's corner; the pendulum gets a square
        # box centred in what remains, so its full reach is always visible.
        cell = grid[0, column].get_position(fig)
        fig.text(cell.x0, cell.y1 - 0.02, "●", color=COLORS[lib], fontsize=18, va="center")
        fig.text(cell.x0 + 0.022, cell.y1 - 0.02, LABELS[lib], va="center",
                 color=TEXT, fontsize=16, fontweight="bold")
        fig.text(cell.x0 + 0.022, cell.y1 - 0.065, f"import {lib} as sp", va="center",
                 color=TEXT_2, fontsize=12, family="monospace")
        status = fig.text(cell.x0 + 0.022, cell.y1 - 0.115, "", va="center",
                          color=TEXT, fontsize=15)
        ax = fig.add_axes([cell.x0 + 0.2*cell.width, cell.y0, 0.8*cell.width, cell.height])
        ax.set_facecolor(SURFACE)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")
        trail, = ax.plot([], [], color=COLORS[lib], lw=1.5, alpha=0.35)
        arm, = ax.plot([], [], "-o", color=COLORS[lib], lw=3, ms=8, mec=SURFACE, mew=2)
        ax.plot([0], [0], "s", color=TEXT_2, ms=6)
        artists[lib] = (trail, arm, status)

    draw_scaling(fig.add_subplot(grid[1, :]), data, highlight=n)
    footer(fig, data)
    trail_frames = fps

    def update(frame):
        video_time = frame/fps
        clock = video_time*speed
        for lib, (trail, arm, status) in artists.items():
            t = ready[lib]
            if t is not None and clock >= t:
                index = min(int((video_time - t/speed)*fps), len(frames) - 1)
                tip = [p[-1] for p in frames[max(0, index - trail_frames):index + 1]]
                arm.set_data(frames[index][:, 0], frames[index][:, 1])
                arm.set_alpha(1)
                trail.set_data([p[0] for p in tip], [p[1] for p in tip])
                status.set_text(f"prêt en {format_seconds(t)}")
                status.set_color(TEXT)
            else:
                arm.set_data(frames[0][:, 0], frames[0][:, 1])
                arm.set_alpha(0.3)
                trail.set_data([], [])
                if t is None and clock >= timeout:
                    status.set_text(f"abandon après {timeout:.0f} s")
                else:
                    status.set_text(f"calcul en cours... {clock:.1f} s")
                status.set_color(TEXT_2)
        return [a for group in artists.values() for a in group] + [clock_note]

    count = int(length*fps)
    animation = FuncAnimation(fig, update, frames=count, blit=False)
    if str(out).endswith(".gif"):
        writer = PillowWriter(fps=fps)
    else:
        writer = FFMpegWriter(fps=fps, bitrate=6000, extra_args=["-pix_fmt", "yuv420p"])
    animation.save(out, writer=writer, savefig_kwargs={"facecolor": SURFACE},
                   progress_callback=lambda i, total: print(f"\rframe {i + 1}/{total}", end=""))
    print(f"\nwrote {out}")

    still = plt.figure(figsize=(10, 5), dpi=150, facecolor=SURFACE)
    draw_scaling(still.add_axes([0.08, 0.15, 0.86, 0.75]), data, highlight=n)
    footer(still, data)
    still_path = Path(out).with_name(Path(out).stem + "_scaling.png")
    still.savefig(still_path, facecolor=SURFACE)
    print(f"wrote {still_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--worker", nargs=2, metavar=("LIB", "N"), help=argparse.SUPPRESS)
    parser.add_argument("--max-n", type=int, default=10, help="largest chain benchmarked")
    parser.add_argument("--anim-n", type=int, default=6, help="chain shown in the video")
    parser.add_argument("--timeout", type=float, default=120, help="seconds per measurement")
    parser.add_argument("--results", default="pendulum_race.json",
                        help="benchmark cache, reused when it exists")
    parser.add_argument("--out", default="pendulum_race.mp4", help=".mp4 or .gif")
    parser.add_argument("--fps", type=int, default=30)
    options = parser.parse_args()
    if options.worker:
        worker(options.worker[0], int(options.worker[1]))
        return
    if Path(options.results).exists():
        data = load(options.results)
    else:
        data = bench(options.max_n, options.timeout)
        Path(options.results).write_text(json.dumps(data, indent=1))
    report(data)
    render(data, options.anim_n, options.out, options.fps)


if __name__ == "__main__":
    main()
