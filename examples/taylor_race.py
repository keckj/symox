#!/usr/bin/env python
"""Taylor series race: symox vs SymPy, same code, one import apart.

Each library expands f(x) = exp(sin x) cos x / (1 + x**2/3) with
series(f, x, 0, n) for n = 2, 4, 6, ... in a fresh subprocess, for a fixed
time budget. The video shows, side by side, the highest order each library
has reached so far hugging f inside its radius of convergence (sqrt(3), set
by the complex poles at +-i sqrt(3)) and diverging outside it.

    python examples/taylor_race.py --budget 90 --out taylor_race.mp4

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
import threading
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
RADIUS = math.sqrt(3)
X_RANGE = (-3.2, 3.2)
Y_RANGE = (-1.6, 2.6)


def target(sp, x):
    return sp.exp(sp.sin(x))*sp.cos(x)/(1 + x**2/3)


def target_numeric(xs):
    import numpy as np
    return np.exp(np.sin(xs))*np.cos(xs)/(1 + xs**2/3)


def worker(lib, max_order):
    sp = importlib.import_module(lib)
    clear_cache = importlib.import_module("sympy.core.cache").clear_cache if lib == "sympy" else None
    print(json.dumps({"version": sp.__version__}), flush=True)
    x = sp.symbols("x")
    elapsed = 0.0
    for order in range(2, max_order + 1, 2):
        if clear_cache:
            clear_cache()
        start = time.perf_counter()
        polynomial = sp.series(target(sp, x), x, 0, order).removeO()
        elapsed += time.perf_counter() - start
        coeffs = [float(polynomial.coeff(x, k)) for k in range(order)]
        print(json.dumps({"order": order, "elapsed": elapsed, "coeffs": coeffs}), flush=True)


def run_worker(lib, max_order, budget):
    """Rows streamed by a worker until it finishes or the budget runs out."""
    process = subprocess.Popen([sys.executable, __file__, "--worker", lib, str(max_order)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    lines = []

    def pump():
        for line in process.stdout:
            lines.append(json.loads(line))

    reader = threading.Thread(target=pump, daemon=True)
    reader.start()
    reader.join(budget)
    finished = not reader.is_alive()
    process.kill()
    process.wait()
    reader.join(1)
    if finished and process.returncode:
        message = process.stderr.read().strip().splitlines()
        raise SystemExit(f"{lib} worker failed: {message[-1] if message else process.returncode}")
    version = lines[0]["version"] if lines else "?"
    return {"rows": lines[1:], "finished": finished}, version


def cpu_name():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def bench(max_order, budget):
    results, versions = {}, {}
    for lib in LIBS:
        results[lib], versions[lib] = run_worker(lib, max_order, budget)
        rows = results[lib]["rows"]
        last = rows[-1] if rows else None
        state = "terminé" if results[lib]["finished"] else f"arrêté après {budget:.0f} s"
        print(f"{LABELS[lib]:>6}: ordre {last['order'] if last else 0} "
              f"en {last['elapsed'] if last else 0:.3f} s ({state})")
    return {"results": results, "versions": versions, "cpu": cpu_name(), "budget": budget}


def report(data):
    rows = {lib: {row["order"]: row for row in data["results"][lib]["rows"]} for lib in LIBS}
    print(f"\n{'ordre':>5} {'symox (s)':>10} {'SymPy (s)':>10} {'speedup':>9} {'max rel diff':>13}")
    for order in sorted(rows["symox"]):
        fast, slow = rows["symox"][order], rows["sympy"].get(order)
        if slow is None:
            if order > max(rows["sympy"], default=0) + 2:
                continue
            print(f"{order:>5} {fast['elapsed']:>10.4f} {'—':>10}")
            continue
        diff = max(abs(a - b)/max(abs(a), abs(b), 1e-300)
                   for a, b in zip(fast["coeffs"], slow["coeffs"]) if a or b)
        print(f"{order:>5} {fast['elapsed']:>10.4f} {slow['elapsed']:>10.4f} "
              f"{slow['elapsed']/fast['elapsed']:>8.0f}x {diff:>13.1e}")


def format_seconds(seconds):
    if seconds < 1:
        return f"{seconds*1000:.1f} ms" if seconds < 0.01 else f"{seconds*1000:.0f} ms"
    return f"{seconds:.1f} s"


def latest(rows, clock):
    reached = None
    for row in rows:
        if row["elapsed"] > clock:
            break
        reached = row
    return reached


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=TEXT_2, labelsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def draw_progress(ax, data):
    """Order reached against elapsed time; returns the moving cursor line."""
    budget = data["budget"]
    style_axes(ax)
    top = 0
    for lib in LIBS:
        result = data["results"][lib]
        rows = result["rows"]
        if not rows:
            continue
        times = [row["elapsed"] for row in rows]
        orders = [row["order"] for row in rows]
        end = times[-1] if result["finished"] else budget
        ax.step(times + [end], orders + [orders[-1]], where="post", color=COLORS[lib], lw=2,
                label=LABELS[lib], zorder=3)
        ax.plot(times[-1:], orders[-1:], "o", color=COLORS[lib], ms=8, mec=SURFACE, mew=2, zorder=4)
        text = f"{LABELS[lib]} : ordre {orders[-1]} en {format_seconds(times[-1])}"
        if not result["finished"]:
            ax.plot([budget], [orders[-1]], "x", color=COLORS[lib], ms=10, mew=2, zorder=4)
            text = f"{LABELS[lib]} : ordre {orders[-1]}, l'ordre {orders[-1] + 2} dépasse {budget:.0f} s"
        ax.annotate(text, (times[-1], orders[-1]), xytext=(0, 10), textcoords="offset points",
                    ha="right", va="bottom", color=TEXT, fontsize=11)
        top = max(top, orders[-1])
    first = min(row["elapsed"] for lib in LIBS for row in data["results"][lib]["rows"][:1])
    ax.set_xscale("log")
    ax.set_xlim(first/2, budget*3)
    ax.set_ylim(0, top*1.3)
    ax.set_xlabel("temps de calcul cumulé (s, échelle log)", color=TEXT_2, fontsize=11)
    ax.set_ylabel("ordre atteint", color=TEXT_2, fontsize=11)
    ax.grid(True, color=GRID, lw=0.8)
    legend = ax.legend(loc="upper left", frameon=False, fontsize=11)
    for text in legend.get_texts():
        text.set_color(TEXT)
    return ax.axvline(first/2, color=TEXT_2, lw=1, alpha=0.6, zorder=2)


def footer(fig, data):
    versions = data["versions"]
    fig.text(0.5, 0.012,
             f"symox {versions.get('symox', '?')} · SymPy {versions.get('sympy', '?')} · {data['cpu']}"
             " · même script, un processus par bibliothèque, cache SymPy vidé à chaque ordre",
             ha="center", color=TEXT_2, fontsize=9)


def render(data, out, fps, hold=6.0, target_length=25.0, history=4):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
    from numpy.polynomial import polynomial

    results, budget = data["results"], data["budget"]
    last = max(results[lib]["rows"][-1]["elapsed"] if results[lib]["finished"] else budget
               for lib in LIBS)
    # The clocks run in real time unless SymPy would make the video too long.
    speed = max(1.0, last/target_length)
    length = last/speed + hold

    xs = np.linspace(*X_RANGE, 800)
    curves = {lib: [polynomial.polyval(xs, row["coeffs"]) for row in results[lib]["rows"]] for lib in LIBS}

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100, facecolor=SURFACE)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.9, 1], left=0.07, right=0.97,
                            top=0.84, bottom=0.1, hspace=0.3, wspace=0.12)
    fig.text(0.07, 0.94, "Série de Taylor de exp(sin x)·cos x / (1 + x²/3)",
             color=TEXT, fontsize=20, fontweight="bold")
    fig.text(0.07, 0.895, "series(f, x, 0, n) pour n = 2, 4, 6, … : même code, un import d'écart",
             color=TEXT_2, fontsize=12)
    fig.text(0.97, 0.94, f"chronos accélérés ×{speed:.1f}" if speed > 1.05 else "",
             ha="right", color=TEXT_2, fontsize=11)
    clock_text = fig.text(0.97, 0.895, "", ha="right", color=TEXT, fontsize=13, family="monospace")

    artists = {}
    for column, lib in enumerate(LIBS):
        ax = fig.add_subplot(grid[0, column])
        style_axes(ax)
        ax.axvspan(-RADIUS, RADIUS, color=GRID, alpha=0.45, lw=0, zorder=0)
        ax.axhline(0, color=GRID, lw=0.8, zorder=1)
        ax.plot(xs, target_numeric(xs), color=TEXT, lw=2.5, zorder=4)
        ax.set_xlim(*X_RANGE)
        ax.set_ylim(*Y_RANGE)
        ax.set_xticks([-3, -RADIUS, 0, RADIUS, 3], ["−3", "−√3", "0", "√3", "3"])
        ax.text(0.02, 0.95, LABELS[lib], transform=ax.transAxes, color=TEXT, fontsize=16,
                fontweight="bold", va="top")
        ax.text(0.02, 0.84, f"import {lib} as sp", transform=ax.transAxes, color=TEXT_2,
                fontsize=11, family="monospace", va="top")
        status = ax.text(0.98, 0.95, "", transform=ax.transAxes, color=TEXT, fontsize=16,
                         ha="right", va="top")
        detail = ax.text(0.98, 0.84, "", transform=ax.transAxes, color=TEXT_2, fontsize=11,
                         ha="right", va="top")
        ghosts = [ax.plot([], [], color=COLORS[lib], lw=1.2, alpha=0.18*(k + 1), zorder=2)[0]
                  for k in range(history)]
        current, = ax.plot([], [], color=COLORS[lib], lw=2.2, zorder=5)
        artists[lib] = (current, ghosts, status, detail)
        if column == 0:
            ax.text(0, Y_RANGE[0] + 0.12, "rayon de convergence : |x| < √3", ha="center",
                    color=TEXT_2, fontsize=10)

    cursor = draw_progress(fig.add_subplot(grid[1, :]), data)
    first = cursor.get_xdata()[0]
    footer(fig, data)

    def update(frame):
        clock = min(frame/fps*speed, last)
        clock_text.set_text(f"t = {clock:5.1f} s")
        cursor.set_xdata([max(clock, first)]*2)
        for lib, (current, ghosts, status, detail) in artists.items():
            rows = results[lib]["rows"]
            reached = latest(rows, clock)
            index = rows.index(reached) if reached else -1
            current.set_data(xs, curves[lib][index]) if reached else current.set_data([], [])
            for k, ghost in enumerate(ghosts):
                j = index - history + k
                ghost.set_data((xs, curves[lib][j]) if reached and j >= 0 else ([], []))
            if reached:
                status.set_text(f"ordre {reached['order']}")
                detail.set_text(f"atteint en {format_seconds(reached['elapsed'])}")
            else:
                status.set_text("calcul en cours…")
                detail.set_text("")
            done = results[lib]["finished"] and index == len(rows) - 1
            if not done and not results[lib]["finished"] and clock >= budget:
                detail.set_text(f"ordre {rows[-1]['order'] + 2 if rows else 2} : abandon après {budget:.0f} s")
        return []

    count = int(length*fps)
    animation = FuncAnimation(fig, update, frames=count, blit=False)
    if str(out).endswith(".gif"):
        writer = PillowWriter(fps=fps)
    else:
        writer = FFMpegWriter(fps=fps, bitrate=6000, extra_args=["-pix_fmt", "yuv420p"])
    animation.save(out, writer=writer, savefig_kwargs={"facecolor": SURFACE},
                   progress_callback=lambda i, total: print(f"\rframe {i + 1}/{total}", end=""))
    print(f"\nwrote {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--worker", nargs=2, metavar=("LIB", "ORDER"), help=argparse.SUPPRESS)
    parser.add_argument("--max-order", type=int, default=200, help="highest order requested")
    parser.add_argument("--budget", type=float, default=90, help="seconds per library")
    parser.add_argument("--results", default="taylor_race.json",
                        help="benchmark cache, reused when it exists")
    parser.add_argument("--out", default="taylor_race.mp4", help=".mp4 or .gif")
    parser.add_argument("--fps", type=int, default=30)
    options = parser.parse_args()
    if options.worker:
        worker(options.worker[0], int(options.worker[1]))
        return
    if Path(options.results).exists():
        data = json.loads(Path(options.results).read_text())
    else:
        data = bench(options.max_order, options.budget)
        Path(options.results).write_text(json.dumps(data))
    report(data)
    render(data, options.out, options.fps)


if __name__ == "__main__":
    main()
