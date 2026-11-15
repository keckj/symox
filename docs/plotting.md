# Plotting

[README](../README.md)

Plot functions sample native Symox expressions with NumPy. Rendering uses the
optional Matplotlib dependency:

```bash
python -m pip install 'symox[plotting]'
```

```python
from symox import symbols, sin, plot, plot3d, PlotGrid

x, y = symbols('x y')
curve = plot(sin(x), (x, -6, 6), n=500, show=False)
surface = plot3d(x*y, (x, -2, 2), (y, -2, 2), n=(40, 50), show=False)
figure = PlotGrid(1, 2, curve, surface, show=False, size=(10, 4))
figure.save('plots.png')
figure.close()
```

`plot`, `plot_parametric`, `plot3d`, `plot3d_parametric_line`,
`plot3d_parametric_surface`, `plot_contour` and `plot_implicit` are available from
`symox` and `symox.plotting`. The latter also exports `plot_backends`,
`PygletPlot`, `PlotGrid` and `textplot`.

## Samples and composition

Ranges are `(symbol, start, end)` with finite real endpoints; as in SymPy, a
degenerate range with equal endpoints is sampled at that single value. Missing
ranges are inferred from free symbols, with endpoints -10 and 10. A scalar
plot can contain several expressions, or tuples containing an expression,
a local range and an optional label. Parametric curves take two or three
coordinate expressions; parametric surfaces take three coordinates and two
ranges. Vectorized Python numerical callables are supported with explicit ranges.
Expressions must use heads supported by the native numerical compiler.

`n` controls the sample count, with a tuple specifying one count per surface
parameter. `n1` and `n2` override individual counts. Constants broadcast to the
grid. Nonfinite and nonreal samples become NaN. `detect_poles=True` also breaks
large opposite-sign jumps in scalar lines; this is a heuristic, not a proof
that every discontinuity has been found.

A plot contains series accessible through `p[index]`, iteration, `append` and
`extend`. `get_data()` returns numeric arrays: `(x,y)` for lines, `(x,y,t)` or
`(x,y,z,t)` for parametric curves, and `(X,Y,Z)` for surfaces and contours.
Implicit series return `(X,Y,Z,mode)` where mode is `contour` for a zero set or
`contourf` for an inequality mask. Rectangular surface arrays have shape
`(n2,n1)`, using NumPy's `meshgrid(indexing='xy')` convention.

`show=False` allows data access without importing Matplotlib; NumPy itself is
imported on first use, so `import symox` does not load it. `save` supports
Matplotlib formats such as PNG, SVG and PDF. `close` releases the figure.
`title`, axis labels and limits, `legend`, `line_color`, `surface_color`, `size`
and per-series `rendering_kw` control rendering. As in SymPy's Matplotlib
backend, `markers`, `annotations`, `rectangles` and `fill` take Matplotlib
keyword dictionaries for decorations, as used by the beam, cable, truss and
arch drawings. `params` supplies static parameter substitutions for top-level
expressions.

## Scope

Sampling is uniform and deterministic. Adaptive interval subdivision,
interactive parameter widgets, complex-domain plots and animations are not
implemented. The compatibility options `adaptive` and `depth` do not activate
adaptive sampling; series report `adaptive=False` and `sampling='uniform'`.
Implicit equality and inequality plots use sampled contours, so thin or
isolated components can be missed. `And`, `Or`, `Not` and `Xor` combinations
of inequalities are drawn as filled regions (`contourf`) from the Boolean
combination of each inequality's sampled mask; a combination containing an
equality raises `NotImplementedError`, as in SymPy.

`PygletPlot` is a legacy entry point for these Matplotlib figures, selecting a
surface when two ranges are supplied. It does not create a Pyglet/OpenGL
window. `plot_backends` provides `matplotlib`, `default` and `text` backends.
The text backend supports scalar 2D lines; `textplot(expr,a,b,W=55,H=21)` prints
an ASCII graph directly without Matplotlib.
