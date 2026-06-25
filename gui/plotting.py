"""Pure drawing functions for the GUI plots.

Each function takes a matplotlib ``Axes`` plus data and renders into it. They do
not own Tk canvas state — callers are responsible for axes projection (2D vs 3D)
and for calling ``canvas.draw()`` afterwards. This keeps the drawing logic
unit-testable without a Tk root.
"""

import statistics

import matplotlib

_SUBSCRIPT_DIGITS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def subscript(n: int) -> str:
    return str(n).translate(_SUBSCRIPT_DIGITS)


def draw_placeholder(ax):
    """Render the 'run an algorithm' placeholder into a single axes."""
    ax.clear()
    ax.set_facecolor("#f7f7f7")
    ax.text(0.5, 0.5, "Run an algorithm to see results",
            ha="center", va="center", color="#aaaaaa", fontsize=11,
            transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])


def draw_pareto_cartesian(ax, result):
    """Scatter the Pareto front (2D or 3D, depending on axes projection)."""
    ax.clear()
    if not result:
        return

    nobjs = len(result[0])
    color = "#2577c8"

    if nobjs >= 3:
        ax.scatter([r[0] for r in result],
                   [r[1] for r in result],
                   [r[2] for r in result],
                   s=14, c=color, alpha=0.75, depthshade=True)
        ax.set_xlabel("f₁", labelpad=4)
        ax.set_ylabel("f₂", labelpad=4)
        ax.set_zlabel("f₃", labelpad=4)
        if nobjs > 3:
            ax.set_title(f"Pareto Front (f₁–f₃ shown, {nobjs} total)  ·  {len(result)} solutions")
        else:
            ax.set_title(f"Pareto Front  ·  {len(result)} solutions")
    else:
        ax.scatter([r[0] for r in result],
                   [r[1] for r in result],
                   s=16, c=color, alpha=0.75)
        ax.set_xlabel(f"f{subscript(1)}")
        ax.set_ylabel(f"f{subscript(2)}")
        ax.set_title(f"Pareto Front  ·  {len(result)} solutions")
        ax.grid(True, alpha=0.25)


def draw_pareto_parallel(ax, result, current_nobjs):
    """Parallel-axes view of the Pareto front, one polyline per solution."""
    ax.clear()
    if not result:
        return

    # Downsample very large populations for live-update responsiveness.
    if len(result) > 500:
        step = max(1, len(result) // 500)
        rows = result[::step]
    else:
        rows = result

    nobjs = len(rows[0])
    axis_labels = [f"f{subscript(i + 1)}" for i in range(nobjs)]
    draw_parallel_axes(
        ax, rows, axis_labels,
        color_by=[r[0] for r in rows],
        cmap="viridis",
        title=f"Pareto Front (parallel)  ·  {len(result)} solutions",
    )


def draw_progress_cartesian(ax, nfe_hist, size_hist):
    """Line plot of non-dominated set size against NFE."""
    ax.clear()
    ax.plot(nfe_hist, size_hist, color="#2577c8", linewidth=1.8)
    ax.fill_between(nfe_hist, size_hist, alpha=0.12, color="#2577c8")
    ax.set_xlabel("NFE")
    ax.set_ylabel("Non-dominated solutions")
    ax.set_title("Pareto Front Size")
    ax.grid(True, alpha=0.25)


def draw_progress_parallel(ax, last_result, nfe_hist):
    """Parallel-axes view of per-objective min/median/max of the population."""
    ax.clear()
    if not last_result:
        return

    result = last_result
    nobjs = len(result[0])
    cols = [[r[j] for r in result] for j in range(nobjs)]
    mins = [min(c) for c in cols]
    meds = [statistics.median(c) for c in cols]
    maxs = [max(c) for c in cols]
    axis_labels = [f"f{subscript(i + 1)}" for i in range(nobjs)]
    colors = ["#2577c8", "#1e7d1e", "#b87a00"]

    nfe_txt = f"{nfe_hist[-1]:,}" if nfe_hist else "—"
    draw_parallel_axes(
        ax, [mins, meds, maxs], axis_labels,
        line_colors=colors,
        alpha=0.9, lw=1.6,
        title=f"Progress (parallel)  ·  NFE {nfe_txt}",
    )
    ax.legend(["min", "median", "max"], loc="upper right")


def draw_parallel_axes(ax, rows, axis_labels, *, color_by=None,
                       line_colors=None, cmap="viridis", alpha=0.5,
                       lw=0.8, title=None):
    """Parallel-coordinates primitive: normalize each axis to [0, 1] and draw
    one polyline per row across ``len(axis_labels)`` vertical axes."""
    rows = [list(r) for r in rows]
    if not rows:
        return
    n_axes = len(axis_labels)
    cols = [[r[j] for r in rows] for j in range(n_axes)]
    col_min = [min(c) for c in cols]
    col_max = [max(c) for c in cols]
    ranges = [(hi - lo) if (hi - lo) > 1e-12 else 1e-12
              for lo, hi in zip(col_min, col_max)]

    # Gridlines, ticks, limits.
    for j in range(n_axes):
        ax.axvline(j, color="#cccccc", lw=0.8, zorder=0)
    ax.set_xticks(range(n_axes))
    ax.set_xticklabels(axis_labels)
    ax.set_xlim(-0.2, n_axes - 1 + 0.2)
    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks([])

    # Left-axis min/max annotations (original scale).
    ax.text(-0.12, 0.0, f"{col_min[0]:.3g}", ha="right", va="center",
            fontsize=8, color="#666666")
    ax.text(-0.12, 1.0, f"{col_max[0]:.3g}", ha="right", va="center",
            fontsize=8, color="#666666")

    # Resolve per-line colors.
    if line_colors is not None:
        per_line = list(line_colors)
    elif color_by is not None:
        cb = list(color_by)
        cb_lo, cb_hi = min(cb), max(cb)
        cb_rng = (cb_hi - cb_lo) if (cb_hi - cb_lo) > 1e-12 else 1e-12
        cmap_obj = matplotlib.colormaps[cmap]
        per_line = [cmap_obj((v - cb_lo) / cb_rng) for v in cb]
    else:
        per_line = ["#2577c8"] * len(rows)

    xs = list(range(n_axes))
    for row, col in zip(rows, per_line):
        normalized = [(row[j] - col_min[j]) / ranges[j] for j in range(n_axes)]
        ax.plot(xs, normalized, color=col, alpha=alpha, lw=lw)

    if title:
        ax.set_title(title)
