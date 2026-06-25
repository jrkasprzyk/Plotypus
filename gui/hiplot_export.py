"""HiPlot export helpers.

Heavy/optional dependencies (``parasolpy``, ``hiplot``, ``pandas``) are imported
lazily inside the functions so that ``import gui.hiplot_export`` succeeds even
when they are not installed.
"""

import pathlib
import tempfile
import webbrowser


def hiplot_available() -> bool:
    try:
        import parasolpy  # noqa: F401
        import hiplot  # noqa: F401
        return True
    except Exception:
        return False


def open_in_hiplot(last_result):
    """Render the latest Pareto front in HiPlot and open it in the browser.

    Assumes minimization for all objectives (consistent with Plotypus's
    benchmark registry). If a future GUI exposes maximization problems, flip the
    corresponding entries in ``obj_directions`` below.

    Raises on failure; the caller is responsible for surfacing the error.
    """
    import pandas as pd
    import parasolpy

    nobjs = len(last_result[0])
    cols = [f"f{i + 1}" for i in range(nobjs)]
    df = pd.DataFrame(last_result, columns=cols)
    exp = parasolpy.parallel_plot_hp(
        df,
        obj_names=cols,
        obj_directions=["minimize"] * nobjs,
    )
    tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
    tmp.close()
    exp.to_html(tmp.name)
    webbrowser.open(pathlib.Path(tmp.name).resolve().as_uri())
