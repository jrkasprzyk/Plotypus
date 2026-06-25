#!/usr/bin/env python3
"""Plotypus Explorer — GUI for multiobjective optimization."""

import queue
import tkinter as tk
import webbrowser
from tkinter import ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from gui import plotting
from gui.hiplot_export import hiplot_available, open_in_hiplot
from gui.registries import ALGORITHMS, PROBLEMS, make_problem
from gui.tooltip import Tooltip
from gui.worker import Worker


# ── Main application ──────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Plotypus Explorer")
        self.geometry("1150x700")
        self.minsize(800, 500)

        style = ttk.Style()
        style.theme_use("vista")

        self._worker: Worker | None = None
        self._queue: queue.Queue = queue.Queue()
        self._nfe_hist: list[int] = []
        self._size_hist: list[int] = []
        self._max_nfe = 10000
        self._current_nobjs = 2
        self._last_result: list[list[float]] = []
        self._pareto_view = tk.StringVar(value="Cartesian")
        self._progress_view = tk.StringVar(value="Cartesian")

        self._build_ui()
        self._build_menu()
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_menu(self):
        menubar = tk.Menu(self)
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="Documentation",
                              command=lambda: webbrowser.open(
                                  "https://plotypus.readthedocs.io"))
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)
        self.config(menu=menubar)

    def _show_about(self):
        import plotypus
        messagebox.showinfo(
            "About Plotypus Explorer",
            f"Plotypus Explorer\n"
            f"Version {plotypus.__version__}\n\n"
            f"A GUI for multiobjective optimization\n"
            f"built on the Plotypus framework.\n\n"
            f"https://plotypus.readthedocs.io"
        )

    def _build_ui(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        self._build_sidebar()
        self._build_plots()

    def _build_sidebar(self):
        side = ttk.Frame(self, padding=(12, 12, 8, 12))
        side.grid(row=0, column=0, sticky="ns")

        r = 0
        ttk.Label(side, text="Plotypus Explorer", font=("Segoe UI", 12, "bold")).grid(
            row=r, column=0, columnspan=2, pady=(0, 14)); r += 1

        # ── Problem ──
        ttk.Label(side, text="Problem", font=("Segoe UI", 9, "bold")).grid(
            row=r, column=0, columnspan=2, sticky="w", pady=(0, 2)); r += 1

        ttk.Label(side, text="Problem").grid(row=r, column=0, sticky="w")
        self._problem_var = tk.StringVar(value="DTLZ2")
        prob_cb = ttk.Combobox(side, textvariable=self._problem_var,
                               values=list(PROBLEMS), state="readonly", width=10)
        prob_cb.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        prob_cb.bind("<<ComboboxSelected>>", self._on_problem_change)
        Tooltip(prob_cb, "Benchmark problem to optimize. DTLZ/WFG support variable objective counts; ZDT/UF/CF are fixed at 2.")

        ttk.Label(side, text="Objectives").grid(row=r, column=0, sticky="w")
        self._nobjs_var = tk.IntVar(value=2)
        self._nobjs_spin = ttk.Spinbox(side, from_=2, to=8,
                                       textvariable=self._nobjs_var, width=5)
        self._nobjs_spin.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        Tooltip(self._nobjs_spin, "Number of objectives for the selected problem (disabled for problems with a fixed objective count).")

        ttk.Separator(side, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8); r += 1

        # ── Algorithm ──
        ttk.Label(side, text="Algorithm", font=("Segoe UI", 9, "bold")).grid(
            row=r, column=0, columnspan=2, sticky="w", pady=(0, 2)); r += 1

        ttk.Label(side, text="Algorithm").grid(row=r, column=0, sticky="w")
        self._algo_var = tk.StringVar(value="NSGA-II")
        algo_cb = ttk.Combobox(side, textvariable=self._algo_var,
                               values=list(ALGORITHMS), state="readonly", width=10)
        algo_cb.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        Tooltip(algo_cb, "Multiobjective evolutionary algorithm to run.")

        ttk.Label(side, text="Pop size").grid(row=r, column=0, sticky="w")
        self._pop_var = tk.IntVar(value=100)
        pop_spin = ttk.Spinbox(side, from_=10, to=2000, increment=10,
                               textvariable=self._pop_var, width=7)
        pop_spin.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        Tooltip(pop_spin, "Number of candidate solutions maintained each generation. Larger values improve coverage but increase runtime.")

        ttk.Label(side, text="Max NFE").grid(row=r, column=0, sticky="w")
        self._nfe_var = tk.IntVar(value=10000)
        nfe_entry = ttk.Entry(side, textvariable=self._nfe_var, width=8)
        nfe_entry.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        Tooltip(nfe_entry, "Maximum number of function evaluations (NFE). The algorithm stops after this many calls to the problem's objective functions.")

        ttk.Label(side, text="Update every").grid(row=r, column=0, sticky="w")
        self._freq_var = tk.IntVar(value=500)
        freq_spin = ttk.Spinbox(side, from_=50, to=5000, increment=50,
                                textvariable=self._freq_var, width=7)
        freq_spin.grid(row=r, column=1, pady=2, padx=(6, 0), sticky="e"); r += 1
        Tooltip(freq_spin, "Refresh the plot every N function evaluations. Lower values give smoother live updates; higher values are faster.")

        ttk.Separator(side, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8); r += 1

        # ── Controls ──
        self._run_btn = ttk.Button(side, text="▶  Run", command=self._run, width=16)
        self._run_btn.grid(row=r, column=0, columnspan=2, pady=2); r += 1

        self._stop_btn = ttk.Button(side, text="■  Stop", command=self._stop,
                                    state="disabled", width=16)
        self._stop_btn.grid(row=r, column=0, columnspan=2, pady=2); r += 1

        ttk.Separator(side, orient="horizontal").grid(
            row=r, column=0, columnspan=2, sticky="ew", pady=8); r += 1

        # ── Progress ──
        ttk.Label(side, text="Progress", font=("Segoe UI", 9, "bold")).grid(
            row=r, column=0, columnspan=2, sticky="w", pady=(0, 4)); r += 1

        self._progress = ttk.Progressbar(side, length=180, mode="determinate")
        self._progress.grid(row=r, column=0, columnspan=2, pady=2); r += 1

        self._nfe_label = ttk.Label(side, text="NFE: —", font=("Consolas", 9))
        self._nfe_label.grid(row=r, column=0, columnspan=2, pady=(2, 0)); r += 1

        self._status_label = ttk.Label(side, text="Ready", foreground="#666666")
        self._status_label.grid(row=r, column=0, columnspan=2, pady=(2, 0)); r += 1

    def _build_plots(self):
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        nb = ttk.Notebook(right)
        nb.grid(sticky="nsew")

        # Tab 1 — Pareto front
        t1 = ttk.Frame(nb)
        nb.add(t1, text="  Pareto Front  ")
        t1.rowconfigure(1, weight=1)
        t1.columnconfigure(0, weight=1)

        bar1 = ttk.Frame(t1)
        bar1.grid(row=0, column=0, sticky="ew", padx=4, pady=(4, 0))
        ttk.Label(bar1, text="View:").pack(side="left")
        pareto_combo = ttk.Combobox(bar1, textvariable=self._pareto_view,
                                    values=["Cartesian", "Parallel axes"],
                                    state="readonly", width=14)
        pareto_combo.pack(side="left", padx=(4, 0))
        pareto_combo.bind("<<ComboboxSelected>>", self._on_pareto_view_change)
        Tooltip(pareto_combo, "Switch between Cartesian (scatter/line) and parallel-axes (one polyline per solution across normalized objective axes) views.")

        self._hiplot_btn = ttk.Button(bar1, text="Open in HiPlot…",
                                      command=self._open_in_hiplot)
        self._hiplot_btn.pack(side="left", padx=(8, 0))
        if not hiplot_available():
            self._hiplot_btn.configure(state="disabled")
            Tooltip(self._hiplot_btn,
                    "Install `parasolpy` and `hiplot` (pip install parasolpy hiplot) "
                    "to enable interactive HiPlot export.")
        else:
            Tooltip(self._hiplot_btn,
                    "Render the latest Pareto front in HiPlot and open it in your browser.")

        self._fig_pareto = Figure(tight_layout=True)
        self._ax_pareto = self._fig_pareto.add_subplot(111)
        canvas1 = FigureCanvasTkAgg(self._fig_pareto, master=t1)
        canvas1.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        tb1 = NavigationToolbar2Tk(canvas1, t1, pack_toolbar=False)
        tb1.grid(row=2, column=0, sticky="ew")
        self._canvas_pareto = canvas1

        # Tab 2 — Progress
        t2 = ttk.Frame(nb)
        nb.add(t2, text="  Progress  ")
        t2.rowconfigure(1, weight=1)
        t2.columnconfigure(0, weight=1)

        bar2 = ttk.Frame(t2)
        bar2.grid(row=0, column=0, sticky="ew", padx=4, pady=(4, 0))
        ttk.Label(bar2, text="View:").pack(side="left")
        progress_combo = ttk.Combobox(bar2, textvariable=self._progress_view,
                                      values=["Cartesian", "Parallel axes"],
                                      state="readonly", width=14)
        progress_combo.pack(side="left", padx=(4, 0))
        progress_combo.bind("<<ComboboxSelected>>", self._on_progress_view_change)
        Tooltip(progress_combo, "Switch between Cartesian (line) and parallel-axes (per-objective min/median/max of the latest population) views.")

        self._fig_prog = Figure(tight_layout=True)
        self._ax_prog = self._fig_prog.add_subplot(111)
        canvas2 = FigureCanvasTkAgg(self._fig_prog, master=t2)
        canvas2.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        tb2 = NavigationToolbar2Tk(canvas2, t2, pack_toolbar=False)
        tb2.grid(row=2, column=0, sticky="ew")
        self._canvas_prog = canvas2

        self._draw_placeholders()

    # ── Event handlers ────────────────────────────────────────────────────────

    def _on_problem_change(self, _event=None):
        _, accepts_nobjs = PROBLEMS[self._problem_var.get()]
        self._nobjs_spin.configure(state="normal" if accepts_nobjs else "disabled")

    def _run(self):
        prob_name = self._problem_var.get()
        algo_key  = ALGORITHMS[self._algo_var.get()]
        nobjs     = self._nobjs_var.get()
        pop_size  = self._pop_var.get()
        max_nfe   = self._nfe_var.get()
        freq      = self._freq_var.get()

        self._max_nfe = max_nfe
        self._current_nobjs = nobjs
        self._nfe_hist.clear()
        self._size_hist.clear()

        try:
            problem = make_problem(prob_name, nobjs)
        except Exception as e:
            messagebox.showerror("Problem error", str(e))
            return

        self._reset_pareto_axes(nobjs)

        self._run_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal")
        self._progress["value"] = 0
        self._set_status("Running…", "#1a6bbf")

        self._worker = Worker(self._queue, problem, algo_key, pop_size, max_nfe, freq)
        self._worker.start()
        self.after(100, self._poll)

    def _stop(self):
        if self._worker:
            self._worker.request_stop()
        self._set_status("Stopping…", "#b87a00")

    def _on_close(self):
        if self._worker:
            self._worker.request_stop()
        self.destroy()

    # ── Queue polling ─────────────────────────────────────────────────────────

    def _poll(self):
        try:
            while True:
                msg = self._queue.get_nowait()
                t = msg["t"]
                if t == "update":
                    self._apply_update(msg)
                elif t == "done":
                    self._apply_update(msg)
                    self._on_finish("Done  ✓", "#1e7d1e")
                    return
                elif t == "stopped":
                    self._on_finish("Stopped", "#b87a00")
                    return
                elif t == "error":
                    self._on_finish("Error", "red")
                    messagebox.showerror("Algorithm error",
                                         f"{msg['msg']}\n\n{msg.get('tb', '')}")
                    return
        except queue.Empty:
            pass

        if self._worker and self._worker.is_alive():
            self.after(100, self._poll)

    # ── Update helpers ────────────────────────────────────────────────────────

    def _apply_update(self, msg):
        nfe, result = msg["nfe"], msg["result"]
        self._nfe_hist.append(nfe)
        self._size_hist.append(len(result))
        self._last_result = result
        self._progress["value"] = min(100, 100 * nfe / self._max_nfe)
        self._nfe_label.configure(text=f"NFE: {nfe:,} / {self._max_nfe:,}")
        self._redraw_pareto(result)
        self._redraw_progress()

    def _on_finish(self, text, color):
        self._set_status(text, color)
        self._run_btn.configure(state="normal")
        self._stop_btn.configure(state="disabled")

    def _set_status(self, text, color):
        self._status_label.configure(text=text, foreground=color)

    # ── Plot helpers ──────────────────────────────────────────────────────────

    def _draw_placeholders(self):
        for ax, canvas in [(self._ax_pareto, self._canvas_pareto),
                           (self._ax_prog, self._canvas_prog)]:
            plotting.draw_placeholder(ax)
            canvas.draw()

    def _reset_pareto_axes(self, nobjs):
        """Swap the Pareto axes between 2D and 3D projections.

        Stays on ``App`` because it recreates a Tk-managed figure subplot; the
        pure plotting helpers must not assume a projection.
        """
        self._fig_pareto.clear()
        if self._pareto_view.get() == "Parallel axes":
            self._ax_pareto = self._fig_pareto.add_subplot(111)
        elif nobjs >= 3:
            self._ax_pareto = self._fig_pareto.add_subplot(111, projection="3d")
        else:
            self._ax_pareto = self._fig_pareto.add_subplot(111)
        self._ax_prog.clear()
        self._canvas_pareto.draw()
        self._canvas_prog.draw()

    # ── View-toggle handlers ──────────────────────────────────────────────────

    def _on_pareto_view_change(self, _event=None):
        if self._last_result:
            nobjs = len(self._last_result[0])
            is_3d = getattr(self._ax_pareto, "name", None) == "3d"
            want_3d = self._pareto_view.get() != "Parallel axes" and nobjs >= 3
            if want_3d != is_3d:
                self._reset_pareto_axes(nobjs)
        self._redraw_pareto(self._last_result)

    def _on_progress_view_change(self, _event=None):
        self._redraw_progress()

    # ── Pareto dispatcher ─────────────────────────────────────────────────────

    def _redraw_pareto(self, result):
        if self._pareto_view.get() == "Parallel axes":
            self._redraw_pareto_parallel(result)
        else:
            self._redraw_pareto_cartesian(result)

    def _redraw_pareto_cartesian(self, result):
        if not result:
            self._ax_pareto.clear()
            self._canvas_pareto.draw()
            return

        nobjs = len(result[0])
        is_3d = getattr(self._ax_pareto, "name", None) == "3d"
        if (nobjs >= 3) != is_3d:
            self._reset_pareto_axes(nobjs)

        plotting.draw_pareto_cartesian(self._ax_pareto, result)
        self._canvas_pareto.draw()

    def _redraw_pareto_parallel(self, result):
        # Ensure axes are 2D for parallel view.
        if getattr(self._ax_pareto, "name", None) == "3d":
            nobjs = len(result[0]) if result else self._current_nobjs
            self._reset_pareto_axes(nobjs)

        plotting.draw_pareto_parallel(self._ax_pareto, result, self._current_nobjs)
        self._canvas_pareto.draw()

    # ── Progress dispatcher ───────────────────────────────────────────────────

    def _redraw_progress(self):
        if self._progress_view.get() == "Parallel axes":
            self._redraw_progress_parallel()
        else:
            self._redraw_progress_cartesian()

    def _redraw_progress_cartesian(self):
        plotting.draw_progress_cartesian(self._ax_prog, self._nfe_hist, self._size_hist)
        self._canvas_prog.draw()

    def _redraw_progress_parallel(self):
        plotting.draw_progress_parallel(self._ax_prog, self._last_result, self._nfe_hist)
        self._canvas_prog.draw()

    # ── HiPlot export ─────────────────────────────────────────────────────────

    def _open_in_hiplot(self):
        if not self._last_result:
            messagebox.showinfo("HiPlot export", "Run an algorithm first.")
            return
        try:
            open_in_hiplot(self._last_result)
        except Exception as e:
            messagebox.showerror("HiPlot export failed", str(e))


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
