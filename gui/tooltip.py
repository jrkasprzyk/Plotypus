"""Lightweight hover tooltip for tkinter widgets."""

import tkinter as tk


class Tooltip:
    """Lightweight hover tooltip for any tkinter widget."""

    def __init__(self, widget, text):
        self._widget = widget
        self._text = text
        self._tip = None
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def _show(self, _event=None):
        x = self._widget.winfo_rootx() + self._widget.winfo_width() + 4
        y = self._widget.winfo_rooty()
        self._tip = tw = tk.Toplevel(self._widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        tk.Label(tw, text=self._text, justify="left", background="#ffffe0",
                 relief="solid", borderwidth=1, font=("Segoe UI", 8),
                 wraplength=220).pack()

    def _hide(self, _event=None):
        if self._tip:
            self._tip.destroy()
            self._tip = None
