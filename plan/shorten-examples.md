# Plan: Make the CI `examples` job faster

## Goal & Scope
The `examples` job runs every `examples/*.py` end-to-end (~6 min before trims).
Cut its wall-clock time **without degrading examples for human readers**. No
users yet — free to change example NFE.

## Current State (already done)
- MPI portion removed (no `lamboot`/`mpirun`, no `--extras full`/MPI apt libs).
- `examples` gated by a `changes` job (dorny/paths-filter): skipped unless
  `plotypus/`, `examples/`, lockfile/pyproject, or the workflow change.
- `madlibs_moea.py` (interactive) and `parallel_mpi.py` excluded.
- Runtime driver = per-example NFE/seed counts, run **serially**.
- Option B applied to the 7 measured hot spots via `examples/_ci.py`
  (`scaled(full, ci=500)`, gated on the `CI` env var). Cheap examples left
  untouched so the doc-embedded ones stay free of CI plumbing.

## Runtime hot spots (measured locally, full -> CI=true after gating)
- `experimenter.py`: 70s -> 2.4s (`nfe=scaled(10000)`, `seeds=scaled(10, 2)`).
- `portfolio_optimization.py`: 35s -> 3.0s.
- `experimenter_plot.py`: 28s -> 7.9s; `experimenter_parallel.py`: 26s -> 3.5s.
- `tsp.py`: 19s -> 2.3s; `particle_swarm.py`: 12s -> 2.5s
  (`max_iterations = scaled(500, 5)` so OMOPSO's mutation schedule stays
  consistent with the run length).
- `parallel_multiprocess.py`: 10s -> 2.8s.
- Remaining examples: `run(10000)` ≈ 3-5s each — left at full fidelity.

## Options (cheapest -> deepest)

**A. Parallelize CI run (no example edits) — do first.** In
`test-and-publish.yml`, change the serial `xargs -I {}` to `xargs -P 4 -I {}`.
~4x win, zero fidelity loss, reversible. Caveat: interleaved stdout;
`save_results_*.py` write JSON to CWD (distinct names, low collision risk —
verify).

**B. CI quick-mode env gate (touches each script).**
`NFE = 500 if os.getenv("CI") else 10000`. `CI=true` is auto-set by GitHub
Actions — no workflow change. For `experimenter.py` also drop `seeds` (10->2)
under CI. ~20 files; optionally a tiny `examples/_ci.py` helper to keep scripts
clean.

**C. Trim overkill literals outright.** Lower `seeds=10`->3, `100k`->10k
unconditionally. Simplest, but changes what humans see — decide per example.

## Recommendation
**B is done** for the hot spots (see above); they now account for ~25s instead
of ~190s of script time. **A** (xargs -P 4) deliberately skipped for now: with
the hogs gated, the job is dominated by setup + ~15 cheap examples (~1 min
serial), and serial output keeps failure logs readable. Revisit A or C only if
the job is still too slow in practice.

Note: docs literalinclude the experimenter examples, so those three show the
`from _ci import scaled` import. The other 6 doc-embedded examples are
untouched.

## Validation
Push a workflow/examples change so the gated `examples` job runs; confirm green
and note new duration vs ~6 min baseline.
