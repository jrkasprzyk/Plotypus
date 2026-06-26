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

## Runtime hot spots (from grep of NFE/seed counts)
- Most examples: `run(10000)`.
- `experimenter.py`: `nfe=10000, seeds=10` x multiple algos/problems — biggest hog.
- `experimenter_parallel.py`, `experimenter_plot.py`: `nfe=10000`.
- `portfolio_optimization.py`: `run(100_000)`; `tsp.py`: `run(100000)`.
- `parallel_multiprocess.py`: `run(10000)` + `sum(range(100000))` busywork.

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
Do **A** now. If still slow, layer **B** (`CI` is free). Reserve **C** for
examples where big NFE isn't needed for the demo output.

## Validation
Push a workflow/examples change so the gated `examples` job runs; confirm green
and note new duration vs ~6 min baseline.
