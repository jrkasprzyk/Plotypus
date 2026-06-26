# Copyright 2026 Joseph Kasprzyk
#
# This file is part of plotypus, a fork of Platypus -- a Python module for
# designing and using evolutionary algorithms (EAs) and multiobjective
# evolutionary algorithms (MOEAs).  Platypus is Copyright 2015-2024 David
# Hadka.
#
# plotypus is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# plotypus is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with plotypus.  If not, see <http://www.gnu.org/licenses/>.

"""Equality-constrained DTLZ test problems (Cuate et al., 2020).

This family augments the classic DTLZ problems with equality constraints,
restricting the Pareto-optimal set to one or more (hyper)spherical rings in
the position-variable subspace.  Two constraint patterns are provided:

* ``Eq1`` problems impose a single ring constraint (radius ``r = 0.4``).
* ``Eq2`` problems impose two intersecting ring constraints (radius
  ``r = 0.5``), reducing the feasible front to the points where the rings
  meet.

The objective functions are unchanged from the underlying DTLZ / inverted
DTLZ formulations; only the constraints differ.

Reference
---------
Oliver Cuate, Lourdes Uribe, Adriana Lara, Oliver Schuetze, "Dataset on a
Benchmark for Equality Constrained Multi-objective Optimization", Data in
Brief, Volume 29, 2020, 105130.  https://doi.org/10.1016/j.dib.2020.105130
"""

import functools
import math
import operator
import random

from ..core import Problem, Solution
from ..types import Real

# Tolerance band that turns the analytical equality h(x) = 0 into a thin
# feasible shell, matching the original MATLAB benchmark (PopCon ... - 1e-4).
EPSILON = 1e-4


def _g_multimodal(distance_vars):
    """g function used by DTLZ1 / DTLZ3 (highly multimodal)."""
    k = len(distance_vars)
    return 100.0 * (k + sum(math.pow(x - 0.5, 2.0) -
                            math.cos(20.0 * math.pi * (x - 0.5))
                            for x in distance_vars))


def _g_quadratic(distance_vars):
    """g function used by DTLZ2 / DTLZ4 (single, smooth basin)."""
    return sum(math.pow(x - 0.5, 2.0) for x in distance_vars)


def _shape_dtlz1(position_vars, g, nobjs):
    """Linear (simplex) DTLZ1-style objective shape."""
    f = [0.5 * (1.0 + g)] * nobjs
    for i in range(nobjs):
        f[i] *= functools.reduce(operator.mul,
                                 position_vars[:nobjs - i - 1], 1)
        if i > 0:
            f[i] *= 1.0 - position_vars[nobjs - i - 1]
    return f


def _shape_dtlz2(position_vars, g, nobjs):
    """Spherical DTLZ2/3/4-style objective shape."""
    f = [1.0 + g] * nobjs
    for i in range(nobjs):
        f[i] *= functools.reduce(
            operator.mul,
            [math.cos(0.5 * math.pi * x)
             for x in position_vars[:nobjs - i - 1]], 1)
        if i > 0:
            f[i] *= math.sin(0.5 * math.pi * position_vars[nobjs - i - 1])
    return f


def _eq1_constraints(position_vars):
    """Single ring constraint, centre 0.5, radius 0.4 (Eq1 variants)."""
    r = 0.4
    xx = [x - 0.5 for x in position_vars]
    return [abs(sum(v * v for v in xx) - r * r) - EPSILON]


def _eq2_constraints(position_vars):
    """Two intersecting ring constraints, radius 0.5 (Eq2 variants)."""
    r = 0.5
    xx = [x - 0.5 for x in position_vars]
    yy = list(xx)
    yy[-1] -= r
    return [abs(sum(v * v for v in xx) - r * r) - EPSILON,
            abs(sum(v * v for v in yy) - r * r) - EPSILON]


class _EqDTLZ(Problem):
    """Common scaffolding for the equality-constrained DTLZ problems.

    Subclasses define :attr:`_default_k` (number of distance variables), the
    constraint pattern via :meth:`_constraints`, and the objective via
    :meth:`_objectives`.
    """

    #: number of distance ("g") variables when ``nvars`` is not supplied
    _default_k = 10
    #: number of equality constraints (1 for Eq1, 2 for Eq2)
    _nconstrs = 1

    def __init__(self, nobjs=3, nvars=None):
        if nvars is None:
            nvars = nobjs + self._default_k - 1
        super().__init__(nvars, nobjs, self._nconstrs)
        self.types[:] = Real(0, 1)
        self.constraints[:] = "<=0"

    def _split(self, solution):
        position = solution.variables[:self.nobjs - 1]
        distance = solution.variables[self.nobjs - 1:]
        return position, distance

    def _objectives(self, position, distance):
        raise NotImplementedError

    def _constraints(self, position):
        raise NotImplementedError

    def evaluate(self, solution):
        position, distance = self._split(solution)
        solution.objectives[:] = self._objectives(position, distance)
        solution.constraints[:] = self._constraints(position)

    def random(self):
        # Reference solution on the unconstrained DTLZ front (distance
        # variables at 0.5).  Note this point need not satisfy the equality
        # constraints; it mirrors the DTLZ family's convenience generator.
        solution = Solution(self)
        solution.variables[:self.nobjs - 1] = \
            [random.uniform(0.0, 1.0) for _ in range(self.nobjs - 1)]
        solution.variables[self.nobjs - 1:] = 0.5
        solution.evaluate()
        return solution


# ---------------------------------------------------------------------------
# Eq1 family: single ring constraint (r = 0.4)
# ---------------------------------------------------------------------------

class Eq1DTLZ1(_EqDTLZ):
    _default_k = 5
    _nconstrs = 1

    def _objectives(self, position, distance):
        return _shape_dtlz1(position, _g_multimodal(distance), self.nobjs)

    def _constraints(self, position):
        return _eq1_constraints(position)


class Eq1DTLZ2(_EqDTLZ):
    _default_k = 10
    _nconstrs = 1

    def _objectives(self, position, distance):
        return _shape_dtlz2(position, _g_quadratic(distance), self.nobjs)

    def _constraints(self, position):
        return _eq1_constraints(position)


class Eq1DTLZ3(_EqDTLZ):
    _default_k = 10
    _nconstrs = 1

    def _objectives(self, position, distance):
        return _shape_dtlz2(position, _g_multimodal(distance), self.nobjs)

    def _constraints(self, position):
        return _eq1_constraints(position)


class Eq1DTLZ4(_EqDTLZ):
    _default_k = 10
    _nconstrs = 1

    def __init__(self, nobjs=3, nvars=None, alpha=100.0):
        super().__init__(nobjs, nvars)
        self.alpha = alpha

    def _objectives(self, position, distance):
        biased = [math.pow(x, self.alpha) for x in position]
        return _shape_dtlz2(biased, _g_quadratic(distance), self.nobjs)

    def _constraints(self, position):
        return _eq1_constraints(position)


class Eq1IDTLZ1(_EqDTLZ):
    _default_k = 5
    _nconstrs = 1

    def _objectives(self, position, distance):
        g = _g_multimodal(distance)
        base = _shape_dtlz1(position, g, self.nobjs)
        return [(1.0 + g) / 2.0 - base[i] for i in range(self.nobjs)]

    def _constraints(self, position):
        return _eq1_constraints(position)


class Eq1IDTLZ2(_EqDTLZ):
    _default_k = 10
    _nconstrs = 1

    def _objectives(self, position, distance):
        g = _g_quadratic(distance)
        base = _shape_dtlz2(position, g, self.nobjs)
        return [(1.0 + g) - base[i] for i in range(self.nobjs)]

    def _constraints(self, position):
        return _eq1_constraints(position)


# ---------------------------------------------------------------------------
# Eq2 family: two intersecting ring constraints (r = 0.5)
# ---------------------------------------------------------------------------

class Eq2DTLZ1(_EqDTLZ):
    _default_k = 5
    _nconstrs = 2

    def _objectives(self, position, distance):
        return _shape_dtlz1(position, _g_multimodal(distance), self.nobjs)

    def _constraints(self, position):
        return _eq2_constraints(position)


class Eq2DTLZ2(_EqDTLZ):
    _default_k = 10
    _nconstrs = 2

    def _objectives(self, position, distance):
        return _shape_dtlz2(position, _g_quadratic(distance), self.nobjs)

    def _constraints(self, position):
        return _eq2_constraints(position)


class Eq2DTLZ3(_EqDTLZ):
    _default_k = 10
    _nconstrs = 2

    def _objectives(self, position, distance):
        return _shape_dtlz2(position, _g_multimodal(distance), self.nobjs)

    def _constraints(self, position):
        return _eq2_constraints(position)


class Eq2DTLZ4(_EqDTLZ):
    _default_k = 10
    _nconstrs = 2

    def __init__(self, nobjs=3, nvars=None, alpha=100.0):
        super().__init__(nobjs, nvars)
        self.alpha = alpha

    def _objectives(self, position, distance):
        biased = [math.pow(x, self.alpha) for x in position]
        return _shape_dtlz2(biased, _g_quadratic(distance), self.nobjs)

    def _constraints(self, position):
        return _eq2_constraints(position)


class Eq2IDTLZ1(_EqDTLZ):
    _default_k = 5
    _nconstrs = 2

    def _objectives(self, position, distance):
        g = _g_multimodal(distance)
        base = _shape_dtlz1(position, g, self.nobjs)
        return [(1.0 + g) / 2.0 - base[i] for i in range(self.nobjs)]

    def _constraints(self, position):
        return _eq2_constraints(position)


class Eq2IDTLZ2(_EqDTLZ):
    _default_k = 10
    _nconstrs = 2

    def _objectives(self, position, distance):
        g = _g_quadratic(distance)
        base = _shape_dtlz2(position, g, self.nobjs)
        return [(1.0 + g) - base[i] for i in range(self.nobjs)]

    def _constraints(self, position):
        return _eq2_constraints(position)
