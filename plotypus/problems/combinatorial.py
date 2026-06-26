# Copyright 2015-2024 David Hadka
#
# This file is part of Platypus, a Python module for designing and using
# evolutionary algorithms (EAs) and multiobjective evolutionary algorithms
# (MOEAs).
#
# Platypus is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Platypus is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with plotypus.  If not, see <http://www.gnu.org/licenses/>.

import math

from ..core import Constraint, Direction, Problem
from ..types import Binary, Permutation, Real


class Knapsack(Problem):

    def __init__(self, weights, profits, capacity):
        super().__init__(1, 1, 1)
        self.types[0] = Binary(len(weights))
        self.directions[0] = Direction.MAXIMIZE
        self.constraints[0] = Constraint("<=", capacity)
        self.weights = weights
        self.profits = profits

    def evaluate(self, solution):
        selection = solution.variables[0]
        n = len(self.weights)
        solution.objectives[0] = sum(self.profits[i] for i in range(n) if selection[i])
        solution.constraints[0] = sum(self.weights[i] for i in range(n) if selection[i])


class TSP(Problem):

    def __init__(self, cities):
        super().__init__(1, 1)
        self.types[0] = Permutation(range(len(cities)))
        self.directions[0] = Direction.MINIMIZE
        self.cities = cities

    def evaluate(self, solution):
        tour = solution.variables[0]
        n = len(self.cities)
        solution.objectives[0] = sum(
            self._dist(self.cities[tour[i]], self.cities[tour[(i + 1) % n]])
            for i in range(n)
        )

    @staticmethod
    def _dist(a, b):
        return round(math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2))


class PortfolioOptimization(Problem):
    """Minimize portfolio variance and concentration; maximize expected return."""

    def __init__(self, assets, returns, vols, corr):
        n = len(assets)
        super().__init__(n, 3)
        self.types[:] = Real(0, 1)
        self.directions[1] = Direction.MAXIMIZE
        self.assets = assets
        self.returns = returns
        self.cov = [[corr[i][j] * vols[i] * vols[j] for j in range(n)] for i in range(n)]

    def evaluate(self, solution):
        raw = solution.variables[:]
        total = sum(raw)
        n = self.nvars
        w = [r / total for r in raw] if total > 1e-9 else [1/n] * n
        variance = sum(w[i] * w[j] * self.cov[i][j] for i in range(n) for j in range(n))
        ret = sum(w[i] * self.returns[i] for i in range(n))
        hhi = sum(wi**2 for wi in w)
        solution.objectives[:] = [variance, ret, hhi]
