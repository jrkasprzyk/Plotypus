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

import functools
import math
import operator
import random

from ..core import Problem, Solution
from ..types import Real


class DTLZ1(Problem):

    def __init__(self, nobjs=2):
        super().__init__(nobjs+4, nobjs)
        self.types[:] = Real(0, 1)

    def evaluate(self, solution):
        k = self.nvars - self.nobjs + 1
        g = 100.0 * (k + sum([math.pow(x - 0.5, 2.0) - math.cos(20.0 * math.pi * (x - 0.5)) for x in solution.variables[self.nvars-k:]]))
        f = [0.5 * (1.0 + g)]*self.nobjs

        for i in range(self.nobjs):
            f[i] *= functools.reduce(operator.mul,
                                     [x for x in solution.variables[:self.nobjs-i-1]],
                                     1)

            if i > 0:
                f[i] *= 1 - solution.variables[self.nobjs-i-1]

        solution.objectives[:] = f

    def random(self):
        solution = Solution(self)
        solution.variables[:self.nobjs-1] = [random.uniform(0.0, 1.0) for _ in range(self.nobjs-1)]
        solution.variables[self.nobjs-1:] = 0.5
        solution.evaluate()
        return solution

class DTLZ2(Problem):

    def __init__(self, nobjs=2, nvars=None):
        super().__init__(nobjs+9 if nvars is None else nvars, nobjs)
        self.types[:] = Real(0, 1)

    def evaluate(self, solution):
        k = self.nvars - self.nobjs + 1
        g = sum([math.pow(x - 0.5, 2.0) for x in solution.variables[self.nvars-k:]])
        f = [1.0+g]*self.nobjs

        for i in range(self.nobjs):
            f[i] *= functools.reduce(operator.mul,
                                     [math.cos(0.5 * math.pi * x) for x in solution.variables[:self.nobjs-i-1]],
                                     1)

            if i > 0:
                f[i] *= math.sin(0.5 * math.pi * solution.variables[self.nobjs-i-1])

        solution.objectives[:] = f

    def random(self):
        solution = Solution(self)
        solution.variables[:self.nobjs-1] = [random.uniform(0.0, 1.0) for _ in range(self.nobjs-1)]
        solution.variables[self.nobjs-1:] = 0.5
        solution.evaluate()
        return solution

class DTLZ3(Problem):

    def __init__(self, nobjs=2, nvars=None):
        super().__init__(nobjs+9 if nvars is None else nvars, nobjs)
        self.types[:] = Real(0, 1)

    def evaluate(self, solution):
        k = self.nvars - self.nobjs + 1
        g = 100.0 * (k + sum([math.pow(x - 0.5, 2.0) - math.cos(20.0 * math.pi * (x - 0.5)) for x in solution.variables[self.nvars-k:]]))
        f = [1.0+g]*self.nobjs

        for i in range(self.nobjs):
            f[i] *= functools.reduce(operator.mul,
                                     [math.cos(0.5 * math.pi * x) for x in solution.variables[:self.nobjs-i-1]],
                                     1)

            if i > 0:
                f[i] *= math.sin(0.5 * math.pi * solution.variables[self.nobjs-i-1])

        solution.objectives[:] = f

    def random(self):
        solution = Solution(self)
        solution.variables[:self.nobjs-1] = [random.uniform(0.0, 1.0) for _ in range(self.nobjs-1)]
        solution.variables[self.nobjs-1:] = 0.5
        solution.evaluate()
        return solution

class DTLZ4(Problem):

    def __init__(self, nobjs=2, alpha=100.0):
        super().__init__(nobjs+9, nobjs)
        self.types[:] = Real(0, 1)
        self.alpha = alpha

    def evaluate(self, solution):
        k = self.nvars - self.nobjs + 1
        g = sum([math.pow(x - 0.5, 2.0) for x in solution.variables[self.nvars-k:]])
        f = [1.0+g]*self.nobjs

        for i in range(self.nobjs):
            f[i] *= functools.reduce(operator.mul,
                                     [math.cos(0.5 * math.pi * math.pow(x, self.alpha)) for x in solution.variables[:self.nobjs-i-1]],
                                     1)

            if i > 0:
                f[i] *= math.sin(0.5 * math.pi * math.pow(solution.variables[self.nobjs-i-1], self.alpha))

        solution.objectives[:] = f

    def random(self):
        solution = Solution(self)
        solution.variables[:self.nobjs-1] = [random.uniform(0.0, 1.0) for _ in range(self.nobjs-1)]
        solution.variables[self.nobjs-1:] = 0.5
        solution.evaluate()
        return solution

class DTLZ7(Problem):

    def __init__(self, nobjs=2):
        super().__init__(nobjs+19, nobjs)
        self.types[:] = Real(0, 1)

    def evaluate(self, solution):
        k = self.nvars - self.nobjs + 1
        g = 1.0 + (9.0 * sum(solution.variables[self.nvars-k:])) / k
        h = self.nobjs - sum([x / (1.0 + g) * (1.0 + math.sin(3.0 * math.pi * x)) for x in solution.variables[:self.nobjs-1]])

        solution.objectives[:self.nobjs-1] = solution.variables[:self.nobjs-1]
        solution.objectives[-1] = (1.0 + g) * h

    def random(self):
        solution = Solution(self)
        solution.variables[:self.nobjs-1] = [random.uniform(0.0, 1.0) for _ in range(self.nobjs-1)]
        solution.variables[self.nobjs-1:] = 0.0
        solution.evaluate()
        return solution
