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
import random
from abc import ABCMeta

from ..core import Problem, Solution
from ..types import Real
from ._wfg_functions import (
    _normalize_z,
    _r_sum,
    _WFG1_t1,
    _WFG1_t2,
    _WFG1_t3,
    _WFG1_t4,
    _WFG1_shape,
    _WFG2_t2,
    _WFG2_t3,
    _WFG2_shape,
    _WFG3_shape,
    _WFG4_t1,
    _WFG4_shape,
    _WFG5_t1,
    _WFG6_t2,
    _WFG7_t1,
    _WFG8_t1,
    _WFG9_t1,
    _WFG9_t2,
)


class WFG(Problem, metaclass=ABCMeta):

    def __init__(self, k, l, m):
        super().__init__(k+l, m)
        self.k = k
        self.l = l
        self.m = m
        self.types[:] = [Real(0.0, 2.0*(i+1)) for i in range(k+l)]

class WFG1(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG1_t1(y, self.k)
        y = _WFG1_t2(y, self.k)
        y = _WFG1_t3(y)
        y = _WFG1_t4(y, self.k, self.m)
        y = _WFG1_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [math.pow(random.uniform(0.0, 1.0), 50.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG2(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG1_t1(y, self.k)
        y = _WFG2_t2(y, self.k)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG2_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG3(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG1_t1(y, self.k)
        y = _WFG2_t2(y, self.k)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG3_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG4(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG4_t1(y)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG5(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG5_t1(y)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG6(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG1_t1(y, self.k)
        y = _WFG6_t2(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG7(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG7_t1(y, self.k)
        y = _WFG1_t1(y, self.k)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        solution = Solution(self)
        solution.variables[:self.k] = [random.uniform(0.0, 1.0) for _ in range(self.k)]
        solution.variables[self.k:] = 0.35
        solution.variables[:] = [solution.variables[i] * 2.0 * (i+1) for i in range(self.nvars)]
        self.evaluate(solution)
        return solution

class WFG8(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG8_t1(y, self.k)
        y = _WFG1_t1(y, self.k)
        y = _WFG2_t3(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        result = [random.uniform(0.0, 1.0) for _ in range(self.k)] + [0.0]*self.l

        for i in range(self.k, self.nvars):
            w = [1.0]*(self.nvars)
            u = _r_sum(result, w)
            tmp1 = abs(math.floor(0.5 - u) + 0.98 / 49.98)
            tmp2 = 0.02 + 49.98 * (0.98 / 49.98 - (1.0 - 2.0*u) * tmp1)
            result[i] = math.pow(0.35, math.pow(tmp2, -1.0))

        result = [result[i] * 2.0 * (i+1) for i in range(self.nvars)]

        solution = Solution(self)
        solution.variables[:] = result
        self.evaluate(solution)
        return solution

class WFG9(WFG):

    def __init__(self, nobjs=2):
        super().__init__(nobjs-1, 10, nobjs)

    def evaluate(self, solution):
        y = _normalize_z(solution.variables[:])
        y = _WFG9_t1(y)
        y = _WFG9_t2(y, self.k)
        y = _WFG6_t2(y, self.k, self.m)
        y = _WFG4_shape(y)
        solution.objectives[:] = y

    def random(self):
        result = [random.uniform(0.0, 1.0) for _ in range(self.k)] + [0.0]*(self.l-1) + [0.35]

        for i in range(self.nvars-2, self.k-1, -1):
            result_sub = result[i+1:self.nvars]
            w = [1.0]*len(result_sub)
            tmp1 = _r_sum(result_sub, w)
            result[i] = math.pow(0.35, math.pow(0.02 + 1.96 * tmp1, -1.0))

        result = [result[i] * 2.0 * (i+1) for i in range(self.nvars)]

        solution = Solution(self)
        solution.variables[:] = result
        self.evaluate(solution)
        return solution
