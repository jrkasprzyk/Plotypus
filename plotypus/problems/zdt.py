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
from abc import ABCMeta

from ..core import Problem
from ..types import Binary, Real


class ZDT(Problem, metaclass=ABCMeta):

    def __init__(self, nvars):
        super().__init__(nvars, 2)
        self.types[:] = Real(0, 1)

class ZDT1(ZDT):

    def __init__(self):
        super().__init__(30)

    def evaluate(self, solution):
        x = solution.variables[:]
        g = (9.0 / (self.nvars - 1.0))*sum(x[1:]) + 1.0
        h = 1.0 - math.sqrt(x[0] / g)
        solution.objectives[:] = [x[0], g*h]

class ZDT2(ZDT):

    def __init__(self):
        super().__init__(30)

    def evaluate(self, solution):
        x = solution.variables[:]
        g = (9.0 / (self.nvars - 1.0))*sum(x[1:]) + 1.0
        h = 1.0 - math.pow(x[0] / g, 2.0)
        solution.objectives[:] = [x[0], g*h]

class ZDT3(ZDT):

    def __init__(self):
        super().__init__(30)

    def evaluate(self, solution):
        x = solution.variables[:]
        g = (9.0 / (self.nvars - 1.0))*sum(x[1:]) + 1.0
        h = 1.0 - math.sqrt(x[0]/g) - (x[0]/g)*math.sin(10.0*math.pi*x[0])
        solution.objectives[:] = [x[0], g*h]

class ZDT4(ZDT):

    def __init__(self):
        super().__init__(10)

    def evaluate(self, solution):
        x = solution.variables[:]
        g = 1.0 + 10.0*(self.nvars-1) + sum([math.pow(x[i], 2.0) - 10.0*math.cos(4.0*math.pi*x[i]) for i in range(1, self.nvars)])
        h = 1.0 - math.sqrt(x[0] / g)
        solution.objectives[:] = [x[0], g*h]

class ZDT5(ZDT):

    def __init__(self):
        super().__init__(11)
        self.types[0] = Binary(30)
        self.types[1:] = Binary(5)

    def evaluate(self, solution):
        f = 1.0 + sum(solution.variables[0])
        g = sum([2+sum(v) if sum(v) < 5 else 1 for v in solution.variables[1:]])
        h = 1.0 / f
        solution.objectives[:] = [f, g*h]

class ZDT6(ZDT):

    def __init__(self):
        super().__init__(10)

    def evaluate(self, solution):
        x = solution.variables[:]
        f = 1.0 - math.exp(-4.0*x[0])*math.pow(math.sin(6.0*math.pi*x[0]), 6.0)
        g = 1.0 + 9.0*math.pow(sum(x[1:]) / (self.nvars-1.0), 0.25)
        h = 1.0 - math.pow(f / g, 2.0)
        solution.objectives[:] = [f, g*h]
