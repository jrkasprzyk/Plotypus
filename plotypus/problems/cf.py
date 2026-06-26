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

from ..core import Problem
from ..types import Real


class CF1(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 1)
        self.types[:] = Real(0, 1)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        sum1 = 0.0
        sum2 = 0.0
        N = 10.0
        a = 1.0

        for j in range(2, self.nvars+1):
            yj = x[j-1] - math.pow(x[0], 0.5 * (1.0 + 3.0 * (j - 2.0) / (self.nvars - 2.0)))

            if j % 2 == 1:
                sum1 += yj * yj
                count1 += 1
            else:
                sum2 += yj * yj
                count2 += 1

        f1 = x[0] + 2.0 * sum1 / count1
        f2 = 1.0 - x[0] + 2.0 * sum2 / count2
        solution.objectives[:] = [f1, f2]
        solution.constraints[0] = f1 + f2 - a * abs(math.sin(N * math.pi * (f1 - f2 + 1.0))) - 1.0

class CF2(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 1)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-1, 1)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        sum1 = 0.0
        sum2 = 0.0
        N = 2.0
        a = 1.0

        for j in range(2, self.nvars+1):
            if j % 2 == 1:
                yj = x[j-1] - math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum1 += yj * yj
                count1 += 1
            else:
                yj = x[j-1] - math.cos(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum2 += yj * yj
                count2 += 1

        f1 = x[0] + 2.0 * sum1 / count1
        f2 = 1.0 - math.sqrt(x[0]) + 2.0 * sum2 / count2
        t = f2 + math.sqrt(f1) - a * math.sin(N * math.pi * (math.sqrt(f1) - f2 + 1.0)) - 1.0
        solution.objectives[:] = [f1, f2]
        solution.constraints[0] = (1 if t >= 0 else -1) * abs(t) / (1.0 + math.exp(4.0 * abs(t)))

class CF3(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 1)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        sum1 = 0.0
        sum2 = 0.0
        prod1 = 1.0
        prod2 = 1.0
        N = 2.0
        a = 1.0

        for j in range(2, self.nvars+1):
            yj = x[j-1] - math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)
            pj = math.cos(20.0 * yj * math.pi / math.sqrt(j))

            if j % 2 == 1:
                sum1 += yj * yj
                prod1 *= pj
                count1 += 1
            else:
                sum2 += yj * yj
                prod2 *= pj
                count2 += 1

        f1 = x[0] + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / count1
        f2 = 1.0 - x[0]**2 + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / count2
        solution.objectives[:] = [f1, f2]
        solution.constraints[0] = f2 + f1**2 - a * math.sin(N * math.pi * (f1**2 - f2 + 1.0)) - 1.0

class CF4(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 1)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        sum1 = 0.0
        sum2 = 0.0

        for j in range(2, self.nvars+1):
            yj = x[j-1] - math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)

            if j % 2 == 1:
                sum1 += yj**2
            else:
                if j == 2:
                    sum2 += abs(yj) if yj < 1.5 - 0.75 * math.sqrt(2.0) else 0.125 + (yj - 1)**2
                else:
                    sum2 += yj**2

        f1 = x[0] + sum1
        f2 = 1.0 - x[0] + sum2
        t = x[1] - math.sin(6.0*x[0]*math.pi + 2.0*math.pi/self.nvars) - 0.5*x[0] + 0.25
        solution.objectives[:] = [f1, f2]
        solution.constraints[0] = (1 if t >= 0 else -1) * abs(t) / (1.0 + math.exp(4.0 * abs(t)))

class CF5(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 1)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        sum1 = 0.0
        sum2 = 0.0

        for j in range(2, self.nvars+1):
            if j % 2 == 1:
                yj = x[j-1] - 0.8*x[0]*math.cos(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum1 += 2.0*yj**2 - math.cos(4.0*math.pi*yj) + 1.0
            else:
                yj = x[j-1] - 0.8*x[0]*math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)

                if j == 2:
                    sum2 += abs(yj) if yj < 1.5 - 0.75*math.sqrt(2.0) else 0.125 + (yj - 1)**2
                else:
                    sum2 += 2.0*yj**2 - math.cos(4.0*math.pi*yj) + 1.0

        f1 = x[0] + sum1
        f2 = 1.0 - x[0] + sum2
        solution.objectives[:] = [f1, f2]
        solution.constraints[0] = x[1] - 0.8*x[0]*math.sin(6.0*x[0]*math.pi + 2.0*math.pi/self.nvars) - 0.5*x[0] + 0.25

class CF6(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 2)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        sum1 = 0.0
        sum2 = 0.0

        for j in range(2, self.nvars+1):
            if j % 2 == 1:
                yj = x[j-1] - 0.8*x[0]*math.cos(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum1 += yj**2
            else:
                yj = x[j-1] - 0.8*x[0]*math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum2 += yj**2

        f1 = x[0] + sum1
        f2 = (1.0 - x[0])**2 + sum2
        c1 = x[1] - 0.8*x[0]*math.sin(6.0*x[0]*math.pi + 2.0*math.pi/self.nvars) - \
            (1 if (x[0]-0.5)*(1.0-x[0]) >= 0 else -1) * math.sqrt(abs((x[0]-0.5)*(1.0-x[0])))
        c2 = x[3] - 0.8*x[0]*math.sin(6.0*x[0]*math.pi + 4.0*math.pi/self.nvars) - \
            (1 if 0.25 * math.sqrt(1.0-x[0]) - 0.5*(1.0-x[0]) >= 0 else -1)*math.sqrt(abs(0.25 * math.sqrt(1-x[0]) - 0.5*(1.0-x[0])))
        solution.objectives[:] = [f1, f2]
        solution.constraints[:] = [c1, c2]

class CF7(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 2, 2)
        self.types[0] = Real(0, 1)
        self.types[1:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        sum1 = 0.0
        sum2 = 0.0

        for j in range(2, self.nvars+1):
            if j % 2 == 1:
                yj = x[j-1] - math.cos(6.0*math.pi*x[0] + j*math.pi/self.nvars)
                sum1 += 2.0*yj**2 - math.cos(4.0*math.pi*yj) + 1.0
            else:
                yj = x[j-1] - math.sin(6.0*math.pi*x[0] + j*math.pi/self.nvars)

                if j == 2 or j == 4:
                    sum2 += yj**2
                else:
                    sum2 += 2.0*yj**2 - math.cos(4.0*math.pi*yj) + 1.0

        f1 = x[0] + sum1
        f2 = (1.0 - x[0])**2 + sum2
        c1 = x[1] - math.sin(6.0*x[0]*math.pi + 2.0*math.pi/self.nvars) - \
            (1 if (x[0]-0.5)*(1.0-x[0]) >= 0 else -1) * math.sqrt(abs((x[0]-0.5)*(1.0-x[0])))
        c2 = x[3] - math.sin(6.0*x[0]*math.pi + 4.0*math.pi/self.nvars) - \
            (1 if 0.25 * math.sqrt(1.0-x[0]) - 0.5*(1.0-x[0]) >= 0 else -1)*math.sqrt(abs(0.25 * math.sqrt(1-x[0]) - 0.5*(1.0-x[0])))
        solution.objectives[:] = [f1, f2]
        solution.constraints[:] = [c1, c2]

class CF8(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 3, 1)
        self.types[0:2] = Real(0, 1)
        self.types[2:] = Real(-4, 4)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        count3 = 0
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        N = 2.0
        a = 4.0

        for j in range(3, self.nvars+1):
            yj = x[j-1] - 2.0*x[1]*math.sin(2.0*math.pi*x[0] + j*math.pi/self.nvars)

            if j % 3 == 1:
                sum1 += yj**2
                count1 += 1
            elif j % 3 == 2:
                sum2 += yj**2
                count2 += 1
            else:
                sum3 += yj**2
                count3 += 1

        f1 = math.cos(0.5*math.pi*x[0]) * math.cos(0.5*math.pi*x[1]) + 2.0*sum1/count1
        f2 = math.cos(0.5*math.pi*x[0]) * math.sin(0.5*math.pi*x[1]) + 2.0*sum2/count2
        f3 = math.sin(0.5*math.pi*x[0]) + 2.0*sum3/count3
        c1 = (f1**2 + f2**2) / (1.0 - f3**2) - a*abs(math.sin(N*math.pi*((f1**2 - f2**2) / (1.0 - f3**2) + 1.0))) - 1.0
        solution.objectives[:] = [f1, f2, f3]
        solution.constraints[:] = [c1]

class CF9(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 3, 1)
        self.types[0:2] = Real(0, 1)
        self.types[2:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        count3 = 0
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        N = 2.0
        a = 3.0

        for j in range(3, self.nvars+1):
            yj = x[j-1] - 2.0*x[1]*math.sin(2.0*math.pi*x[0] + j*math.pi/self.nvars)

            if j % 3 == 1:
                sum1 += yj**2
                count1 += 1
            elif j % 3 == 2:
                sum2 += yj**2
                count2 += 1
            else:
                sum3 += yj**2
                count3 += 1

        f1 = math.cos(0.5*math.pi*x[0]) * math.cos(0.5*math.pi*x[1]) + 2.0*sum1/count1
        f2 = math.cos(0.5*math.pi*x[0]) * math.sin(0.5*math.pi*x[1]) + 2.0*sum2/count2
        f3 = math.sin(0.5*math.pi*x[0]) + 2.0*sum3/count3
        c1 = (f1**2 + f2**2) / (1.0 - f3**2) - a*math.sin(N*math.pi*((f1**2 - f2**2) / (1.0 - f3**2) + 1.0)) - 1.0
        solution.objectives[:] = [f1, f2, f3]
        solution.constraints[:] = [c1]

class CF10(Problem):

    def __init__(self, nvars=10):
        super().__init__(nvars, 3, 1)
        self.types[0:2] = Real(0, 1)
        self.types[2:] = Real(-2, 2)
        self.constraints[:] = ">=0"

    def evaluate(self, solution):
        x = solution.variables[:]
        count1 = 0
        count2 = 0
        count3 = 0
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        N = 2.0
        a = 1.0

        for j in range(3, self.nvars+1):
            yj = x[j-1] - 2.0*x[1]*math.sin(2.0*math.pi*x[0] + j*math.pi/self.nvars)
            hj = 4.0*yj**2 - math.cos(8.0*math.pi*yj) + 1.0

            if j % 3 == 1:
                sum1 += hj
                count1 += 1
            elif j % 3 == 2:
                sum2 += hj
                count2 += 1
            else:
                sum3 += hj
                count3 += 1

        f1 = math.cos(0.5*math.pi*x[0]) * math.cos(0.5*math.pi*x[1]) + 2.0*sum1/count1
        f2 = math.cos(0.5*math.pi*x[0]) * math.sin(0.5*math.pi*x[1]) + 2.0*sum2/count2
        f3 = math.sin(0.5*math.pi*x[0]) + 2.0*sum3/count3
        c1 = (f1**2 + f2**2) / (1.0 - f3**2) - a*math.sin(N*math.pi*((f1**2 - f2**2) / (1.0 - f3**2) + 1.0)) - 1.0
        solution.objectives[:] = [f1, f2, f3]
        solution.constraints[:] = [c1]
