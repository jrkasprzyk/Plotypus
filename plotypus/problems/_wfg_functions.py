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

from .._math import EPSILON


def _normalize_z(z):
    return [z[i] / (2.0 * (i+1)) for i in range(len(z))]

def _correct_to_01(a):
    if a <= 0.0 and a >= -EPSILON:
        return 0.0
    elif a >= 1.0 and a <= 1.0 + EPSILON:
        return 1.0
    else:
        return a

def _vector_in_01(x):
    return all([a >= 0.0 and a <= 1.0 for a in x])

def _s_linear(y, A):
    return _correct_to_01(abs(y - A) / abs(math.floor(A - y) + A))

def _s_multi(y, A, B, C):
    tmp1 = abs(y - C) / (2.0 * (math.floor(C - y) + C))
    tmp2 = (4.0 * A + 2.0) * math.pi * (0.5 - tmp1)
    return _correct_to_01((1.0 + math.cos(tmp2) + 4.0 * B * math.pow(tmp1, 2.0)) / (B + 2.0))

def _s_decept(y, A, B, C):
    tmp1 = math.floor(y - A + B) * (1.0 - C + (A - B) / B) / (A - B)
    tmp2 = math.floor(A + B - y) * (1.0 - C + (1.0 - A - B) / B) / (1.0 - A - B)
    return _correct_to_01(1.0 + (abs(y - A) - B) * (tmp1 + tmp2 + 1.0 / B))

def _b_flat(y, A, B, C):
    return _correct_to_01(A +
                          min(0.0, math.floor(y - B)) * A * (B - y) / B -
                          min(0.0, math.floor(C - y)) * (1.0 - A) * (y - C))

def _b_poly(y, alpha):
    return _correct_to_01(math.pow(y, alpha))

def _b_param(y, u, A, B, C):
    return _correct_to_01(math.pow(y, B + (C-B) * (A - (1.0 - 2.0*u) * abs(math.floor(0.5 - u) + A))))

def _subvector(v, head, tail):
    return [v[i] for i in range(head, tail)]

def _r_sum(y, w):
    numerator = sum([w[i]*y[i] for i in range(len(y))])
    denominator = sum([w[i] for i in range(len(y))])
    return _correct_to_01(numerator / denominator)

def _r_nonsep(y, A):
    numerator = sum([y[j] + sum([abs(y[j] - y[(j + k + 1) % len(y)]) for k in range(A-1)]) for j in range(len(y))])
    tmp = math.ceil(A / 2.0)
    denominator = len(y) * tmp * (1.0 + 2.0*A - 2.0*tmp) / A
    return _correct_to_01(numerator / denominator)

def _WFG1_t1(y, k):
    return y[:k] + list(map(functools.partial(_s_linear, A=0.35), y[k:]))

def _WFG1_t2(y, k):
    return y[:k] + list(map(functools.partial(_b_flat, A=0.8, B=0.75, C=0.85), y[k:]))

def _WFG1_t3(y):
    return list(map(functools.partial(_b_poly, alpha=0.02), y))

def _WFG1_t4(y, k, M):
    w = [2.0*(i+1) for i in range(len(y))]
    t = []

    for i in range(M-1):
        head = i * k // (M-1)
        tail = (i+1) * k // (M-1)
        y_sub = _subvector(y, head, tail)
        w_sub = _subvector(w, head, tail)
        t.append(_r_sum(y_sub, w_sub))

    y_sub = _subvector(y, k, len(y))
    w_sub = _subvector(w, k, len(y))
    t.append(_r_sum(y_sub, w_sub))
    return t

def _WFG2_t2(y, k):
    l = len(y) - k
    t = y[:k]

    for i in range(k+1, k+(l//2)+1):
        head = k + 2 * (i - k) - 2
        tail = k + 2 * (i - k)
        t.append(_r_nonsep(_subvector(y, head, tail), 2))

    return t

def _WFG2_t3(y, k, M):
    w = [1.0]*len(y)
    t = []

    for i in range(M-1):
        head = i * k // (M-1)
        tail = (i+1) * k // (M-1)
        y_sub = _subvector(y, head, tail)
        w_sub = _subvector(w, head, tail)
        t.append(_r_sum(y_sub, w_sub))

    y_sub = _subvector(y, k, len(y))
    w_sub = _subvector(w, k, len(y))
    t.append(_r_sum(y_sub, w_sub))
    return t

def _WFG4_t1(y):
    return list(map(functools.partial(_s_multi, A=30, B=10, C=0.35), y))

def _WFG5_t1(y):
    return list(map(functools.partial(_s_decept, A=0.35, B=0.001, C=0.05), y))

def _WFG6_t2(y, k, M):
    t = []

    for i in range(M-1):
        head = i * k // (M-1)
        tail = (i+1) * k // (M-1)
        y_sub = _subvector(y, head, tail)
        t.append(_r_nonsep(y_sub, k // (M-1)))

    y_sub = _subvector(y, k, len(y))
    t.append(_r_nonsep(y_sub, len(y)-k))
    return t

def _WFG7_t1(y, k):
    w = [1.0]*len(y)
    t = []

    for i in range(k):
        y_sub = _subvector(y, i+1, len(y))
        w_sub = _subvector(w, i+1, len(y))
        u = _r_sum(y_sub, w_sub)
        t.append(_b_param(y[i], u, 0.98 / 49.98, 0.02, 50))

    for i in range(k, len(y)):
        t.append(y[i])

    return t

def _WFG8_t1(y, k):
    w = [1.0]*len(y)
    t = y[:k]

    for i in range(k, len(y)):
        y_sub = _subvector(y, 0, i)
        w_sub = _subvector(w, 0, i)
        u = _r_sum(y_sub, w_sub)
        t.append(_b_param(y[i], u, 0.98 / 49.98, 0.02, 50))

    return t

def _WFG9_t1(y):
    w = [1.0]*len(y)
    t = []

    for i in range(len(y)-1):
        y_sub = _subvector(y, i + 1, len(y))
        w_sub = _subvector(w, i + 1, len(y))
        u = _r_sum(y_sub, w_sub)
        t.append(_b_param(y[i], u, 0.98 / 49.98, 0.02, 50))

    t.append(y[-1])
    return t

def _WFG9_t2(y, k):
    return list(map(functools.partial(_s_decept, A=0.35, B=0.001, C=0.05), y[:k])) + list(map(functools.partial(_s_multi, A=30, B=95, C=0.35), y[k:]))

def _create_A(M, degenerate):
    if degenerate:
        return [1.0 if i == 0 else 0.0 for i in range(M-1)]
    else:
        return [1.0]*(M-1)

def _calculate_x(t_p, A):
    return [max(t_p[-1], A[i]) * (t_p[i] - 0.5) + 0.5 for i in range(len(t_p)-1)] + [t_p[-1]]

def _convex(x, m):
    result = functools.reduce(operator.mul,
                              [1.0 - math.cos(x[i-1] * math.pi / 2.0) for i in range(1, len(x)-m+1)],
                              1.0)

    if m != 1:
        result *= 1.0 - math.sin(x[len(x)-m] * math.pi / 2.0)

    return _correct_to_01(result)

def _concave(x, m):
    result = functools.reduce(operator.mul,
                              [math.sin(x[i-1] * math.pi / 2.0) for i in range(1, len(x)-m+1)],
                              1.0)

    if m != 1:
        result *= math.cos(x[len(x)-m] * math.pi / 2.0)

    return _correct_to_01(result)

def _linear(x, m):
    result = functools.reduce(operator.mul, x[:len(x)-m], 1.0)

    if m != 1:
        result *= 1.0 - x[len(x)-m]

    return _correct_to_01(result)

def _mixed(x, A, alpha):
    tmp = 2.0 * A * math.pi
    return _correct_to_01(math.pow(1.0 - x[0] - math.cos(tmp * x[0] + math.pi / 2.0) / tmp, alpha))

def _disc(x, A, alpha, beta):
    tmp = A * math.pow(x[0], beta) * math.pi
    return _correct_to_01(1.0 - math.pow(x[0], alpha) * math.pow(math.cos(tmp), 2.0))

def _calculate_f(D, x, h, S):
    return [D * x[-1] + S[i]*h[i] for i in range(len(h))]

def _WFG_calculate_f(x, h):
    S = [m * 2.0 for m in range(1, len(h)+1)]
    return _calculate_f(1.0, x, h, S)

def _WFG1_shape(t_p):
    A = _create_A(len(t_p), False)
    x = _calculate_x(t_p, A)
    h = [_convex(x, m) for m in range(1, len(t_p))] + [_mixed(x, 5, 1.0)]
    return _WFG_calculate_f(x, h)

def _WFG2_shape(t_p):
    A = _create_A(len(t_p), False)
    x = _calculate_x(t_p, A)
    h = [_convex(x, m) for m in range(1, len(t_p))] + [_disc(x, 5, 1.0, 1.0)]
    return _WFG_calculate_f(x, h)

def _WFG3_shape(t_p):
    A = _create_A(len(t_p), True)
    x = _calculate_x(t_p, A)
    h = [_linear(x, m) for m in range(1, len(t_p)+1)]
    return _WFG_calculate_f(x, h)

def _WFG4_shape(t_p):
    A = _create_A(len(t_p), False)
    x = _calculate_x(t_p, A)
    h = [_concave(x, m) for m in range(1, len(t_p)+1)]
    return _WFG_calculate_f(x, h)
