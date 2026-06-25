"""Problem and algorithm registries plus their factory functions."""

import math

from plotypus import (
    EpsMOEA, EpsNSGAII, GDE3, IBEA, MOEAD, NSGAII, NSGAIII, OMOPSO, SMPSO, SPEA2,
    DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ7,
    ZDT1, ZDT2, ZDT3, ZDT4, ZDT5, ZDT6,
    WFG1, WFG2, WFG3, WFG4, WFG5, WFG6, WFG7, WFG8, WFG9,
    UF1, UF2, UF3, UF4, UF5, UF6, UF7, UF8, UF9, UF10, UF11, UF12, UF13,
    CF1, CF2, CF3, CF4, CF5, CF6, CF7, CF8, CF9, CF10,
    Schaffer, Belegundu,
)

# (class, accepts nobjs kwarg)
PROBLEMS = {
    "DTLZ1": (DTLZ1, True),
    "DTLZ2": (DTLZ2, True),
    "DTLZ3": (DTLZ3, True),
    "DTLZ4": (DTLZ4, True),
    "DTLZ7": (DTLZ7, True),
    "ZDT1":  (ZDT1,  False),
    "ZDT2":  (ZDT2,  False),
    "ZDT3":  (ZDT3,  False),
    "ZDT4":  (ZDT4,  False),
    "ZDT5":  (ZDT5,  False),
    "ZDT6":  (ZDT6,  False),
    "WFG1":  (WFG1,  True),
    "WFG2":  (WFG2,  True),
    "WFG3":  (WFG3,  True),
    "WFG4":  (WFG4,  True),
    "WFG5":  (WFG5,  True),
    "WFG6":  (WFG6,  True),
    "WFG7":  (WFG7,  True),
    "WFG8":  (WFG8,  True),
    "WFG9":  (WFG9,  True),
    "UF1":   (UF1,   False),
    "UF2":   (UF2,   False),
    "UF3":   (UF3,   False),
    "UF4":   (UF4,   False),
    "UF5":   (UF5,   False),
    "UF6":   (UF6,   False),
    "UF7":   (UF7,   False),
    "UF8":   (UF8,   False),
    "UF9":   (UF9,   False),
    "UF10":  (UF10,  False),
    "UF11":  (UF11,  False),
    "UF12":  (UF12,  False),
    "UF13":  (UF13,  False),
    "CF1":   (CF1,   False),
    "CF2":   (CF2,   False),
    "CF3":   (CF3,   False),
    "CF4":   (CF4,   False),
    "CF5":   (CF5,   False),
    "CF6":   (CF6,   False),
    "CF7":   (CF7,   False),
    "CF8":   (CF8,   False),
    "CF9":   (CF9,   False),
    "CF10":  (CF10,  False),
    "Schaffer":   (Schaffer,   False),
    "Belegundu":  (Belegundu,  False),
}

# display name -> key used in make_algorithm()
ALGORITHMS = {
    "NSGA-II":   "nsgaii",
    "NSGA-III":  "nsgaiii",
    "MOEA/D":    "moead",
    "IBEA":      "ibea",
    "SPEA2":     "spea2",
    "GDE3":      "gde3",
    "OMOPSO":    "omopso",
    "SMPSO":     "smpso",
    "ε-MOEA":    "epsmoea",
    "ε-NSGA-II": "epsnsgaii",
}


def make_problem(name, nobjs):
    cls, accepts_nobjs = PROBLEMS[name]
    return cls(nobjs=nobjs) if accepts_nobjs else cls()


def _nsgaiii_divisions(nobjs, target_pop):
    div = 1
    while math.comb(nobjs + div - 1, div) < target_pop:
        div += 1
    return div


def make_algorithm(key, problem, pop_size):
    if key == "nsgaii":
        return NSGAII(problem, population_size=pop_size)
    if key == "nsgaiii":
        div = _nsgaiii_divisions(problem.nobjs, pop_size)
        return NSGAIII(problem, divisions_outer=div)
    if key == "moead":
        return MOEAD(problem, population_size=pop_size)
    if key == "ibea":
        return IBEA(problem, population_size=pop_size)
    if key == "spea2":
        return SPEA2(problem, population_size=pop_size)
    if key == "gde3":
        return GDE3(problem, population_size=pop_size)
    if key == "omopso":
        return OMOPSO(problem, epsilons=[0.05], swarm_size=pop_size)
    if key == "smpso":
        return SMPSO(problem, swarm_size=pop_size)
    if key == "epsmoea":
        return EpsMOEA(problem, epsilons=[0.05], population_size=pop_size)
    if key == "epsnsgaii":
        return EpsNSGAII(problem, epsilons=[0.05], population_size=pop_size)
    raise ValueError(f"Unknown algorithm: {key}")
