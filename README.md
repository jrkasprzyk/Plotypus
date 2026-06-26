# Plotypus Multiobjective Optimization Library

Fork of [Platypus](https://github.com/Project-Platypus/Platypus) with a focus on visualization and interactive exploration.

### What is Plotypus?

Plotypus inherits rich evolutionary computing infrastructure from Platypus, which is a framework for evolutionary computing in Python with a focus on multiobjective evolutionary algorithms (MOEAs).

### What's different from the original Platypus library?

Plotypus is an experimental repo that adds functionality to the original Platypus including:

- **Interactive GUI** — Select an algorithm and test problem, set parameters, and visualize the Pareto front evolve.
- **Expanded problem library** — We refactored some of the examples within Platypus (the ZDT, WFG, UF, and CF problem families) to be accessible directly from the package alongside the existing DTLZ suite. We also added new families not in Platypus, such as the equality-constrained DTLZ benchmark (`Eq1DTLZ1`–`Eq2IDTLZ2`, Cuate et al. 2020).
- **Human-in-the-loop example** — `examples/madlibs_moea.py` is a simple proof-of-concept where a human is asked to evaluate the quality of mad-libs, and the MOEA uses this as an objective function.

### Installation

To install the latest development version of Plotypus, first install [Poetry](https://python-poetry.org/docs/#installation) if not already installed:

```
curl -sSL https://install.python-poetry.org | python3 -
```

Then clone and install Plotypus:

```
git clone https://github.com/jrkasprzyk/Plotypus.git
cd Plotypus
poetry install
```

### Examples

#### Running the GUI

```
poetry run plotypus-gui
```

Requires `matplotlib` and `tkinter` (tkinter ships with most Python distributions).

#### Define a problem and optimize with NSGAII

Optimizing a simple biobjective problem with a single real-valued decision variable:

```python
from plotypus import NSGAII, Problem, Real

def schaffer(x):
    return [x[0]**2, (x[0]-2)**2]

problem = Problem(1, 2)
problem.types[:] = Real(-10, 10)
problem.function = schaffer

algorithm = NSGAII(problem)
algorithm.run(10000)
```

### Citation

If you use this software in your work, please cite the original Platypus library as follows (APA style):

> Hadka, D. (2024). Platypus: A Framework for Evolutionary Computing in Python (Version 1.4.1) [Computer software].  Retrieved from https<span>://</span>github.com/Project-Platypus/Platypus.

### License

Plotypus is released under the GNU General Public License.
