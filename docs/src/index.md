# StructureSolver.jl

StructureSolver.jl is a Julia package for solving 1D (spherically symmetric) relativistic stellar structure models in Damour–Esposito–Farese (DEF) scalar-tensor gravity, with helper tools for neutron-star equations of state (EoS).

It is aimed at research workflows: single solves, parameter scans, families of models, and sensitivity/derivative extraction.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/AlexBatrakov/StructureSolver.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path="/path/to/StructureSolver.jl")
```

## What to read next

- [Quickstart](quickstart.md): minimal examples for EoS and a single solve.
- [Tutorial: M–R curve](tutorial_mr.md): end-to-end family scan over central pressure.
- [Notation & outputs](notation.md): what the main `model.quantities` and `model.derivatives` keys mean.
- [EoS](eos.md): equation-of-state interfaces and unit conventions.
- [Simulations](simulations.md): regimes and simulation drivers.
- [Examples](examples.md): copy-paste snippets.
- [Troubleshooting](troubleshooting.md): common issues (PyPlot, docs deploy, etc.).
- [API](api.md): auto-generated API reference from docstrings.
