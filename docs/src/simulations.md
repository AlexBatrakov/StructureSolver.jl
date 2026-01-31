# Simulations

StructureSolver.jl separates *physical models* from *regimes* (how parameters are fixed or shot) and *simulation drivers* (single run, families, grids).

## Models

The main model provided by the package is `DEFp_Model`, which combines

- a coupling function (e.g. `DEF1_CouplingFunction`, `DEF3_CouplingFunction`), and
- an EoS (e.g. `PWP_EoS`, `Table_EoS`).

## Regimes

- `Simple_DirectRegime`: direct (non-shooting) regime where you fix inner/external parameters.
- `Simple_ShootingRegime`: single-model shooting regime (adjusts one or more inner parameters).
- `ShootingRegime`: general N-D shooting regime used by `GeneralSimulation`.

## Drivers

- `SingleSimulation`: run one model.
- `FamilySimulation`: scan one parameter (commonly central pressure) and store arrays of quantities/derivatives.
- `GeneralSimulation`: multi-dimensional scans with optional shooting.

## Optional plotting

The package does not depend on `PyPlot` by default. Some helpers (e.g. `plot_radial_structure`) import `PyPlot` only when called.
