# StructureSolver.jl

[![CI](https://github.com/AlexBatrakov/StructureSolver.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/AlexBatrakov/StructureSolver.jl/actions/workflows/ci.yml)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlexBatrakov.github.io/StructureSolver.jl/dev)
[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlexBatrakov.github.io/StructureSolver.jl/stable)

StructureSolver is a Julia package for solving 1D (spherically symmetric) relativistic stellar structure models in Damour–Esposito–Farese (DEF) scalar-tensor gravity, with practical tools for neutron-star equations of state (EoS).

It targets research and engineering-style workflows: single solves, parameter scans (families), N-D grids, and extraction of derived quantities and sensitivities.

## Highlights

- **Physics scope**: relativistic stellar structure in DEF scalar–tensor gravity (TOV-like, with slow-rotation quantities).
- **EoS toolbox**: piecewise polytropes and tabulated EoS + unit-consistent conversion helpers.
- **Workflows**: single runs, one-parameter families (e.g. central pressure scans), and N-D grids with shooting.
- **Outputs**: `model.quantities` (mass, radius, scalarization, etc.) and `model.derivatives` (response coefficients).

## Quick links

- Docs (dev): https://AlexBatrakov.github.io/StructureSolver.jl/dev
- Tutorial (M–R curve): https://AlexBatrakov.github.io/StructureSolver.jl/dev/tutorial_mr/
- Notation & outputs: https://AlexBatrakov.github.io/StructureSolver.jl/dev/notation/
- API reference: https://AlexBatrakov.github.io/StructureSolver.jl/dev/api/

## Quality signals

- CI runs tests on Julia **1.10** and **1.11**.
- Docs build is automated via GitHub Actions (Documenter).
- Coverage is generated in CI and can be uploaded to Codecov (non-blocking; see `RELEASING.md`).
- Plotting is **optional** (no hard `PyPlot`/`matplotlib` dependency at package load).

## What’s inside (at a glance)

- **EoS implementations**: `PWP_EoS`, `Polytropic_EoS`, `Table_EoS`
- **EoS utilities**: `get_pressure`, `get_density`, `get_number_density`, `get_energy_density`, `get_specific_enthalpy`, `get_sound_velocity`, `get_internal_energy`, `find_max_pressure`
- **DEF scalar-tensor model**: `DEFp_Model`, coupling functions `DEF1_CouplingFunction`, `DEF3_CouplingFunction`
- **Simulation drivers**: `SingleSimulation`, `FamilySimulation`, `GeneralSimulation`
- **Regimes**: `Simple_DirectRegime`, `Simple_ShootingRegime`, `ShootingRegime`
- **Utilities/constants**: `LogRange`, `c`, `G`, `mp`, `M_sun`

## Installation

If the package is not in your registries, install directly from Git:

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/AlexBatrakov/StructureSolver.jl")
```

For local development:

```julia
julia> using Pkg
julia> Pkg.develop(path="/path/to/StructureSolver.jl")
```

## Documentation

- Dev docs: https://AlexBatrakov.github.io/StructureSolver.jl/dev
- Stable docs (after tagging a release): https://AlexBatrakov.github.io/StructureSolver.jl/stable

Build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Releases

- The docs site publishes `/stable` after you push a tag like `v0.1.2`.
- See `RELEASING.md` for the recommended release flow (TagBot) and dependency/compat automation (CompatHelper).

## 30-second demo (EoS + single solve)

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=1e4, dtmax=1.0, reltol=1e-10, abstol=1e-10)
regime = Simple_DirectRegime(
  Dict(:φc => 0.0, :pc => 1e35),
  Dict(:α0 => 0.0, :β0 => 0.0),
)

sim = SingleSimulation(model, regime, int_params)
calculate!(sim)

mass_Msun = sim.model.quantities[:mA] / M_sun
radius_km = sim.model.quantities[:R] / 1e5
@show mass_Msun radius_km
```

<details>
<summary><strong>Quick start (EoS)</strong></summary>

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)

rho = 1e15  # g/cm^3
p = get_pressure(eos, rho, from=:density)          # dyn/cm^2
eps = get_energy_density(eos, rho, from=:density)  # erg/cm^3
cs = get_sound_velocity(eos, rho, from=:density)   # cm/s
```

</details>

<details>
<summary><strong>Solve a single star (direct regime)</strong></summary>

Direct regime means you specify the inner parameters (e.g. central values) directly.

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=1e4, dtmax=1.0, reltol=1e-10, abstol=1e-10)

inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)

sim = SingleSimulation(model, regime, int_params)
calculate!(sim)

sim.model.quantities    # computed global quantities (Dict)
sim.model.derivatives   # computed derivatives/sensitivities (Dict)
```

</details>

<details>
<summary><strong>Solve with shooting (match boundary conditions)</strong></summary>

Shooting regime adjusts one or more "shooting" inner parameters to satisfy target boundary conditions on output quantities.

Example: shoot on `:φc` to enforce `:bc_φ∞ == 0.0`.

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF3_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=2e3, dtmax=1.0, reltol=1e-10, abstol=1e-10)

inparams_fixed = Dict(:pc => 5e35)
exparams = Dict(:α0 => -1e-3, :β0 => -5.0)
exparams_symbolic = Dict(:low_eos => :SLy, :high_eos => :SLy)
inparams_shooting = Dict(:φc => 1e-2)   # initial guess
quantities_fixed = Dict(:bc_φ∞ => 0.0)  # target condition

regime = Simple_ShootingRegime(inparams_fixed, exparams, exparams_symbolic, inparams_shooting, quantities_fixed)
sim = SingleSimulation(model, regime, int_params)
calculate!(sim)
```

</details>

<details>
<summary><strong>Families (scan central pressure)</strong></summary>

`FamilySimulation` sweeps one (or several) parameters and stores the resulting `model.inparams`, `model.quantities`, and `model.derivatives` in N-D arrays.

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=1e4, dtmax=1.0, reltol=1e-12, abstol=1e-12)
regime = Simple_DirectRegime(Dict(:φc => 0.0, :pc => 1e35), Dict(:α0 => 0.0, :β0 => 0.0))

family_params = Dict(:pc => LogRange(1e34, 1e36, 100))
sim = FamilySimulation(model, regime, family_params, int_params)
calculate!(sim)

mass_Msun = sim.family.quantities[:mA] ./ M_sun
radius_km = sim.family.quantities[:R] ./ 1e5
```

</details>

<details>
<summary><strong>Plotting (optional)</strong></summary>

This package does **not** depend on `PyPlot` by default (so CI and headless runs do not require `matplotlib`).

- Some helpers (e.g. `plot_radial_structure`) will attempt to `import PyPlot` only when called.
- Some scripts in `examples/` use `PyPlot` explicitly.

If you want plotting:

```julia
using Pkg
Pkg.add("PyPlot")
```

If `PyPlot` fails to import `matplotlib` (common on fresh machines/CI), follow the `PyCall` guidance in the error message (either install `matplotlib` for your system Python, or reconfigure PyCall to use Conda).

</details>

<details>
<summary><strong>N-D grids (GeneralSimulation)</strong></summary>

`GeneralSimulation` uses `ShootingRegime` where each entry in the dicts can be either:

- a scalar (fixed value), or
- a vector (grid dimension).

```julia
using StructureSolver

model_type = DEFp_Model{Float64,DEF3_CouplingFunction,PWP_EoS}
int_params = IntParams(maxiters=1e4, dtmax=1.0, reltol=1e-10, abstol=1e-10)

in_fixed = Dict(:pc => LogRange(1e34, 1e36, 50))
in_shoot = Dict(:φc => 0.1)
q_fixed = Dict(:bc_φ∞ => 0.0)
ex = Dict(
    :α0 => -1e-5,
    :β0 => [-6.0, -5.0, -4.0],
    :low_eos => :SLy,
    :high_eos => :WFF1,
)

regime = ShootingRegime(in_fixed, in_shoot, q_fixed, ex)
sim = GeneralSimulation{Float64}(model_type, regime, int_params)
calculate!(sim)

# N-D arrays with shape `sim.data.dims`
mA = sim.data.quantities[:mA]
αA = sim.data.derivatives[:αA]
```

</details>

## Testing

From the repository root:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

To generate coverage locally:

```bash
julia --project --code-coverage=user -e 'using Pkg; Pkg.test()'
```

## License

MIT. See `LICENSE`.
