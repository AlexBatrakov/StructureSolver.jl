# Equations of State (EoS)

StructureSolver.jl provides multiple EoS backends behind a common interface.

## Implementations

- `PWP_EoS`: piecewise polytrope (with built-in named low/high segments).
- `Polytropic_EoS`: simple polytrope.
- `Table_EoS`: tabulated EoS (default stable implementation).
- `Table_EoS_Hermite`: experimental Hermite-based tabulated EoS (not the default).

## Core interface

The main helpers are

- `get_pressure`, `get_density`, `get_number_density`
- `get_energy_density`, `get_specific_enthalpy`, `get_internal_energy`
- `get_sound_velocity`

Most methods accept a `from` keyword (e.g. `from=:density` or `from=:pressure`).

Example (pressure from density):

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)

ρ = 1e15
p = get_pressure(eos, ρ, from=:density)
```

Example (density from pressure):

```julia
using StructureSolver

eos = Polytropic_EoS(K=1e5, Γ=2.0)

p = 1e34
ρ = get_density(eos, p, from=:pressure)
```

## Units

Where relevant, the package uses CGS conventions:

- density $\rho$ in g/cm^3
- pressure $p$ in dyn/cm^2
- energy density $\epsilon$ in erg/cm^3
- sound speed $c_s$ in cm/s

If an EoS supports additional unit systems, helpers typically expose a `units` keyword.

Tip: when you are mixing different EoS backends in one project, keep an eye on the `units` keyword and always document what your input arrays represent.

## Tabulated EoS

`Table_EoS(:NAME)` loads a built-in tabulated EoS by name (see the package data). This is the default stable implementation.

`Table_EoS_Hermite` is experimental and not used by default.

## Maximum pressure helper

`find_max_pressure(eos)` can be useful for validating tabulated or piecewise models and setting scan ranges.
