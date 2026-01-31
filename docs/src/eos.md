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

## Units

Where relevant, the package uses CGS conventions (e.g. density in g/cm^3, pressure in dyn/cm^2, energy density in erg/cm^3, speed in cm/s). If an EoS supports additional unit systems, the helpers typically expose a `units` keyword.

## Maximum pressure helper

`find_max_pressure(eos)` can be useful for validating tabulated or piecewise models and setting scan ranges.
