# Notation & outputs

This page explains the most common outputs produced by the solver and how to interpret them.

Throughout the docs we use the term **model** to refer to the object returned by `solve_structure(...)` (or produced by the simulation helpers). That model stores derived results in two places:

- `model.quantities`: scalar “headline” results (mass, radius, etc.)
- `model.derivatives`: scalar sensitivities / response coefficients (useful for parameter estimation and post-processing)

The exact set of keys can depend on the gravity regime/model. The lists below describe the current **DEFp** implementation (see `DEFpGravity`).

## Conventions

- Units:
  - Internal integration is done in CGS-like units for microphysics, with standard conversions for plotting/reporting.
  - Unless stated otherwise, masses are in grams (or in solar masses when explicitly converted), radii are in cm (or km when explicitly converted), frequencies are in 1/s.
- The symbol `φ` denotes the scalar field in scalar–tensor gravity. `φc` is its central value.
- The subscript `A` denotes the star (object) evaluated at the surface / asymptotic matching.

Most results are stored with `Symbol` keys (some of them use Unicode). You can always inspect what your particular run produced:

```julia
sort!(collect(keys(model.quantities)))
sort!(collect(keys(model.derivatives)))
```

## `model.quantities`

The DEFp model currently defines the following keys.

### Global star properties

- `:R`
  - Stellar radius (surface where pressure drops to ~0).
  - Typical range: ~10–15 km for neutron stars.

- `:mA`
  - Gravitational (ADM) mass of the star.
  - Typical range: ~1–2 solar masses.

- `:m̃A`
  - “Jordan-frame” mass (a frame-dependent mass used in scalar–tensor conventions).
  - Useful for comparisons depending on which frame your observable is defined in.

### Scalarization-related

- `:φ∞`
  - Asymptotic value of the scalar field far from the star.

- `:αA`
  - Effective scalar coupling of the star (sometimes called the scalar charge per unit mass, depending on conventions).
  - Interpreted as how strongly the star sources / responds to the scalar field.

- `:bc_φ∞`
  - A boundary-condition residual for `φ∞`.
  - Near zero indicates the shooting/matching for the asymptotic scalar field has converged.

### Rotation-related

- `:Ω`
  - Rotation angular frequency used for slow-rotation quantities.

- `:JA`
  - Angular momentum.

- `:IA`
  - Moment of inertia.

### Log-transformed quantities

Some logs are stored to make derivatives and interpolation more stable:

- `:lnIA` — `log(IA)`
- `:lnmA` — `log(mA)`
- `:lnm̃A` — `log(m̃A)`

### Central conditions / EoS-related

- `:pc`
  - Central pressure.

- `:φc`
  - Central scalar field.

- `:cs_c`
  - Speed of sound at the center (as returned by the EoS).

## `model.derivatives`

The DEFp model also provides a small set of response coefficients.

- `:dφ∞_dφc`
  - Sensitivity of the asymptotic scalar `φ∞` to the central scalar `φc`.

- `:αA`, `:βA`, `:kA`
  - Scalarization-related response coefficients.
  - These are typically defined from derivatives of mass / scalar charge w.r.t. `φ∞` (exact conventions depend on the model).

## Sanity checks

A few quick checks that often catch configuration/unit mistakes early:

- `:R` is of order 1e6 cm (10 km) for neutron stars.
- `:mA` is of order 1e33 g (1 solar mass).
- `:bc_φ∞` should be close to 0 for converged solutions.
- If you enable slow rotation, `:IA` should be positive and scale roughly like `mA * R^2`.

## Where these are defined

- DEFp quantities and derivatives are assembled in the DEFp gravity implementation.
- If you want the authoritative list for your current version, inspect the `quantities_names` and derivative rules in the source.
