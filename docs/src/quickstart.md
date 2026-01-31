# Quickstart

## EoS: basic conversions

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)

rho = 1e15  # g/cm^3
p   = get_pressure(eos, rho, from=:density)          # dyn/cm^2
ε   = get_energy_density(eos, rho, from=:density)    # erg/cm^3
cs  = get_sound_velocity(eos, rho, from=:density)    # cm/s
```

## Solve a single star (direct regime)

Direct regime means you specify the inner parameters (e.g. central values) directly.

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=2e3, dtmax=10.0, reltol=1e-8, abstol=1e-8)

inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)

sim = SingleSimulation(model, regime, int_params)
calculate!(sim)

sim.model.quantities
sim.model.derivatives
```

Tip: for quick sanity checks you can inspect typical keys like `:mA`, `:m̃A`, `:R`, `:αA`, `:φ∞`.
