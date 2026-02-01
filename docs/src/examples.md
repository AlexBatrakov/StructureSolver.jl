# Examples

This page collects short, copy-pastable examples. For longer scripts, see the repository `examples/` folder.

## Mass–radius point (single solve)

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=2e3, dtmax=10.0, reltol=1e-8, abstol=1e-8)
regime = Simple_DirectRegime(Dict(:φc => 0.0, :pc => 1e35), Dict(:α0 => 0.0, :β0 => 0.0))

sim = SingleSimulation(model, regime, int_params)
calculate!(sim)

m_Msun = sim.model.quantities[:mA] / M_sun
R_km = sim.model.quantities[:R] / 1e5
(m_Msun, R_km)
```

## Mass–radius curve (family scan)

```julia
using StructureSolver

eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)

int_params = IntParams(maxiters=2e3, dtmax=10.0, reltol=1e-8, abstol=1e-8)
regime = Simple_DirectRegime(Dict(:φc => 0.0, :pc => 1e35), Dict(:α0 => 0.0, :β0 => 0.0))

family_params = Dict(:pc => LogRange(1e34, 1e36, 40))
sim = FamilySimulation(model, regime, family_params, int_params)
calculate!(sim)

m_Msun = sim.family.quantities[:mA] ./ M_sun
R_km = sim.family.quantities[:R] ./ 1e5
(m_Msun, R_km)
```

## Plotting (optional)

The package does not depend on `PyPlot` by default. If you want plotting helpers:

```julia
using Pkg
Pkg.add("PyPlot")
```

Then you can use plotting helpers (if present in your workflow), or plot your own arrays with your preferred backend.
