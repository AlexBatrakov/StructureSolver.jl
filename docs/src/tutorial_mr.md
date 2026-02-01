# Tutorial: Mass–radius (M–R) curve

This tutorial walks through a typical “research workflow” task: compute a mass–radius curve for a chosen EoS by scanning the central pressure $p_c$.

The goal is to help you answer these questions quickly:

- Which parameters do I need to set?
- Which solver knobs matter?
- Where do mass and radius appear in the outputs?

## 1) Pick a model and an EoS

For a first run, stick to GR-like settings ($\alpha_0=\beta_0=0$) and a built-in EoS name.

```julia
using StructureSolver

# Equation of state
# (names depend on what the package ships; SLy is a good default)
eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)

# Coupling function + model
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
```

## 2) Choose a regime (direct vs shooting)

For an M–R curve, you typically want **direct regime**:

- you fix inner parameters directly (here: $p_c$ and $\varphi_c$)
- you fix external parameters (here: $\alpha_0$ and $\beta_0$)

```julia
inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)   # pc will be overwritten by the scan
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)
```

If you need to enforce boundary conditions (e.g. drive `:bc_φ∞` to zero by adjusting `:φc`), switch to a shooting regime. For the M–R curve tutorial we keep it simple.

## 3) Solver controls (IntParams)

`IntParams` controls the ODE solve (tolerances, step size limits, max iterations).

```julia
int_params = IntParams(
    maxiters=2e3,
    dtmax=10.0,
    reltol=1e-8,
    abstol=1e-8,
)
```

For production scans you will often tighten tolerances and increase `maxiters`.

## 4) Run a family scan over central pressure

A family scan is simply a grid over one parameter (here: `:pc`).

```julia
family_params = Dict(:pc => LogRange(1e34, 1e36, 40))

sim = FamilySimulation(model, regime, family_params, int_params)
calculate!(sim)
```

## 5) Extract mass and radius arrays

After `calculate!(sim)`, the results are stored in `sim.family.quantities` as arrays.

```julia
mA = sim.family.quantities[:mA]   # gravitational mass (CGS)
R  = sim.family.quantities[:R]    # radius (cm)

m_Msun = mA ./ M_sun
R_km = R ./ 1e5
```

At this point you can plot `(R_km, m_Msun)` or save it to a file.

## 6) Optional: save as CSV-like table

To keep dependencies minimal, here is a very simple tab-separated export:

```julia
open("mr.tsv", "w") do io
    println(io, "R_km\tm_Msun")
    for i in eachindex(R_km)
        println(io, "$(R_km[i])\t$(m_Msun[i])")
    end
end
```

## 7) Optional plotting

The package does not depend on plotting backends by default. You can plot with whatever you like.

If you prefer `PyPlot`:

```julia
using Pkg
Pkg.add("PyPlot")

using PyPlot
plot(R_km, m_Msun, ".-")
xlabel("R [km]")
ylabel("M [M⊙]")
grid(true)
```

## Sanity checks

A few quick checks you can do when debugging a scan:

- `all(isfinite, m_Msun)` and `all(isfinite, R_km)`
- `maximum(m_Msun)` should be a reasonable neutron-star mass scale (order unity in $M_\odot$)
- If you see many `NaN`s, try increasing `maxiters` and/or tightening tolerances.
