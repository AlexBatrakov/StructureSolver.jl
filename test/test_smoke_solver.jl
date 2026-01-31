using Test
using StructureSolver

@testset "smoke: single solve" begin
    eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
    cf = DEF1_CouplingFunction()
    model = DEFp_Model{Float64}(cf, eos)

    int_params = IntParams(maxiters=2e3, dtmax=10.0, reltol=1e-8, abstol=1e-8)

    inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
    exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
    regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)

    sim = SingleSimulation(model, regime, int_params)
    calculate!(sim)

    @test !isempty(sim.model.quantities)
    @test any(haskey(sim.model.quantities, k) for k in (:mA, :m̃A, :R, :r_s))
end
