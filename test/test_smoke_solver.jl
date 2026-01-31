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

    expected_quantity_keys = (
        :R,
        :mA,
        :m̃A,
        :αA,
        :φ∞,
        :IA,
        :pc,
        :φc,
        :bc_φ∞,
        :cs_c,
    )

    for k in expected_quantity_keys
        @test haskey(sim.model.quantities, k)
        @test isfinite(sim.model.quantities[k])
    end

    @test sim.model.quantities[:R] > 0
    @test sim.model.quantities[:mA] > 0
    @test sim.model.quantities[:m̃A] > 0
    @test sim.model.quantities[:IA] > 0
    @test sim.model.quantities[:pc] > 0

    @test !isempty(sim.model.derivatives)
    @test haskey(sim.model.derivatives, :dφ∞_dφc)
    @test isfinite(sim.model.derivatives[:dφ∞_dφc])
end
