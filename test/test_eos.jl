using Test
using StructureSolver

@testset "EoS" begin
    @testset "PWP_EoS basic conversions" begin
        eos = PWP_EoS(low_eos=:SLy, high_eos=:SLy)

        ρ = 1e15
        p = get_pressure(eos, ρ; from=:density)
        @test isfinite(p)
        @test p > 0

        ρ_back = get_density(eos, p; from=:pressure)
        @test isapprox(ρ_back, ρ; rtol=1e-6)

        ε = get_energy_density(eos, ρ; from=:density)
        @test isfinite(ε)
        @test ε > 0

        cs = get_sound_velocity(eos, ρ; from=:density)
        @test isfinite(cs)
        @test 0 < cs <= 10c
    end

    @testset "Polytropic_EoS inversion" begin
        eos = Polytropic_EoS(1.0, 2.0)
        ρ = 3.0
        p = get_pressure(eos, ρ; from=:density)
        @test p ≈ ρ^2
        @test get_density(eos, p; from=:pressure) ≈ ρ
    end

    @testset "Table_EoS constructor by name" begin
        eos = Table_EoS(:SLy4)

        ρ = 1e15
        p = get_pressure(eos, ρ; from=:density)
        @test isfinite(p)
        @test p > 0

        ε = get_energy_density(eos, ρ; from=:density)
        @test isfinite(ε)
        @test ε > 0
    end
end

