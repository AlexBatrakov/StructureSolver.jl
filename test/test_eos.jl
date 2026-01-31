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

        ε_from_p = get_energy_density(eos, p; from=:pressure)
        @test isfinite(ε_from_p)
        @test isapprox(ε_from_p, ε; rtol=1e-6)

        pmax = find_max_pressure(eos)
        @test isfinite(pmax)
        @test 0 < pmax <= 1e38
    end

    @testset "Polytropic_EoS inversion" begin
        eos = Polytropic_EoS(1.0, 2.0)
        ρ = 3.0
        p = get_pressure(eos, ρ; from=:density)
        @test p ≈ ρ^2
        @test get_density(eos, p; from=:pressure) ≈ ρ

        @test_throws ErrorException get_density(eos, 1.0; from=:unknown)
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

        ρ_back = get_density(eos, p; from=:pressure)
        @test isfinite(ρ_back)
        @test ρ_back > 0
        @test isapprox(ρ_back, ρ; rtol=1e-2)

        ε_from_p = get_energy_density(eos, p; from=:pressure)
        @test isfinite(ε_from_p)
        @test ε_from_p > 0

        @test_throws KeyError Table_EoS(:__definitely_not_a_real_eos__)
    end

    @testset "Invalid arguments" begin
        eos_pwp = PWP_EoS(low_eos=:SLy, high_eos=:SLy)
        @test_throws ErrorException get_density(eos_pwp, 1.0; from=:unknown)
        @test_throws ErrorException get_pressure(eos_pwp, 1e15; from=:density, units=:unknown)

        eos_tab = Table_EoS(:SLy4)
        @test_throws ErrorException get_density(eos_tab, 1.0; from=:unknown)
    end
end

