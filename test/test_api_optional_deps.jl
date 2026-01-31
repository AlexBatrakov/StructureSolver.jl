using Test
using StructureSolver

function _is_pkg_loaded(pkgname::AbstractString)
    return any(k.name == pkgname for k in keys(Base.loaded_modules))
end

@testset "API & optional deps" begin
    # Regression: package load must not implicitly pull plotting stack.
    @test !_is_pkg_loaded("PyPlot")

    # Public-ish API sanity checks (names that should exist in the module).
    @test isdefined(StructureSolver, :Table_EoS)
    @test isdefined(StructureSolver, :PWP_EoS)
    @test isdefined(StructureSolver, :Polytropic_EoS)
    @test isdefined(StructureSolver, :DEFp_Model)
    @test isdefined(StructureSolver, :calculate!)
    @test isdefined(StructureSolver, :calculate_quantities)

    # Keep the experimental Hermite implementation available (not necessarily exported).
    @test isdefined(StructureSolver, :Table_EoS_Hermite)
end
