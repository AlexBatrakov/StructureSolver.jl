using Test
using StructureSolver

@testset "utils" begin
    @test LogRange(1.0, 10.0, 4) ≈ [1.0, 10.0^(1/3), 10.0^(2/3), 10.0]

    xs = LogRange(1e-3, 1e3, 20)
    @test length(xs) == 20
    @test isapprox(first(xs), 1e-3; rtol=0, atol=0)
    @test isapprox(last(xs), 1e3; rtol=0, atol=0)
    @test all(diff(xs) .> 0)
end
