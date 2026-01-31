using Test
using StructureSolver

@testset "StructureSolver.jl" begin
	include("test_utils.jl")
	include("test_eos.jl")
	include("test_smoke_solver.jl")
end