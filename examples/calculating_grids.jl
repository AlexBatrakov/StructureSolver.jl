using Distributed
using Revise
using PyPlot
using ProgressMeter

addprocs(8)

@everywhere using JLD
@everywhere using DelimitedFiles
@everywhere using StructureSolver


#high_eos_names = collect(keys(StructureSolver.high_pwp3_eos_base))
#high_eos_names = [:WFF1, :MPA1, :NL3, :Rs, :GM1, :SKb, :SLY2, :DDME2, :BSP, :BSR6Y, :GM1Y, :SK272, :SkI6, :H4, :TM1, :SLY9, :DD2, :SLY4, :NL3Yss, :ENG, :BSR2, :SkMP, :BSk22, :SKa, :NL3Y, :FSUGarnet, :NL3ωY, :IOPB, :DDHd, :DDME2Y, :FSU2, :SkI4, :BSk24, :APR3, :SLy, :NL3ωYss, :SkI3, :BSR6, :NL3ω, :SLY230a, :MS1, :BSk23, :WFF2, :MS1b, :Model1, :APR4, :SK255, :BSk20, :SkI5, :SkI2, :BSk25, :BSk26, :BSk21, :DD2Y]
#high_eos_names = [:SLy, :APR3, :APR4, :WFF1, :WFF2, :ENG, :MPA1, :MS1, :MS1b, :H4, :ALF2]
high_eos_names = [:BSk22]

@everywhere function calculate_grid(high_eos::Symbol)
	model_type = DEFp_Model{Float64,DEF1_CouplingFunction,PWP_EoS}
	int_params = IntParams(maxiters=1e4, dtmax=1.0, reltol=1e-11, abstol=1e-11)
	pc_grid = vcat(LogRange(1e34, 5e34, 1+20)[1:end-1], LogRange(5e34, 5e35, 1+80), LogRange(5e35, 1e36, 1+20)[2:end])
	inparams_fixed = Dict(:pc => pc_grid)
	imparams_shooting = Dict(:φc => 4e-1)
	quantities_fixed = Dict(:bc_φ∞ => 0.0)
	α0_grid = -LogRange(1e-1, 1e-5, 1+5*20)
#	β0_grid = collect(LinRange(-6.0, +6.0, N_β0))
	β0_grid = vcat(collect(LinRange(-6.0, -5.0, 1+20))[1:end-1],collect(LinRange(-5.0, -4.0, 1+50))[1:end-1], collect(LinRange(-4.0, +10.0, 1+20*14)))
	N_α0, N_β0, N_pc = length(α0_grid), length(β0_grid), length(pc_grid)
	N_α0, N_β0, N_pc = length(α0_grid), length(β0_grid), length(pc_grid)
	exparams = Dict(:α0 => α0_grid, :β0 => β0_grid, :low_eos => :SLy, :high_eos => high_eos)
	regime = ShootingRegime{Float64}(inparams_fixed, imparams_shooting, quantities_fixed, exparams)
	simulation = GeneralSimulation{Float64}(model_type, regime, int_params)
	calculate!(simulation)

	qs = simulation.data.quantities
	dr = simulation.data.derivatives
	gr = merge(simulation.regime.exparams_grid, simulation.regime.inparams_fixed_grid)

	save("$high_eos/$high_eos.jld", "simulation", simulation)
	function write_3Dgrid(grid_name, grid)
		open("$high_eos/$grid_name.dat", "w") do io
			for i in 1:N_α0, j in 1:N_β0
				writedlm(io, transpose(vcat(i, j, round.(grid[:,i,j], sigdigits=10))))
			end
		end
	end

	function write_1Dgrid(grid_name, grid)
		writedlm("$high_eos/$grid_name.dat", grid)
	end

	write_3Dgrid(:massA, qs[:mA] ./ M_sun)
	write_3Dgrid(:alphaA, dr[:αA])
	write_3Dgrid(:betaA, dr[:βA])
	write_3Dgrid(:kA, -dr[:kA])

	write_1Dgrid(:pc, gr[:pc])
	write_1Dgrid(:alpha0, gr[:α0])
	write_1Dgrid(:beta0, gr[:β0])

	writedlm("$high_eos/dims.dat", permutedims(hcat([simulation.data.params_names...], [simulation.data.dims...])))

	return nothing
end

function save_from_jld(high_eos::Symbol)
	simulation = load("$high_eos/$high_eos.jld", "simulation")
	N_pc, N_α0, N_β0 = simulation.data.dims
	qs = simulation.data.quantities
	dr = simulation.data.derivatives
	gr = merge(simulation.regime.exparams_grid, simulation.regime.inparams_fixed_grid)

	function write_3Dgrid(grid_name, grid)
		open("$high_eos/$grid_name.dat", "w") do io
			for i in 1:N_α0, j in 1:N_β0
				writedlm(io, transpose(vcat(i, j, round.(grid[:,i,j], sigdigits=10))))
			end
		end
	end

	function write_1Dgrid(grid_name, grid)
		writedlm("$high_eos/$grid_name.dat", grid)
	end

	mkpath("$high_eos")
	write_3Dgrid(:massA, qs[:mA] ./ M_sun)
	write_3Dgrid(:alphaA, dr[:αA])
	write_3Dgrid(:betaA, dr[:βA])
	write_3Dgrid(:kA, -dr[:kA])

	write_1Dgrid(:pc, gr[:pc])
	write_1Dgrid(:alpha0, gr[:α0])
	write_1Dgrid(:beta0, gr[:β0])

	writedlm("$high_eos/dims.dat", permutedims(hcat([simulation.data.params_names...], [simulation.data.dims...])))

	return nothing
end

@showprogress pmap(calculate_grid, high_eos_names)