module StructureSolver
	
	using ProgressMeter
	using Distributed
	using DifferentialEquations
	#using ODEInterfaceDiffEq
	#using BoundaryValueDiffEq
	#using OrdinaryDiffEq
	using JLD
	#using Dierckx
	using NLsolve
	#using Dates
	using Printf
	using Interpolations
	using ForwardDiff
	using Roots

	include("physical_constants.jl")
	include("eos.jl")
	include("eos_data.jl")
	include("utils.jl")
	include("model_interface.jl")
	include("regime_interface.jl")
	include("simulation_interface.jl")
	include("DEFp_gravity.jl")

	export PWP_EoS, Polytropic_EoS, Table_EoS
	export get_pressure, get_number_density, get_density,
        get_energy_density, get_specific_enthalpy, get_sound_velocity,
        get_internal_energy, find_max_pressure

	export DEFp_Model, DEF1_CouplingFunction, DEF3_CouplingFunction,
        Simple_DirectRegime, SingleSimulation, IntParams, FamilySimulation,
        GridSimulation
	export calculate_old!, calculate!, internal_physics, A_fun, V_fun,
        α_fun, β_fun, φ∞_const, LogRange
	export M_sun, high_pwp3_eos_names, mp, c, G

	export ode_system!, internal_physics, calculate_quantities

	export Simple_ShootingRegime, ShootingRegime, GeneralSimulation

	"""
	    calculate_old!(args...)

	Backward-compatible alias for `calculate!`.
	"""
	calculate_old!(args...) = calculate!(args...)

end # module
  