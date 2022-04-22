abstract type AbstractSimulation end

#fields

#methods

calculate(simulation::AbstractSimulation) = error("It is not defined how to calculate simulation for $(typeof(simulation))")


#realisations
#-----------------------------------------------------

struct IntParams
	maxiters::Float64
	dtmax::Float64
	reltol::Float64
	abstol::Float64
end

IntParams(;maxiters::Float64, dtmax::Float64, reltol::Float64, abstol::Float64) = IntParams(maxiters, dtmax, reltol, abstol)


struct SingleSimulation{T <: AbstractModel, P <: AbstractRegime} <: AbstractSimulation
	model::T
	regime::P
	int_params::IntParams
end

function calculate!(simulation::SingleSimulation{T, P}) where {T <: AbstractModel, P <: AbstractRegime}
	model = simulation.model
	regime = simulation.regime
	int_params = simulation.int_params

	calculate_model!(model, regime, int_params)
	return simulation
end

#-----------------------------------------------------

struct Family{T <: AbstractModel}
	sample_model::T
	family_params::Dict{Symbol,Vector{Float64}}
#	family_params_names::Vector{Symbol}
#	family_params_arr::Array{NTuple{N,Float64} where N}
	dims::NTuple{N,Int64} where N
	inparams::Dict{Symbol,Array{Float64}}
#	exparams::Dict{Symbol,Any}

	quantities::Dict{Symbol,Array{Float64}}
	derivatives::Dict{Symbol,Array{Float64}}
	family_quantities::Dict{Symbol,Float64}
	function Family{T}(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}) where {T <: AbstractModel}

	#	family_params_names = collect(keys(family_params))
	#	family_params_arr = collect(Iterators.product(values(family_params)...))
	#	dims = size(family_params_arr)
		N = length(family_params)
		dims = Tuple(length(family_params[key]) for key in keys(family_params))

		inparams = Dict{Symbol, Array{Float64,N}}()
		quantities = Dict{Symbol, Array{Float64,N}}()
		derivatives = Dict{Symbol, Array{Float64,N}}()
		family_quantities = Dict{Symbol, Float64}()
		return new(sample_model, family_params, dims, inparams, quantities, derivatives, family_quantities)
	end
end

Family(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}) where {T <: AbstractModel} = Family{T}(sample_model, family_params)

function save_result!(family::Family{T}, I::CartesianIndex) where {T <: AbstractModel}
	for field_name in (:inparams, :quantities, :derivatives)
		family_field = getfield(family, field_name)
		model_field = getfield(family.sample_model, field_name)
		for name in keys(model_field)
			if !haskey(family_field, name)
				family_field[name] = Array{Float64}(undef, family.dims)
			end
			family_field[name][I] = model_field[name]
		end
	end
end

struct FamilySimulation{T <: AbstractModel, P <: AbstractRegime} <: AbstractSimulation
	family::Family{T}
	regime::P
	int_params::IntParams
	function FamilySimulation{T,P}(sample_model::T, regime::P, family_params::Dict{Symbol,Vector{Float64}}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime}
		family = Family(sample_model, family_params)
		return new(family, regime, int_params)
	end
end

FamilySimulation(sample_model::T, regime::P, family_params::Dict{Symbol,Vector{Float64}}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime} = FamilySimulation{T,P}(sample_model, regime, family_params, int_params)

function Base.show(io::IO, simulation::FamilySimulation)
	println(io, "Family Simulation")
	println(io, "Regime:  ", simulation.regime)
	println(io, "Family:  ", simulation.family)
	println(io, "Integration parameters:  ", simulation.int_params)
	return nothing
end

function write_to_regime(regime::Simple_DirectRegime, name, value)
	if haskey(regime.inparams_fixed, name)
		regime.inparams_fixed[name] = value
	else
		error("inappropriate parameter for Simple_DirectRegime")
	end
end

function write_to_regime(regime::Simple_ShootingRegime, name, value)
	if haskey(regime.inparams_fixed, name)
		regime.inparams_fixed[name] = value
	elseif haskey(regime.inparams_shooting, name)
		regime.inparams_shooting[name] = value
	elseif haskey(regime.quantities_fixed, name)
		regime.quantities_fixed[name] = value
	else
		error("inappropriate parameter $name for Simple_DirectRegime")
	end
end

function calculate!(simulation::FamilySimulation{T, P}) where {T <: AbstractModel, P <: AbstractRegime}
	family = simulation.family
	regime = simulation.regime
	int_params = simulation.int_params
	for I in CartesianIndices(family.dims)
		for (i, name) in enumerate(keys(family.family_params))
			value = family.family_params[name][I[i]]
			write_to_regime(regime, name, value)
			#regime.quantities_fixed[name] = family.family_params[name][I[i]]
		end

#		regime.inparams_shooting[:φc] = 0.1
#		println(regime)

#		if haskey(family.sample_model.derivatives, :dφc_dpc)
#			dφc_dpc = family.sample_model.derivatives[:dφc_dpc]
#			regime.inparams_shooting[:φc] += dφc_dpc * (family.family_params[:pc][I] - family.family_params[:pc][I-1])
#			println(regime.inparams_shooting[:φc])
#		end
		calculate_model!(family.sample_model, regime, int_params)
		save_result!(family, I)
	end
	return simulation
end

#-----------------------------------------------------

struct Grid_old{T <: AbstractModel}
	sample_model::T
	family_params::Dict{Symbol,Vector{Float64}}
	family_dims::NTuple{N,Int64} where N
	grid_params::Dict{Symbol,Q} where {Q <: Vector}
	grid_dims::NTuple{N,Int64} where N

	inparams::Dict{Symbol,Array{Float64}}
	exparams::Dict{Symbol,Q} where {Q <: Array}
	quantities::Dict{Symbol,Array{Float64}}
	derivatives::Dict{Symbol,Array{Float64}}
	family_quantities::Dict{Symbol,Array{Float64}}
	function Grid_old{T}(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}) where {T <: AbstractModel}
		N = length(family_params) + length(grid_params)

		family_dims = Tuple(length(family_params[key]) for key in keys(family_params))
		grid_dims = Tuple(length(grid_params[key]) for key in keys(grid_params))

		inparams = Dict{Symbol, Array{Float64,N}}()
		quantities = Dict{Symbol, Array{Float64,N}}()
		derivatives = Dict{Symbol, Array{Float64,N}}()
		family_quantities = Dict{Symbol, Float64}()
		return new(sample_model, family_params, family_dims, grid_params, grid_dims, inparams, quantities, derivatives, family_quantities)
	end
end

Grid_old(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}) where {T <: AbstractModel} = Grid{T}(sample_model, family_params, grid_params)

struct GridSimulation_old{T <: AbstractModel, P <: AbstractRegime} <: AbstractSimulation
	regime::P
	grid::Grid_old{T}
	int_params::IntParams
	function GridSimulation_old{T,P}(sample_model::T, regime::P, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime}
		grid = Grid_old(sample_model, family_params, grid_params)
		return new(regime, grid, int_params)
	end
end

GridSimulation_old(sample_model::T, regime::P, family_params::Dict{Symbol,Vector{Float64}},  grid_params::Dict{Symbol,Q} where {Q <: Vector}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime} = GridSimulation_old{T,P}(sample_model, regime, family_params, grid_params, int_params)

function calculate!(simulation::GridSimulation_old{T, P}) where {T <: AbstractModel, P <: AbstractRegime}
	grid = simulation.grid
	regime = simulation.regime
	int_params = simulation.int_params
	family_exparams = copy(grid.sample_model.exparams)
	for I_grid in CartesianIndices(grid.grid_dims)
		for (i, name) in enumerate(keys(grid.grid_params))
			family_exparams[name] = grid.grid_params[name][I_grid[i]]
		end
		family_sample_model = typeof(grid.sample_model)(family_exparams)
		family_simulation = FamilySimulation(family_sample_model, regime, grid.family_params, int_params)
		println(family_simulation.family.sample_model)
		println(regime)
		@time calculate!(family_simulation)
		save_result!(grid, family_simulation.family, I_grid)

			qs = family_simulation.family.quantities
			plot(qs[:pc][:], -qs[:αA][:], label="$I_grid")
			xscale("log")
	end
	return simulation
end

function save_result!(grid::Grid_old{T}, family::Family{T}, I_grid::CartesianIndex) where {T <: AbstractModel}
	for field_name in (:inparams, :quantities, :derivatives)
		grid_field = getfield(grid, field_name)
		family_field = getfield(family, field_name)
		model_field = getfield(family.sample_model, field_name)
		Is_family = CartesianIndices(grid.family_dims)
#		full_dims = Tuple([grid.family_dims..., grid.grid_dims...])
		for name in keys(model_field)
			if !haskey(grid_field, name)
				grid_field[name] = Array{Float64}(undef, grid.family_dims..., grid.grid_dims...)
			end
			grid_field[name][Is_family, I_grid] .= family_field[name]
		end
	end
end

#-----------------------------------------------------

struct Grid{T <: AbstractModel}
	sample_model::T
	family_params::Dict{Symbol,Vector{Float64}}
	family_dims::NTuple{N,Int64} where N
	grid_params::Dict{Symbol,Q} where {Q <: Vector}
	grid_dims::NTuple{N,Int64} where N

	inparams::Dict{Symbol,Array{Float64}}
	exparams::Dict{Symbol,Q} where {Q <: Array}
	quantities::Dict{Symbol,Array{Float64}}
	derivatives::Dict{Symbol,Array{Float64}}
	family_quantities::Dict{Symbol,Array{Float64}}
	function Grid{T}(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}) where {T <: AbstractModel}
		N = length(family_params) + length(grid_params)

		family_dims = Tuple(length(family_params[key]) for key in keys(family_params))
		grid_dims = Tuple(length(grid_params[key]) for key in keys(grid_params))

		inparams = Dict{Symbol, Array{Float64,N}}()
		quantities = Dict{Symbol, Array{Float64,N}}()
		derivatives = Dict{Symbol, Array{Float64,N}}()
		family_quantities = Dict{Symbol, Float64}()
		return new(sample_model, family_params, family_dims, grid_params, grid_dims, inparams, quantities, derivatives, family_quantities)
	end
end

Grid(sample_model::T, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}) where {T <: AbstractModel} = Grid{T}(sample_model, family_params, grid_params)

struct GridSimulation{T <: AbstractModel, P <: AbstractRegime} <: AbstractSimulation
	grid::Grid{T}
	sample_regime::P
	int_params::IntParams
	function GridSimulation{T,P}(sample_model::T, sample_regime::P, family_params::Dict{Symbol,Vector{Float64}}, grid_params::Dict{Symbol,Q} where {Q <: Vector}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime}
		grid = Grid(sample_model, family_params, grid_params)
		return new(grid, sample_regime, int_params)
	end
end

GridSimulation(sample_model::T, sample_regime::P, family_params::Dict{Symbol,Vector{Float64}},  grid_params::Dict{Symbol,Q} where {Q <: Vector}, int_params::IntParams) where {T <: AbstractModel, P <: AbstractRegime} = GridSimulation{T,P}(sample_model, sample_regime, family_params, grid_params, int_params)

function GridSimulation(family_constructor, family_params::Dict{Symbol,Vector{Float64}},  grid_params::Dict{Symbol,Q} where {Q <: Vector}, int_params::IntParams)
		I0 = ones(Int,length(grid_params))
		sample_model, sample_regime = family_constructor(grid_params, I0)
	return GridSimulation(sample_model, sample_regime, family_params,  grid_params, int_params)
end

function calculate!(simulation::GridSimulation{T, P}, family_constructor) where {T <: AbstractModel, P <: AbstractRegime}
	grid = simulation.grid
	regime = simulation.sample_regime
	int_params = simulation.int_params
	family_exparams = copy(grid.sample_model.exparams)
	for I_grid in CartesianIndices(grid.grid_dims)

		family_sample_model, family_regime = family_constructor(grid.grid_params, I_grid)
		family_simulation = FamilySimulation(family_sample_model, family_regime, grid.family_params, int_params)
		println(family_sample_model)
		println(family_regime)
		@time calculate!(family_simulation)
		save_result!(grid, family_simulation.family, I_grid)

			qs = family_simulation.family.quantities
			plot(qs[:pc][:], -qs[:αA][:], label="$I_grid")
			xscale("log")
	end
	return simulation
end

function save_result!(grid::Grid{T}, family::Family{T}, I_grid::CartesianIndex) where {T <: AbstractModel}
	for field_name in (:inparams, :quantities, :derivatives)
		grid_field = getfield(grid, field_name)
		family_field = getfield(family, field_name)
		model_field = getfield(family.sample_model, field_name)
		Is_family = CartesianIndices(grid.family_dims)
#		full_dims = Tuple([grid.family_dims..., grid.grid_dims...])
		for name in keys(model_field)
			if !haskey(grid_field, name)
				grid_field[name] = Array{Float64}(undef, grid.family_dims..., grid.grid_dims...)
			end
			grid_field[name][Is_family, I_grid] .= family_field[name]
		end
	end
end

function calculate_parallel!(simulation::GridSimulation, family_constructor) end


#-----------------------------------------------------

struct SimulationData{T <: Real, N}
	params_names::NTuple{N,Symbol}
	dims::NTuple{N,Int64}
	quantities::Dict{Symbol,Array{T,N}}
	derivatives::Dict{Symbol,Array{T,N}}
end

function Base.show(io::IO, data::SimulationData)
	println(io, "Names of grid parameters:")
	println(io, data.params_names)
	println(io, "Sizes of grid parameters:")
	println(io, data.dims)
	println(io, "Calculated quantities:")
	print_dict(io, data.quantities)
	println(io, "Calculated derivatives:")
	print_dict(io, data.derivatives)
	return nothing
end

struct GeneralSimulation{T1 <: Real, T2 <: AbstractRegime, N} <: AbstractSimulation
	model_type::DataType
	regime::T2
	data::SimulationData{T1, N}
	int_params::IntParams
	function GeneralSimulation{T1}(model_type::DataType, regime::T2, int_params::IntParams) where {T1 <: Real, T2 <: AbstractRegime}
		params_names = get_params_names(regime)
		params_dims = get_params_dims(regime)
		N = length(params_names)
		quantities = Dict{Symbol,Array{T1,N}}()
		derivatives = Dict{Symbol,Array{T1,N}}()
		data = SimulationData(params_names, params_dims, quantities, derivatives)
		return new{T1,T2,N}(model_type, regime, data, int_params)
	end
end

GeneralSimulation(model_type::DataType, regime::T2, int_params::IntParams) where {T2 <: AbstractRegime} = GeneralSimulation{Float64}(model_type, regime, int_params)

function Base.show(io::IO, simulation::GeneralSimulation)
	println(io, "General Simulation")
	println(io, "Model type:")
	println(io, simulation.model_type)
	println(io, "Regime:")
	println(io, simulation.regime)
	println(io, "Caluclated data:")
	println(io, simulation.data)
	return nothing
end


struct ShootingRegime{T <: Real} <: AbstractShootingRegime
	inparams_fixed::Dict{Symbol,T}
	inparams_fixed_grid::Dict{Symbol,Vector{T}}
	inparams_shooting::Dict{Symbol,T}
	quantities_fixed::Dict{Symbol,T}
	quantities_fixed_grid::Dict{Symbol,Vector{T}}
	exparams::Dict{Symbol,T}
	exparams_grid::Dict{Symbol,Vector{T}}
	exparams_symbolic::Dict{Symbol,Symbol}
	exparams_symbolic_grid::Dict{Symbol,Vector{Symbol}}
	function ShootingRegime{T}(inparams_fixed_in::Dict{Symbol}, inparams_shooting_in::Dict{Symbol}, quantities_fixed_in::Dict{Symbol}, exparams_in::Dict{Symbol}) where  {T <: Real}
		inparams_fixed = Dict{Symbol,T}()
		inparams_fixed_grid = Dict{Symbol,Vector{T}}()
		inparams_shooting = Dict{Symbol,T}()
		quantities_fixed = Dict{Symbol,T}()
		quantities_fixed_grid = Dict{Symbol,Vector{T}}()
		exparams = Dict{Symbol,T}()
		exparams_grid = Dict{Symbol,Vector{T}}()
		exparams_symbolic = Dict{Symbol,Symbol}()
		exparams_symbolic_grid = Dict{Symbol,Vector{Symbol}}()

		for key in keys(inparams_fixed_in)
			if typeof(inparams_fixed_in[key]) <: Real
				inparams_fixed[key] = inparams_fixed_in[key]
			elseif typeof(inparams_fixed_in[key]) <: Array{T} where {T <: Real}
				inparams_fixed_grid[key] = copy(inparams_fixed_in[key])
			else
				error("Unexpected type for $key in inparams_fixed_in")
			end
		end

		inparams_shooting = copy(inparams_shooting_in)

		for key in keys(quantities_fixed_in)
			if typeof(quantities_fixed_in[key]) <: Real
				quantities_fixed[key] = quantities_fixed_in[key]
			elseif typeof(quantities_fixed_in[key]) <: Array{T} where {T <: Real}
				quantities_fixed_grid[key] = copy(quantities_fixed_in[key])
			else
				error("Unexpected type for $key in quantities_fixed_in")
			end
		end

		for key in keys(exparams_in)
			if typeof(exparams_in[key]) <: Real
				exparams[key] = exparams_in[key]
			elseif typeof(exparams_in[key]) == Symbol
				exparams_symbolic[key] = exparams_in[key]
			elseif typeof(exparams_in[key]) <: Array{T} where {T <: Real}
				exparams_grid[key] = copy(exparams_in[key])
			elseif typeof(exparams_in[key]) <: Array{Symbol}
				exparams_symbolic_grid[key] = copy(exparams_in[key])
			else
				error("Unexpected type for $key in exparams_in")
			end
		end

		return new(
			inparams_fixed,
			inparams_fixed_grid,
			inparams_shooting,
			quantities_fixed,
			quantities_fixed_grid,
			exparams,
			exparams_grid,
			exparams_symbolic,
			exparams_symbolic_grid)
	end
end

ShootingRegime(inparams_fixed_in::Dict{Symbol}, inparams_shooting_in::Dict{Symbol}, quantities_fixed_in::Dict{Symbol}, exparams_in::Dict{Symbol}) = ShootingRegime{Float64}(inparams_fixed_in, inparams_shooting_in, quantities_fixed_in, exparams_in)

function get_params_names(regime::ShootingRegime)
	inparams_fixed_grid_names = collect(keys(regime.inparams_fixed_grid))
	quantities_fixed_grid_names = collect(keys(regime.quantities_fixed_grid))
	exparams_grid_names = collect(keys(regime.exparams_grid))
	exparams_symbolic_grid_names = collect(keys(regime.exparams_symbolic_grid))
	params_names_arr = vcat(inparams_fixed_grid_names, quantities_fixed_grid_names, exparams_grid_names, exparams_symbolic_grid_names)
	return Tuple(params_names_arr)
end

function get_params_dims(regime::ShootingRegime)
	inparams_fixed_grid_dims = length.(values(regime.inparams_fixed_grid))
	quantities_fixed_grid_dims = length.(values(regime.quantities_fixed_grid))
	exparams_grid_dims = length.(values(regime.exparams_grid))
	exparams_symbolic_grid_dims = length.(values(regime.exparams_symbolic_grid))
	params_dims = vcat(inparams_fixed_grid_dims, quantities_fixed_grid_dims, exparams_grid_dims, exparams_symbolic_grid_dims)
	return Tuple(params_dims)
end

function Base.show(io::IO, regime::ShootingRegime)
	println(io, "Shooting Regime")
	println(io, "Inner parameters fixed")
	print_dict(io, regime.inparams_fixed)
	print_dict(io, regime.inparams_fixed_grid)
	println(io, "Inner parameters shooting:")
	print_dict(io, regime.inparams_shooting)
	println(io, "Quantities fixed:")
	print_dict(io, regime.quantities_fixed)
	print_dict(io, regime.quantities_fixed_grid)
	println(io, "Exteriour parameters:")
	print_dict(io, regime.exparams)
	print_dict(io, regime.exparams_grid)
	print_dict(io, regime.exparams_symbolic)
	print_dict(io, regime.exparams_symbolic_grid)
	return nothing
end

function calculate!(simulation::GeneralSimulation)
	
	regime = simulation.regime
	int_params = simulation.int_params

	sample_model = simulation.model_type()
	sample_regime = Simple_ShootingRegime(regime.inparams_fixed, regime.exparams, regime.exparams_symbolic, regime.inparams_shooting, regime.quantities_fixed)
	sample_direct_regime = Simple_DirectRegime(regime.inparams_fixed, regime.exparams)

	#get_indices(grid) = CartesianIndices(Tuple(length.(values(grid))))
	#@time ex_sym_indices = get_indices(regime.exparams_symbolic_grid)
	#@time ex_indices = get_indices(regime.exparams_grid)
	#@time qf_indices = get_indices(regime.quantities_fixed_grid)
	#@time in_indices = get_indices(regime.inparams_fixed_grid)

	ex_sym_dim = Tuple(length.(values(regime.exparams_symbolic_grid)))
	ex_dim = Tuple(length.(values(regime.exparams_grid)))
	qf_dim = Tuple(length.(values(regime.quantities_fixed_grid)))
	in_dim = Tuple(length.(values(regime.inparams_fixed_grid)))

	function task()
		for I_ex_sym in CartesianIndices(ex_sym_dim)
			for (i, name) in enumerate(keys(regime.exparams_symbolic_grid))
				sample_regime.exparams_symbolic[name] = regime.exparams_symbolic_grid[name][I_ex_sym[i]]
			end
			for I_ex in CartesianIndices(ex_dim)
				for (i, name) in enumerate(keys(regime.exparams_grid))
					sample_regime.exparams[name] = regime.exparams_grid[name][I_ex[i]]
				end
				if is_the_same_model(sample_model, regime.exparams, regime.exparams_symbolic) == false
        			sample_model = typeof(sample_model)(regime.exparams, regime.exparams_symbolic)
    			end
    			sample_model = update_model!(sample_model, sample_regime.exparams, sample_regime.exparams_symbolic)
				#println(sample_regime)
				#println(I_ex_sym, " ", I_ex)
				for I_qf in CartesianIndices(qf_dim)
					for (i, name) in enumerate(keys(regime.quantities_fixed_grid))
						sample_regime.quantities_fixed[name] = regime.quantities_fixed_grid[name][I_qf[i]]
					end
					for I_in in CartesianIndices(in_dim)
						for (i, name) in enumerate(keys(regime.inparams_fixed_grid))
							sample_regime.inparams_fixed[name] = regime.inparams_fixed_grid[name][I_in[i]]
						end

						calculate_model!(sample_model, sample_regime, int_params)
						
						#println("φc", sample_model.quantities[:φc])
						#println("bc_φ∞", sample_model.quantities[:bc_φ∞])
						save_to_data!(simulation, sample_model, (I_in, I_qf, I_ex, I_ex_sym))
					end
				end
			end
		end 
	end
	task()
	return simulation
end

function save_to_data!(simulation, model, I)
	data = simulation.data
	for field_name in (:quantities, :derivatives)
		data_field = getfield(data, field_name)
		model_field = getfield(model, field_name)
		for name in keys(model_field)
			if !haskey(data_field, name)
				data_field[name] = valtype(data_field)(undef, simulation.data.dims)
			end
			data_field[name][I...] = model_field[name]
		end
	end
end