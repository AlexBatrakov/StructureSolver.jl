 #AbstractRegime
abstract type AbstractRegime end

#necessary fields for AbstractRegime

#necessary methods for AbstractRegime
initiate_model!(model::AbstractModel, regime::AbstractRegime) = error("Function for setting initial conditions is not defined for $(typeof(regime))")

#realisations of AbstractRegime


#DirectRegime
abstract type AbstractDirectRegime <: AbstractRegime end

#necessary fields for DirectRegime
inparams_fixed(regime::AbstractDirectRegime) = regime.inparams_fixed
exparams(regime::AbstractDirectRegime) = regime.exparams

#necessary methods for DirectRegime

function Base.show(io::IO, regime::AbstractDirectRegime)
	println(io, "Direct Regime")
	println(io, "Inner parameters fixed:")
	print_dict(io, regime.inparams_fixed)
	println(io, "External parameters fixed:")
	print_dict(io, regime.exparams)
	return nothing
end

#realisations of DirectRegime

#Simple_DirectRegime
struct Simple_DirectRegime <: AbstractDirectRegime
	inparams_fixed::Dict{Symbol, Float64}
	exparams::Dict{Symbol, Float64}
end

function calculate_model!(model::AbstractModel, regime::Simple_DirectRegime, int_params)
	initiate_model!(model, regime)
	rs_problem = get_rs_problem(model)
	rs_solution = calculate_rs_solution(rs_problem, model, int_params)
		#println(length(rs_solution.t), " ", rs_solution.t[end])
		#plot(rs_solution[5,1:end-1], rs_solution[5,2:end] .- rs_solution[5,1:end-1],".")
		#show()
		#plot(rs_solution.t, rs_solution[5,:], ".")
		#plot(rs_solution[5,:], rs_solution[1,:] ./ rs_solution[1,end], ".")
		#show()
		#readline()
	calculate_quantities(rs_solution, model, save_to_model=true)
	calculate_derivatives!(rs_problem, model, int_params)
	return nothing
end

function initiate_model!(model::AbstractModel, regime::Simple_DirectRegime)
	for (key, value) in regime.inparams_fixed
		model.inparams[key] =  value
	end
	for (key, value) in regime.exparams
		model.exparams[key] =  value
	end
	return nothing
end

function get_rs_problem(model)
	rs_vars = get_rs_vars(model)
	rs_span = model.rs_span
	rs_problem = ODEProblem(ode_system!, rs_vars, rs_span, model)
	return rs_problem
end

function calculate_rs_solution(rs_problem, model, int_params)
	maxiters = int_params.maxiters
	dtmax = int_params.dtmax
	reltol = int_params.reltol
	abstol = int_params.abstol
	solver = Vern7()#Vern8()#Rosenbrock23()#Vern8()#Tsit5()
	rs_solution = solve(rs_problem, solver, callback=surface_callback(model), maxiters=maxiters, dtmax=dtmax, reltol=reltol, abstol=abstol, save_everystep=false, dtmin=eps(), force_dtmin=true)
	return rs_solution
end

function calculate_derivatives!(rs_problem, model, int_params) 
	function get_quantities(inparams_values)
		calc_type = eltype(inparams_values)
		model_dual = convert_model(calc_type, model)
		for (i,key) in enumerate(keys(model.inparams))
			model_dual.inparams[key] = inparams_values[i]
		end
		rs_problem_dual = get_rs_problem(model_dual)
		rs_solution_dual = calculate_rs_solution(rs_problem_dual, model_dual, int_params)
		calculate_quantities(rs_solution_dual, model_dual, save_to_model=true)
		return collect(values(model_dual.quantities))
	end

	inparams_values = collect(values(model.inparams))
	DQs_values = ForwardDiff.jacobian(get_quantities, inparams_values)
	DQs = Dict{Symbol,typeof(DQs_values[1,:])}()
	for (i,key) in enumerate(keys(model.quantities))
		DQs[key] = DQs_values[i,:]
	end
	for der_name in keys(model.der_rules)
		A, B, C, N = model.der_rules[der_name]
		dA_dx, dA_dy = DQs[A]
		dB_dx, dB_dy = DQs[B]
		dC_dx, dC_dy = DQs[C]
		model.derivatives[der_name] = (dA_dx * dC_dy - dA_dy * dC_dx) / (dB_dx * dC_dy - dB_dy * dC_dx)
	end
	for (i,key) in enumerate(keys(model.inparams))
		model.inparams[key] = inparams_values[i]
	end
#	println(DQs)
end

#ShootingRegime

abstract type AbstractShootingRegime <: AbstractRegime end

#necessary fields of ShootingRegime

inparams_fixed(regime::AbstractShootingRegime) = regime.inparams_fixed
exparams(regime::AbstractShootingRegime) = regime.exparams
exparams_symbolic(regime::AbstractShootingRegime) = regime.exparams
inparams_shooting(regime::AbstractShootingRegime) = regime.inparams_shooting
quantities_fixed(regime::AbstractShootingRegime) = regime.quantities_fixed

#necessary methods of ShootingRegime
 
boundary_conditions!(model::AbstractModel, regime::AbstractShootingRegime) = error("Bondary conditions are not defined for $(typeof(regime))")

#realisations of ShootingRegime

#Simple_ShootingRegime
struct Simple_ShootingRegime <: AbstractShootingRegime
	inparams_fixed::Dict{Symbol, Float64}
	exparams::Dict{Symbol, Float64}
	exparams_symbolic::Dict{Symbol, Symbol}
	inparams_shooting::Dict{Symbol, Float64}
	quantities_fixed::Dict{Symbol, Float64}
	residuals::Vector{Float64}
	function Simple_ShootingRegime(inparams_fixed::Dict{Symbol, Float64}, exparams::Dict{Symbol, Float64}, exparams_symbolic::Dict{Symbol, Symbol}, inparams_shooting::Dict{Symbol, Float64}, quantities_fixed::Dict{Symbol, Float64})
		residuals = Vector(undef, length(inparams_shooting))
		return new(inparams_fixed, exparams, exparams_symbolic, inparams_shooting, quantities_fixed, residuals)
	end
end

function initiate_model!(model::AbstractModel, regime::Simple_ShootingRegime)
	for (key, value) in regime.inparams_fixed
		model.inparams[key] =  value
	end
	for (key, value) in regime.inparams_shooting
		model.inparams[key] =  value
	end
	for (key, value) in regime.exparams
		model.exparams[key] =  value
	end
	for (key, value) in regime.exparams_symbolic
		model.exparams_symbolic[key] =  value
	end
end

function calculate_rs_problem!(rs_problem, model, int_params; quantities::Bool, derivatives::Bool)
	rs_solution = calculate_rs_solution(rs_problem, model, int_params)
	if quantities == true
		calculate_quantities(rs_solution, model, save_to_model=true)
	end
	if derivatives == true
		calculate_derivatives!(rs_problem, model, int_params)
	end
end

function calculate_model!(model::AbstractModel, regime::Simple_ShootingRegime, int_params)
	initiate_model!(model, regime)
	rs_vars = get_rs_vars(model)
	rs_problem_init = get_rs_problem(model)

	enumerated_shooting_keys = enumerate(keys(regime.inparams_shooting))

	function update_inparams_shooting()
		if model.inparams[:φc] > 0
			model.inparams[:φc] *= 2.0
		else 
			model.inparams[:φc] = - model.inparams[:φc]
		end
	end

	function is_good_guess()
		return model.derivatives[:dφ∞_dφc] < 0.0
	end

	calculate_rs_problem!(rs_problem_init, model, int_params, quantities=true, derivatives=true)

	while is_good_guess()
		update_inparams_shooting()
		update_rs_vars!(rs_vars, model)
		rs_problem = remake(rs_problem_init; u0=rs_vars)
		calculate_rs_problem!(rs_problem, model, int_params, quantities=true, derivatives=true)
	end

	for (i,key) in enumerated_shooting_keys
		regime.inparams_shooting[key] = model.inparams[key]
	end

	inparams_shooting_values = collect(values(regime.inparams_shooting))

	#inparams_shooting_values = atanh.(2.0 .* inparams_shooting_values .- 1.0)

	function get_residuals!(residuals, inparams_shooting_values)
		#println(inparams_shooting_values)
		for (i,key) in enumerated_shooting_keys
			#model.inparams[key] = 0.5 * (tanh(inparams_shooting_values[i]) + 1.0)
			model.inparams[key] = inparams_shooting_values[i]
		end
		update_rs_vars!(rs_vars, model)
		rs_problem = remake(rs_problem_init; u0=rs_vars)
		calculate_rs_problem!(rs_problem, model, int_params, quantities=true, derivatives=false)
		residuals .= boundary_conditions!(model, regime)
#		println(residuals)
		return residuals
	end

	
	shooting_solution = nlsolve(get_residuals!, inparams_shooting_values, ftol=int_params.reltol*1e1, iterations=100, show_trace=false)
	inparams_shooting_values .= shooting_solution.zero

	for (i,key) in enumerated_shooting_keys
		#regime.inparams_shooting[key] = 0.5 * (tanh(inparams_shooting_values[i]) + 1.0)
		regime.inparams_shooting[key] = inparams_shooting_values[i]
		model.inparams[key] = inparams_shooting_values[i]
	end

	converged = shooting_solution.f_converged
	if converged == false && false
		println(shooting_solution)
		println(regime)
	end

	update_rs_vars!(rs_vars, model)
	rs_problem = remake(rs_problem_init; u0=rs_vars)
	calculate_rs_problem!(rs_problem, model, int_params, quantities=true, derivatives=true)
	return nothing
end

function norm_factor(q, q_fixed)
	if q_fixed == 0.0
		return 1.0
	else
		return q_fixed
	end
end

function boundary_conditions!(model::AbstractModel, regime::Simple_ShootingRegime)
	for (i, (key, value_fixed)) in enumerate(regime.quantities_fixed)
		value_calc = model.quantities[key]
		regime.residuals[i] = (value_calc - value_fixed) / norm_factor(value_calc, value_fixed)
	end
	return regime.residuals
end

