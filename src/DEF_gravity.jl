#realisation for DEF gravity
#old

#coupling functions

abstract type AbstractCouplingFunction end

A_fun(cf::AbstractCouplingFunction, φ::Real) = error("A-function is not defined for $(typeof(cf))")
V_fun(cf::AbstractCouplingFunction, φ::Real) = error("V-function is not defined for $(typeof(cf))")
α_fun(cf::AbstractCouplingFunction, φ::Real) = error("α-function is not defined for $(typeof(cf))")
β_fun(cf::AbstractCouplingFunction, φ::Real) = error("β-function is not defined for $(typeof(cf))")
φ∞_const(cf::AbstractCouplingFunction) = error("φ∞-constant is not defined for $(typeof(cf))")

struct DEF1_CouplingFunction <: AbstractCouplingFunction
	α0::Float64
	β0::Float64
	function DEF1_CouplingFunction(;α0, β0)
		return new(α0, β0)
	end
end

DEF1_CouplingFunction(exparams::Dict{Symbol}) = DEF1_CouplingFunction(α0 = exparams[:α0], β0 = exparams[:β0])

function Base.show(io::IO, cf::DEF1_CouplingFunction)
	println(io, "DEF1 coupling function")
	println(io, "  :α0", " => ", cf.α0)
	println(io, "  :β0", " => ", cf.β0)
	println(io, "  :φ∞", " => ", φ∞_const(cf))
	return nothing
end

A_fun(cf::DEF1_CouplingFunction, φ::Real) = exp(cf.α0 * φ + 0.5 * cf.β0 * φ^2)
V_fun(cf::DEF1_CouplingFunction, φ::Real) = cf.α0 * φ + 0.5 * cf.β0 * φ^2
α_fun(cf::DEF1_CouplingFunction, φ::Real) = cf.α0 + cf.β0 * φ
β_fun(cf::DEF1_CouplingFunction, φ::Real) = cf.β0
φ∞_const(cf::DEF1_CouplingFunction) = 0.0
get_params(cf::DEF1_CouplingFunction) = Dict(:α0 => cf.α0, :β0 => cf.β0)

struct DEF3_CouplingFunction <: AbstractCouplingFunction
	α0::Float64
	β0::Float64
	function DEF3_CouplingFunction(;α0, β0)
		return new(α0, β0)
	end
end

DEF3_CouplingFunction(exparams::Dict{Symbol}) = DEF3_CouplingFunction(α0 = exparams[:α0], β0 = exparams[:β0])

function Base.show(io::IO, cf::DEF3_CouplingFunction)
	println(io, "DEF3 coupling function")
	println(io, "  :α0", " => ", cf.α0)
	println(io, "  :β0", " => ", cf.β0)
	println(io, "  :φ∞", " => ", φ∞_const(cf))
	return nothing
end

A_fun(cf::DEF3_CouplingFunction, φ::Real) = exp(0.5 * cf.β0 * φ^2)
V_fun(cf::DEF3_CouplingFunction, φ::Real) = 0.5 * cf.β0 * φ^2
α_fun(cf::DEF3_CouplingFunction, φ::Real) = cf.β0 * φ
β_fun(cf::DEF3_CouplingFunction, φ::Real) = cf.β0
φ∞_const(cf::DEF3_CouplingFunction) = cf.β0 != 0.0 ? cf.α0 / cf.β0 : 0.0
get_params(cf::DEF3_CouplingFunction) = Dict(:α0 => cf.α0, :β0 => cf.β0)


#DEF_Model{T}

struct DEF_Model{T <: AbstractCouplingFunction, P <: AbstractEoS} <: AbstractModel
	inparams::Dict{Symbol,Float64}
	exparams::Dict{Symbol,Q} where Q
	cf::T
	eos::P
	quantities::Dict{Symbol,Float64}
	derivatives::Dict{Symbol,Float64}
	der_rules::Dict{Symbol,Tuple{Symbol,Symbol,Symbol,Int64}}
	n_odes::Int64
	rspan::Tuple{Float64,Float64}
	vars_init::Vector{Float64}
end

	function DEF_Model{T, P}(cf::T, eos::P) where {T <: AbstractCouplingFunction, P <: AbstractEoS}
		inparams = Dict{Symbol,Float64}(:φc => NaN, :pc => NaN)
		exparams = merge(get_params(cf), get_params(eos))
		quantities = Dict{Symbol,Float64}()
		derivatives = Dict{Symbol,Float64}()
		der_rules = Dict(
			:αA => (:lnmA, :φ∞, :m̃A, 1),
			:βA => (:αA, :φ∞, :m̃A, 1),
			:kA => (:lnIA, :φ∞, :m̃A, 1),
			:dm_dR => (:m̃A, :R, :φc, 1)
			)
		n_odes = 8
		rspan = (0.0, 1e10)
		vars_init = Vector{Float64}(undef, n_odes)
		return DEF_Model(inparams, exparams, cf, eos, quantities, derivatives, der_rules, n_odes, rspan, vars_init)
	end

DEF_Model(cf::T, eos::P) where {T <: AbstractCouplingFunction, P <: AbstractEoS} = DEF_Model{T, P}(cf, eos)

function DEF_Model{T, P}(exparams::Dict{Symbol}) where {T <: AbstractCouplingFunction, P <: AbstractEoS}
	cf = T(exparams)
	eos = P(exparams)
	return DEF_Model(cf, eos)
end

#DEF_Model(sample_model::DEF_Model, exparams) = DEF_Model(typeof(model.cf)(exparams), typeof(model.eos)(exparams))

similar(model::DEF_Model) = DEF_Model(model.cf, model.eos)

function Base.show(io::IO, model::DEF_Model)
	println(io, "DEF gravity model")
	println(io, "Number of ODEs:  ", model.n_odes)
	println(io, "Inner parameters:")
	print_dict(io, model.inparams)
	print(io, model.cf)
	println(io, "Equation of state:\n  ", model.eos)
	println(io, "Calculated quantities")
	print_dict(io, model.quantities)
	if !isempty(model.derivatives)
		println(io, "Calculated derivatives")
		print_dict(io, model.derivatives)
		println(io, "Convergence checks:")
		println(io, "  First order derivative relative αA / αA_der - 1: ", model.quantities[:αA] / model.derivatives[:αA] - 1.0)
		println(io, "  First order derivative absolute αA - αA_der: ", model.quantities[:αA] - model.derivatives[:αA])
		println(io, "  Surface accuracy Dν_s_theory / Dν_s_real - 1: ", model.quantities[:Dν_s] / model.quantities[:Dν_s_check] - 1.0)
	end
	return nothing
end

#function set_vars_init!(model::DEF_Model)
#	model.vars_init[1] = 0.0
#	model.vars_init[2] = 0.0
#	model.vars_init[3] = model.inparams[:φc]
#	model.vars_init[4] = 0.0
#	model.vars_init[5] = model.inparams[:pc]
#	model.vars_init[6] = 0.0
#	model.vars_init[7] = 1.0
#	model.vars_init[8] = 0.0
#	return nothing
#end

function get_rs_vars(inparams::Dict, model::DEF_Model)
	rs_vars = zeros(valtype(inparams),model.n_odes)
	rs_vars[3] = inparams[:φc]
	rs_vars[5] = inparams[:pc]
	rs_vars[7] = 1.0
	return rs_vars
end

function update_rs_vars!(rs_vars::Vector, inparams::Dict, model::DEF_Model)
	if valtype(rs_vars) != valtype(inparams)
		error("The type of rs_vars is inconsistent with the type of inparams")
	end
	rs_vars[3] = inparams[:φc]
	rs_vars[5] = inparams[:pc]
	return nothing
end

function ode_system!(DVars::Vector, Vars::Vector, model::DEF_Model, r::Float64)
	μ, ν, φ, ψ, p̃, M̃, ω, ϖ = Vars

    if p̃ < 0 || (r - 2.0*μ) < 0
        DVars[1:8] .= 0.0
        return DVars
    end

    A, α, ñ, ε̃ = internal_physics(Vars, model, r)
    κ = 4π * G / c^4

    DVars[1] = Dμ = κ * r^2 * A^4 * ε̃ + 0.5 * r * (r - 2.0*μ) * ψ^2
    DVars[2] = Dν = 2κ * r^2 * A^4 * p̃ / (r - 2.0*μ) + r * ψ^2 + 2.0 * μ / (r * (r - 2.0*μ))
    DVars[3] = Dφ = ψ
    DVars[4] = Dψ = κ * r * A^4 / (r - 2.0*μ) * (α * (ε̃ - 3.0*p̃) + r * ψ * (ε̃ - p̃)) - 2.0 * (r - μ) * ψ / (r * (r - 2.0*μ))
#    DVars[5] = Dp̃ = - (ε̃ + p̃) * (4π * G / c^4 * r^2 * A^4 * p̃ / (r - 2.0*μ) + 0.5 * r * ψ^2 + μ / (r * (r - 2.0*μ)) + α * ψ)
	DVars[5] = Dp̃ = - (ε̃ + p̃) * (0.5*Dν + α * ψ)
    DVars[6] = DM̃ = 4π * mp * ñ * A^3 * r^2 * (1.0 - 2.0 * μ / r)^-0.5
    DVars[7] = Dω = ϖ
    DVars[8] = Dϖ = κ * r^2 / (r - 2.0*μ) * A^4 * (ε̃ + p̃) * (ϖ + 4ω / r) + (ψ^2 * r - 4.0 / r) * ϖ

    if r == 0.0
        DVars[1:8] .= 0.0
        DVars[4] = Dψ = 1/3 * κ * A^4 * α * (ε̃ - 3.0*p̃)
        DVars[8] = Dϖ = 4/5 * κ * A^4 * (ε̃ + p̃) * ω
    end

#	println("r=",rnd(r))
#    println("μ=",rnd(μ) ," ν=", rnd(ν)," φ=", rnd(φ)," ψ=", rnd(ψ)," p̃=", rnd(p̃)," M̃=", rnd(M̃))
#
#    println("Dμ=",rnd(DVars[1]) ," Dν=", rnd(DVars[2])," Dφ=", rnd(DVars[3])," Dψ=", rnd(DVars[4])," Dp̃=", rnd(DVars[5])," DM̃=", rnd(DVars[6]))

    return DVars
end

function internal_physics(Vars::Vector, model::DEF_Model, r::Real)
    μ, ν, φ, ψ, p̃, M̃, ω, ϖ = Vars

    A = A_fun(model.cf, φ)
    α = α_fun(model.cf, φ)

#    A = A < 1e-30 ? 1e-30 : A
#    A = A > 1e30 ? 1e30 : A
#    p̃ = p̃ > 0.0 ? p̃ : 0.0
#    p̃ = p̃ < 1e50 ? p̃ : 1e50

    n = get_number_density(model.eos, p̃; from=:pressure, units=:cgs)
    ε = get_energy_density(model.eos, p̃; from=:pressure, units=:cgs)

    return A, α, n, ε
end

function calculate_quantities_old!(rs_solution, model::DEF_Model)
	quantities = model.quantities

    μ_c, ν_c, φ_c, ψ_c, p̃_c, M̃_c, ω_c, ϖ_c = Vars_c = rs_solution[1]
    r_c = rs_solution.t[1]

    quantities[:p̃_c] = p̃_c

    A_c, α_c, ñ_c, ε̃_c = internal_physics(Vars_c, model, r_c)

    μ_s, ν_s, φ_s, ψ_s, p̃_s, M̃_s, ω_s, ϖ_s = Vars_s = rs_solution[end]


	quantities[:ψ_s] = ψ_s
    quantities[:μ_s] = μ_s / (4π * G / c^2)

    r_s = rs_solution.t[end]

    A_s, α_s, ñ_s, ε̃_s = internal_physics(Vars_s, model, r_s) 

    quantities[:R]    = R    = r_s
    quantities[:Dν_s] = Dν_s = R * ψ_s^2 + 2.0 * μ_s / (R * (R - 2.0*μ_s))
    quantities[:Dν_s_check] = Dν_s_check = 8π * G / c^4 * r_s^2 * A_s^4 * p̃_s / (r_s - 2.0*μ_s) + r_s * ψ_s^2 + 2.0 * μ_s / (r_s * (r_s - 2.0*μ_s))
    quantities[:αA]   = αA   = 2.0 * ψ_s / Dν_s
    quantities[:Q1]   = Q1   = sqrt(1.0 + αA^2)
    quantities[:Q2]   = Q2   = sqrt(1.0 - 2.0 * μ_s / R)
    quantities[:ν̂_s]  = ν̂_s  = - 2.0 / Q1 * atanh(Q1 / (1.0 + 2.0 / (R * Dν_s)))
    quantities[:φ0]   = φ0   = φ_s - 0.5 * αA * ν̂_s
    quantities[:φ∞]   = φ∞   = φ_s - 0.5 * αA * ν̂_s
    quantities[:mA]   = mA   = c^2 / G * 0.5 * Dν_s * R^2 * Q2 * exp(0.5 * ν̂_s)
    quantities[:m̃A]   = m̃A   = M̃_s
    quantities[:JA]   = JA   = c^2 / G * 1/6 * ϖ_s * R^4 * Q2 * exp(-0.5 * ν̂_s)
    quantities[:Ω]    = Ω    = ω_s - c^4 / G^2 * 3.0 * JA / (4.0 * mA^3 * (3 - αA^2)) * (exp(2ν̂_s) - 1.0 + 4G * mA / (R * c^2) * exp(ν̂_s) * (2G * mA / (R * c^2) + exp(0.5 * ν̂_s) * cosh(0.5 * Q1 * ν̂_s)))
    quantities[:IA]   = IA   = JA / Ω

    quantities[:lnIA] = lnIA = IA > 0 ? log(IA) : NaN
    quantities[:lnmA] = lnmA = log(mA)
    quantities[:lnm̃A] = lnm̃A = log(m̃A)

    quantities[:pc] = pc = p̃_c
    quantities[:φc] = φc = φ_c  

    quantities[:p_s] = p̃_s 
    quantities[:Dp_s] = Dp_s =  - (ε̃_s + p̃_s) * (0.5*Dν_s_check + α_s * ψ_s)

    quantities[:bc_φ∞] = bc_φ∞ = φ∞ - φ∞_const(model.cf)


#    φ_max = maximum(rs_solution[3,:])
#    φ_min = minimum(rs_solution[3,:])
#    φ_span = φ_max - φ_min

#    quantities_names   = [:R, :Dν_s, :αA, :Q1, :Q2, :ν̂_s, :φ0, :mA, :m̃A, :p̃_c, :φ_c, :φ_max, :φ_span]
#    quantities_values = [R, Dν_s,  αA,  Q1,  Q2,  ν̂_s,  φ0,  mA,  m̃A, p̃_c, φ_c, φ_max, φ_span]
#    
#    for (i, quantity) in enumerate(quantities_names)
#    	model.quantities[quantity] = quantities_values[i]
#    end

    return nothing
end

function calculate_quantities(rs_solution, model::DEF_Model; save_to_model::Bool)

    μ_c, ν_c, φ_c, ψ_c, p̃_c, M̃_c, ω_c, ϖ_c = Vars_c = rs_solution[1]
    r_c = rs_solution.t[1]

    A_c, α_c, ñ_c, ε̃_c = internal_physics(Vars_c, model, r_c)

    μ_s, ν_s, φ_s, ψ_s, p̃_s, M̃_s, ω_s, ϖ_s = Vars_s = rs_solution[end]
    r_s = rs_solution.t[end]

    A_s, α_s, ñ_s, ε̃_s = internal_physics(Vars_s, model, r_s) 

    R    = r_s
    Dν_s = R * ψ_s^2 + 2.0 * μ_s / (R * (R - 2.0*μ_s))
    Dν_s_check = 8π * G / c^4 * r_s^2 * A_s^4 * p̃_s / (r_s - 2.0*μ_s) + r_s * ψ_s^2 + 2.0 * μ_s / (r_s * (r_s - 2.0*μ_s))
    αA   = 2.0 * ψ_s / Dν_s
    Q1   = sqrt(1.0 + αA^2)
    Q2   = sqrt(1.0 - 2.0 * μ_s / R)
    ν̂_s  = - 2.0 / Q1 * atanh(Q1 / (1.0 + 2.0 / (R * Dν_s)))
    φ0   = φ_s - 0.5 * αA * ν̂_s
    φ∞   = φ_s - 0.5 * αA * ν̂_s
    mA   = c^2 / G * 0.5 * Dν_s * R^2 * Q2 * exp(0.5 * ν̂_s)
    m̃A   = M̃_s
    JA   = c^2 / G * 1/6 * ϖ_s * R^4 * Q2 * exp(-0.5 * ν̂_s)
    Ω    = ω_s - c^4 / G^2 * 3.0 * JA / (4.0 * mA^3 * (3 - αA^2)) * (exp(2ν̂_s) - 1.0 + 4G * mA / (R * c^2) * exp(ν̂_s) * (2G * mA / (R * c^2) + exp(0.5 * ν̂_s) * cosh(0.5 * Q1 * ν̂_s)))
    IA   = JA / Ω

    lnIA = IA > 0 ? log(IA) : NaN
    lnmA = log(mA)
    lnm̃A = log(m̃A)

    pc = p̃_c
    φc = φ_c  

    Dp_s =  - (ε̃_s + p̃_s) * (0.5*Dν_s_check + α_s * ψ_s)

    bc_φ∞ = φ∞ - φ∞_const(model.cf)

    quantities_values = (R, Dν_s, Dν_s_check, αA, Q1, Q2, ν̂_s, φ0, φ∞, mA, m̃A, JA, Ω, IA, lnIA, lnmA, lnm̃A, pc, φc, bc_φ∞)
    quantities_names  = (:R, :Dν_s, :Dν_s_check, :αA, :Q1, :Q2, :ν̂_s, :φ0, :φ∞, :mA, :m̃A, :JA, :Ω, :IA, :lnIA, :lnmA, :lnm̃A, :pc, :φc, :bc_φ∞)

    if save_to_model == true
    	for i in 1:length(quantities_names)
    		model.quantities[quantities_names[i]] = quantities_values[i]
    	end
    end

    return quantities_values, quantities_names
end

function surface_callback(model::DEF_Model)
	function surface_condition(Vars, r, integrator)
    	p̃ = Vars[5]
		p̃_min = 1e15
		return p̃ < p̃_min
	end
	affect!(integrator) = terminate!(integrator)
	return DiscreteCallback(surface_condition, affect!)
end









