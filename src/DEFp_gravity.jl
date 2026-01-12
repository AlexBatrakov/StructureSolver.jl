#realisation for DEF gravity

#coupling functions

abstract type AbstractCouplingFunction end

raw"""
    A_fun(cf, exparams, φ)

Coupling function $A(\varphi)$ used in DEF-type scalar-tensor gravity.
"""
A_fun(cf::AbstractCouplingFunction, exparams::Dict, φ::Real) = error("A-function is not defined for $(typeof(cf))")

"""
    V_fun(cf, exparams, φ)

Convenience potential-like function (implementation-specific) associated with the coupling.
"""
V_fun(cf::AbstractCouplingFunction, exparams::Dict, φ::Real) = error("V-function is not defined for $(typeof(cf))")

raw"""
    α_fun(cf, exparams, φ)

Coupling strength $\alpha(\varphi) = d\ln A / d\varphi$.
"""
α_fun(cf::AbstractCouplingFunction, exparams::Dict, φ::Real) = error("α-function is not defined for $(typeof(cf))")

raw"""
    β_fun(cf, exparams, φ)

Second coupling derivative (implementation-dependent; often $\beta(\varphi) = d\alpha/d\varphi$).
"""
β_fun(cf::AbstractCouplingFunction, exparams::Dict, φ::Real) = error("β-function is not defined for $(typeof(cf))")

"""
    φ∞_const(cf, exparams)

Asymptotic scalar-field value used as a reference for boundary conditions.
"""
φ∞_const(cf::AbstractCouplingFunction, exparams::Dict) = error("φ∞-constant is not defined for $(typeof(cf))")

"""DEF1 coupling function.

Typically parameterized by external parameters `:α0` and `:β0`.
"""
struct DEF1_CouplingFunction <: AbstractCouplingFunction end

A_fun(cf::DEF1_CouplingFunction, exparams::Dict, φ::Real) = exp(exparams[:α0] * φ + 0.5 * exparams[:β0] * φ^2)
V_fun(cf::DEF1_CouplingFunction, exparams::Dict, φ::Real) = exparams[:α0] * φ + 0.5 * exparams[:β0] * φ^2
α_fun(cf::DEF1_CouplingFunction, exparams::Dict, φ::Real) = exparams[:α0] + exparams[:β0] * φ
β_fun(cf::DEF1_CouplingFunction, exparams::Dict, φ::Real) = exparams[:β0]
φ∞_const(cf::DEF1_CouplingFunction, exparams::Dict) = 0.0
get_params(cf::DEF1_CouplingFunction) = Dict(:α0 => NaN, :β0 => NaN)

"""DEF3 coupling function.

Typically parameterized by external parameters `:α0` and `:β0`.
"""
struct DEF3_CouplingFunction <: AbstractCouplingFunction end

A_fun(cf::DEF3_CouplingFunction, exparams::Dict, φ::Real) = exp(0.5 * exparams[:β0] * φ^2)
V_fun(cf::DEF3_CouplingFunction, exparams::Dict, φ::Real) = 0.5 * exparams[:β0] * φ^2
α_fun(cf::DEF3_CouplingFunction, exparams::Dict, φ::Real) = exparams[:β0] * φ
β_fun(cf::DEF3_CouplingFunction, exparams::Dict, φ::Real) = exparams[:β0]
φ∞_const(cf::DEF3_CouplingFunction, exparams::Dict) = exparams[:β0] != 0.0 ? exparams[:α0] / exparams[:β0] : 0.0
get_params(cf::DEF3_CouplingFunction) = Dict(:α0 => NaN, :β0 => NaN)


"""
    DEFp_Model{T}(cf, eos)

Damour-Esposito-Farese scalar-tensor stellar-structure model.

`cf` is a coupling-function instance (e.g. `DEF1_CouplingFunction()`), and `eos` is an equation of
state (e.g. `PWP_EoS(...)`). The model stores

- inner parameters `inparams` (e.g. `:pc`, `:φc`),
- external parameters `exparams` / `exparams_symbolic`,
- computed global `quantities` and `derivatives`.
"""
struct DEFp_Model{T1 <: Real, T2 <: AbstractCouplingFunction, T3 <: AbstractEoS} <: AbstractModel
    inparams::Dict{Symbol,T1}
    exparams::Dict{Symbol,T1}
    exparams_symbolic::Dict{Symbol,Symbol}
    cf::T2
    eos::T3
    quantities::Dict{Symbol,T1}
    derivatives::Dict{Symbol,T1}
    der_rules::Dict{Symbol,Tuple{Symbol,Symbol,Symbol,Int64}}
    n_odes::Int64
    rs_span::Tuple{Float64,Float64}
end
function DEFp_Model{T1}(cf::T2, eos::T3) where {T1 <: Real, T2 <: AbstractCouplingFunction, T3 <: AbstractEoS}
    inparams = Dict{Symbol,T1}(:φc => NaN, :pc => NaN)
    exparams_all = merge(get_params(cf), get_params(eos))
    exparams = Dict{Symbol,T1}()
    exparams_symbolic = Dict{Symbol,Symbol}()
    for (exparam_name, exparam_value) in exparams_all
        if typeof(exparam_value) <: Real
            exparams[exparam_name] = exparam_value
        else
            exparams_symbolic[exparam_name] = exparam_value
        end
    end
    quantities = Dict{Symbol,T1}()
    derivatives = Dict{Symbol,T1}()
    der_rules = Dict(
        :αA => (:lnmA, :φ∞, :m̃A, 1),
        :βA => (:αA, :φ∞, :m̃A, 1),
        :kA => (:lnIA, :φ∞, :m̃A, 1),
        :dφ∞_dφc => (:φ∞, :φc, :pc, 1)
        #:dlnφc_dlnpc => (:lnφc, :lnpc, :φ∞, 1),
        #:dφc_dpc => (:φc, :pc, :φ∞, 1)
        #:dm_dR => (:mA, :R, :φ∞, 1),
        #:dR_dp => (:R, :pc, :φ∞, 1),
        )
    n_odes = 8
    rs_span = (1.0, 1e-13)
    return DEFp_Model(inparams, exparams, exparams_symbolic, cf, eos, quantities, derivatives, der_rules, n_odes, rs_span)
end

function DEFp_Model{T1, T2, T3}(exparams::Dict{Symbol}, exparams_symbolic::Dict{Symbol,Symbol}) where {T1 <: Real, T2 <: AbstractCouplingFunction, T3 <: AbstractEoS}
    cf = T2()
    eos = T3(exparams, exparams_symbolic)
    return DEFp_Model{T1}(cf, eos)
end

function DEFp_Model{T1, T2, T3}() where {T1 <: Real, T2 <: AbstractCouplingFunction, T3 <: AbstractEoS}
    cf = T2()
    eos = T3()
    return DEFp_Model{T1}(cf, eos)
end

function is_the_same_model(model::DEFp_Model, exparams::Dict{Symbol}, exparams_symbolic::Dict{Symbol,Symbol})
    eos_params = get_params(model.eos)
    answer = true
    for (key, value) in eos_params
        if typeof(key) == Symbol
            value_new = exparams_symbolic[key]
        else
            value_new = exparams[key]
        end
        if value_new != value
            answer = false
        end
    end
    return answer
end

function update_model!(model::DEFp_Model, exparams::Dict{Symbol}, exparams_symbolic::Dict{Symbol,Symbol})
#    if is_the_same_model(model, exparams, exparams_symbolic)
#        model = model
#    else
#        model = typeof(model)(exparams, exparams_symbolic)
#    end
    for (key, value) in exparams
        model.exparams[key] = value
    end
    for (key, value) in exparams_symbolic
        model.exparams_symbolic[key] = value
    end
    return model
end

function convert_model(T1::DataType, model::DEFp_Model)
    model_new = DEFp_Model{T1}(model.cf, model.eos)
    for (key, value) in model.inparams
        model_new.inparams[key] = value
    end
    for (key, value) in model.exparams
        model_new.exparams[key] = value
    end
    return model_new
end

#DEFp_Model{T1}(cf::T2, eos::T3) where {T1 <: Real, T2 <: AbstractCouplingFunction, T3 <: AbstractEoS} = DEFp_Model{T1, T2, T3}(cf, eos)

#function DEFp_Model{T, P}(exparams::Dict{Symbol}) where {T <: AbstractCouplingFunction, P <: AbstractEoS}
#   cf = T(exparams)
#   eos = P(exparams)
#   return DEFp_Model(cf, eos)
#end

#DEFp_Model(sample_model::DEFp_Model, exparams) = DEFp_Model(typeof(model.cf)(exparams), typeof(model.eos)(exparams))

#similar(model::DEFp_Model) = DEFp_Model(model.cf, model.eos)

function Base.show(io::IO, model::DEFp_Model)
    println(io, "DEF gravity model")
    println(io, "Numerical type:  ", typeof(model).parameters[1])
    println(io, "Number of ODEs:  ", model.n_odes)
    println(io, "Inner parameters:")
    print_dict(io, model.inparams)
    println(io, "External parameters:")
    print_dict(io, model.exparams)
    println(io, "External parameters symbolic:")
    print_dict(io, model.exparams_symbolic)
    print(io, "Coupling function:\n  ", model.cf)
    println(io, "\nEquation of state:\n  ", model.eos)
    println(io, "Calculated quantities")
    print_dict(io, model.quantities)
    if !isempty(model.derivatives)
        println(io, "Calculated derivatives")
        print_dict(io, model.derivatives)
        println(io, "Convergence checks:")
        println(io, "  First order derivative relative αA / αA_der - 1: ", model.quantities[:αA] / model.derivatives[:αA] - 1.0)
        println(io, "  First order derivative absolute αA - αA_der: ", model.quantities[:αA] - model.derivatives[:αA])
#        println(io, "  Surface accuracy Dν_s_theory / Dν_s_real - 1: ", model.quantities[:Dν_s] / model.quantities[:Dν_s_check] - 1.0)
    end
    return nothing
end

function get_rs_vars(model::DEFp_Model)
    inparams_type = typeof(model).parameters[1]
    rs_vars = zeros(inparams_type, model.n_odes)
    update_rs_vars!(rs_vars, model)
    return rs_vars
end

function update_rs_vars!(rs_vars::Vector, model::DEFp_Model)
    inparams_type = typeof(model).parameters[1]
    if eltype(rs_vars) != typeof(model.inparams[:φc])
        error("The type of rs_vars is inconsistent with the type of inparams")
    end
    p̃_c = model.inparams[:pc]
    r_min = 1e-4
    A_c, α_c, ñ_c, ε̃_c = internal_physics(rs_vars, model, 1.0)
    rs_vars[1] = μ_c = 0.0
    rs_vars[2] = ν_c = 0.0
    rs_vars[3] = φ_c = model.inparams[:φc]
    rs_vars[4] = ψ_c = 1/3 * r_min * 4π * G / c^4 * A_c^4 * α_c * (ε̃_c - 3.0*p̃_c)
    rs_vars[5] = r_c = r_min
    rs_vars[6] = M̃_c = 0.0
    rs_vars[7] = ω_c = 1.0
    rs_vars[8] = ϖ_c = 4/5 * r_min * 4π * G / c^4 * A_c^4 * (ε̃_c + p̃_c) * ω_c
    return nothing
end

"""
    ode_system!(D_rs_vars::Vector, rs_vars::Vector, model::DEFp_Model, q::Real)

Right-hand side for the stellar structure ODE system in DEF scalar-tensor gravity.

`q` is a dimensionless pressure coordinate used internally by the solver (with physical pressure
`p̃ = q * pc`, where `pc = model.inparams[:pc]`).
"""
function ode_system!(D_rs_vars::Vector, rs_vars::Vector, model::DEFp_Model, q::Real)
    μ, ν, φ, ψ, r, M̃, ω, ϖ = rs_vars
    p̃_c = model.inparams[:pc]
    p̃ = q * p̃_c

    if (r - 2.0*μ) < 0.0 || 1.0 - 2.0 * μ / r < 0.0
        D_rs_vars[1:8] .= 0.0
        return D_rs_vars
    end

    A, α, ñ, ε̃ = internal_physics(rs_vars, model, q)

    κ = 4π * G / c^4

    Dμ_Dr = κ * r^2 * A^4 * ε̃ + 0.5 * r * (r - 2.0*μ) * ψ^2
    Dν_Dr = 2κ * r^2 * A^4 * p̃ / (r - 2.0*μ) + r * ψ^2 + 2.0 * μ / (r * (r - 2.0*μ))
    Dφ_Dr = ψ
    Dψ_Dr = κ * r * A^4 / (r - 2.0*μ) * (α * (ε̃ - 3.0*p̃) + r * ψ * (ε̃ - p̃)) - 2.0 * (r - μ) * ψ / (r * (r - 2.0*μ))
    Dp̃_Dr = - (ε̃ + p̃) * (0.5*Dν_Dr + α * ψ)
    DM̃_Dr = 4π * mp * ñ * A^3 * r^2 * (1.0 - 2.0 * μ / r)^-0.5
    Dω_Dr = ϖ
    Dϖ_Dr = κ * r^2 / (r - 2.0*μ) * A^4 * (ε̃ + p̃) * (ϖ + 4ω / r) + (ψ^2 * r - 4.0 / r) * ϖ

    Dr_Dq = p̃_c / Dp̃_Dr

    D_rs_vars[1] = Dμ_Dr * Dr_Dq
    D_rs_vars[2] = Dν_Dr * Dr_Dq
    D_rs_vars[3] = Dφ_Dr * Dr_Dq
    D_rs_vars[4] = Dψ_Dr * Dr_Dq
    D_rs_vars[5] = Dr_Dq
    D_rs_vars[6] = DM̃_Dr * Dr_Dq
    D_rs_vars[7] = Dω_Dr * Dr_Dq
    D_rs_vars[8] = Dϖ_Dr * Dr_Dq


    if Dp̃_Dr > 0.0
        D_rs_vars[1:8] .= 0.0
        return D_rs_vars
    end

#    println(q, " ", A, " ", (r - 2.0*μ), " ", Dr_Dq)

    return D_rs_vars
end

"""
    internal_physics(rs_vars::Vector, model::DEFp_Model, q::Real)

Compute local matter quantities used by the ODE system.

Returns `(A, α, ñ, ε̃)` evaluated at the current state.
"""
function internal_physics(rs_vars::Vector, model::DEFp_Model, q::Real)
    μ, ν, φ, ψ, r, M̃, ω, ϖ = rs_vars
    p̃_c = model.inparams[:pc]
#    p̃_c = 1e35
    p̃ = q * p̃_c

    A = A_fun(model.cf, model.exparams, φ)
    α = α_fun(model.cf, model.exparams, φ)
#    A = 1.0 
#    α = 0.0

#    A = A < 1e-30 ? 1e-30 : A
#    A = A > 1e30 ? 1e30 : A
#    p̃ = p̃ > 0.0 ? p̃ : 0.0
#    p̃ = p̃ < 1e50 ? p̃ : 1e50

    n = get_number_density(model.eos, p̃; from=:pressure, units=:cgs)
    ε = get_energy_density(model.eos, p̃; from=:pressure, units=:cgs)

    return A, α, n, ε
end

"""
    calculate_quantities(rs_solution, model::DEFp_Model; save_to_model::Bool)

Post-process an ODE solution and compute global quantities/diagnostics.

If `save_to_model=true`, results are written into `model.quantities`.
"""
function calculate_quantities(rs_solution, model::DEFp_Model; save_to_model::Bool)

    μ_c, ν_c, φ_c, ψ_c, r_c, M̃_c, ω_c, ϖ_c = rs_vars_c = rs_solution[1]
    q_c = rs_solution.t[1]
    p̃_c = model.inparams[:pc]

    A_c, α_c, ñ_c, ε̃_c = internal_physics(rs_vars_c, model, q_c)

    μ_s, ν_s, φ_s, ψ_s, r_s, M̃_s, ω_s, ϖ_s = rs_vars_s = rs_solution[end]
    q_s = rs_solution.t[end]
    p̃_s = q_s * p̃_c


    A_s, α_s, ñ_s, ε̃_s = internal_physics(rs_vars_s, model, q_s) 

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
    lnmA = mA > 0 ? log(mA) : NaN
    lnm̃A = m̃A > 0 ? log(m̃A) : NaN

    pc = p̃_c
    φc = φ_c  
    lnpc = pc > 0 ? log(pc) : NaN
    lnφc = φc > 0 ? log(φc) : NaN
    ps = p̃_s

    Dp_s =  - (ε̃_s + p̃_s) * (0.5*Dν_s_check + α_s * ψ_s)
    test = r_s - 2.0*μ_s
    μ_s = μ_s

#    bc_φ∞ = abs(φ∞ - φ∞_const(model.cf, model.exparams)) + abs(φc > 0.0 ? 0.0 : φc)
    bc_φ∞ = φ∞ - φ∞_const(model.cf, model.exparams)

    cs_c = get_sound_velocity(model.eos, p̃_c; from=:pressure, units=:cgs)

    #quantities_values = (R, Dν_s, Dν_s_check, αA, Q1, Q2, ν̂_s, φ0, φ∞, mA, m̃A, JA, Ω, IA, lnIA, lnmA, lnm̃A, pc, φc, bc_φ∞, ps, lnpc, lnφc, Dp_s, test, μ_s)
    #quantities_names  = (:R, :Dν_s, :Dν_s_check, :αA, :Q1, :Q2, :ν̂_s, :φ0, :φ∞, :mA, :m̃A, :JA, :Ω, :IA, :lnIA, :lnmA, :lnm̃A, :pc, :φc, :bc_φ∞, :ps, :lnpc, :lnφc, :Dp_s, :test, :μ_s)
    quantities_values = (R, αA, φ∞, mA, m̃A, JA, Ω, IA, lnIA, lnmA, lnm̃A, pc, φc, bc_φ∞, cs_c)
    quantities_names  = (:R, :αA, :φ∞, :mA, :m̃A, :JA, :Ω, :IA, :lnIA, :lnmA, :lnm̃A, :pc, :φc, :bc_φ∞, :cs_c)

    if save_to_model == true
        for i in 1:length(quantities_names)
            model.quantities[quantities_names[i]] = quantities_values[i]
        end
    end

    return quantities_values, quantities_names
end

surface_callback(model::DEFp_Model) = nothing









