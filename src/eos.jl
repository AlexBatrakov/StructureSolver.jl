abstract type AbstractEoS end

raw"""
    get_density(eos, val; from, units=:cgs)

Convert an input `val` to mass density $\rho$.

`from` specifies what `val` represents (e.g. `:pressure`, `:density`, `:number_density`).
`units` controls unit conventions; most implementations support `:cgs` and sometimes `:natural`.
"""
get_density(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Density method is not implemented for $(typeof(eos))")

raw"""
    get_number_density(eos, val; from, units=:cgs)

Convert an input `val` to baryon number density $n$.

`from` specifies what `val` represents.
"""
get_number_density(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Number density method is not implemented for $(typeof(eos))")

"""
    get_pressure(eos, val; from, units=:cgs)

Return pressure corresponding to `val`.

`from` specifies what `val` represents. In CGS, pressure is typically returned in dyn/cm^2.
"""
get_pressure(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Pressure method is not implemented for $(typeof(eos))")

"""
    get_energy_density(eos, val; from, units=:cgs)

Return energy density corresponding to `val`.

In CGS, this is typically returned in erg/cm^3.
"""
get_energy_density(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Energy density method is not implemented for $(typeof(eos))")

"""
    get_specific_enthalpy(eos, val; from, units=:cgs)

Return specific enthalpy (implementation-dependent normalization) corresponding to `val`.
"""
get_specific_enthalpy(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Specific enthalpy method is not implemented for $(typeof(eos))")

"""
    get_sound_velocity(eos, val; from, units=:cgs)

Return the (adiabatic) sound speed corresponding to `val`.

In CGS, returned in cm/s.
"""
get_sound_velocity(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Sound velocity method is not implemented for $(typeof(eos))")

"""
    get_internal_energy(eos, val; from, units=:cgs)

Return internal energy per unit mass (implementation-dependent normalization) corresponding to `val`.
"""
get_internal_energy(eos::AbstractEoS, val::Real; from::Symbol, units=:cgs) = error("Internal energy method is not implemented for $(typeof(eos))")

get_pressure(eos::AbstractEoS, val_arr::Vector{T}; from::Symbol, units=:cgs) where {T<:Real} = [get_pressure(eos, val, from=from, units=units) for val in val_arr]
get_energy_density(eos::AbstractEoS, val_arr::Vector{T}; from::Symbol, units=:cgs) where {T<:Real} = [get_energy_density(eos, val, from=from, units=units) for val in val_arr]
get_sound_velocity(eos::AbstractEoS, val_arr::Vector{T}; from::Symbol, units=:cgs) where {T<:Real} = [get_sound_velocity(eos, val, from=from, units=units) for val in val_arr]

#const get_funs = Dict(
#:pressure => :get_pressure,
#:number_density => :get_number_density,
#:density => :get_density,
#:energy_density => :get_energy_density,
#:specific_enthalpy => :get_specific_enthalpy,
#:sound_velocity => :get_sound_velocity,
#:internal_energy => :get_internal_energy
#)

#function call_eos(eos::AbstractEoS, val::Float64; what::Symbol, from::Symbol, units=:cgs)
#    return getfield(ST_gravity, get_funs[what])(eos, val, from=from, units=units)
#end

"""Piecewise-polytropic equation of state.

`PWP_EoS` supports several constructors, including a convenient named form:

    PWP_EoS(low_eos=:SLy, high_eos=:MPA1)

The EoS provides `get_density`, `get_pressure`, `get_energy_density`, etc.
"""
struct PWP_EoS <: AbstractEoS
    ρ::Vector{Float64}
    Γ::Vector{Float64}
    K::Vector{Float64}
    p::Vector{Float64}
    a::Vector{Float64}
    N::Int64
    params::Dict{Symbol,Any}

    function PWP_EoS(p1::Float64, ρ::Vector{Float64}, Γ::Vector{Float64}, params::Dict{Symbol}=Dict{Symbol,Any}())
        N = length(Γ)
        K = Vector{Float64}(undef, N)
        p = Vector{Float64}(undef, N)
        a = Vector{Float64}(undef, N)
        ϵ = Vector{Float64}(undef, N)

        p[1] = p1
        K[1] = p[1] / (c^2*ρ[1]^Γ[1])
        a[1] = 0.0

        for i in 2:N
            K[i] = K[i-1] * ρ[i-1] ^ (Γ[i-1] - Γ[i])
            p[i] = K[i] * c^2 * ρ[i] ^ Γ[i]
            a[i] = a[i-1] + (K[i-1] / (Γ[i-1] - 1.0) * ρ[i-1] ^ (Γ[i-1] - 1.0) - K[i] / (Γ[i] - 1.0) * ρ[i-1] ^ (Γ[i] - 1.0))
        end
        return new(ρ, Γ, K, p, a, N, params)
    end
end
function PWP_EoS(ρ::Array{Float64,1}, Γ::Array{Float64,1}, K::Array{Float64,1}, params::Dict{Symbol}=Dict{Symbol,Any}())
    p1 = K[1] * c^2 * ρ[1] ^ Γ[1]
    return PWP_EoS(p1, ρ, Γ, params)
end

function Base.show(io::IO, eos::PWP_EoS) 
    println(io, "Piecewise Polytropic EoS with ", eos.N, " pieces")
    print_dict(io, eos.params)
    return nothing
end

get_params(eos::PWP_EoS) = eos.params

find_piece(val, arr) = findfirst(x -> x > val, arr)
#find_piece(val, arr) = searchsortedfirst(arr, val)

function factor(units::Symbol)
    if units == :cgs
        return c^2
    elseif units == :natural
        return 1.0
    else 
        error("Units type $units is not defined")
    end
end

function get_density(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs)
    if from == :pressure
        p_val = val
        i = find_piece(p_val, eos.p)
        return (p_val / (eos.K[i] * factor(units))) ^ (1.0 / eos.Γ[i])
    elseif from == :density
        return val
    elseif from == :number_density
        return val * mp
    else
        error("There is no method from $from to ρ")
    end
end

function get_number_density(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return ρ_val / mp
end

function get_pressure(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs) #in dyn/cm^2
    ρ_val = get_density(eos, val, from=from, units=units)
    i = find_piece(ρ_val, eos.ρ)
    return eos.K[i] * ρ_val^eos.Γ[i] * factor(units)
end

function get_energy_density(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs) #with c^2
    ρ_val = get_density(eos, val, from=from, units=units)
    i = find_piece(ρ_val, eos.ρ)
    return ((1.0 + eos.a[i]) * ρ_val + eos.K[i] / (eos.Γ[i] - 1.0) * ρ_val^eos.Γ[i]) * factor(units)
end

function get_specific_enthalpy(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    i = find_piece(ρ_val, eos.ρ)
    return (1.0 + eos.a[i] + eos.Γ[i] / (eos.Γ[i] - 1.0) * eos.K[i] * eos.ρ[i] ^ (eos.Γ[i] - 1.0)) * factor(units)
end

function get_sound_velocity(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs) #cm/s
    ρ_val = get_density(eos, val, from=from, units=units)
    p_val = get_pressure(eos, val, from=from, units=units)
    ε_val = get_energy_density(eos, val, from=from, units=units)
    i = find_piece(ρ_val, eos.ρ)
    return sqrt(eos.Γ[i] * p_val / (ε_val + p_val)) * c
end

function get_internal_energy(eos::PWP_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    i = find_piece(ρ_val, eos.ρ)
    return (eos.a[i] + eos.K[i] / (eos.Γ[i] - 1.0) * ρ_val^(eos.Γ[i] - 1.0)) * factor(units)
end

function combine_eos(eos_soft::PWP_EoS, eos_hard::PWP_EoS)
        for i in eos_soft.N:-1:0
            if i > 0
                ρ_m = (eos_hard.K[1] / eos_soft.K[i]) ^ (1.0 / ( eos_soft.Γ[i] - eos_hard.Γ[1]))

                if i > 1 && ρ_m > eos_soft.ρ[i-1]
                    p1 = eos_soft.p[1]
                    ρ_comb = vcat(eos_soft.ρ[1:i-1], ρ_m, eos_hard.ρ)
                    Γ_comb = vcat(eos_soft.Γ[1:i], eos_hard.Γ)
                    params = merge(eos_soft.params, eos_hard.params)
                    eos_comb = PWP_EoS(p1, ρ_comb, Γ_comb, params)
                    return eos_comb
                elseif i == 1
                    p1 = get_pressure(eos_soft, ρ_m, from=:density)
                    ρ_comb = vcat(ρ_m, eos_hard.ρ)
                    Γ_comb = vcat(eos_soft.Γ[1], eos_hard.Γ)
                    params = merge(eos_soft.params, eos_hard.params)
                    eos_comb = PWP_EoS(p1, ρ_comb, Γ_comb, params)
                    return eos_comb
                else
                    continue
                end
            else
                return eos_hard
            end
        end
end

#function combine_EoS(eos_soft::PWP_EoS, eos_hard::PWP_EoS)
#               ρ_m = (eos_hard.K[1] / eos_soft.K[end]) ^ (1.0 / ( eos_soft.Γ[end] - eos_hard.Γ[1]))
#               p1 = eos_soft.p[1]
#               ρ_comb = vcat(eos_soft.ρ[1:end-1], ρ_m, eos_hard.ρ)
#               Γ_comb = vcat(eos_soft.Γ, eos_hard.Γ)
#               eos_comb = PWP_EoS(p1, ρ_comb, Γ_comb)
#           return eos_comb
#end

function PWP_EoS(exparams::Dict{Symbol}, exparams_symbolic::Dict{Symbol,Symbol})
    if haskey(exparams_symbolic, :low_eos) && haskey(exparams_symbolic, :high_eos)
        return PWP_EoS(low_eos = exparams_symbolic[:low_eos], high_eos = exparams_symbolic[:high_eos])
    elseif haskey(exparams_symbolic, :low_eos) && haskey(exparams, :log10p1) && haskey(exparams, :Γ1) && haskey(exparams, :Γ2) &&  haskey(exparams, :Γ3)
        return PWP_EoS(low_eos = exparams_symbolic[:low_eos], exparams[:log10p1], exparams[:Γ1], exparams[:Γ2], exparams[:Γ3])
    else
        error("Unable to define PWP_EoS with existing parameters")
    end
end


function PWP_EoS(log10p1::Float64, Γ1::Float64, Γ2::Float64, Γ3::Float64; low_eos::Symbol=:SLy)
    p1 = exp10(log10p1)
    Γ = [Γ1, Γ2, Γ3]
    ρ = [exp10(14.7), exp10(15), Inf]
    low_params = Dict(:low_eos => low_eos)
    low_eos = PWP_EoS(low_pwp3_eos_base[low_eos]..., low_params)
    high_params = Dict(:log10p1 => log10p1, :Γ1 => Γ1, :Γ2 => Γ2, :Γ3 => Γ3)
    high_eos = PWP_EoS(p1, ρ, Γ, high_params)
    return combine_eos(low_eos, high_eos)
end

function PWP_EoS(;low_eos::Symbol=:SLy, high_eos::Symbol=:SLy)
    low_params = Dict(:low_eos => low_eos)
    high_params = Dict(:high_eos => high_eos)
    low_eos = PWP_EoS(low_pwp3_eos_base[low_eos]..., low_params)
    if haskey(high_pwp3_eos_base, high_eos)
        log10p1, Γ1, Γ2, Γ3 = high_pwp3_eos_base[high_eos]
        p1 = exp10(log10p1)
        Γ = [Γ1, Γ2, Γ3]
        ρ = [exp10(14.7), exp10(15), Inf]
        high_eos = PWP_EoS(p1, ρ, Γ, high_params)
    elseif haskey(high_pwp_eos_base, high_eos)
        high_eos = PWP_EoS(high_pwp_eos_base[high_eos]..., high_params)
    else 
        error("Unable to find high eos $(high_eos)")
    end
    return combine_eos(low_eos, high_eos)
end

#PWP_EoS(;eos_name::Symbol) = PWP_EoS(low_eos=:SLy, high_eos=eos_name)

const low_pwp3_eos_base = Dict(
    :SLy => ([2.44034e7, 3.78358e11, 2.62780e12, Inf],
        [1.58425, 1.28733, 0.62223, 1.35692],
        [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8])
)

const high_pwp3_eos_base = Dict(
    :PAL6      => [34.380, 2.227, 2.189, 2.159], #Read2009
    :SLy       => [34.384, 3.005, 2.988, 2.851],
    :APR1      => [33.943, 2.442, 3.256, 2.908],
    :APR2      => [34.126, 2.643, 3.014, 2.945],
    :APR3      => [34.392, 3.166, 3.573, 3.281],
    :APR4      => [34.269, 2.830, 3.445, 3.348],
    :FPS       => [34.283, 2.985, 2.863, 2.600],
    :WFF1      => [34.031, 2.519, 3.791, 3.660],
    :WFF2      => [34.233, 2.888, 3.475, 3.517],
    :WFF3      => [34.283, 3.329, 2.952, 2.589],
    :BBB2      => [34.331, 3.418, 2.835, 2.832],
    :BPAL12    => [34.358, 2.209, 2.201, 2.176],
    :ENG       => [34.437, 3.514, 3.130, 3.168],
    :MPA1      => [34.495, 3.446, 3.572, 2.887],
    :MS1       => [34.858, 3.224, 3.033, 1.325],
    :MS2       => [34.605, 2.447, 2.184, 1.855],
    :MS1b      => [34.855, 3.456, 3.011, 1.425],
    :PS        => [34.671, 2.216, 1.640, 2.365],
    :GS1       => [34.504, 2.350, 1.267, 2.421],
    :GS2       => [34.642, 2.519, 1.571, 2.314],
    :BGN1H1    => [34.623, 3.258, 1.472, 2.464],
    :GNH3      => [34.648, 2.664, 2.194, 2.304],
    :H1        => [34.564, 2.595, 1.845, 1.897],
    :H2        => [34.617, 2.775, 1.855, 1.858],
    :H3        => [34.646, 2.787, 1.951, 1.901],
    :H4        => [34.669, 2.909, 2.246, 2.144],
    :H5        => [34.609, 2.793, 1.974, 1.915],
    :H6a       => [34.593, 2.637, 2.121, 2.064],
    :H7        => [34.559, 2.621, 2.048, 2.006],
    :PCL2      => [34.507, 2.554, 1.880, 1.977],
    :ALF1      => [34.055, 2.013, 3.389, 2.033],
    :ALF2      => [34.616, 4.070, 2.411, 1.890],
    :ALF3      => [34.283, 2.883, 2.653, 1.952],
    :ALF4      => [34.314, 3.009, 3.438, 1.803],
    :BCPM      => [34.385, 2.784, 2.920, 2.687], #Kumar2020
    :BKA20     => [34.599, 2.811, 2.451, 1.930],
    :BSk20     => [34.377, 3.141, 3.196, 3.042],
    :BSk21     => [34.539, 3.456, 3.073, 2.657],
    :BSk22     => [34.593, 3.147, 2.865, 2.668],
    :BSk23     => [34.571, 3.285, 2.954, 2.659],
    :BSk24     => [34.540, 3.457, 3.072, 2.656],
    :BSk25     => [34.525, 3.747, 3.067, 2.417],
    :BSk26     => [34.381, 3.141, 3.193, 3.040],
    :BSP       => [34.556, 3.204, 2.637, 2.218],
    :BSR2      => [34.661, 3.310, 2.951, 2.271],
    :BSR2Y     => [34.676, 3.378, 2.216, 1.892],
    :BSR6      => [34.664, 3.028, 3.046, 2.224],
    :BSR6Y     => [34.678, 3.075, 2.257, 1.915],
    :DD2       => [34.638, 3.414, 3.097, 2.322],
    :DD2Y      => [34.660, 3.523, 2.427, 2.004],
    :DDHd      => [34.597, 3.573, 2.649, 2.346],
    :DDME2     => [34.665, 3.639, 3.137, 2.259],
    :DDME2Y    => [34.679, 3.723, 2.376, 2.081],
    :FSU2      => [34.655, 2.675, 2.477, 1.830],
    :FSUGarnet => [34.624, 3.538, 2.556, 1.825],
    :G3        => [34.516, 3.115, 2.735, 2.194],
    :GM1       => [34.679, 2.937, 2.815, 2.438],
    :GM1Y      => [34.702, 3.032, 2.716, 2.013],
    :IOPB      => [34.640, 3.253, 2.664, 1.786],
    :KDE0v1    => [34.366, 2.791, 2.897, 2.779],
    :Model1    => [34.601, 3.247, 2.560, 1.830],
#   :MPA1      => [34.477, 3.441, 3.580, 2.884],
    :MS1b      => [34.845, 3.410, 3.030, 1.467],
    :NL3       => [34.847, 3.246, 3.098, 1.298],
    :NL3ω      => [34.821, 3.974, 3.127, 1.552],
    :NL3ωY     => [34.809, 3.922, 2.264, 2.166],
    :NL3ωYss   => [34.805, 3.913, 1.895, 2.106],
    :NL3Y      => [34.810, 3.092, 2.222, 2.214],
    :NL3Yss    => [34.802, 3.062, 1.766, 2.051],
    :Rs        => [34.555, 2.674, 2.670, 2.670],
    :SINPA     => [34.593, 3.321, 2.563, 1.839],
    :SK255     => [34.549, 2.623, 2.758, 2.703],
    :SK272     => [34.574, 2.730, 2.848, 2.766],
    :SKa       => [34.546, 2.810, 2.873, 2.783],
    :SKb       => [34.507, 3.143, 2.909, 2.808],
    :SkI2      => [34.613, 2.658, 2.588, 2.649],
    :SkI3      => [34.632, 2.824, 2.676, 2.697],
    :SkI4      => [34.507, 3.111, 2.909, 2.734],
    :SkI5      => [34.663, 2.587, 2.572, 2.718],
    :SkI6      => [34.519, 3.107, 2.918, 2.734],
    :SkMP      => [34.508, 2.782, 2.777, 2.729],
    :SKOp      => [34.451, 2.672, 2.712, 2.635],
    :SLY230a   => [34.399, 3.150, 3.082, 2.789],
    :SLY2      => [34.392, 2.959, 2.984, 2.829],
    :SLY4      => [34.380, 2.979, 2.999, 2.849],
    :SLY9      => [34.493, 2.992, 2.936, 2.750],
    :TM1       => [34.701, 2.754, 2.472, 1.870],
    :PhTr1     => [34.495, 3.446, 0.000, 4.000],
    :PhTr2     => [34.495, 3.446, 0.000, 2.887]
)

"""
    high_pwp3_eos_names

Names (symbols) of built-in high-density piecewise-polytrope parameter sets available for `PWP_EoS`.
"""
const high_pwp3_eos_names = collect(keys(high_pwp3_eos_base))

const high_pwp_eos_base = Dict(
    :ACB4 => [63.178e6*1e39*1.6021772e-12, [0.3174,0.5344,0.7500,Inf]*1e39*mp, [4.921, 0.0, 4.000, 2.800]],
    :test => [exp10(34.495), [exp10(14.7), exp10(15), exp10(15.3), Inf], [3.446, 3.572, 2.887, 2.887]]
)

const high_pwp_eos_names = collect(keys(high_pwp_eos_base))


"""
    find_max_pressure(eos::PWP_EoS; c_max=c)

Return an estimate of the maximum pressure (in CGS) for which the sound speed does not exceed `c_max`.

This is useful for building central-pressure grids that avoid superluminal sound speeds.
"""
function find_max_pressure(eos::PWP_EoS; c_max=c)
    cs_38 = get_sound_velocity(eos, 1e38, from=:pressure)
    if cs_38 > c_max
        p_max = find_zero(p -> get_sound_velocity(eos, p, from=:pressure) - c_max, (1e30,1e38), Bisection())
    else
        p_max = 1e38
    end
    return p_max
end

raw"""
    Polytropic_EoS(K, Γ)

Simple polytropic equation of state with $p = K\,\rho^{\Gamma}$.
"""
struct Polytropic_EoS <: AbstractEoS
    K::Float64
    Γ::Float64
    function Polytropic_EoS(K::Float64, Γ::Float64)
        return new(K, Γ)
    end
end

Polytropic_EoS(exparams::Dict{Symbol}, exparams_symbolic::Dict{Symbol,Symbol}) = Polytropic_EoS(exparams[:K], exparams[:Γ])

Base.show(io::IO, eos::Polytropic_EoS) = print(io, "Polytropic EoS with K = ", eos.K, " and Γ = ", eos.Γ)

get_params(eos::Polytropic_EoS) = Dict(:K => eos.K, :Γ => eos.Γ)

function get_density(eos::Polytropic_EoS, val::Real; from::Symbol, units=:cgs)
    if from == :pressure
        p_val = val
        return (p_val / eos.K) ^ (1.0 / eos.Γ)
    elseif from == :density
        return val
    elseif from == :number_density
        return val * mp
    else
        error("There is no method from $from to ρ")
    end
end

function get_number_density(eos::Polytropic_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return ρ_val / mp
end

function get_pressure(eos::Polytropic_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return eos.K * ρ_val ^ eos.Γ
end

function get_energy_density(eos::Polytropic_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return ((ρ_val * c^2 + eos.K / (eos.Γ - 1.0) * (ρ_val) ^ eos.Γ))
end


"""
        Table_EoS(name::Symbol)

Tabulated equation of state using simple linear interpolation in log-space.

This is the "known-working" table EoS implementation.
The input tables are taken from `table_eos_base` in `src/eos_data.jl`.

Notes:
- `table_eos_base` uses number density `n` in fm^-3, so this constructor converts `n -> ρ` via
    `ρ = n * 10^39 * m_p` (in g/cm^3).
- The Hermite-based experimental implementation is kept as `Table_EoS_Hermite`.
"""
struct Table_EoS <: AbstractEoS
    N::Int64
    lgρlgε::Any
    lgρlgp::Any
    lgplgρ::Any
    lgplgε::Any
    name::Symbol

    function Table_EoS(name::Symbol)
        table = sort(table_eos_base[name], dims=1)
        lgρ = log10.(table[:,1] .* 1e39 .* mp)
        lgε = log10.(table[:,2] .* c^2)
        lgp = log10.(table[:,3])
        lgρlgε = linear_interpolation(lgρ, lgε; extrapolation_bc = Line())
        lgρlgp = linear_interpolation(lgρ, lgp; extrapolation_bc = Line())
        lgplgρ = linear_interpolation(lgp, lgρ; extrapolation_bc = Line())
        lgplgε = linear_interpolation(lgp, lgε; extrapolation_bc = Line())
        N = length(lgρ)
        return new(N, lgρlgε, lgρlgp, lgplgρ, lgplgε, name)
    end
end

Table_EoS(eos_params::Dict{Symbol}) = Table_EoS(eos_params[:eos])
get_params(eos::Table_EoS) = Dict(:eos => eos.name)

Base.show(io::IO, eos::Table_EoS) = print(io, "Tabulated EoS ", eos.name, " with ", eos.N, " points")

function get_density(eos::Table_EoS, val::Real; from::Symbol, units=:cgs)
    if from == :pressure
        p_val = val
        return exp10(eos.lgplgρ(log10(p_val)))
    elseif from == :density
        return val
    elseif from == :number_density
        return val * mp
    else
        error("There is no method from $from to ρ")
    end
end

function get_number_density(eos::Table_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return ρ_val / mp
end

function get_pressure(eos::Table_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    return exp10(eos.lgρlgp(log10(ρ_val))) * factor(units)
end

function get_energy_density(eos::Table_EoS, val::Real; from::Symbol, units=:cgs)
    if from == :pressure
        p_val = val
        return exp10(eos.lgplgε(log10(p_val))) * factor(units)
    else
        ρ_val = get_density(eos, val, from=from, units=units)
        return exp10(eos.lgρlgε(log10(ρ_val))) * factor(units)
    end
end

function get_sound_velocity(eos::Table_EoS, val::Real; from::Symbol, units=:cgs)
    ρ_val = get_density(eos, val, from=from, units=units)
    ρ_val <= 0 && error("Density must be positive")

    # Finite-difference estimate of c_s^2 = dp/dε.
    δ = 1e-4
    ρ1 = ρ_val * (1 - δ)
    ρ2 = ρ_val * (1 + δ)

    p1 = get_pressure(eos, ρ1, from=:density, units=units)
    p2 = get_pressure(eos, ρ2, from=:density, units=units)
    ε1 = get_energy_density(eos, ρ1, from=:density, units=units)
    ε2 = get_energy_density(eos, ρ2, from=:density, units=units)

    denom = (ε2 - ε1)
    denom == 0 && error("Unable to estimate dp/dε (ε is locally flat)")
    cs2 = (p2 - p1) / denom
    cs2 < 0 && error("Estimated dp/dε < 0 (table is not thermodynamically consistent near ρ=$(ρ_val))")
    return sqrt(cs2) * c
end

function smooth(η)
    η_sm = similar(η)
    η_sm[1] = 0.5*η[1] + 0.5*η[2]
    for j in 2:length(η)-1
         η_sm[j] = 0.25*η[j-1] + 0.5*η[j] + 0.25*η[j+1]
    end
    η_sm[end] = 0.5*η[end-1] + 0.5*η[end]
    return η_sm
end

"""
    Table_EoS_Hermite(n, ε, p; name=:unknown)

Experimental Hermite-interpolated table EoS.

This implementation is intentionally *not* the default `Table_EoS` because it is
unfinished / not validated for the bundled tables.
"""
struct Table_EoS_Hermite <: AbstractEoS
    N::Int64
    n::Vector{Float64}
    ε::Vector{Float64}
    p::Vector{Float64}
    ε_n_l::Vector{Float64}
    ε_n_r::Vector{Float64}
    ε_nn_l::Vector{Float64}
    ε_nn_r::Vector{Float64}
    name::Symbol

    function Table_EoS_Hermite(n::Vector, ε::Vector, p::Vector; name::Symbol=:unknown)
        N = length(n)
#        dE_dn = get_derivative(E, n)
        dε_dn = (ε .+ p) ./ n
        for i in 1:5
            dε_dn = smooth(dε_dn)
        end
        d2ε_dn2 = get_derivative(p, n) ./ n 
        ε_n_l = dε_dn[1:N-1] .* (n[2:N] - n[1:N-1])
        ε_n_r = dε_dn[2:N] .* (n[2:N] - n[1:N-1])
        ε_nn_l = d2ε_dn2[1:N-1] .* (n[2:N] - n[1:N-1]).^2
        ε_nn_r = d2ε_dn2[2:N] .* (n[2:N] - n[1:N-1]).^2

        return new(N, n, ε, p, ε_n_l, ε_n_r, ε_nn_l, ε_nn_r, name)
    end
end

function Table_EoS_Hermite(name::Symbol)
    table = table_eos_base[name]
    n = vec(table[:, 1])
    ε = vec(table[:, 2])
    p = vec(table[:, 3])
    return Table_EoS_Hermite(n, ε, p; name=name)
end
function get_derivative0(Q::Vector, x::Vector)
    N = length(Q)
    dQ_dx = similar(Q)
    for i in 2:N-1
        dQ_dx[i] = (Q[i+1] - Q[i]) / (x[i+1] - x[i])
    end
    dQ_dx[N] = dQ_dx[N-1]
    return dQ_dx
end

function get_derivative(Q::Vector, x::Vector)
    N = length(Q)
    dQ_dx = similar(Q)
    h1 = x[2] - x[1]
    h2 = x[3] - x[2]
    dQ_dx[1] = (2.0 * h1 * h2 * (-Q[1] + Q[2]) + h2^2 * (-Q[1] + Q[2]) + h1^2 * (Q[2] - Q[3])) / (h1 * h2 * (h1 + h2))
    for i in 2:N-1
        h1 = x[i] - x[i-1]
        h2 = x[i+1] - x[i]
        dQ_dx[i] = (h2^2 * (-Q[i-1] + Q[i]) + h1^2 * (-Q[i] + Q[i+1])) / (h1*h2*(h1 + h2)) 
    end
    h1 = x[end-1] - x[end-2]
    h2 = x[end] - x[end-1]
    dQ_dx[N] = (h2^2 * (Q[end-2] - Q[end-1]) + h1^2 * (-Q[end-1] + Q[end]) + 2.0 * h1 * h2 * (-Q[end-1] + Q[end])) / (h1 * h2 * (h1 + h2))
    return dQ_dx
end

function get_derivative2(Q::Vector, x::Vector)
    N = length(Q)
    dQ_dx = similar(Q)
    dQ_dx[1] = (-3*Q[1] + 4*Q[2] - Q[3]) / (-3*x[1] + 4*x[2] - x[3])
    for i in 2:N-1
        dQ_dx[i] = (Q[i+1] - Q[i-1]) / (x[i+1] - x[i-1])
    end
    dQ_dx[N] = (Q[N-2] - 4*Q[N-1] + 3*Q[N]) / (x[N-2] - 4*x[N-1] + 3*x[N])
    return dQ_dx
end

ψ0(z) = -6z^5 + 15z^4 - 10z^3 + 1
ψ1(z) = -3z^5 + 8z^4 - 6z^3 + z
ψ2(z) = 0.5*(-z^5 +3z^4 - 3z^3 + z^2)

dψ0(z) = -30z^4 + 60z^3 - 30z^2
dψ1(z) = -15z^4 + 32z^3 - 18z^2 + 1
dψ2(z) = 0.5*(-5z^4 +12z^3 - 9z^2 + 2z)

H0(z) = 1 - 3z^2 + 2z^3
H1(z) = z - 2z^2 + z^3
H2(z) = -z^2 + z^3
H3(z) = 3z^2 - 2z^3 

dH0(z) = - 6z + 6z^2
dH1(z) = 1 - 4z + 3z^2
dH2(z) = -2z + 3z^2
dH3(z) = 6z - 6z^2 


function get_energy_density(eos::Table_EoS_Hermite, val::Real; from::Symbol, units=:cgs)
    if from == :number_density
        n_val = val
        i = searchsortedfirst(eos.n, n_val) - 1
        if i == 0
            error("out of range")
        end
        x = (n_val - eos.n[i]) / (eos.n[i+1] - eos.n[i])
#        ε_val = (eos.ε[i]*ψ0(x) + eos.ε[i+1]*ψ0(1-x) + eos.ε_n_l[i]*ψ1(x) - eos.ε_n_r[i]*ψ1(1-x) + eos.ε_nn_l[i]*ψ2(x) + eos.ε_nn_r[i]*ψ2(1-x))
        ε_val = (eos.ε[i]*H0(x) + eos.ε[i+1]*H3(x) + eos.ε_n_l[i]*H1(x) + eos.ε_n_r[i]*H2(x))
        return ε_val * factor(units)
    else
        error("Table_EoS_Hermite only implements from=:number_density")
    end
end

function get_de_dn(eos, n_val)
    i = searchsortedfirst(eos.n, n_val) - 1
    x = (n_val - eos.n[i]) / (eos.n[i+1] - eos.n[i])
    de_dn_val = (eos.ε[i]*dH0(x) + eos.ε[i+1]*dH3(x) + eos.ε_n_l[i]*dH1(x) + eos.ε_n_r[i]*dH2(x)) / (eos.n[i+1] - eos.n[i])
#    de_dn_val = (eos.ε[i] * dψ0(x) - eos.ε[i+1]*dψ0(1-x) + eos.ε_n_l[i]*dψ1(x) + eos.ε_n_r[i]*dψ1(1-x) + eos.ε_nn_l[i]*dψ2(x) - eos.ε_nn_r[i]*dψ2(1-x)) / (eos.n[i+1] - eos.n[i])
    return de_dn_val
end

function get_pressure(eos::Table_EoS_Hermite, val::Real; from::Symbol, units=:cgs)
    if from == :number_density
        n_val = val
        i = searchsortedfirst(eos.n, n_val) - 1
        if i == 0
            error("out of range")
        end
        x = (n_val - eos.n[i]) / (eos.n[i+1] - eos.n[i])
        ε_val = get_energy_density(eos, val; from=from, units=units)
#        p_val = (eos.ε[i] * dψ0(x) - eos.ε[i+1]*dψ0(1-x) + eos.ε_n_l[i]*dψ1(x) + eos.ε_n_r[i]*dψ1(1-x) + eos.ε_nn_l[i]*dψ2(x) - eos.ε_nn_r[i]*dψ2(1-x)) / (eos.n[i+1] - eos.n[i]) * n_val - ε_val
        p_val = (eos.ε[i]*dH0(x) + eos.ε[i+1]*dH3(x) + eos.ε_n_l[i]*dH1(x) + eos.ε_n_r[i]*dH2(x)) / (eos.n[i+1] - eos.n[i]) * n_val - ε_val
        return p_val * factor(units)
    else
        error("Table_EoS_Hermite only implements from=:number_density")
    end
end

struct Table_EoSnew2 <: AbstractEoS
    N::Int64
    n::Vector{Float64}
    ε::Vector{Float64}
    p::Vector{Float64}
    K::Vector{Float64}
    Γ::Vector{Float64}
    a::Vector{Float64}
    name::Symbol

    function Table_EoSnew2(n::Vector, ε::Vector, p::Vector)
        N = length(n)
        K_arr = similar(n)
        Γ_arr = similar(n)
        a_arr = similar(n)
        for i in 1:N-1
            ρ0, ρ1 = (n[i], n[i+1]) .* mp
            p0, p1 = p[i], p[i+1]
            e0, e1 = ε[i], ε[i+1]
            Γ0 = log(p1/p0) / log(ρ1/ρ0)
            K0 = p0 / ρ0 ^ Γ0            
            a0 = (e0 - K0 * ρ0 ^ Γ0 / (Γ0 - 1)) / (ρ0*c^2)
            K_arr[i], Γ_arr[i], a_arr[i] = K0, Γ0, a0
        end
        name = :unknown
        return new(N, n, ε, p, K_arr, Γ_arr, a_arr, name)
    end
end


