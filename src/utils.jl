rnd(U) = @sprintf("%0.4g", U)

"""
    LogRange(L, R, N)

Return `N` logarithmically spaced points from `L` to `R` (inclusive).

Equivalent to `10 .^ LinRange(log10(L), log10(R), N)`.
"""
LogRange(L, R, N) = 10.0 .^ LinRange(log10(L), log10(R), N)

    function print_dict(io::IO, dict::Dict)
        for key in keys(dict)
            println(io, "  :", key, " => ", dict[key])
        end
    end

function plot_radial_structure(RS_sol, model)
    α0, β0  = model.exparams[:α0], model.exparams[:β0]
    Qs = model.quantities

    μ = RS_sol[1,:]
    ν = RS_sol[2,:]
    φ = RS_sol[3,:]
    ψ = RS_sol[4,:]
    p̃ = RS_sol[5,:]
    M̃ = RS_sol[6,:]

    r = RS_sol.t

    plot(r ./ km, μ ./ μ[end], label="μ / μ_s")
    plot(r ./ km, ν ./ ν[end], label="ν / ν_s")
    plot(r ./ km, φ, label="φ")
    plot(r ./ km, ψ, label="ψ")
    plot(r ./ km, p̃ ./ p̃[1], label="p̃ / p̃_c")
    plot(r ./ km, M̃ ./ M̃[end], label="M̃ / M̃_s")
    plot(r ./ km, 2.0 * μ ./ r, label="2μ/r")

    title("α0 = $(rnd(α0)), β0 = $(rnd(β0))\nm̃A = $(rnd(Qs[:m̃A]/M_sun)) Ms, R = $(rnd(Qs[:R]/km)) km, p̃_c = $(rnd(Qs[:p̃_c])) dyn/cm^3")
    xlabel("r, km")
    yscale("symlog", linthreshy=1e-4)
    xscale("log")
    xlim(left=1e-3)
    legend()
    show()

    return nothing
end

function show_status(res, RS_sol, Qs, Param_Vars)
    μ_c, ν_c, φ_c, ψ_c, p̃_c, M̃_c = Vars_c = RS_sol[1]
    r_c = RS_sol.t[1]
    A_c, α_c, ñ_c, ε̃_c = internal_physics(Vars_c, Param_Vars)
    μ_s, ν_s, φ_s, ψ_s, p̃_s, M̃_s = Vars_s = RS_sol[end]
    r_s = RS_sol.t[end]
    A_s, α_s, ñ_s, ε̃_s = internal_physics(Vars_s, Param_Vars)
    N = length(RS_sol.t)

    println("r_c = $r_c, μ_c = $μ_c, ν_c = $ν_c, φ_c = $φ_c, ψ_c = $ψ_c, p̃_c = $p̃_c, M̃_c = $M̃_c")
    println("r_s = $r_s, μ_s = $μ_s, ν_s = $ν_s, φ_s = $φ_s, ψ_s = $ψ_s, p̃_s = $p̃_s, M̃_s = $M̃_s")
    println("A_c = $A_c, α_c = $α_c, A_s = $A_s, α_s = $α_s")
    println("N = $N, αA = $(Qs[:αA]), φ0 = $(Qs[:φ0]), φ_span = $(Qs[:φ_span])")
    println("residuals = $res")
    println("")
    return res
end

