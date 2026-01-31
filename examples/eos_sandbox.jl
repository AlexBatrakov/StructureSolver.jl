@time using PyPlot
@time using StructureSolver
@info "Initialization is completed"

#LogRange(L, R, N) = 10.0 .^ LinRange(log10(L), log10(R), N)

ρ_arr = LogRange(1e4, 1e16, 2000)

#eos = Table_EoS(:SLy4)
#p_arr = get_pressure(eos, ρ_arr, from=:density)
#plot(ρ_arr, p_arr, label="Tabled SLY4")

for high_name in keys(StructureSolver.high_pwp3_eos_base)
    eos = PWP_EoS(low_eos=:SLy, high_eos=high_name)
    p_arr = get_pressure(eos, ρ_arr, from=:density)
    plot(ρ_arr, p_arr, label=high_name)
end
xscale("log")
yscale("log")
xlabel("Density, g/cm^3")
ylabel("Pressure, dyn/cm^2")
legend()
show()


eos = Table_EoS(:SLy4)
ε_arr = get_energy_density(eos, ρ_arr, from=:density)
plot(ρ_arr, ε_arr, label="Tabled SLY4")

for high_name in keys(StructureSolver.high_pwp3_eos_base)
    eos = PWP_EoS(low_eos=:SLy, high_eos=high_name)
    ε_arr = get_energy_density(eos, ρ_arr, from=:density)
    plot(ρ_arr, ε_arr, label=high_name)
end
xscale("log")
yscale("log")
xlabel("Density, g/cm^3")
ylabel("Energy density, erg/cm^3")
legend()
show()


for high_name in keys(StructureSolver.high_pwp3_eos_base)
    eos = PWP_EoS(low_eos=:SLy, high_eos=high_name)
    plot(eos.ρ, eos.Γ, ".", label=high_name)
end
xscale("log")
legend()
show()

#-------------------------------------------------------------------------------
#18.01.2022 check eos with phase transition

StructureSolver.high_pwp_eos_base[:test] .= [exp10(34.495), [exp10(14.7), exp10(14.8), exp10(15.0), Inf], [3.446, 3.572, 3.572, 2.887]]
eos_PhTr1 = PWP_EoS(low_eos=:SLy, high_eos=:test)
eos_PhTr2 = PWP_EoS(34.6, 4.0, 0.000, 3.0; low_eos=:SLy)
eos_MPA1  = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)
eos_WFF1  = PWP_EoS(low_eos=:SLy, high_eos=:WFF1)

#eos_PhTr2 = PWP_EoS(34.495, [exp10(14.7), exp10(15),exp10(15.3), Inf], [3.446, 3.572, 2.887, 2.887]; low_eos=:SLy)
#eos_PhTr2 = PWP_EoS(log10(63.178e6*1e39*1.6021772e-12), [0.3174,0.5344,0.7500,Inf]*1e39*mp, [4.921, 0.0, 4.000, 2.800]; low_eos=:SLy)
#log10(63.178e6*1e39*1.6021772e-12)

ρ_arr = LogRange(1e4, 1e16, 2000)
p_arr_PhTr1 = get_pressure(eos_PhTr1, ρ_arr, from=:density)
p_arr_PhTr2 = get_pressure(eos_PhTr2, ρ_arr, from=:density)
p_arr_MPA1  = get_pressure(eos_MPA1, ρ_arr, from=:density)

plot(ρ_arr, p_arr_PhTr1, label="PhTr1")
plot(ρ_arr, p_arr_PhTr2, label="PhTr2")
plot(ρ_arr, p_arr_MPA1, label="MPA1")
xscale("log")
yscale("log")
xlabel("Density, g/cm^3")
ylabel("Pressure, dyn/cm^2")
legend()
show()
axvline(exp10(14.7), color="grey")
axvline(exp10(15.0), color="grey")
axhline(exp10(34.495), color="grey")

ε_arr_PhTr1 = get_energy_density(eos_PhTr1, ρ_arr, from=:density)
ε_arr_MPA1  = get_energy_density(eos_MPA1, ρ_arr, from=:density)

plot(ρ_arr, ε_arr_PhTr1, label="PhTr1")
plot(ρ_arr, ε_arr_MPA1, label="MPA1")
xscale("log")
yscale("log")
xlabel("Density, g/cm^3")
ylabel("Energy density, erg/cm^3")
legend()
show()

plot(p_arr_PhTr1, ε_arr_PhTr1, label="PhTr1")
plot(p_arr_MPA1, ε_arr_MPA1, label="MPA1")
xscale("log")
yscale("log")
xlabel("Pressure, dyn/cm^2")
ylabel("Energy density, erg/cm^3")
legend()
show()

cs_arr_PhTr1 = get_sound_velocity(eos_PhTr1, ρ_arr, from=:density)
cs_arr_MPA1  = get_sound_velocity(eos_MPA1, ρ_arr, from=:density)
cs_arr_WFF1  = get_sound_velocity(eos_WFF1, ρ_arr, from=:density)

plot(ρ_arr, cs_arr_PhTr1 ./c, label="PhTr1")
plot(ρ_arr, cs_arr_MPA1 ./c, label="MPA1")
plot(ρ_arr, cs_arr_WFF1 ./c, label="WFF1")
xscale("log")
yscale("log")
xlabel("Density, g/cm^3")
ylabel("Sound velocity, c")
legend()
show()

plot(p_arr_PhTr1, cs_arr_PhTr1 ./c, label="PhTr1")
plot(p_arr_MPA1, cs_arr_MPA1 ./c, label="MPA1")
xscale("log")
yscale("log")
xlabel("Pressure, dyn/cm^2")
ylabel("Sound velocity, c")
legend()
show()

get_pressure(eos_PhTr1, exp10(14.7), from=:density)
get_pressure(eos_PhTr1, exp10(15.0), from=:density)
get_density(eos_PhTr1, 3.1260793671239247e34, from=:pressure)

StructureSolver.high_pwp3_eos_base[:PhTr1] .= [34.85, 4.5, 0.000, 5.0] #partially shielded
StructureSolver.high_pwp3_eos_base[:PhTr1] .= [34.9, 4.5, 0.000, 5.0]  #fully shielded
StructureSolver.high_pwp3_eos_base[:PhTr1] .= [35.0, 4.5, 0.000, 5.0]  #over shielded
StructureSolver.high_pwp3_eos_base[:PhTr1] .= [34.8, 3.0, 1.05, 3.0]

#check GR case

#find_max_pressure(eos)

eos = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
int_params = IntParams(1e4, 1.0, 1e-12, 1e-12)
inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
family_params = Dict(:pc => LogRange(1e34, 1e37, 1000))
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)
simulation = FamilySimulation(model, regime, family_params, int_params)
calculate!(simulation)

plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun, label="MPA1 GR")
xlabel("Radius, km")
ylabel("Mass, M_sun")

eos = PWP_EoS(low_eos=:SLy, high_eos=:WFF1)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
int_params = IntParams(1e4, 1.0, 1e-12, 1e-12)
inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
family_params = Dict(:pc => LogRange(1e34, 1e37, 1000))
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)
simulation = FamilySimulation(model, regime, family_params, int_params)
calculate!(simulation)

plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun, label="WFF1 GR")
xlabel("Radius, km")
ylabel("Mass, M_sun")

#[34.031, 2.519, 3.791, 3.660] ENG
#StructureSolver.high_pwp_eos_base[:test] .=[exp10(34.495), [exp10(14.7), exp10(15), exp10(15.3), Inf], [3.446, 3.572, 2.887, 2.887]]
StructureSolver.high_pwp_eos_base[:test] .= [exp10(34.0), [exp10(14.7), exp10(15.0), exp10(15.1), Inf], [2.50, 3.8, 0.0, 4.0]]

eos = PWP_EoS(low_eos=:SLy, high_eos=:test)
cf = DEF1_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
int_params = IntParams(1e4, 1.0, 1e-12, 1e-12)
inparams_fixed = Dict(:φc => 0.0, :pc => 1e35)
exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
family_params = Dict(:pc => LogRange(1e34, 1e37, 1000))
regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)
simulation = FamilySimulation(model, regime, family_params, int_params)
calculate!(simulation)

plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun, label="PhTr1 GR")

#check DEF case

eos = PWP_EoS(low_eos=:SLy, high_eos=:MPA1)
cf = DEF3_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
int_params = IntParams(1e4, 1.0, 1e-10, 1e-10)
inparams_fixed = Dict(:pc => 1e35)
exparams = Dict(:α0 => -1e-4, :β0 => -4.5)
exparams_symbolic = Dict(:low_eos => :SLy, :high_eos => :MPA1)
imparams_shooting = Dict(:φc => 0.1)
quantities_fixed = Dict(:bc_φ∞ => 0.0)
family_params = Dict(:pc => LogRange(1e34, 1e36, 100))
regime = Simple_ShootingRegime(inparams_fixed, exparams, exparams_symbolic, imparams_shooting, quantities_fixed)
simulation = FamilySimulation(model, regime, family_params, int_params)
calculate!(simulation)

plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun, label="MPA1 DEF")

eos = PWP_EoS(low_eos=:SLy, high_eos=:PhTr1)
cf = DEF3_CouplingFunction()
model = DEFp_Model{Float64}(cf, eos)
int_params = IntParams(1e4, 1.0, 1e-10, 1e-10)
inparams_fixed = Dict(:pc => 1e35)
exparams = Dict(:α0 => -1e-4, :β0 => -4.5)
exparams_symbolic = Dict(:low_eos => :SLy, :high_eos => :test)
imparams_shooting = Dict(:φc => 0.1)
quantities_fixed = Dict(:bc_φ∞ => 0.0)
family_params = Dict(:pc => LogRange(1e34, 1e36, 100))
regime = Simple_ShootingRegime(inparams_fixed, exparams, exparams_symbolic, imparams_shooting, quantities_fixed)
simulation = FamilySimulation(model, regime, family_params, int_params)
calculate!(simulation)

plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun, label="PhTr1 DEF")

plot(simulation.family.quantities[:mA]./M_sun, -simulation.family.quantities[:αA])
xlabel("Mass, M_sun")
ylabel(L"\log_{10}(|\alpha_0|)")
yscale("log")

plot(simulation.family.quantities[:pc], -simulation.family.quantities[:αA])
xscale("log")
yscale("log")

plot(simulation.family.quantities[:mA], simulation.family.quantities[:cs_c]./c)
xscale("log")
axhline(1.0, color="grey")

Γ1_arr = LinRange(2.0, 5.0, 4)
Γ2_arr = [1.005] #vcat(0.0, 0.5, LinRange(2.0, 5.0, 7))
Γ3_arr = LinRange(2.0, 5.0, 4)
log10p1_arr = LinRange(34.6, 34.95, 8)

for log10p1 in log10p1_arr, Γ1 in Γ1_arr, Γ2 in Γ2_arr, Γ3 in Γ3_arr
#    println([log10p1, Γ1, Γ2, Γ3])
#    StructureSolver.high_pwp3_eos_base[:PhTr1] .= [log10p1, Γ1, Γ2, Γ3]
    StructureSolver.high_pwp_eos_base[:test] .= [exp10(34.0), [exp10(14.7), exp10(15.0), exp10(15.1), Inf], [2.50, 3.8, 0.0, 4.0]]

    eos_name = :test
    eos = PWP_EoS(low_eos=:SLy, high_eos=eos_name)
    cf = DEF1_CouplingFunction()
    model = DEFp_Model{Float64}(cf, eos)
    int_params = IntParams(1e4, 1.0, 1e-9, 1e-9)
    inparams_fixed = Dict(:φc => 0.0, :pc => 1e34)
    exparams_fixed = Dict(:α0 => 0.0, :β0 => 0.0)
    family_params = Dict(:pc => LogRange(1e34, find_max_pressure(eos, c_max=1.0*c), 100))
    # find_max_pressure(eos, c_max=1.2*c)
    regime = Simple_DirectRegime(inparams_fixed, exparams_fixed)
    simulation = FamilySimulation(model, regime, family_params, int_params)
    calculate!(simulation)
    mA, i = findmax(simulation.family.quantities[:mA]./M_sun)
    R = simulation.family.quantities[:R][i]
    if 2.0 < mA < 3.0 && R < 17.5e5

        inparams_fixed = Dict(:pc => 1e35)
        exparams = Dict(:α0 => -1e-4, :β0 => -4.5)
        exparams_symbolic = Dict(:low_eos => :SLy, :high_eos => eos_name)
        imparams_shooting = Dict(:φc => 0.1)
        quantities_fixed = Dict(:bc_φ∞ => 0.0)
        family_params = Dict(:pc => LogRange(1e34, find_max_pressure(eos, c_max=1.0*c), 100))
        regime = Simple_ShootingRegime(inparams_fixed, exparams, exparams_symbolic, imparams_shooting, quantities_fixed)
        simulation = FamilySimulation(model, regime, family_params, int_params)
        calculate!(simulation)
#
#        plot(simulation.family.quantities[:R]./1e5, simulation.family.quantities[:mA]./M_sun)
#        xlabel("Radius, km")
#        ylabel("Mass, M_sun")
#        xlim(0.0,20.0)
        plot(simulation.family.quantities[:mA]./M_sun, -simulation.family.quantities[:αA])
        xlabel("Mass, M_sun")
        ylabel(L"\log_{10}(|\alpha_0|)")
        yscale("log")

    end
end

