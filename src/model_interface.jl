abstract type AbstractModel end

#fields
inparams(model::AbstractModel) = model.inparams
eos(model::AbstractModel) = model.eos
quantities(model::AbstractModel) = model.quantities
n_odes(model::AbstractModel) = model.n_odes
rspan(model::AbstractModel) = model.rspan
vars_init(model::AbstractModel) = model.vars_init


#methods
#ode_system!(DVars::Vector{Float64}, Vars::Vector{Float64}, model::AbstractModel, r::Float64) = error("Ordinary differential equation system is not defined for $(typeof(model))")
#internal_physics(Vars::Vector{Float64}, model::AbstractModel, r::Float64) = error("Internal physics method is not implemented for $(typeof(model))")
#calculate_quantities!(model::AbstractModel) = error("Method for calculating quantities is not implemented for $(typeof(model))")
#surface_callback(model::AbstractModel) = error("Callback for surface conditions is not implemented for $(typeof(model))")
#boundary_conditions(model::AbstractModel) = error("Boundary conditions are not implemented for $(typeof(model))")
