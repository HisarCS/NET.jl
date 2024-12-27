module NET

# ---- imports ----
using DifferentialEquations
using Plots
using LinearAlgebra
using Optim
using Random
using Logging

# ---- include ----
include("driven_stohcastic_systems.jl")
include("entropy_optimization.jl.jl")
include("entropy_production_rate.jl")
include("fluctuation_theorems.jl")
include("fokker_planck.jl")
include("langevin_dynamics.jl")
include("linear_onsager.jl")
include("local_entropy_production.jl")
include("noise_induced_bifuractions.jl")
include("nonlinear_onsager.jl")
include("non_equilibrium_potential.jl")
include("non_linear_onsager.jl")
include("path_integral_formulation.jl")
include("stochastic_onsager_relations.jl")
include("time_dependent_onsager.jl")
include("time_dependent_entropy.jl")

### ---- exports ----

export compute_stochastic_fluxes, validate_dimensions_stochastic, visualize_stochastic_fluxes
export compute_entropy, entropy_gradient, gradient_descent
export compute_entropy_production_rate, validate_dimensions_entropy
export thermal_noise, mean, jarzynski_equality, crooks_theorem
export solve_fokker_planck
export langevin_dynamics
export compute_fluxesi validate_dimensions_lo, visualize_fluxes
export compute_local_entropy_production, validate_dimensions_entropy_le, visualize_local_entropy_production_heatmap
export detect_bifurcation
export  calculate_non_equilibrium_potential
export compute_total_fluxes, validate_dimensions_multinonlinear
export calculate_action
export analyze_time_dependent_entropy_production
export compute_time_dependent_fluxes, validate_dimensions_time_dependent, visualize_time_dependent_fluxes

end