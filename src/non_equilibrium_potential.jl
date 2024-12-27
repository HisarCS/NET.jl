using Logging
using Plots

"""
    calculate_non_equilibrium_potential(grid, drift, diffusion)

Calculates the non-equilibrium potential landscape for a given drift and diffusion.

# Arguments
- grid::Vector{Float64}: Spatial grid for the potential.
- drift::Function: Drift function, `drift(x)`.
- diffusion::Float64: Diffusion coefficient.

# Returns
- potential::Vector{Float64}: Non-equilibrium potential landscape.
"""
function calculate_non_equilibrium_potential(grid::Vector{Float64}, drift::Function, diffusion::Float64)
    potential = zeros(length(grid))
    for i in 2:length(grid)
        dx = grid[i] - grid[i-1]
        potential[i] = potential[i-1] + drift(grid[i]) * dx / diffusion
    end
    info("Non-equilibrium potential landscape calculation complete for grid of size $(length(grid)).")
    return potential
end

"""
    visualize_potential(grid::Vector{Float64}, potential::Vector{Float64}, title_name::String)

Plots the calculated non-equilibrium potential landscape.

# Arguments
- grid::Vector{Float64}: Spatial grid for the potential.
- potential::Vector{Float64}: Non-equilibrium potential values.
- title_name::String: Title for the plot.
"""
function visualize_potential(grid::Vector{Float64}, potential::Vector{Float64}, title_name::String)
    plot(
        grid,
        potential,
        xlabel="Grid Points",
        ylabel="Potential",
        title=title_name,
        legend=false,
        color=:blue
    )
end
