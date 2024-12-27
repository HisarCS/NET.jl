"""
    solve_fokker_planck(x_range, t_range, D, drift)

Manually solves the Fokker-Planck equation for a given drift field and diffusion coefficient
using an explicit finite difference method.

# Arguments
- `x_range::Vector{Float64}`: Spatial range for the solution.
- `t_range::Vector{Float64}`: Time range for the solution.
- `D::Float64`: Diffusion coefficient.
- `drift::Function`: Drift function, `drift(x)`.

# Returns
- `solution::Matrix{Float64}`: Probability distribution over space and time.
"""
function solve_fokker_planck(x_range::Vector{Float64}, t_range::Vector{Float64}, D::Float64, drift::Function)
    dx = x_range[2] - x_range[1]
    dt = t_range[2] - t_range[1]
    if dt > dx^2 / (2 * D)
        error("The time step size is too large for stability. Choose a smaller dt.")
    end
    P = zeros(length(x_range), length(t_range))
    for i in 1:length(x_range)
        P[i, 1] = exp(-x_range[i]^2 / 2) / sqrt(2π)
    end
    for t in 1:(length(t_range) - 1)
        for i in 2:(length(x_range) - 1)
            drift_term = -drift(x_range[i]) * (P[i+1, t] - P[i-1, t]) / (2 * dx)
            diff_term = D * (P[i+1, t] - 2 * P[i, t] + P[i-1, t]) / dx^2
            P[i, t+1] = P[i, t] + dt * (drift_term + diff_term)
        end
        P[1, t+1] = P[1, t]
        P[end, t+1] = P[end, t]
    end
    return P
end

"""
    compute_steady_state(x_range, D, drift)

Manually computes the steady-state solution to the Fokker-Planck equation
by solving a tridiagonal linear system.

# Arguments
- `x_range::Vector{Float64}`: Spatial range for the solution.
- `D::Float64`: Diffusion coefficient.
- `drift::Function`: Drift function, `drift(x)`.

# Returns
- `steady_state::Vector{Float64}`: Steady-state probability distribution.
"""
function compute_steady_state(x_range::Vector{Float64}, D::Float64, drift::Function)
    n = length(x_range)
    dx = x_range[2] - x_range[1]
    A = zeros(n, n)
    b = zeros(n)
    for i in 2:(n-1)
        drift_coeff = drift(x_range[i]) / (2 * dx)
        diffusion_coeff = D / dx^2
        A[i, i-1] = -diffusion_coeff - drift_coeff
        A[i, i] = 2 * diffusion_coeff
        A[i, i+1] = -diffusion_coeff + drift_coeff
    end
    A[1, 1] = 1
    A[end, end] = 1
    b[end] = 1
    for i in 1:(n-1)
        pivot = A[i, i]
        for j in (i+1):n
            factor = A[j, i] / pivot
            A[j, i:n] .= A[j, i:n] .- factor .* A[i, i:n]
            b[j] -= factor * b[i]
        end
    end
    steady_state = zeros(n)
    for i in n:-1:1
        steady_state[i] = (b[i] - sum(A[i, i+1:n] .* steady_state[i+1:n])) / A[i, i]
    end
    steady_state ./= sum(steady_state * dx)
    return steady_state
end