"""
    calculate_action(path, drift, diffusion)

Calculates the action functional for a given path in a stochastic system.

# Arguments
- path::Vector{Float64}: The trajectory of the system.
- drift::Function: Drift function, `drift(x)`.
- diffusion::Float64: Diffusion coefficient.

# Returns
- action::Float64: The action value.
"""
function calculate_action(path::Vector{Float64}, drift::Function, diffusion::Float64)
    action = 0.0
    for i in 1:(length(path) - 1)
        dx = path[i+1] - path[i]
        drift_term = (drift(path[i])^2 / (2 * diffusion)) - drift(path[i]) * dx / diffusion
        action += (dx^2 / (2 * diffusion)) + drift_term
    end
    return action
end

"""
    minimize_action(initial_path, drift, diffusion, max_iters, tol)

Minimizes the action functional for a given path using manual optimization.

# Arguments
- initial_path::Vector{Float64}: Initial guess for the path.
- drift::Function: Drift function, `drift(x)`.
- diffusion::Float64: Diffusion coefficient.
- max_iters::Int: Maximum number of iterations for optimization.
- tol::Float64: Tolerance for convergence.

# Returns
- optimal_path::Vector{Float64}: Path that minimizes the action functional.
"""
function minimize_action(
    initial_path::Vector{Float64},
    drift::Function,
    diffusion::Float64,
    max_iters::Int,
    tol::Float64
)
    path = copy(initial_path)
    action_prev = calculate_action(path, drift, diffusion)
    n = length(path)

    for iter in 1:max_iters
        grad = zeros(n)
        for i in 2:(n-1)
            grad[i] = (path[i+1] - 2 * path[i] + path[i-1]) / diffusion -
                      (drift(path[i]) - (path[i+1] - path[i]) / diffusion) *
                      (path[i+1] - path[i]) / diffusion
        end

        step_size = 0.01
        path[2:(n-1)] -= step_size * grad[2:(n-1)]

        action_curr = calculate_action(path, drift, diffusion)
        if abs(action_curr - action_prev) < tol
            break
        end
        action_prev = action_curr
    end
    return path
end
