"""
    analyze_time_dependent_entropy_production(J_func::Function, X_func::Function, t_span::Tuple{Float64, Float64}, dt::Float64) -> Vector{Float64}

Analyzes time-dependent entropy production based on time-varying fluxes and forces.
- J_func, X_func: Functions that return J and X matrices as a function of time.
- t_span: Tuple specifying the start and end time for the analysis.
- dt: Time step for the manual ODE integration.
Returns a vector of entropy production values over time.
"""
function analyze_time_dependent_entropy_production(
    J_func::Function,
    X_func::Function,
    t_span::Tuple{Float64, Float64},
    dt::Float64
)
    t_start, t_end = t_span
    time_points = collect(t_start:dt:t_end)
    σ_values = zeros(length(time_points))
    σ_values[1] = 0.0  # Initial entropy production

    for i in 2:length(time_points)
        t = time_points[i]
        J = J_func(t)
        X = X_func(t)
        dσ = sum(J .* X) * dt
        σ_values[i] = σ_values[i-1] + dσ
    end

    return time_points, σ_values
end

"""
    visualize_entropy_production(time_points::Vector{Float64}, σ_values::Vector{Float64}, title_name::String)

Creates a plot to visualize entropy production over time.
- time_points: Vector of time points.
- σ_values: Vector of entropy production values corresponding to the time points.
"""
function visualize_entropy_production(
    time_points::Vector{Float64},
    σ_values::Vector{Float64},
    title_name::String
)
    plot(
        time_points,
        σ_values,
        xlabel="Time",
        ylabel="Entropy Production",
        title=title_name,
        legend=false,
        color=:blue
    )
end
