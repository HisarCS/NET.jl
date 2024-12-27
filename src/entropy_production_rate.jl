"""
    compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                         J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64) -> Float64

Calculates the rate of irreversible entropy production by comparing entropy production at two different states over time Δt.
- J1, X1: 2D arrays of fluxes and thermodynamic forces at the initial state.
- J2, X2: 2D arrays of fluxes and thermodynamic forces at the final state.
- Δt: Time interval between the two states.
Returns: The rate of entropy production in the system per unit time.
"""
function compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                              J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64)
    validate_dimensions_entropy(J1, X1)
    validate_dimensions_entropy(J2, X2)

    rows1, cols1 = size(J1)
    rows2, cols2 = size(J2)

    entropy_initial = 0.0
    for i in 1:rows1
        for j in 1:cols1
            entropy_initial += J1[i, j] * X1[i, j]
        end
    end

    entropy_final = 0.0
    for i in 1:rows2
        for j in 1:cols2
            entropy_final += J2[i, j] * X2[i, j]
        end
    end

    Δentropy = entropy_final - entropy_initial
    entropy_rate = Δentropy / Δt

    log_info("Entropy rate calculation complete: rate = $(entropy_rate) per unit time.")
    return entropy_rate
end

"""
    validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})

Ensures that the input dimensions of the fluxes (J) and forces (X) match.
"""
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    rows_J, cols_J = size(J)
    rows_X, cols_X = size(X)
    if rows_J != rows_X || cols_J != cols_X
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
end
