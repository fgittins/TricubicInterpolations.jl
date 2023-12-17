"Calculate vector of coefficients for interpolation."
function calculate_coefficients!(α, i, j, k, X, Y, Z, F, ∂F∂X, ∂F∂Y, ∂F∂Z,
                                 ∂²F∂X∂Y, ∂²F∂X∂Z, ∂²F∂Y∂Z, ∂³F∂X∂Y∂Z)
    B⁻¹ = SA[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -3 3 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           2 -2 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           9 -9 -9 0 9 0 0 0 6 3 -6 0 -3 0 0 0 6 -6 3 0 -3 0 0 0 0 0 0 0 0 0 0 0 4 2 2 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -6 6 6 0 -6 0 0 0 -3 -3 3 0 3 0 0 0 -4 4 -2 0 2 0 0 0 0 0 0 0 0 0 0 0 -2 -2 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -6 6 6 0 -6 0 0 0 -4 -2 4 0 2 0 0 0 -3 3 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 -2 -1 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           4 -4 -4 0 4 0 0 0 2 2 -2 0 -2 0 0 0 2 -2 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 3 0 0 0 0 0 0 -2 -1 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -2 0 0 0 0 0 0 1 1 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 -9 -9 0 9 0 0 0 0 0 0 0 0 0 0 0 6 3 -6 0 -3 0 0 0 6 -6 3 0 -3 0 0 0 4 2 2 0 1 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 6 0 -6 0 0 0 0 0 0 0 0 0 0 0 -3 -3 3 0 3 0 0 0 -4 4 -2 0 2 0 0 0 -2 -2 -1 0 -1 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 6 0 -6 0 0 0 0 0 0 0 0 0 0 0 -4 -2 4 0 2 0 0 0 -3 3 -3 0 3 0 0 0 -2 -1 -2 0 -1 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 -4 -4 0 4 0 0 0 0 0 0 0 0 0 0 0 2 2 -2 0 -2 0 0 0 2 -2 2 0 -2 0 0 0 1 1 1 0 1 0 0 0
           -3 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 -3 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           9 -9 0 -9 0 9 0 0 6 3 0 -6 0 -3 0 0 0 0 0 0 0 0 0 0 6 -6 0 3 0 -3 0 0 0 0 0 0 0 0 0 0 4 2 0 2 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -6 6 0 6 0 -6 0 0 -3 -3 0 3 0 3 0 0 0 0 0 0 0 0 0 0 -4 4 0 -2 0 2 0 0 0 0 0 0 0 0 0 0 -2 -2 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 0 -1 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 -9 0 -9 0 9 0 0 0 0 0 0 0 0 0 0 6 3 0 -6 0 -3 0 0 0 0 0 0 0 0 0 0 6 -6 0 3 0 -3 0 0 4 2 0 2 0 1 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 0 -3 -3 0 3 0 3 0 0 0 0 0 0 0 0 0 0 -4 4 0 -2 0 2 0 0 -2 -2 0 -1 0 -1 0 0
           9 0 -9 -9 0 0 9 0 0 0 0 0 0 0 0 0 6 0 3 -6 0 0 -3 0 6 0 -6 3 0 0 -3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 2 2 0 0 1 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 9 0 -9 -9 0 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 0 3 -6 0 0 -3 0 6 0 -6 3 0 0 -3 0 0 0 0 0 0 0 0 0 4 0 2 2 0 0 1 0
           -27 27 27 27 -27 -27 -27 27 -18 -9 18 18 9 9 -18 -9 -30 30 -9 30 9 -30 9 -9 -18 18 18 -9 -18 9 9 -9 -18 -12 -6 18 -3 12 6 3 -12 -6 12 -6 6 -3 6 3 -20 20 -6 -10 6 10 -3 3 -12 -8 -4 -6 -2 -4 -2 -1
           18 -18 -18 -18 18 18 18 -18 9 9 -9 -9 -9 -9 9 9 24 -24 6 -24 -6 24 -6 6 12 -12 -12 6 12 -6 -6 6 12 12 3 -12 3 -12 -3 -3 6 6 -6 3 -6 3 -3 -3 16 -16 4 8 -4 -8 2 -2 8 8 2 4 2 4 1 1
           -6 0 6 6 0 0 -6 0 0 0 0 0 0 0 0 0 -3 0 -3 3 0 0 3 0 -4 0 4 -2 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -2 -1 0 0 -1 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 -6 0 6 6 0 0 -6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3 0 -3 3 0 0 3 0 -4 0 4 -2 0 0 2 0 0 0 0 0 0 0 0 0 -2 0 -2 -1 0 0 -1 0
           18 -18 -18 -18 18 18 18 -18 12 6 -12 -12 -6 -6 12 6 21 -21 9 -21 -9 21 -9 9 12 -12 -12 6 12 -6 -6 6 12 9 6 -12 3 -9 -6 -3 8 4 -8 4 -4 2 -4 -2 14 -14 6 7 -6 -7 3 -3 8 6 4 4 2 3 2 1
           -12 12 12 12 -12 -12 -12 12 -6 -6 6 6 6 6 -6 -6 -18 18 -6 18 6 -18 6 -6 -8 8 8 -4 -8 4 4 -4 -9 -9 -3 9 -3 9 3 3 -4 -4 4 -2 4 -2 2 2 -12 12 -4 -6 4 6 -2 2 -6 -6 -2 -3 -2 -3 -1 -1
           2 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 2 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           -6 6 0 6 0 -6 0 0 -4 -2 0 4 0 2 0 0 0 0 0 0 0 0 0 0 -3 3 0 -3 0 3 0 0 0 0 0 0 0 0 0 0 -2 -1 0 -2 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           4 -4 0 -4 0 4 0 0 2 2 0 -2 0 -2 0 0 0 0 0 0 0 0 0 0 2 -2 0 2 0 -2 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6 6 0 6 0 -6 0 0 0 0 0 0 0 0 0 0 -4 -2 0 4 0 2 0 0 0 0 0 0 0 0 0 0 -3 3 0 -3 0 3 0 0 -2 -1 0 -2 0 -1 0 0
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 -4 0 -4 0 4 0 0 0 0 0 0 0 0 0 0 2 2 0 -2 0 -2 0 0 0 0 0 0 0 0 0 0 2 -2 0 2 0 -2 0 0 1 1 0 1 0 1 0 0
           -6 0 6 6 0 0 -6 0 0 0 0 0 0 0 0 0 -4 0 -2 4 0 0 2 0 -3 0 3 -3 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 0 -1 -2 0 0 -1 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 -6 0 6 6 0 0 -6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4 0 -2 4 0 0 2 0 -3 0 3 -3 0 0 3 0 0 0 0 0 0 0 0 0 -2 0 -1 -2 0 0 -1 0
           18 -18 -18 -18 18 18 18 -18 12 6 -12 -12 -6 -6 12 6 24 -24 6 -24 -6 24 -6 6 9 -9 -9 9 9 -9 -9 9 14 10 4 -14 2 -10 -4 -2 6 3 -6 6 -3 3 -6 -3 14 -14 3 10 -3 -10 3 -3 8 6 2 6 1 4 2 1
           -12 12 12 12 -12 -12 -12 12 -6 -6 6 6 6 6 -6 -6 -20 20 -4 20 4 -20 4 -4 -6 6 6 -6 -6 6 6 -6 -10 -10 -2 10 -2 10 2 2 -3 -3 3 -3 3 -3 3 3 -12 12 -2 -8 2 8 -2 2 -6 -6 -1 -4 -1 -4 -1 -1
           4 0 -4 -4 0 0 4 0 0 0 0 0 0 0 0 0 2 0 2 -2 0 0 -2 0 2 0 -2 2 0 0 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0
           0 0 0 0 0 0 0 0 4 0 -4 -4 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 2 -2 0 0 -2 0 2 0 -2 2 0 0 -2 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 0
           -12 12 12 12 -12 -12 -12 12 -8 -4 8 8 4 4 -8 -4 -18 18 -6 18 6 -18 6 -6 -6 6 6 -6 -6 6 6 -6 -10 -8 -4 10 -2 8 4 2 -4 -2 4 -4 2 -2 4 2 -11 11 -3 -7 3 7 -3 3 -6 -5 -2 -4 -1 -3 -2 -1
           8 -8 -8 -8 8 8 8 -8 4 4 -4 -4 -4 -4 4 4 16 -16 4 -16 -4 16 -4 4 4 -4 -4 4 4 -4 -4 4 8 8 2 -8 2 -8 -2 -2 2 2 -2 2 -2 2 -2 -2 10 -10 2 6 -2 -6 2 -2 5 5 1 3 1 3 1 1]

    b = SA[F[i, j, k],
         F[i + 1, j, k],
         F[i, j + 1, k],
         F[i, j, k + 1],
         F[i + 1, j + 1, k],
         F[i + 1, j, k + 1],
         F[i, j + 1, k + 1],
         F[i + 1, j + 1, k + 1],
         #####
         (X[i + 1] - X[i])*∂F∂X[i, j, k],
         (X[i + 1] - X[i])*∂F∂X[i + 1, j, k],
         (X[i + 1] - X[i])*∂F∂X[i, j + 1, k],
         (X[i + 1] - X[i])*∂F∂X[i, j, k + 1],
         (X[i + 1] - X[i])*∂F∂X[i + 1, j + 1, k],
         (X[i + 1] - X[i])*∂F∂X[i + 1, j, k + 1],
         (X[i + 1] - X[i])*∂F∂X[i, j + 1, k + 1],
         (X[i + 1] - X[i])*∂F∂X[i + 1, j + 1, k + 1],
         #####
         (Y[j + 1] - Y[j])*∂F∂Y[i, j, k],
         (Y[j + 1] - Y[j])*∂F∂Y[i + 1, j, k],
         (Y[j + 1] - Y[j])*∂F∂Y[i, j + 1, k],
         (Y[j + 1] - Y[j])*∂F∂Y[i, j, k + 1],
         (Y[j + 1] - Y[j])*∂F∂Y[i + 1, j + 1, k],
         (Y[j + 1] - Y[j])*∂F∂Y[i + 1, j, k + 1],
         (Y[j + 1] - Y[j])*∂F∂Y[i, j + 1, k + 1],
         (Y[j + 1] - Y[j])*∂F∂Y[i + 1, j + 1, k + 1],
         #####
         (Z[k + 1] - Z[k])*∂F∂Z[i, j, k],
         (Z[k + 1] - Z[k])*∂F∂Z[i + 1, j, k],
         (Z[k + 1] - Z[k])*∂F∂Z[i, j + 1, k],
         (Z[k + 1] - Z[k])*∂F∂Z[i, j, k + 1],
         (Z[k + 1] - Z[k])*∂F∂Z[i + 1, j + 1, k],
         (Z[k + 1] - Z[k])*∂F∂Z[i + 1, j, k + 1],
         (Z[k + 1] - Z[k])*∂F∂Z[i, j + 1, k + 1],
         (Z[k + 1] - Z[k])*∂F∂Z[i + 1, j + 1, k + 1],
         #####
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i, j, k],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i + 1, j, k],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i, j + 1, k],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i, j, k + 1],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i + 1, j + 1, k],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i + 1, j, k + 1],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i, j + 1, k + 1],
         (X[i + 1] - X[i])*(Y[j + 1] - Y[j])*∂²F∂X∂Y[i + 1, j + 1, k + 1],
         #####
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i, j, k],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i + 1, j, k],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i, j + 1, k],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i, j, k + 1],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i + 1, j + 1, k],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i + 1, j, k + 1],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i, j + 1, k + 1],
         (X[i + 1] - X[i])*(Z[k + 1] - Z[k])*∂²F∂X∂Z[i + 1, j + 1, k + 1],
         #####
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i, j, k],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i + 1, j, k],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i, j + 1, k],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i, j, k + 1],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i + 1, j + 1, k],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i + 1, j, k + 1],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i, j + 1, k + 1],
         (Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])*∂²F∂Y∂Z[i + 1, j + 1, k + 1],
         #####
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i, j, k]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i + 1, j, k]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i, j + 1, k]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i, j, k + 1]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i + 1, j + 1, k]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i + 1, j, k + 1]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i, j + 1, k + 1]),
         ((X[i + 1] - X[i])*(Y[j + 1] - Y[j])*(Z[k + 1] - Z[k])
          *∂³F∂X∂Y∂Z[i + 1, j + 1, k + 1])]
    
    α .= B⁻¹*b
end
