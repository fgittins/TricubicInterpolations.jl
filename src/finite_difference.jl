"Five-point finite difference for partial derivative with respect to `x`."
function build_∂F∂X(F, X)
    ∂F∂X = zero(F)
    ∂F∂X[1, :, :] = ((-25*F[1, :, :] + 48*F[2, :, :] - 36*F[3, :, :]
                      + 16*F[4, :, :] - 3*F[5, :, :])
                     / (-25*X[1] + 48*X[2] - 36*X[3] + 16*X[4] - 3*X[5]))
    ∂F∂X[2, :, :] = ((-3*F[1, :, :] - 10*F[2, :, :] + 18*F[3, :, :]
                      - 6*F[4, :, :] + F[5, :, :])
                     / (-3*X[1] - 10*X[2] + 18*X[3] - 6*X[4] + X[5]))
    for i = 3:length(X) - 2
        ∂F∂X[i, :, :] = ((F[i - 2, :, :] - 8*F[i - 1, :, :] + 8*F[i + 1, :, :]
                          - F[i + 2, :, :])
                         / (X[i - 2] - 8*X[i - 1] + 8*X[i + 1] - X[i + 2]))
    end
    ∂F∂X[end - 1, :, :] = ((3*F[end - 4, :, :] + 10*F[end - 3, :, :]
                            - 18*F[end - 2, :, :] + 6*F[end - 1, :, :]
                            - F[end, :, :]) 
                           / (3*X[end - 4] + 10*X[end - 3] - 18*X[end - 2] 
                              + 6*X[end - 1] - X[end]))
    ∂F∂X[end, :, :] = ((25*F[end - 4, :, :] - 48*F[end - 3, :, :]
                        + 36*F[end - 2, :, :] - 16*F[end - 1, :, :]
                        + 3*F[end, :, :])
                       / (25*X[end - 4] - 48*X[end - 3] + 36*X[end - 2]
                          - 16*X[end - 1] + 3*X[end]))
    ∂F∂X
end

"Five-point finite difference for partial derivative with respect to `y`."
function build_∂F∂Y(F, Y)
    ∂F∂Y = zero(F)
    ∂F∂Y[:, 1, :] = ((-25*F[:, 1, :] + 48*F[:, 2, :] - 36*F[:, 3, :]
                      + 16*F[:, 4, :] - 3*F[:, 5, :])
                     / (-25*Y[1] + 48*Y[2] - 36*Y[3] + 16*Y[4] - 3*Y[5]))
    ∂F∂Y[:, 2, :] = ((-3*F[:, 1, :] - 10*F[:, 2, :] + 18*F[:, 3, :]
                      - 6*F[:, 4, :] + F[:, 5, :])
                     / (-3*Y[1] - 10*Y[2] + 18*Y[3] - 6*Y[4] + Y[5]))
    for j = 3:length(Y) - 2
        ∂F∂Y[:, j, :] = ((F[:, j - 2, :] - 8*F[:, j - 1, :] + 8*F[:, j + 1, :]
                          - F[:, j + 2, :])
                         / (Y[j - 2] - 8*Y[j - 1] + 8*Y[j + 1] - Y[j + 2]))
    end
    ∂F∂Y[:, end - 1, :] = ((3*F[:, end - 4, :] + 10*F[:, end - 3, :]
                            - 18*F[:, end - 2, :] + 6*F[:, end - 1, :]
                            - F[:, end, :]) 
                           / (3*Y[end - 4] + 10*Y[end - 3] - 18*Y[end - 2] 
                              + 6*Y[end - 1] - Y[end]))
    ∂F∂Y[:, end, :] = ((25*F[:, end - 4, :] - 48*F[:, end - 3, :]
                        + 36*F[:, end - 2, :] - 16*F[:, end - 1, :]
                        + 3*F[:, end, :])
                       / (25*Y[end - 4] - 48*Y[end - 3] + 36*Y[end - 2]
                          - 16*Y[end - 1] + 3*Y[end]))
    ∂F∂Y
end

"Five-point finite difference for partial derivative with respect to `z`."
function build_∂F∂Z(F, Z)
    ∂F∂Z = zero(F)
    ∂F∂Z[:, :, 1] = ((-25*F[:, :, 1] + 48*F[:, :, 2] - 36*F[:, :, 3]
                      + 16*F[:, :, 4] - 3*F[:, :, 5])
                     / (-25*Z[1] + 48*Z[2] - 36*Z[3] + 16*Z[4] - 3*Z[5]))
    ∂F∂Z[:, :, 2] = ((-3*F[:, :, 1] - 10*F[:, :, 2] + 18*F[:, :, 3]
                      - 6*F[:, :, 4] + F[:, :, 5])
                     / (-3*Z[1] - 10*Z[2] + 18*Z[3] - 6*Z[4] + Z[5]))
    for k = 3:length(Z) - 2
        ∂F∂Z[:, :, k] = ((F[:, :, k - 2] - 8*F[:, :, k - 1] + 8*F[:, :, k + 1]
                          - F[:, :, k + 2])
                         / (Z[k - 2] - 8*Z[k - 1] + 8*Z[k + 1] - Z[k + 2]))
    end
    ∂F∂Z[:, :, end - 1] = ((3*F[:, :, end - 4] + 10*F[:, :, end - 3]
                            - 18*F[:, :, end - 2] + 6*F[:, :, end - 1]
                            - F[:, :, end]) 
                           / (3*Z[end - 4] + 10*Z[end - 3] - 18*Z[end - 2] 
                              + 6*Z[end - 1] - Z[end]))
    ∂F∂Z[:, :, end] = ((25*F[:, :, end - 4] - 48*F[:, :, end - 3]
                        + 36*F[:, :, end - 2] - 16*F[:, :, end - 1]
                        + 3*F[:, :, end])
                       / (25*Z[end - 4] - 48*Z[end - 3] + 36*Z[end - 2]
                          - 16*Z[end - 1] + 3*Z[end]))
    ∂F∂Z
end
