"Five-point finite difference for partial derivative with respect to `x`."
function calculate_∂F∂X!(∂F∂X, F, X)
    dims = size(∂F∂X)
    for k = 1:dims[3], j = 1:dims[2]
        ∂F∂X[1, j, k] = ((-25*F[1, j, k] + 48*F[2, j, k] - 36*F[3, j, k]
                          + 16*F[4, j, k] - 3*F[5, j, k])
                         / (-25*X[1] + 48*X[2] - 36*X[3] + 16*X[4] - 3*X[5]))
        ∂F∂X[2, j, k] = ((-3*F[1, j, k] - 10*F[2, j, k] + 18*F[3, j, k]
                          - 6*F[4, j, k] + F[5, j, k])
                         / (-3*X[1] - 10*X[2] + 18*X[3] - 6*X[4] + X[5]))
        for i = 3:dims[1] - 2
            ∂F∂X[i, j, k] = ((F[i - 2, j, k] - 8*F[i - 1, j, k]
                              + 8*F[i + 1, j, k] - F[i + 2, j, k])
                             / (X[i - 2] - 8*X[i - 1] + 8*X[i + 1] - X[i + 2]))
        end
        ∂F∂X[end - 1, j, k] = ((3*F[end - 4, j, k] + 10*F[end - 3, j, k]
                                - 18*F[end - 2, j, k] + 6*F[end - 1, j, k]
                                - F[end, j, k]) 
                               / (3*X[end - 4] + 10*X[end - 3] - 18*X[end - 2]
                                  + 6*X[end - 1] - X[end]))
        ∂F∂X[end, j, k] = ((25*F[end - 4, j, k] - 48*F[end - 3, j, k]
                            + 36*F[end - 2, j, k] - 16*F[end - 1, j, k]
                            + 3*F[end, j, k])
                           / (25*X[end - 4] - 48*X[end - 3] + 36*X[end - 2]
                              - 16*X[end - 1] + 3*X[end]))
    end
   ∂F∂X
end

"Five-point finite difference for partial derivative with respect to `y`."
function calculate_∂F∂Y!(∂F∂Y, F, Y)
    dims = size(∂F∂Y)
    for k = 1:dims[3], i = 1:dims[1]
        ∂F∂Y[i, 1, k] = ((-25*F[i, 1, k] + 48*F[i, 2, k] - 36*F[i, 3, k]
                          + 16*F[i, 4, k] - 3*F[i, 5, k])
                         / (-25*Y[1] + 48*Y[2] - 36*Y[3] + 16*Y[4] - 3*Y[5]))
        ∂F∂Y[i, 2, k] = ((-3*F[i, 1, k] - 10*F[i, 2, k] + 18*F[i, 3, k]
                          - 6*F[i, 4, k] + F[i, 5, k])
                         / (-3*Y[1] - 10*Y[2] + 18*Y[3] - 6*Y[4] + Y[5]))
        for j = 3:dims[2] - 2
            ∂F∂Y[i, j, k] = ((F[i, j - 2, k] - 8*F[i, j - 1, k]
                              + 8*F[i, j + 1, k] - F[i, j + 2, k])
                             / (Y[j - 2] - 8*Y[j - 1] + 8*Y[j + 1] - Y[j + 2]))
        end
        ∂F∂Y[i, end - 1, k] = ((3*F[i, end - 4, k] + 10*F[i, end - 3, k]
                                - 18*F[i, end - 2, k] + 6*F[i, end - 1, k]
                                - F[i, end, k]) 
                               / (3*Y[end - 4] + 10*Y[end - 3] - 18*Y[end - 2]
                                  + 6*Y[end - 1] - Y[end]))
        ∂F∂Y[i, end, k] = ((25*F[i, end - 4, k] - 48*F[i, end - 3, k]
                            + 36*F[i, end - 2, k] - 16*F[i, end - 1, k]
                            + 3*F[i, end, k])
                           / (25*Y[end - 4] - 48*Y[end - 3] + 36*Y[end - 2]
                              - 16*Y[end - 1] + 3*Y[end]))
    end
    ∂F∂Y
end

"Five-point finite difference for partial derivative with respect to `z`."
function calculate_∂F∂Z!(∂F∂Z, F, Z)
    dims = size(∂F∂Z)
    for j = 1:dims[2], i = 1:dims[1]
        ∂F∂Z[i, j, 1] = ((-25*F[i, j, 1] + 48*F[i, j, 2] - 36*F[i, j, 3]
                          + 16*F[i, j, 4] - 3*F[i, j, 5])
                         / (-25*Z[1] + 48*Z[2] - 36*Z[3] + 16*Z[4] - 3*Z[5]))
        ∂F∂Z[i, j, 2] = ((-3*F[i, j, 1] - 10*F[i, j, 2] + 18*F[i, j, 3]
                          - 6*F[i, j, 4] + F[i, j, 5])
                         / (-3*Z[1] - 10*Z[2] + 18*Z[3] - 6*Z[4] + Z[5]))
        for k = 3:dims[3] - 2
            ∂F∂Z[i, j, k] = ((F[i, j, k - 2] - 8*F[i, j, k - 1]
                              + 8*F[i, j, k + 1] - F[i, j, k + 2])
                             / (Z[k - 2] - 8*Z[k - 1] + 8*Z[k + 1] - Z[k + 2]))
        end
        ∂F∂Z[i, j, end - 1] = ((3*F[i, j, end - 4] + 10*F[i, j, end - 3]
                                - 18*F[i, j, end - 2] + 6*F[i, j, end - 1]
                                - F[i, j, end]) 
                               / (3*Z[end - 4] + 10*Z[end - 3] - 18*Z[end - 2]
                                  + 6*Z[end - 1] - Z[end]))
        ∂F∂Z[i, j, end] = ((25*F[i, j, end - 4] - 48*F[i, j, end - 3]
                            + 36*F[i, j, end - 2] - 16*F[i, j, end - 1]
                            + 3*F[i, j, end])
                           / (25*Z[end - 4] - 48*Z[end - 3] + 36*Z[end - 2]
                              - 16*Z[end - 1] + 3*Z[end]))
    end
    ∂F∂Z
end
