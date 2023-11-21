mutable struct Tricubic{V₁<:AbstractVector, V₂<:AbstractVector,
                        V₃<:AbstractVector, A<:AbstractArray}
    const X::V₁
    const Y::V₂
    const Z::V₃
    const F::A
    const ∂F∂X::A
    const ∂F∂Y::A
    const ∂F∂Z::A
    const ∂²F∂X∂Y::A
    const ∂²F∂X∂Z::A
    const ∂²F∂Y∂Z::A
    const ∂³F∂X∂Y∂Z::A
    initialised::Bool
    Xᵢ::Float64
    Xᵢ₊₁::Float64
    Yⱼ::Float64
    Yⱼ₊₁::Float64
    Zₖ::Float64
    Zₖ₊₁::Float64
    α::Vector{Float64}

    function Tricubic(X, Y, Z, F)
        ∂F∂X = build_∂F∂X(F, X)
        ∂F∂Y = build_∂F∂Y(F, Y)
        ∂F∂Z = build_∂F∂Z(F, Z)
        ∂²F∂X∂Y = build_∂F∂X(∂F∂Y, X)
        ∂²F∂X∂Z = build_∂F∂X(∂F∂Z, X)
        ∂²F∂Y∂Z = build_∂F∂Y(∂F∂Z, Y)
        ∂³F∂X∂Y∂Z = build_∂F∂X(∂²F∂Y∂Z, X)
        initialised = false
        Xᵢ = 0.0
        Xᵢ₊₁ = 0.0
        Yⱼ = 0.0
        Yⱼ₊₁ = 0.0
        Zₖ = 0.0
        Zₖ₊₁ = 0.0
        α = zeros(64)
        new{typeof(X), typeof(Y), typeof(Z), typeof(F)}(
                X, Y, Z, F, ∂F∂X, ∂F∂Y, ∂F∂Z,
                ∂²F∂X∂Y, ∂²F∂X∂Z, ∂²F∂Y∂Z, ∂³F∂X∂Y∂Z,
                initialised, Xᵢ, Xᵢ₊₁, Yⱼ, Yⱼ₊₁, Zₖ, Zₖ₊₁, α)
    end
end

"Calculate aspects of cube at `(x, y, z)` for interpolation."
function calculate_cube(tricubic::Tricubic, x, y, z)
    i = j = k = 1
    for I = 1:length(tricubic.X)
        if tricubic.X[I] ≥ x
            i = I
            break
        end
    end
    for J = 1:length(tricubic.Y)
        if tricubic.Y[J] ≥ y
            j = J
            break
        end
    end
    for K = 1:length(tricubic.Z)
        if tricubic.Z[K] ≥ z
            k = K
            break
        end
    end
    if i ≠ 1
        i -= 1
    end
    if j ≠ 1
        j -= 1
    end
    if k ≠ 1
        k -= 1
    end

    tricubic.α = calculate_coefficients(
            i, j, k, tricubic.X, tricubic.Y, tricubic.Z,
            tricubic.F, tricubic.∂F∂X, tricubic.∂F∂Y, tricubic.∂F∂Z,
            tricubic.∂²F∂X∂Y, tricubic.∂²F∂X∂Z, tricubic.∂²F∂Y∂Z,
            tricubic.∂³F∂X∂Y∂Z)

    tricubic.initialised = true

    tricubic.Xᵢ, tricubic.Xᵢ₊₁ = tricubic.X[i], tricubic.X[i + 1]
    tricubic.Yⱼ, tricubic.Yⱼ₊₁ = tricubic.Y[j], tricubic.Y[j + 1]
    tricubic.Zₖ, tricubic.Zₖ₊₁ = tricubic.Z[k], tricubic.Z[k + 1]
end

"Tricubic interpolator."
function (tricubic::Tricubic)(x, y, z)
    if (!tricubic.initialised || !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁
                                   && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
                                   && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁))
        calculate_cube(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1, η, η^2, η^3)
    ζarray = (1, ζ, ζ^2, ζ^3)

    f = 0.0
    for c = 1:4
        for d = 1:4
            f += ((tricubic.α[1 + 4*c + 16*d - 20]
                   + ξ*(tricubic.α[2 + 4*c + 16*d - 20]
                        + ξ*(tricubic.α[3 + 4*c + 16*d - 20]
                             + ξ*tricubic.α[4 + 4*c + 16*d - 20])))
                  *ηarray[c]*ζarray[d])
        end
    end
    f
end

function partial_derivative_x(tricubic::Tricubic, x, y, z)
    if (!tricubic.initialised || !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁
                                   && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
                                   && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁))
        calculate_cube(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ξarray = (1, ξ, ξ^2)
    ζarray = (1, ζ, ζ^2, ζ^3)

    ∂f∂x = 0.0
    for a = 2:4
        for d = 1:4
            ∂f∂x += ((tricubic.α[a + 4 + 16*d - 20]
                      + η*(tricubic.α[a + 8 + 16*d - 20]
                           + η*(tricubic.α[a + 12 + 16*d - 20]
                                + η*tricubic.α[a + 16 + 16*d - 20])))
                     *(a - 1)*ξarray[a - 1]/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
                     *ζarray[d])
        end
    end
    ∂f∂x
end

function partial_derivative_y(tricubic::Tricubic, x, y, z)
    if (!tricubic.initialised || !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁
                                   && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
                                   && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁))
        calculate_cube(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1, η, η^2)
    ζarray = (1, ζ, ζ^2, ζ^3)

    ∂f∂y = 0.0
    for c = 2:4
        for d = 1:4
            ∂f∂y += ((tricubic.α[1 + 4*c + 16*d - 20]
                      + ξ*(tricubic.α[2 + 4*c + 16*d - 20]
                           + ξ*(tricubic.α[3 + 4*c + 16*d - 20]
                                + ξ*tricubic.α[4 + 4*c + 16*d - 20])))
                     *(c - 1)*ηarray[c - 1]/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
                     *ζarray[d])
        end
    end
    ∂f∂y
end

function partial_derivative_z(tricubic::Tricubic, x, y, z)
    if (!tricubic.initialised || !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁
                                   && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
                                   && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁))
        calculate_cube(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1, η, η^2, η^3)
    ζarray = (1, ζ, ζ^2)

    ∂f∂z = 0.0
    for c = 1:4
        for d = 2:4
            ∂f∂z += ((tricubic.α[1 + 4*c + 16*d - 20]
                      + ξ*(tricubic.α[2 + 4*c + 16*d - 20]
                           + ξ*(tricubic.α[3 + 4*c + 16*d - 20]
                                + ξ*tricubic.α[4 + 4*c + 16*d - 20])))
                     *ηarray[c]
                     *(d - 1)*ζarray[d - 1]/(tricubic.Zₖ₊₁ - tricubic.Zₖ))
        end
    end
    ∂f∂z
end

function Base.show(io::IO, mime::MIME"text/plain", tricubic::Tricubic)
    println(io, "Tricubic interpolator")
    print(io, "X: ")
    show(io, mime, tricubic.X)
    println(io)
    print(io, "Y: ")
    show(io, mime, tricubic.Y)
    println(io)
    print(io, "Z: ")
    show(io, mime, tricubic.Z)
    println(io)
    i, j, k = size(tricubic.F)
    print(io, "F: $(i)x$(j)x$(k) ", typeof(tricubic.F))
end
