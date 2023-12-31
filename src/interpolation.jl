mutable struct Tricubic{V₁ <: AbstractVector, V₂ <: AbstractVector,
                        V₃ <: AbstractVector, A <: AbstractArray}
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
    const Xmin::Float64
    const Xmax::Float64
    const Ymin::Float64
    const Ymax::Float64
    const Zmin::Float64
    const Zmax::Float64
    Xᵢ::Float64
    Xᵢ₊₁::Float64
    Yⱼ::Float64
    Yⱼ₊₁::Float64
    Zₖ::Float64
    Zₖ₊₁::Float64
    α::Vector{Float64}

    function Tricubic(X, Y, Z, F)
        dims = size(F)
        ∂F∂X = zeros(dims)
        calculate_∂F∂X!(∂F∂X, F, X)
        ∂F∂Y = zeros(dims)
        calculate_∂F∂Y!(∂F∂Y, F, Y)
        ∂F∂Z = zeros(dims)
        calculate_∂F∂Z!(∂F∂Z, F, Z)
        ∂²F∂X∂Y = zeros(dims)
        calculate_∂F∂X!(∂²F∂X∂Y, ∂F∂Y, X)
        ∂²F∂X∂Z = zeros(dims)
        calculate_∂F∂X!(∂²F∂X∂Z, ∂F∂Z, X)
        ∂²F∂Y∂Z = zeros(dims)
        calculate_∂F∂Y!(∂²F∂Y∂Z, ∂F∂Z, Y)
        ∂³F∂X∂Y∂Z = zeros(dims)
        calculate_∂F∂X!(∂³F∂X∂Y∂Z, ∂²F∂Y∂Z, X)
        Xmin = minimum(X)
        Xmax = maximum(X)
        Ymin = minimum(Y)
        Ymax = maximum(Y)
        Zmin = minimum(Z)
        Zmax = maximum(Z)
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
                Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                Xᵢ, Xᵢ₊₁, Yⱼ, Yⱼ₊₁, Zₖ, Zₖ₊₁, α)
    end
end

"Calculate aspects of cube at `(x, y, z)` for interpolation."
function calculate_cube!(tricubic::Tricubic, x, y, z)
    if (x < tricubic.Xmin || x > tricubic.Xmax
        || y < tricubic.Ymin || y > tricubic.Ymax
        || z < tricubic.Zmin || z > tricubic.Zmax)
        throw(ArgumentError("(x, y, z) inputs are outside grid"))
    end

    i = findfirst(a -> a ≥ x, tricubic.X)
    j = findfirst(a -> a ≥ y, tricubic.Y)
    k = findfirst(a -> a ≥ z, tricubic.Z)
    if i ≠ 1
        i -= 1
    end
    if j ≠ 1
        j -= 1
    end
    if k ≠ 1
        k -= 1
    end

    calculate_coefficients!(
            tricubic.α, i, j, k, tricubic.X, tricubic.Y, tricubic.Z,
            tricubic.F, tricubic.∂F∂X, tricubic.∂F∂Y, tricubic.∂F∂Z,
            tricubic.∂²F∂X∂Y, tricubic.∂²F∂X∂Z, tricubic.∂²F∂Y∂Z,
            tricubic.∂³F∂X∂Y∂Z)

    tricubic.Xᵢ, tricubic.Xᵢ₊₁ = tricubic.X[i], tricubic.X[i + 1]
    tricubic.Yⱼ, tricubic.Yⱼ₊₁ = tricubic.Y[j], tricubic.Y[j + 1]
    tricubic.Zₖ, tricubic.Zₖ₊₁ = tricubic.Z[k], tricubic.Z[k + 1]

    (tricubic.α, tricubic.Xᵢ, tricubic.Xᵢ₊₁, tricubic.Yⱼ, tricubic.Yⱼ₊₁,
     tricubic.Zₖ, tricubic.Zₖ₊₁)
end

"Tricubic interpolator."
function (tricubic::Tricubic)(x, y, z)
    if !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁ && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
         && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁)
        calculate_cube!(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1.0, η, η^2, η^3)
    ζarray = (1.0, ζ, ζ^2, ζ^3)

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
    if !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁ && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
         && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁)
        calculate_cube!(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ξarray = (1.0, ξ, ξ^2)
    ζarray = (1.0, ζ, ζ^2, ζ^3)

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
    if !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁ && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
         && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁)
        calculate_cube!(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1.0, η, η^2)
    ζarray = (1.0, ζ, ζ^2, ζ^3)

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
    if !(tricubic.Xᵢ ≤ x < tricubic.Xᵢ₊₁ && tricubic.Yⱼ ≤ y < tricubic.Yⱼ₊₁
         && tricubic.Zₖ ≤ z < tricubic.Zₖ₊₁)
        calculate_cube!(tricubic, x, y, z)
    end

    ξ = (x - tricubic.Xᵢ)/(tricubic.Xᵢ₊₁ - tricubic.Xᵢ)
    η = (y - tricubic.Yⱼ)/(tricubic.Yⱼ₊₁ - tricubic.Yⱼ)
    ζ = (z - tricubic.Zₖ)/(tricubic.Zₖ₊₁ - tricubic.Zₖ)

    ηarray = (1.0, η, η^2, η^3)
    ζarray = (1.0, ζ, ζ^2)

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
