@testset "Constant" begin
    n = 5
    X = Y = Z = LinRange(0, 1, n)
    val = rand()
    F = fill(val, (n, n, n))

    t = Tricubic(X, Y, Z, F)

    Xtest = Ytest = Ztest = LinRange(0, 1, 4*n)
    for x in Xtest
        for y in Ytest
            for z in Ztest
                @test t(x, y, z) ≈ val atol=1e-12
                @test partial_derivative_x(t, x, y, z) ≈ 0 atol=1e-12
                @test partial_derivative_y(t, x, y, z) ≈ 0 atol=1e-12
                @test partial_derivative_z(t, x, y, z) ≈ 0 atol=1e-12
            end
        end
    end
end

@testset "Constant uneven grid" begin
    n = 10
    X = LinRange(0, 1, n - 1)
    Y = LinRange(0, 1, n)
    Z = LinRange(0, 1, n + 1)
    val = rand()
    F = fill(val, (n - 1, n, n + 1))

    t = Tricubic(X, Y, Z, F)

    Xtest = Ytest = Ztest = LinRange(0, 1, 4*n)
    for x in Xtest
        for y in Ytest
            for z in Ztest
                @test t(x, y, z) ≈ val atol=1e-12
                @test partial_derivative_x(t, x, y, z) ≈ 0 atol=1e-12
                @test partial_derivative_y(t, x, y, z) ≈ 0 atol=1e-12
                @test partial_derivative_z(t, x, y, z) ≈ 0 atol=1e-12
            end
        end
    end
end

@testset "Cube corners" begin
    n = 11
    X = Y = Z = LinRange(0, 1, n)
    F = rand(Float64, (n, n, n))

    t = Tricubic(X, Y, Z, F)

    for i = 1:length(X)
        for j = 1:length(Y)
            for k = 1:length(Z)
                @test t(X[i], Y[j], Z[k]) ≈ F[i, j, k] atol=1e-12
            end
        end
    end
end

@testset "Linear" begin
    f(x, y, z) = 1 + 2*x + 3*y - 4*z

    n = 11
    X = Y = Z = LinRange(-5, 5, n)
    F = [f(x, y, z) for x=X, y=Y, z=Z]

    t = Tricubic(X, Y, Z, F)

    Xtest = Ytest = Ztest = LinRange(0, 1, 3*n)
    for x in Xtest
        for y in Ytest
            for z in Ztest
                @test t(x, y, z) ≈ f(x, y, z) atol=1e-10
                @test partial_derivative_x(t, x, y, z) ≈ 2 atol=1e-10
                @test partial_derivative_y(t, x, y, z) ≈ 3 atol=1e-10
                @test partial_derivative_z(t, x, y, z) ≈ -4 atol=1e-10
            end
        end
    end
end

@testset "Cubic" begin
    f(x, y, z) = - 10 + 0.1*x - 2*x^2 + x^3 - 5*y^2

    n = 11
    X = Y = Z = LinRange(-5, 5, n)
    F = [f(x, y, z) for x=X, y=Y, z=Z]

    t = Tricubic(X, Y, Z, F)

    Xtest = Ytest = Ztest = LinRange(-2, -0.5, 3*n)
    for x in Xtest
        for y in Ytest
            for z in Ztest
                @test t(x, y, z) ≈ f(x, y, z) atol=1e-10
                @test partial_derivative_x(t, x, y, z) ≈ 0.1 - 4*x + 3*x^2 atol=1e-10
                @test partial_derivative_y(t, x, y, z) ≈ -10*y atol=1e-10
                @test partial_derivative_z(t, x, y, z) ≈ 0 atol=1e-10
            end
        end
    end
end

@testset "Ricker wavelet-like" begin
    f(x, y, z) = x*exp(- x^2 - y^2 - z^2)

    n = 51
    X = Y = Z = LinRange(-2, 2, n)
    F = [f(x, y, z) for x=X, y=Y, z=Z]

    t = Tricubic(X, Y, Z, F)

    x, y, z = rand(Float64, 3)
    @test t(x, y, z) ≈ f(x, y, z) atol=1e-4
    @test partial_derivative_x(t, x, y, z) ≈ (1 - 2*x^2)*exp(- x^2 - y^2 - z^2) atol=1e-4
    @test partial_derivative_y(t, x, y, z) ≈ -2*x*y*exp(- x^2 - y^2 - z^2) atol=1e-4
    @test partial_derivative_z(t, x, y, z) ≈ -2*x*z*exp(- x^2 - y^2 - z^2) atol=1e-4
end
