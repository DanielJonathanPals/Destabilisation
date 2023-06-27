using Polynomials

@testset "Time Scales" begin

    f(x) = x + 10*x^2 + 3*x^3 + x^5
    @test coeffs(toPolynomial(f, deg = 5)) ≈ [0.,1.,10.,3.,0.,1.] atol = 1e-10
   

    Random.seed!(3)
    x = zeros(1,100)
    for i in 4:100
        x[1,i] = 0.5*x[1,i-1] + 0.01*x[1,i-2] + 0.3*x[1,i-3] + 0.02*randn()
    end

    model1 = fitVARmodel(x, p = 3)

    @test timeScale(model1) ≈ 1.170210112086673 atol = 1e-5

    model2 = fitVARmodel(x, p = 1, h = x -> x.^2)

    @test_throws ErrorException timeScale(model2)

    model3 = fitVARmodel(x, p = 1, p_traj = x.^2)

    @test_throws ErrorException timeScale(model3)
end