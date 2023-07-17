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
    ts, err = timeScale(model1, include_error = true)
    @test ts ≈ 1.170210112086673 atol = 1e-5
    @test err ≈ 0.09183570 atol = 1e-5

    model2 = fitVARmodel(x, p = 1, h = x -> x.^2)

    @test_throws ErrorException timeScale(model2)

    model3 = fitVARmodel(x, p = 1, p_traj = x.^2)

    @test_throws ErrorException timeScale(model3)

    @test timeScale(model1.B_hat) ≈ 1.170210112086673 atol = 1e-5
    @test_throws ErrorException timeScale(ones(2,4))

    prog(v,p,r) = [-0.01*v[1]^3 + p[1]*v[1]^2 + 1.01*v[1] + p[2] + 1e-4*r[1]]
    obs(v,p) = v
    v_init = [-1.]
    p_init = [0.,0.]
    DS = DynamicalSystem(prog,obs,v_init,p_init;random_vec_length=1)

    ts, err = timeScale(prog,v_init,p_init,random_vec_length=1)
    @test ts ≈ 1.019521301 atol = 1e-5
    @test err ≈ 0.001792830 atol = 1e-5

    ts, err = timeScale(DS)
    @test ts ≈ 1.0196729 atol = 1e-5
    @test err ≈ 0.00180075 atol = 1e-5
end