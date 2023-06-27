using Random

@testset "Fit VAR model" begin

    p_init = [0.,0.9]
    g(p) = p_init*0.1 + 0.9*p + 0.002*randn(2)
    x_tr = [0.0  4.87634e-11  -0.002317  -0.0043605  -0.00450047  -0.00369875  -0.0013696  0.0003507  0.000639107  0.00263507  0.00162606]
    p_tr = [0.0  -0.002317  -0.00227054  -0.000575291  0.000360702  0.00195049  0.00158475  0.000323396  0.00206192  -0.000745936  0.00197454
            0.9   0.901549   0.902012     0.900167     0.902006     0.897626    0.901026    0.900232     0.896804     0.900162     0.90211]
    h(x,p) = x[1] .* p
    model = fitVARmodel(x_tr,p_traj=p_tr,h=h,p=1)
    model2 = fitVARmodel(x_tr,p_traj=p_tr)

    h2(x) = x.^2
    model3 = fitVARmodel(x_tr,h=h2)
    model4 = fitVARmodel(x_tr)

    Random.seed!(3)
    x = zeros(1,103)
    for i in 4:103
        x[1,i] = 0.5*x[1,i-1] + 0.01*x[1,i-2] + 0.3*x[1,i-3] + 0.002*randn()
    end

    @test model.h([1.],[1.,2.]) == h([1.],[1.,2.])
    @test model.p == 1
    @test model.d_x == 1
    @test model.d_p == 2
    @test model.d_h == 2
    @test model.d_y == 5
    @test model.B_hat ≈ [1.08838e-6  4.77665e-5  0.999997  -1.20926e-6  -0.000806239  0.999947] atol=1e-4
    @test model.Σ_hat_u ≈ reshape([9.012087104866685e-18],1,1) atol=1e-4
    @test model.Σ_tilde_u ≈ reshape([3.604834841946674e-18],1,1) atol=1e-4
    @test model.Σ_β_hat ≈ [  1.28117e-12   3.10134e-10  -1.26739e-12  -1.42389e-12  -3.50158e-10  -3.44751e-10
                                3.10134e-10   1.41405e-7   -3.20914e-10  -3.44746e-10  -1.49325e-7   -1.57129e-7
                            -1.26739e-12  -3.20914e-10   1.87494e-12   1.40867e-12   5.30883e-10   3.56748e-10
                            -1.42389e-12  -3.44746e-10   1.40867e-12   1.58252e-12   3.89277e-10   3.83227e-10
                            -3.50158e-10  -1.49325e-7    5.30883e-10   3.89277e-10   2.98483e-7    1.65932e-7
                            -3.44751e-10  -1.57129e-7    3.56748e-10   3.83227e-10   1.65932e-7    1.74601e-7] atol=1e-4
    @test model.Γ_hat ≈ [  1.0         -0.00126214    3.72491e-5   0.900158    -3.77317e-7   -0.00113679
                            -0.00126214   6.76683e-6   -3.77317e-7  -0.00113679   9.53717e-10   6.09254e-6
                            3.72491e-5  -3.77317e-7    2.22133e-6   3.18072e-5  -2.72648e-9   -3.37771e-7
                            0.900158    -0.00113679    3.18072e-5   0.810288    -3.37771e-7   -0.00102389
                            -3.77317e-7   9.53717e-10  -2.72648e-9  -3.37771e-7   9.89756e-12   8.50793e-10
                            -0.00113679   6.09254e-6   -3.37771e-7  -0.00102389   8.50793e-10   5.48545e-6] atol=1e-4
    @test model.Σ_x1_hat ≈ reshape([1.4419339367786697e-17],1,1) atol=1e-4

    @test_throws ErrorException fitVARmodel(x_tr,p_traj=p_tr[:,1:end-1],h=h,p=1)
    @test_throws ErrorException fitVARmodel(x_tr,p_traj=p_tr,h=h,p=0)

    @test model2.B_hat ≈ [0.00182334  0.899947  0.99925  -0.00202657] atol=1e-4

    @test model3.B_hat ≈ [-1.85643e-5  0.817222  -7.31837] atol=1e-4

    @test model4.B_hat ≈ [-4.38661e-5  0.836412] atol=1e-4

    @test VARorder(x, criterion="AIC",p_max=5) == 3
    @test VARorder(x, criterion="FPE",p_max=5) == 3
    @test VARorder(x, criterion="HQ",p_max=5) == 3
    @test VARorder(x, criterion="SC",p_max=5) == 1
end