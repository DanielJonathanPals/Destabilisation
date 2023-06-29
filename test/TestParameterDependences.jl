using Random

@testset "Test Parameter Dependences" begin

    Random.seed!(142)
    prog(v,p) = [p[1] + p[2]*v[1] + 1e-4*randn()]
    obs(v,p) = v
    v_init = [0.]
    p_init = [0.,0.9]
    DS = DynamicalSystem(prog,obs,v_init,p_init)
    g(p) = p_init*0.1 + 0.9*p + 0.002*randn(2)
    v_tr, p_tr, x_tr = integrateTraj(DS,g,1001,v_init,p_init)
    h(x,p) = x[1] .* p
    model1 = fitVARmodel(v_tr,p_traj=p_tr,h=h,p=1)
    model2 = fitVARmodel(v_tr,p_traj=p_tr[1:1,:],h=h,p=2)
    model3 = fitVARmodel(v_tr,p=3)

    @test getParamLocations(model1, 1) == [2, 4]
    @test getParamLocations(model1, 2) == [3, 5]
    @test getParamLocations(model2, 1) == [2, 3]

    @test_throws ArgumentError getParamLocations(model3, 1)
    @test_throws ArgumentError getParamLocations(model1, 3)

    @test testParamCausality(model1,1,1001)[1] == true
    @test testParamCausality(model1,1,1001)[2] ≈ 0. atol = 1e-6
    @test testParamCausality(model1,2,1001)[1] == true
    @test testParamCausality(model1,2,1001)[2] ≈ 0. atol = 1e-6
end