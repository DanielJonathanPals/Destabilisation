using Random

@testset "Test VAR model" begin

    Random.seed!(123)
    prog(v,p) = [p[1] + p[2]*v[1] + 1e-10*randn()]
    obs(v,p) = v
    v_init = [0.]
    p_init = [0.,0.9]

    DS = DynamicalSystem(prog,obs,v_init,p_init)

    g(p) = p_init*0.1 + 0.9*p + 0.002*randn(2)

    v, p, x = integrateTraj(DS,g,110,v_init,p_init)
    x_train = x[:,1:100]
    x_test = x[:,101:end]
    p_train = p[:,1:100]
    p_test = p[:,101:end]

    h(x,p) = x[1] .* p

    model = fitVARmodel(x_train,p_traj=p_train,h=h,p=1)
    
    h1(x,p) = reshape(Array(x ⊗ p),length(x ⊗ p))

    model2 = VARmodel(h1,4,2,3,6,11,ones(2,45),ones(2,2),ones(2,2),ones(90,90),ones(45,45),ones(2,2))

    @test_throws ErrorException testPredictions(model,ones(2,3),p_traj=ones(3,3))
    @test_throws ErrorException testPredictions(model2,ones(2,3),p_traj=ones(3,3))

    @test testPredictions(model,x_test,p_traj=p_test)[1] == true
    @test testPredictions(model,x_test,p_traj=p_test)[2] ≈ 0.08722288299560255 atol=1e-4

    @test_throws ErrorException LMtest(model,5,ones(2,10),p_traj=ones(3,10))
    @test_throws ErrorException LMtest(model,20,x_test,p_traj=p_test)

    @test LMtest(model,5,x_train,p_traj=p_train)[1] == true
    @test LMtest(model,5,x_train,p_traj=p_train)[2] ≈ 0.24964172685128727 atol=1e-4
end