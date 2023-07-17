using Random

@testset "Integrate" begin

    prog(v,p) = v + 0.01*(p[1].*(v-[p[2],p[3]]))
    prog2(v,p,r) = v + 0.01*(p[1].*(v-[p[2],p[3]])) + 0.1 .* r
    obs(v,p) = [v[1] + v[2]]
    v_init = [0.,0.]
    p_init = [1.,4.,2.]
    DS = DynamicalSystem(prog,obs,v_init,p_init)
    DS2 = DynamicalSystem(prog2,obs,v_init,p_init,random_vec_length=2)
    g(p) = p + [0.,0.01,0.02]

    v_tr = [0.0 -0.04 -0.0805; 0.0 -0.02 -0.0404]
    p_tr = [1.0 1.0 1.0; 4.0 4.01 4.02; 2.0 2.02 2.04]
    x_tr = [0.0 -0.06 -0.1209]

    Random.seed!(321)
    tr = integrateTraj(DS2,g,3)
    Random.seed!(321)
    tr2 = integrateTraj(DS2,g,3,include_reference_traj=true)
    Random.seed!(321)
    tr3 = integrateTraj(DS2,p_tr)
    Random.seed!(321)
    tr4 = integrateTraj(DS2,p_tr,include_reference_traj=true)

    @test integrateTraj(DS,g,3)[1] ≈ v_tr atol = 1e-4	
    @test integrateTraj(DS,g,3)[2] ≈ p_tr atol = 1e-4	
    @test integrateTraj(DS,g,3)[3] ≈ x_tr atol = 1e-4	
    @test integrateTraj(DS,g,3,include_reference_traj=true)[4] ≈ [0.0 -0.04 -0.0804; 0.0 -0.02 -0.0402] atol = 1e-4
    @test integrateTraj(DS,g,3,include_reference_traj=true)[5] ≈ [0.0 -0.06 -0.1206] atol = 1e-4

    @test tr[1] ≈ [0.0 0.065364728 -0.1553659; 0.0 -0.0365064 -0.032768] atol = 1e-4
    @test tr[2] ≈ [1.0 1.0 1.0; 4.0 4.01 4.02; 2.0 2.02 2.04] atol = 1e-4
    @test tr[3] ≈ [0.0 0.02885827281 -0.1881346] atol = 1e-4
    @test tr[4] ≈ [0.0 1.053647 -1.812842; 0.0 -0.16506 0.243027] atol = 1e-4
    @test tr2[4] ≈ [0.0 0.0653647 -0.155265; 0.0 -0.03650 -0.03256] atol = 1e-4
    @test tr2[5] ≈ [0.0 0.028858 -0.187834] atol = 1e-4

    @test_throws ArgumentError integrateTraj(DS,g,0)
    @test_throws ArgumentError integrateTraj(DS,(x,p) -> x,3) 
    @test_throws ArgumentError integrateTraj(DS,(x) -> x[1],3)
    @test_throws ArgumentError integrateTraj(DS,(p) -> [1.,2.],3)

    @test integrateTraj(DS,p_tr)[1] ≈ v_tr atol = 1e-4	
    @test integrateTraj(DS,p_tr)[2] ≈ p_tr atol = 1e-4	
    @test integrateTraj(DS,p_tr)[3] ≈ x_tr atol = 1e-4	
    @test integrateTraj(DS,p_tr,include_reference_traj=true)[4] ≈ [0.0 -0.04 -0.0804; 0.0 -0.02 -0.0402] atol = 1e-4
    @test integrateTraj(DS,p_tr,include_reference_traj=true)[5] ≈ [0.0 -0.06 -0.1206] atol = 1e-4

    @test tr3[1] ≈ [0.0 0.065364728 -0.1553659; 0.0 -0.0365064 -0.032768] atol = 1e-4
    @test tr3[2] ≈ [1.0 1.0 1.0; 4.0 4.01 4.02; 2.0 2.02 2.04] atol = 1e-4
    @test tr3[3] ≈ [0.0 0.02885827281 -0.1881346] atol = 1e-4
    @test tr3[4] ≈ [0.0 1.053647 -1.812842; 0.0 -0.16506 0.243027] atol = 1e-4
    @test tr4[4] ≈ [0.0 0.0653647 -0.155265; 0.0 -0.03650 -0.03256] atol = 1e-4
    @test tr4[5] ≈ [0.0 0.028858 -0.187834] atol = 1e-4

    @test_throws ArgumentError integrateTraj(DS,p_tr[1:2,:])

    @test integrateTraj(DS,3)[1] ≈ [0.0 -0.04 -0.0804
                                    0.0 -0.02 -0.0402] atol = 1e-4	
    @test integrateTraj(DS,3)[2] ≈ [1.0 1.0 1.0
                                    4.0 4.0 4.0
                                    2.0 2.0 2.0] atol = 1e-4	
    @test integrateTraj(DS,3)[3] ≈ [0.0 -0.06 -0.1206] atol = 1e-4	

end