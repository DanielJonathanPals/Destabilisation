@testset "Integrate" begin

    prog(v,p) = v + 0.01*(p[1].*(v-[p[2],p[3]]))
    obs(v,p) = [v[1] + v[2]]
    v_init = [0.,0.]
    p_init = [1.,4.,2.]
    DS = DynamicalSystem(prog,obs,v_init,p_init)
    g(p) = p + [0.,0.01,0.02]

    v_tr = [0.0 -0.04 -0.0805; 0.0 -0.02 -0.0404]
    p_tr = [1.0 1.0 1.0; 4.0 4.01 4.02; 2.0 2.02 2.04]
    x_tr = [0.0 -0.06 -0.1209]

    @test integrateTraj(DS,g,3,v_init,p_init)[1] ≈ v_tr
    @test integrateTraj(DS,g,3,v_init,p_init)[2] ≈ p_tr
    @test integrateTraj(DS,g,3,v_init,p_init)[3] ≈ x_tr

end