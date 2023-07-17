using Random

@testset "EOFs" begin
    Random.seed!(1234)

    prog(v,p,r) = [-0.01*v[1]^3 + p[1]*v[1]^2 + 1.01*v[1] + p[2] + 1e-2*r[1]]
    obs(v,p) = v
    v_init = [-1.]
    p_init = [0.,0.]
    DS = DynamicalSystem(prog,obs,v_init,p_init,random_vec_length=1)

    vals, vecs = EOF(prog, v_init, p_init, random_vec_length=1)
    vals2, vecs2 = EOF(DS)

    @test vals[1] ≈ 0.0027093 atol=1e-5
    @test length(vals) == 1 
    @test vecs == ones(1,1)

    @test vals2[1] ≈ 0.00195096 atol=1e-5
    @test length(vals2) == 1
    @test vecs2 == ones(1,1)

    EOFs = [0. 1. 0.
            1. 0. 0.
            0. 0. 1.]
    eigvals = [3., 1., 2.]
    obs = EOF_to_obs(EOFs, eigvals,2)

    @test obs([4., 5., 1.], [1.]) == [5., 1.]

    @test_throws ArgumentError EOF_to_obs(EOFs, eigvals, 4)
    @test_throws ArgumentError EOF_to_obs(EOFs, eigvals, 0)
    @test_throws ArgumentError EOF_to_obs(ones(2,2), eigvals, 2)
    @test_throws ArgumentError EOF_to_obs(ones(2,3), eigvals, 2)
end