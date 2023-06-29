using Random

@testset "Parameter Series" begin

    Random.seed!(8)
    p_init = [-1.,4.,2.]
    noise = [0.1,0.2,0.3]

    @test_throws ArgumentError parameterSeriesGenerator(0.9, p_init, noise)
    @test_throws ArgumentError parameterSeriesGenerator(3., p_init[1:2], noise)
    @test parameterSeriesGenerator(3., p_init, noise) ≈ [-1.0  -1.04055  -0.996342  -0.879676  -0.976499  -0.803002  -0.831435  -0.887087  -0.978577  -0.862375
                                                        4.0   3.99408   4.17316    3.72233    4.20178    4.14912    3.92293    3.86299    4.12487    3.40655
                                                        2.0   1.48447   1.89611    2.20974    2.26591    1.99513    2.44639    2.28505    2.14657    1.73267] atol = 1e-4
end