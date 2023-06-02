@testset "FormatTests.jl" begin

    @test slice_traj(ones(2, 20), T=10, p=1) == ones(Float64, 2, 11)
    @test_throws ErrorException slice_traj(ones(2, 20))
    
    @test Y(ones(3, 20), T=10, p=10) == ones(Float64, 3, 10)
    @test Y(ones(3, 20), p=3) == ones(Float64, 3, 17)
    @test_throws ErrorException Y(ones(3, 20), T=100, p=3)
    @test_throws ErrorException Y(ones(3, 20), T=10, p=3)
    @test_throws ErrorException Y(ones(3, 3), T=100, p=3)

    @test Z_t(reshape(collect(1:20),2,10), 5, p=2) == reshape([1.,13.,14.,11.,12.],5,1)
    @test_throws ErrorException Z_t(reshape(collect(1:12),2,6), 5, p=2)

    @test Z(reshape(collect(1:20),2,10), p=2) == Z(reshape(collect(1:20),2,10), p=2, T=8)
    @test size(Z(reshape(collect(1:20),2,10), p=2)) == (5,8)
    @test Z(reshape(collect(1:20),2,10), p=2)[:,1] == [1.,3.,4.,1.,2.]
    @test Z(reshape(collect(1:20),2,10), p=2)[:,end] == [1.,17.,18.,15.,16.]
    @test_throws ErrorException Z(reshape(collect(1:20),2,10), p=2, T=7)
    @test_throws ErrorException Z(reshape(collect(1:20),2,10), p=2, T=9)
end