using SparseArrays

@testset "Matracies" begin

    @test slice_traj(ones(2, 20), T=10, p=1) == ones(Float64, 2, 11)
    @test_throws ErrorException slice_traj(ones(2, 20))
    
    @test X(ones(3, 20), T=10, p=10) == ones(Float64, 3, 10)
    @test X(ones(3, 20), p=3) == ones(Float64, 3, 17)
    @test_throws ErrorException X(ones(3, 20), T=100, p=3)
    @test_throws ErrorException X(ones(3, 20), T=10, p=3)
    @test_throws ErrorException X(ones(3, 3), T=100, p=3)

    @test Y_t(reshape(collect(1:20),2,10), 5, p=2) == reshape([1.,13.,14.,11.,12.],5,1)
    @test_throws ErrorException Y_t(reshape(collect(1:12),2,6), 5, p=2)

    @test Y(reshape(collect(1:20),2,10), p=2) == Y(reshape(collect(1:20),2,10), p=2, T=8)
    @test size(Y(reshape(collect(1:20),2,10), p=2)) == (5,8)
    @test Y(reshape(collect(1:20),2,10), p=2)[:,1] == [1.,3.,4.,1.,2.]
    @test Y(reshape(collect(1:20),2,10), p=2)[:,end] == [1.,17.,18.,15.,16.]
    @test_throws ErrorException Y(reshape(collect(1:20),2,10), p=2, T=7)
    @test_throws ErrorException Y(reshape(collect(1:20),2,10), p=2, T=9)

    @test create_y_traj(reshape(collect(1:6),2,3), p_traj = ones(1,3), h = (x,p) -> reshape(Array(x ⊗  p), length(x ⊗  p))) == [1.0 3.0 5.0; 2.0 4.0 6.0; 1.0 1.0 1.0; 1.0 3.0 5.0; 2.0 4.0 6.0]
    @test create_y_traj(reshape(collect(1:6),2,3), p_traj = ones(1,3)) == [1.0 3.0 5.0; 2.0 4.0 6.0; 1.0 1.0 1.0]
    @test create_y_traj(reshape(collect(1:6),2,3), h = x -> x) == [1.0 3.0 5.0; 2.0 4.0 6.0; 1.0 3.0 5.0; 2.0 4.0 6.0]
    @test create_y_traj(reshape(collect(1:6),2,3)) == [1.0 3.0 5.0; 2.0 4.0 6.0]

    @test Array(F_i(3,5)) == [0.0  0.0  0.0  0.0  0.0
                                0.0  0.0  0.0  0.0  0.0
                                0.0  0.0  0.0  0.0  0.0
                                1.0  0.0  0.0  0.0  0.0
                                0.0  1.0  0.0  0.0  0.0]

    @test_throws ErrorException F(5,5)
    @test Array(F(3,4)) == [ 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                            1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                            0.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                            0.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0]

    # test C
    # System with dx = 2, dp = 3, dh = 0 and p = 1
    @test C([4], 2, 5, 1) == sparse([0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0
                                    0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0])
    # System with dx = 2, dp = 3, dh = 6 and p = 1
    @test C([3, 6, 9], 2, 11, 1) == sparse([0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                            0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
                                            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0])
    # System with dx = 1, dp = 3, dh = 0 and p = 2
    @test C([4], 1, 4, 2) == sparse([0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
                                    0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0])

    # System with dx = 1, dp = 3, dh = 3 and p = 2
    @test C([2, 5], 1, 7, 2) == sparse([0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                                        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
                                        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0])
end