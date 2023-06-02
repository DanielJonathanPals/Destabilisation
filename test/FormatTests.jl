@testset "FormatTests.jl" begin
    
    @test_throws ErrorException check_traj([1,2,3])
    @test_throws ErrorException check_traj(1.0)
    @test_throws ErrorException check_traj(["a" 2; 3 5])
    @test_throws ErrorException check_traj(ones(3,5,4))

    @test check_traj(ones(Int64, 2, 3)) === nothing
    @test check_traj(ones(Float64, 3, 1)) === nothing

end