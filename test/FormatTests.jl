@testset "FormatTests.jl" begin
    
    @test_throws ErrorException check_traj([1,2,3])
    @test_throws ErrorException check_traj(1.0)
    @test_throws ErrorException check_traj(["a" 2; 3 5])
    @test_throws ErrorException check_traj(ones(3,5,4))

    @test check_traj(ones(Int64, 2, 3)) === nothing
    @test check_traj(ones(Float64, 3, 1)) === nothing

    @test check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 1, 1) != 0
    @test check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], [(-1.,2.)], 2, 1, 1) != 0

    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 3, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 0, 1)
    @test_throws ErrorException check_DynamicalSystem((x) -> x.^2, (x,p)-> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x[1], (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> [1.], (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (p) -> [x[1]+p[1]], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> x[1]+p[1], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1],5.], [1.,3.], [0.], nothing, 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], [(-1.,2.),(-3.,5.)], 2, 1, 1)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], [(1.,2.)], 2, 1, 1)

    @test_throws ErrorException check_DynamicalSystem((x) -> x.^2, (x,p)-> [x[1]+p[1]], [1.,3.], [0.], nothing)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x[1], (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> [1.], (x,p) -> [x[1]+p[1]], [1.,3.], [0.], nothing)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (p) -> [x[1]+p[1]], [1.,3.], [0.], nothing)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> x[1]+p[1], [1.,3.], [0.], nothing)
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], [(-1.,2.),(-3.,5.)])
    @test_throws ErrorException check_DynamicalSystem((x,p) -> x.^2, (x,p) -> [x[1]+p[1]], [1.,3.], [0.], [(1.,2.)])
end