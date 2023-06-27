using Kronecker

@testset "Format Tests" begin
    
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

    @test check_h((x,p) -> reshape(Array(x ⊗  p), length(x ⊗  p)), ones(3,10), ones(2,10)) === nothing
    @test check_h((x) -> [x[1]], ones(3,10), nothing) === nothing

    @test_throws ErrorException check_h((x) -> x, ones(3,10), ones(2,10))
    @test_throws ErrorException check_h((x,p) -> x[1]+p[1], ones(3,10), ones(2,10))
    @test_throws ErrorException check_h((x) -> x, ones(3,10), ones(2,10))
    @test_throws ErrorException check_h((x,p) -> x.*p[1], ones(3,10), nothing)
    @test_throws ErrorException check_h((x) -> x[1], ones(3,10), nothing)

    @test check_xph(ones(3,10), ones(2,10), (x,p) -> reshape(Array(x ⊗  p), length(x ⊗  p))) === (3,2,6)
    @test check_xph(ones(3,10), ones(2,10), nothing) === (3,2,0)
    @test check_xph(ones(3,10), nothing, nothing) === (3,0,0)
    @test check_xph(ones(3,10), nothing, x-> ones(10)) === (3,0,10)

    @test_throws ErrorException check_xph(ones(3,10), ones(2,9), (x,p) -> ones(3))
end