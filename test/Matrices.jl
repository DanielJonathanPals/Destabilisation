@testset "FormatTests.jl" begin
    
    @test Y(ones(3, 20), T=10) == ones(Float64, 3, 10)
    @test Y(ones(3, 20), p=3) == ones(Float64, 3, 17)
    @test_throws ErrorException Y(ones(3, 3), p=3)

end