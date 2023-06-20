using Destabilisation
using Test

@testset "Destabilisation.jl" begin
    include("FormatTests.jl")
    include("Matrices.jl")
    include("Integrate.jl")
end
