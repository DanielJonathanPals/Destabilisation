using Destabilisation
using Test

@testset "Destabilisation.jl" begin
    include("FormatTests.jl")
    include("Matrices.jl")
    include("Integrate.jl")
    include("FitVARmodel.jl")
    include("TestVARmodel.jl")
    include("VARmodel.jl")
end
