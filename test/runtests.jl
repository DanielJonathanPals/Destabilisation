using Destabilisation
using Test

@testset "Destabilisation" begin
    include("FormatTests.jl")
    include("Matrices.jl")
    include("Integrate.jl")
    include("FitVARmodel.jl")
    include("TestVARmodel.jl")
    include("VARmodel.jl")
    include("TimeScales.jl")
    include("ParameterSeries.jl")
end
