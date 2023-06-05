module Destabilisation

export slice_traj
export Y
export Z_t
export Z
export Id
export J

export check_traj
export check_DynamicalSystem

export DynamicalSystem

export integrateTraj

include("FormatTests.jl")
include("Matrices.jl")
include("Integrate.jl")

end
