module Destabilisation

export slice_traj
export X
export Y_t
export Y
export Id
export J
export create_y_traj

export check_traj
export check_DynamicalSystem
export check_h
export check_xph

export DynamicalSystem

export integrateTraj

export VARmodel
export VARorder

include("FormatTests.jl")
include("FitVARmodel.jl")
include("Integrate.jl")

end
