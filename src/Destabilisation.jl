module Destabilisation

export slice_traj
export X
export Y_t
export Y
export Id
export J
export F_i
export F
export C
export create_y_traj

export check_traj
export check_DynamicalSystem
export check_h
export check_xph

export DynamicalSystem

export integrateTraj

export fitVARmodel
export VARorder
export testPredictions

export VARmodel
export oneStepPred

export testPredictions
export LMtest

export timeScale
export toPolynomial

export parameterSeriesGenerator

export getParamLocations
export testParamCausality

include("FormatTests.jl")
include("VARmodel.jl")
include("ParameterCalibration.jl")

end
