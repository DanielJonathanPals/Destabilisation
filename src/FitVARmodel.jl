# Here we introduce the functions which can be used to fit a VAR model to a trajectory

using .FormatTests

include("Matrices.jl")   


"""
    VARorder(traj::Matrix; p_max::Int64=10, criterion::String="AIC")

This function computes the optimal order of a VAR model for a given trajectory `x_traj` by
    minimizing the AIC, FPE, HQ or SC criterion. The optimal order is returned as an integer.
    The `p_max` argument denotes the maximal order which is considered. If the system is forced
    by varying parameter values the respective time series should be given in the `p_traj` argument
    which must have the same number of datapoints as the `x_traj` argument. The `h` argument allows
    for non-linearities to be included in the model. To this end `h` must be a vector valued function
    which takes two arguments x and p if `p_traj !== nothing` and one argument x if `p_traj === nothing`.
    So in total the optimal VAR order of a model of 
    the form `x_traj_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'` is computed.
"""
function VARorder(x_traj::Matrix;
                    p_traj::Union{Matrix,Nothing}=nothing,
                    h::Union{Function,Nothing}=nothing,
                    p_max::Int64=10, 
                    criterion::String="AIC")
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    (d_x, T) = size(traj)
    if p_traj !== nothing
        check_traj(p_traj)
        p_traj = convert(Matrix{Float64}, p_traj)
        (d_p, T_p) = size(p_traj)
        if T != T_p
            error("The number of datapoints in the trajectory `traj` and the parameter trajectory `p_traj` must be the same")
        end
    end
end


"""
    VARmodel(x_traj::Matrix; p_traj::Union{Matrix,Nothing}=nothing, h::Union(Function,Nothing)=nothing, p::Int64=1)

Fits a VAR(p) model of the form `x_traj_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'` to the given input data.
"""
function VARmodel(x_traj::Matrix;
    p_traj::Union{Matrix,Nothing}=nothing,
    h::Union{Function,Nothing}=nothing,
    p::Int64=1)

    # check input data
    check_traj(x_traj)
    x_traj = convert(Matrix{Float64}, x_traj)
    d_x = size(x_traj)[1]
    if p_traj !== nothing
        check_traj(p_traj)
        p_traj = convert(Matrix{Float64}, p_traj)
        d_p = size(p_traj)[1]
        if size(x_traj)[2] != size(p_traj)[2]
            error("The number of datapoints in the trajectory `traj` and the parameter trajectory `p_traj` must be the same")
        end
    end
    if p < 1
        error("The order of the VAR model must be at least 1")
    end
    check_h(h,x_traj,p_traj)

    # Create y_traj from x_traj, p_traj and h
    y_traj = create_y_traj(x_traj,p_traj=p_traj,h=h)

    # set up the matrices needed for computing the VAR model
    T = size(x_traj)[2] - p
    X_arr = X(x_traj,p,T)
    Y_arr = Y(y_traj,p,T)
end
