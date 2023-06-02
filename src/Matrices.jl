# Here we define some matrices which are often used in order to compute the parameter estimates

using .FormatTests


"""
    Y(traj::Matrix; p::Int64=1, T::Int64=100)

This function turns a trajectory `traj` which is a matrix with K rows into a matrix Y as defined
    in Lütkepohls book in equation (3.2.1). If `T+p` excedes the length of the trajectory, i.e.
    `traj` has less than `T+p` columns, then Y is retruned with `T = size(traj)[2]-p`. If on the 
    other hand the trajectory is longer than `T+p` then Y given by slicing the last `T` columns 
    out of `traj`.
"""
function Y(traj::Matrix; p::Int64=1, T::Int64=100)
    println(traj)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    if size(traj)[2] < p+1
        error("You are currently trying to compute the Matix Y from a trajectory that is too short!")
    end
    if size(traj)[2] < T+p
        @warn "You are currently trying to compute the Matix Y from a trajectory that is too short so instead of T = $T collumns the Y matrix will only have $(size(traj)[2]-p) columns"
        return traj[:,p+1:end]
    else
        return traj[:,end-T+1:end]
    end
end