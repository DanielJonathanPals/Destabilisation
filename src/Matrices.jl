# Here we define some matrices which are often used in order to compute the parameter estimates

using .FormatTests
using LinearAlgebra
using SparseArrays


"""
    slice_traj(data::Matrix; p::Int64=1, T::Int64=100)

Slices out the last `T+p` columns from the `data` argument. Here `p` denotes the number of presamples
and `T` encodes the number of actual samples. 
"""
function slice_traj(data::Matrix; p::Int64=1, T::Int64=100)
    check_traj(data)
    data = convert(Matrix{Float64}, data)
    if size(data)[2] < T+p
        error("The time series given in the input data is to small (length of dataseries = $(size(data)[2])) in order to slice a trajectory with T+p=$(T+p) datapoints out of it")
    else
        return data[:,end-T-p+1:end]
    end
end


"""
    X(traj::Matrix; p::Int64=1, T::Union{Int64,Nothing}=nothing)

This function turns a trajectory `traj`, which is a matrix with K rows, into a matrix X as defined
    in Lütkepohls book in equation (3.2.1) (where it is called Y) by slicing the last `T` colomns from the `traj` argument. 
    If `T===nothing` then `T` is set to the number of
    columns of `traj` minus `p` i.e. the `traj` argument is interpreted in such a way that it
    consists of all the presamples plus all the actual samples. If on the other hand `T` is given
    and `T+p` is not equal to the length of the trajectory, i.e. `traj` has less than `T+p` columns,
    an error is throuwn.
    As always `p` denotes the number of presamples and `T` encodes the number of actual samples. 
"""
function X(traj::Matrix; p::Int64=1, T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    (T === nothing) && (T = size(traj)[2]-p)
    if size(traj)[2] != T+p
        error("Your trajcetory `traj` should contain T+p=$(T+p) columns but it instead contains $(size(traj)[2]) columns")
    end
    return traj[:,p+1:end]
end


"""
    Y_t(traj::Matrix, t; p::Int64=1)

Returns the vector Yₜ as defined e.g. in equation (3.2.1) in Lütkepohls book (where it is called Zₜ). Here it is assumed
    that the trajectory `traj` includes the presample datapoints i.e. yₜ from the book corresponds
    to the column of `traj` with index `t+p`.
    As always `p` denotes the number of presamples.
"""
function Y_t(traj::Matrix, t; p::Int64=1)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    if size(traj)[2] < t+p
        error("The time series encoded in the `traj` argument is not long enought (only includes $(size(traj)[2]) datapoints) as it must include at least T+p=$(t+p) datapoints")
    end
    K = size(traj)[1]
    arr = ones(Float64, K*p+1, 1)
    for i in 0:(p-1)
        arr[2+K*i:1+K*(i+1),1] = traj[:,t+p-i]
    end
    return arr
end


"""
    Y(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)

Returns Y as defined in Lütkepohls book in equation (3.2.1) (where it is called Z). If the number of columns of
    the trajectory `traj` does not equal T+p then an error is throuwn. 
    As always `p` denotes the number of presamples and `T` encodes the number of actual samples. 
"""
function Y(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    (T === nothing) && (T = size(traj)[2]-p)
    if size(traj)[2] != T+p
        error("Your trajcetory `traj` should contain T+p=$(T+p) columns but it instead contains $(size(traj)[2]) columns")
    end
    K = size(traj)[1]
    arr = ones(Float64, K*p+1, T)
    for (idx,t) in enumerate(0:T-1)
        arr[:,idx] = Y_t(traj, t, p=p)
    end
    return arr
end


"""
    Id(d)

Returns a `d` by `d` sparse identity matrix
"""
Id(d) = sparse(Matrix{Float64}(I,d,d))


"""
    J(K::Int64, p::Int64)

Returns the matrix J as defined in equation (2.1.11)
"""
function J(K::Int64, p::Int64)
    arr = spzeros(K,K*p)
    arr[1:K,1:K] = Id(K)
    return arr
end


"""
    create_y_traj(x_traj::Matrix{Float64};p_traj::Union{Matrix{Float64},Nothing}=nothing,h::Union{Function,Nothing}=nothing)

This function creates a trajectory `y_traj` from `x_traj`, `p_traj` and `h` s.t. the i-th column of `y_traj` is given by
    `(x_traj[:,i]', p_traj[:,i]', h(x_traj[:,i],p_traj[:,i])')'`.
"""
function create_y_traj(x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing,h::Union{Function,Nothing}=nothing)
    check_traj(x_traj)
    x_traj = convert(Matrix{Float64}, x_traj)
    d_x = size(x_traj)[1]
    if p_traj === nothing
        if h === nothing
            y_traj = x_traj
        else
            check_h(h,x_traj,p_traj)
            d_h = length(h(x_traj[:,1]))
            y_traj = zeros(Float64, d_h + d_x, size(x_traj)[2])
            y_traj[1:d_x,:] = x_traj
            for i in 1:size(x_traj)[2]
                y_traj[d_x+1:end,i] = h(x_traj[:,i])
            end
        end
    else
        check_traj(p_traj)
        p_traj = convert(Matrix{Float64}, p_traj)
        d_p = size(p_traj)[1]
        if h === nothing
            y_traj = zeros(Float64, d_p + d_x, size(x_traj)[2])
            y_traj[1:d_x,:] = x_traj
            y_traj[d_x+1:end,:] = p_traj
        else
            check_h(h,x_traj,p_traj)
            d_h = length(h(x_traj[:,1],p_traj[:,1]))
            y_traj = zeros(Float64, d_h + d_p + d_x, size(x_traj)[2])
            y_traj[1:d_x,:] = x_traj
            y_traj[d_x+1:d_x+d_p,:] = p_traj
            for i in 1:size(x_traj)[2]
                y_traj[d_x+d_p+1:end,i] = h(x_traj[:,i],p_traj[:,i])
            end
        end
    end
    return y_traj
end