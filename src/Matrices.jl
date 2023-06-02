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
    Y(traj::Matrix; p::Int64=1, T::Union{Int64,Nothing}=nothing)

This function turns a trajectory `traj`, which is a matrix with K rows into a matrix, Y as defined
    in Lütkepohls book in equation (3.2.1) by slicing the last `T` colomns from the `traj` argument. 
    If `T===nothing` then `T` is set to the number of
    columns of `traj` minus `p` i.e. the `traj` argument is interpreted in such a way that it
    consists of all the presamples plus all the actual samples. If on the other hand `T` is given
    and `T+p` is not equal to the length of the trajectory, i.e. `traj` has less than `T+p` columns,
    an error is throuwn.

    As always `p` denotes the number of presamples and `T` encodes the number of actual samples. 
"""
function Y(traj::Matrix; p::Int64=1, T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    (T === nothing) && (T = size(traj)[2]-p)
    if size(traj)[2] != T+p
        error("Your trajcetory `traj` should contain T+p=$(T+p) columns but it instead contains $(size(traj)[2]) columns")
    end
    return traj[:,p+1:end]
end


"""
    Z_t(traj::Matrix, t; p::Int64=1)

Returns the vector Zₜ as defined e.g. in equation (3.2.1) in Lütkepohls book. Here it is assumed
    that the trajectory `traj` includes the presample datapoints i.e. yₜ from the book corresponds
    to the column of `traj` with index `t+p`.

    As always `p` denotes the number of presamples.
"""
function Z_t(traj::Matrix, t; p::Int64=1)
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
    Z(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)

Returns Z as defined in Lütkepohls book in equation (3.2.1). If the number of columns of
the trajectory `traj` does not equal T+p then an error is throuwn. 

As always `p` denotes the number of presamples and `T` encodes the number of actual samples. 
"""
function Z(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)
    check_traj(traj)
    (T === nothing) && (T = size(traj)[2]-p)
    if size(traj)[2] != T+p
        error("Your trajcetory `traj` should contain T+p=$(T+p) columns but it instead contains $(size(traj)[2]) columns")
    end
    K = size(traj)[1]
    arr = ones(Float64, K*p+1, T)
    for (idx,t) in enumerate(0:T-1)
        arr[:,idx] = Z_t(traj, t, p=p)
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