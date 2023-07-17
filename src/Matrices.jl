# Here we define some matrices which are often used in order to compute the parameter estimates

using .FormatTests
using LinearAlgebra
using SparseArrays
using Kronecker

# Simple auxilary function to check whether there exists a natrual number n s.t. x = n*y 
isDivisible(x,y) = (x == (length(collect(0:y:x))-1)*y)

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
    Y(x_traj::Matrix{Float64}, p_traj::Matrix{Float64}, p_init::Vector{Float64} , h::Function, p::Int64 = 1)

Returns the matrix Y that is needed when performing a VAR regression where a reference trajectory is included.
    Here Y is a matrix of size `(dp+dh*p) x (T-p)` where `T = size(x_traj, 2)`, `dp = size(p_traj, 1)` and
    `dh = size(h(x_traj[:,1], p_init), 1)`. The function returns a matrix where the t-th column is given by
    `Yₜ = [pₜ₊ₚ₋₁ - p_init; h(xₜ₊ₚ₋₁, pₜ₊ₚ₋₁ - p_init); ...; h(xₜ, pₜ - p_init)]` where `pₜ` is the t-th column of 
    `p_traj` and `xₜ` is the t-th column of `x_traj`. 
"""
function Y(x_traj::Matrix{Float64}, p_traj::Matrix{Float64}, p_init::Vector{Float64} , h::Function; p::Int64 = 1)
    dx,dp,dh = check_xph(x_traj,p_traj,h)
    T = size(x_traj)[2]
    if length(p_init) != dp
        throw(ArgumentError("The length of the initial parameter vector `p_init` must be equal to the number of parameters in the trajectory `p_traj`"))
    end
    p_init_traj = repeat(p_init,1,T)
    y_traj = create_y_traj(x_traj, p_traj = p_traj - p_init_traj, h = h)[dx+dp+1:end,:]
    Y_arr = ones(Float64, dh*p + dp, T - p)
    Y_arr[1:dp,:] = p_traj[:,p:end-1] - p_init_traj[:,p:end-1]
    Y_arr[dp+1:end,:] = Y(y_traj, p = p)[2:end,:]

    return Y_arr
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
    F_i(i::Int64, T::Int64)

Returns the matrix Fᵢ as defined in equation (4.4.2)
"""
function F_i(i::Int64, T::Int64)
    arr = spzeros(T,T)
    arr[i+1:end,1:T-i] = Id(T-i)
    return arr
end


"""
    F(h::Int64, T::Int64)

Returns the matrix Fᵢ as defined in equation (4.4.2)
"""
function F(h::Int64,T::Int64)
    if h > T-1
        error("The argument `h` must be smaller than T-1")
    end
    arr = spzeros(T,T*h)
    for i in 0:h-1
        arr[:,i*T+1:(i+1)*T] = F_i(i+1,T)
    end
    return arr
end


"""
    C(param_locations::Vector{Int64}, d_x::Int64, d_y::Int64, p::Int64)

C is the matrix as defined in section 3.6.1 in Lütkepohls book. It can be used to test if the observables in the 
    fitted VAR model show a dependence on a given parameter p. To this end the positions of the p dependent entries
    of y = (x',p',h(x,p)')' must be known and encoded in the `param_locations` argument. The argument `d_x` denotes
    the dimension of the x vector, `d_y` the dimension of the y vector and `p` the number of presamples. The function
    then returns the respective matrix C of dimentions `length(param_locations) × d_x * (p * d_y + 1)` as a sparse matrix.
"""
function C(param_locations::Vector{Int64}, d_x::Int64, d_y::Int64, p::Int64)

    # check that the param_locations are valid
    if length(param_locations) > d_y
        error("The number of param_locations must be smaller than d_y")
    end

    n = length(param_locations) * d_x * p
    l = length(param_locations)

    arr = spzeros(n, d_x * (p * d_y + 1))
    for (idx,loc) in enumerate(param_locations), x_pos in 1:d_x, p_pos in 1:p
        arr[x_pos + (idx-1)*d_x + (p_pos-1)*d_x*l, d_x + d_x*(loc-1) + x_pos + (p_pos-1)*d_y] = 1.
    end
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


"""
    cut_traj(traj::Matrix{Float64}, keep_const::Int64, cut_time_steps::Int64)

This function devides a trajectory `traj` into pieces of length `keep_const` and removes `cut_time_steps` time 
    steps from the beginning of each piece. The resulting trajectory is returned as a matrix.
"""
function cut_traj(traj::Matrix{Float64}, keep_const::Int64, cut_time_steps::Int64)
    check_traj(traj)
    traj = convert(Matrix{Float64}, traj)
    n = size(traj)[2]
    if !(isDivisible(n, keep_const))
        throw(ArgumentError("The number of time steps must be divisible by keep_const"))
    end
    if cut_time_steps >= keep_const
        throw(ArgumentError("The number of time steps to cut must be smaller than keep_const"))
    end

    # this is the result of n/keep_const
    quot = length(collect(0:keep_const:n))-1
    new_traj = zeros(Float64, size(traj)[1], n - quot*cut_time_steps)
    for i in 1:quot
        new_traj[:,(i-1)*(keep_const-cut_time_steps)+1:i*(keep_const-cut_time_steps)] = traj[:,(i-1)*keep_const+cut_time_steps+1:i*keep_const]
    end
    return new_traj
end