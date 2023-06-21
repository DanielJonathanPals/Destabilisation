# Here is a collection of functions that can be used to test wether certain objects have the correct format

module FormatTests

export check_traj
export check_DynamicalSystem
export check_h
export check_xph


"""
    check_traj(traj)

Checks if `traj` has the format of a trajectory i.e. if `traj` is a two dimentional matrix with all entries
    given by numbers. Further it is always assumed that each column of trajectory represents a datapoint
    of the time series, i.e. each row denotes the time evolution of a single variable.
"""
function check_traj(traj)
    if !(typeof(traj) <: Matrix)
        error("You are currently treating an object as an trajectory which is not a matrix")
    end
    if any(x -> !(typeof(x) <: Number), traj)
        error("You are currently treating an object as an trajectory which does not only contain numbers")
    end
end

"""
    check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)

Compatibility test for initialisation of an object of `DynamicalSystem`.
"""
function check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)
    # Check if length of x_init is compatible with x_length
    if length(v_init) != v_length
        error("`v_init` is incompatible with `v_length`, as `length(v_init)` = $(length(v_init)) and `v_length` = $v_length")
    end

    # Check if length of p_init is compatible with p_length
    if length(p_init) != p_length
        error("`p_init` is incompatible with `p_length`, as `length(p_init)` = $(length(p_init)) and `x_length` = $p_length")
    end

    # Test if the function 'progressor' has the correct form
    try 
        progressor(v_init,p_init)
    catch
        error("The function `progressor` is not of the correct form. It is supposed to take two vectors representing v and p as inputs and return v of type `Vector{Float64}` for the next timestep")
    end
    v_1 = progressor(v_init,p_init)
    if typeof(v_1) != Vector{Float64}
        error("The output of the function `progressor must be of the type `Vector{Float64}`")
    end
    if length(v_1) != v_length
        error("The output of the function in the argument `progressor` must have length given by `v_length` = $v_length")
    end

    # Test if the function 'observable' has the correct form
    try
        observable(v_init,p_init)
    catch
        error("The function `observable` is not of the correct form. It is supposed to take two vectors representing v and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    x = observable(v_init,p_init)
    if typeof(x) != Vector{Float64}
        error("The output of the function `observable` must be of type Vector{Float64}")
    end
    if length(x) != x_length
        error("The output of the function in the argument `observable` must have length given by `x_length` = $x_length")
    end

    # Check that p_bounds has the correct length
    if p_bounds !== nothing
        if length(p_bounds) != p_length
            error("`p_bounds` does not have the correct length. It is expected to have length given by `p_length` = $p_length")
        end
        for (i,p) in enumerate(p_init)
            if (p < p_bounds[i][1]) || (p > p_bounds[i][2])
                error("All initial parameter values given in `p_init` must respect the restrictions given by `p_bounds`")
            end
        end
    end
end


"""
    check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds)

Compatibility test for initialisation of an object of `DynamicalSystem`.
"""
function check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds)
    # Test if the function 'progressor' has the correct form
    try 
        progressor(v_init,p_init)
    catch
        error("The function `progressor` is not of the correct form. It is supposed to take two vectors representing v and p as inputs and return v of type `Vector{Float64}` for the next timestep")
    end
    v_1 = progressor(v_init,p_init)
    if typeof(v_1) != Vector{Float64}
        error("The output of the function `progressor must be of the type `Vector{Float64}`")
    end
    if length(v_1) != length(v_init)
        error("The output of the function in the argument `progressor` must have length given by `length(v_init)` = $(length(v_init))")
    end

    # Test if the function 'observable' has the correct form
    try
        observable(v_init,p_init)
    catch
        error("The function `observable` is not of the correct form. It is supposed to take two vectors representing v and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    x = observable(v_init,p_init)
    if typeof(x) != Vector{Float64}
        error("The output of the function `observable` must be of type Vector{Float64}")
    end

    # Check that p_bounds has the correct length
    if p_bounds !== nothing
        if length(p_bounds) != length(p_init)
            error("`p_bounds` does not have the correct length. It is expected to have length given by `length(p_init)` = $(length(p_init))")
        end
        for (i,p) in enumerate(p_init)
            if (p < p_bounds[i][1]) || (p > p_bounds[i][2])
                error("All initial parameter values given in `p_init` must respect the restrictions given by `p_bounds`")
            end
        end
    end
end


"""
    check_h(h::Function, x_traj::Union{Matrix{Float64},Nothing}, p_traj::Union{Matrix{Float64},Nothing})

Compatibility test for the function h used for the inclusion of non-linearities in the model. h is supposed
    to be given by a `Vector{Float64}` valued function which takes one or two vector arguments, depending on whether
    a timeseries of parameters is included or not.
"""
function check_h(h::Function, x_traj::Matrix{Float64}, p_traj::Union{Matrix{Float64},Nothing})
    if p_traj !== nothing
        try
            h(x_traj[:,1], p_traj[:,1])
        catch
            error("The function `h` is not of the correct form. It is supposed to take two vectors representing x and p at a given point in time as inputs and return an object of type `Vector{Float64}`")
        end
        x = h(x_traj[:,1], p_traj[:,1])
        if typeof(x) != Vector{Float64}
            error("The output of the function `h` must be of type Vector{Float64}")
        end
    else
        try
            h(x_traj[:,1])
        catch
            error("The function `h` is not of the correct form. It is supposed to take one vector representing x at a given point in time as input and return an object of type `Vector{Float64}`. Note that 'h' only takes one argument if no timeseries of parameters is included")
        end
        x = h(x_traj[:,1])
        if typeof(x) != Vector{Float64}
            error("The output of the function `h` must be of type Vector{Float64}")
        end
    end
end


"""
    check_xph(x_traj::Matrix,p_traj::Union{Matrix,Nothing}=nothing,h::Union{Function,Nothing}=nothing)

Compatibility test for simultaneously checking the trajectory `x_traj`, the parameter trajectory `p_traj` and 
    the function `h`. `x_traj`, `p_traj` and `h` are supposed to satify their individual compatibility tests 
    and in addition, the number of datapoints in `x_traj` and `p_traj` must be the same. 
    Besides checking the compatibility of `x_traj`, `p_traj` and `h`, the function also returns the dimensions
    of the trajectory `x_traj`, the parameter trajectory `p_traj` and the output of `h` at the first datapoint 
    as `d_x`, `d_p` and `d_h`, respectively.
"""
function check_xph(x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing,h::Union{Function,Nothing}=nothing)
    check_traj(x_traj)
    d_x = size(x_traj,1)
    if p_traj !== nothing
        check_traj(p_traj)
        d_p = size(p_traj,1)
        if size(x_traj,2) != size(p_traj,2)
            error("The number of datapoints in the trajectory `traj` and the parameter trajectory `p_traj` must be the same")
        end
    else
        d_p = 0
    end
    if h !== nothing
        check_h(h,x_traj,p_traj)
        if p_traj === nothing
            d_h = size(h(x_traj[:,1]),1)
        else
            d_h = size(h(x_traj[:,1],p_traj[:,1]),1)
        end
    else
        d_h = 0
    end
    return d_x,d_p,d_h
end


end