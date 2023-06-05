# Here is a collection of functions that can be used to test wether certain objects have the correct format

module FormatTests

export check_traj
export check_DynamicalSystem

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
    check_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds,x_length,p_length,obs_length)

Compatibility test for initialisation of an object of `DynamicalSystem`.
"""
function check_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds,x_length,p_length,obs_length)
    # Check if length of x_init is compatible with x_length
    if length(x_init) != x_length
        error("`x_init` is incompatible with `x_length`, as `length(x_init)` = $(length(x_init)) and `x_length` = $x_length")
    end

    # Check if length of p_init is compatible with p_length
    if length(p_init) != p_length
        error("`p_init` is incompatible with `p_length`, as `length(p_init)` = $(length(p_init)) and `x_length` = $p_length")
    end

    # Test if the function 'progressor' has the correct form
    try 
        progressor(x_init,p_init)
    catch
        error("The function `progressor` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return x of type `Vector{Float64}` for the next timestep")
    end
    x_1 = progressor(x_init,p_init)
    if typeof(x_1) != Vector{Float64}
        error("The output of the function `progressor must be of the type `Vector{Float64}`")
    end
    if length(x_1) != x_length
        error("The output of the function in the argument `progressor` must have length given by `x_length` = $x_length")
    end

    # Test if the function 'observable' has the correct form
    try
        observable(x_init,p_init)
    catch
        error("The function `observable` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    obs = observable(x_init,p_init)
    if typeof(obs) != Vector{Float64}
        error("The output of the function `observable` must be of type Vector{Float64}")
    end
    if length(obs) != obs_length
        error("The output of the function in the argument `observable` must have length given by `obs_length` = $obs_length")
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
    check_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds)

Compatibility test for initialisation of an object of `DynamicalSystem`.
"""
function check_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds)
    # Test if the function 'progressor' has the correct form
    try 
        progressor(x_init,p_init)
    catch
        error("The function `progressor` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return x of type `Vector{Float64}` for the next timestep")
    end
    x_1 = progressor(x_init,p_init)
    if typeof(x_1) != Vector{Float64}
        error("The output of the function `progressor must be of the type `Vector{Float64}`")
    end
    if length(x_1) != length(x_init)
        error("The output of the function in the argument `progressor` must have length given by `length(x_init)` = $(length(x_init))")
    end

    # Test if the function 'observable' has the correct form
    try
        observable(x_init,p_init)
    catch
        error("The function `observable` is not of the correct form. It is supposed to take two vectors representing x and p as inputs and return the resulting observables for the current timestep as `Vector{Float64}`")
    end
    obs = observable(x_init,p_init)
    if typeof(obs) != Vector{Float64}
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

end