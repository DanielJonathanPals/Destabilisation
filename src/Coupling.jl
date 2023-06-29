# In this file functions are introduced, which makes the system we want to analyse acceceble for the analysis.

using .FormatTests


"""
    DynamicalSystem

An instance of the struct `DynamicalSystem` contains all the relevant information of the stochastic dynamical 
    system which we want to destabilize.

# Fields
- `progressor::Function`: The Argument `progressor` should contains a function that takes two arguments, each of type
    `Vector{Float64}`, representing the current state variables v and the current parameter values `p` and returns
    the state variables for the next timestep as a `Vector{Float64}`. I.e. this function describes the progression
    of the dynamical system for one time step
- `observable::Function`: This function should take the current state variables v and the current parameter values p,
    both as `Vector{Float64}`, as input and return an output of type `Vector{Float64}` containing all the observables
    of the system which are relevant for further analysis.
- `v_init::Vector{Float64}`: Initial values for the state space variables.
- `p_init::Vector{Float64}`: Initial parameter values.
- `p_bounds::Union{Vector{Tuple},Nothing}`: In case the parameter values are bounded to certain intervalls, the
    bounds of each of these intervalls can be specified in this argument. I.e. the i-th tuple in `p_bounds` 
    represent the interval bounds for the i-th parameter. If `p_bounds = nothing` then it is assumes that there
    are no restictions to the parameter values
- `v_length::Int64`: Number of state space variables
- `p_length::Int64`: Number of parameters
- `x_length::Int64`: Number of observables
"""
mutable struct DynamicalSystem
    progressor::Function
    observable::Function
    v_init::Vector{Float64}
    p_init::Vector{Float64}
    p_bounds::Union{Vector{Tuple},Nothing}
    v_length::Int64
    p_length::Int64
    x_length::Int64
    function DynamicalSystem(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)
        check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)
        new(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)
    end
end


"""
    DynamicalSystem(progressor,observable,v_init,p_init;p_bounds=nothing)

Alternative initialisation of an element of type `DynamicalSystem` where `v_length`, `p_length` and 
    `x_length` are automatically computed.
"""
function DynamicalSystem(progressor,observable,v_init,p_init;p_bounds=nothing)
    check_DynamicalSystem(progressor,observable,v_init,p_init,p_bounds)
    v_length = length(v_init)
    p_length = length(p_init)
    x_length = length(observable(v_init,p_init))
    return DynamicalSystem(progressor,observable,v_init,p_init,p_bounds,v_length,p_length,x_length)
end