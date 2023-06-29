# In this file we define methods that allow us to integrate a given System dependent on the dynamics of the
# parameters

using Kronecker

include("Coupling.jl")


"""
    integrateTraj(DS::DynamicalSystem, g::Function, T::Int64, v_0::Vector{Float64}, p_0::Vector{Float64})

This function iteratively computes the evolution of the state variables v, the parameters p and the
    observables.

# Arguments
- `DS::DynamicalSystem`: Element of the type `DynamicalSystem` which encodes the dynamics of the system under
    considerantion
- `g::Function`: Function determining the one-step evolution of the parameter values. I.e. g must be such that
    it takes one Argument `p::Vector{Float64}` containing the current parameter values and returns the new 
    parameter values again in the form of `Vector{Float64}`
- `T::Int64`: Number of time steps
- `v_0::Vector{Float64}`: Initial values of the state variables
- `p_0::Vector{Float64}`: initial parameter values

# Returns
- `v_arr': Vector of dimension (DS.length_v,T) containing the time evolution of the state space variable
- `p_arr': Vector of dimension (DS.length_p,T) containing the time evolution of the parameters
- `x_arr': Vector of dimension (DS.length_x,T) containing the time evolution of the observables
"""
function integrateTraj(DS::DynamicalSystem, g::Function, T::Int64, v_0::Vector{Float64}, p_0::Vector{Float64})

    # check that the input is valid
    if length(v_0) != DS.v_length
        throw(ArgumentError("The length of v_0 must be equal to the number of state variables."))
    end
    if length(p_0) != DS.p_length
        throw(ArgumentError("The length of p_0 must be equal to the number of parameters."))
    end
    if T <= 1
        throw(ArgumentError("The number of time steps must be larger than 1."))
    end
    try 
        g(p_0)
    catch
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64}."))
    end
    if typeof(g(p_0)) != Vector{Float64}
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64}."))
    end
    if length(g(p_0)) != DS.p_length
        throw(ArgumentError("The function g must take one argument of type Vector{Float64} and return a Vector{Float64} of the same length."))
    end

    p_arr = Array(ones(Float64, 1, T) ⊗ reshape(p_0, DS.p_length, 1))
    for t in 2:T
        p_arr[:,t] = g(p_arr[:,t-1])
    end
    return integrateTraj(DS, v_0, p_arr)
end


"""
    integrateTraj(DS::DynamicalSystem, v_0::Vector{Float64}, p_traj::Matrix{Float64})

This function iteratively computes the evolution of the state variables v and the observables of a system forced by
    a given trajectory of the parameters p.
"""
function integrateTraj(DS::DynamicalSystem, v_0::Vector{Float64}, p_traj::Matrix{Float64})

    # check that the input is valid
    if size(p_traj,1) != DS.p_length
        throw(ArgumentError("The number of rows of p_traj must be equal to the number of parameters."))
    end

    T = size(p_traj,2)
    v_arr = zeros(Float64, DS.v_length, T)
    x_arr = zeros(Float64, DS.x_length, T)
    v_arr[:,1] = v_0
    x_arr[:,1] = DS.observable(v_arr[:,1],p_traj[:,1])
    for t in 2:T
        v_arr[:,t] = DS.progressor(v_arr[:,t-1],p_traj[:,t-1])
        x_arr[:,t] = DS.observable(v_arr[:,t],p_traj[:,t])
    end
    return v_arr, p_traj, x_arr
end


"""
    integrateTraj(DS::DynamicalSystem, T::Int64, v_0::Vector{Float64}, p_0::Vector{Float64})

Integrates the trajectory of a system for T time steps with constant parameters.
"""
function integrateTraj(DS::DynamicalSystem, T::Int64, v_0::Vector{Float64}, p_0::Vector{Float64})
    g(p) = p
    return integrateTraj(DS, g, T, v_0, p_0)
end