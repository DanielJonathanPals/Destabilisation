# In this file we define methods that allow us to integrate a given System dependent on the dynamics of the
# parameters

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
    v_arr = zeros(Float64, DS.v_length, T)
    p_arr = zeros(Float64, DS.p_length, T)
    x_arr = zeros(Float64, DS.x_length, T)
    v_arr[:,1] = v_0
    p_arr[:,1] = p_0
    x_arr[:,1] = DS.observable(v_arr[:,1],p_arr[:,1])
    for t in 2:T
        v_arr[:,t] = DS.progressor(v_arr[:,t-1],p_arr[:,t-1])
        p_arr[:,t] = g(p_arr[:,t-1])
        x_arr[:,t] = DS.observable(v_arr[:,t],p_arr[:,t])
    end
    return v_arr, p_arr, x_arr
end