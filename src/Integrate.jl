# In this file we define methods that allow us to integrate a given System dependent on the dynamics of the
# parameters

include("Coupling.jl")


function integrateTraj(DS::DynamicalSystem, g::Function, T::Int64, x_0::Vector{Float64}, p_0::Vector{Float64})
    x_arr = zeros(Float64, DS.x_length, T)
    p_arr = zeros(Float64, DS.p_length, T)
    x_arr[:,1] = x_0
    p_arr[:,1] = p_0
    for t in 2:T
        x_arr[:,t] = DS.progressor(x_arr[:,t-1],p_arr[:,t-1])
        p_arr[:,t] = g(p_arr[:,t-1])
    end
    return x_arr, p_arr
end