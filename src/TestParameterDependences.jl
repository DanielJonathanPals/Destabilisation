# This file contains methods allowing us to test wheather or not the observables in an estimated VAR model depend on 
# certain parameters or not.

using .VARmodel_module
using Random
using Kronecker
using Distributions

include("Matrices.jl")


"""
    get_param_locations(model::VARmodel, param_idx::Int64)

Returns the indices of the observables that depend on the parameter with index `param_idx` in the form of an array
    which can be used as the argument `param_locations` in the function `C`.
"""
function getParamLocations(model::VARmodel, param_idx::Int64)

    # Check that the input is valid.
    if model.d_p == 0
        throw(ArgumentError("The model does not have any parameters."))
    end
    if param_idx > model.d_p
        throw(ArgumentError("The parameter index is too large."))
    end

    if model.h === nothing
        return [model.d_x + param_idx]
    else
        h = model.h
        locs = zeros(Int64, model.d_y)
        locs[model.d_x + param_idx] = 1
        for i in 1:10
            x = randn(model.d_x)
            p1 = randn(model.d_p)
            p2 = zeros(model.d_p)
            p2 += p1
            p2[param_idx] += randn()
            arr = findall(x -> x != 0, round.(h(x, p1) - h(x, p2), digits = 6))
            locs[arr .+ model.d_x .+ model.d_p] .= 1
        end
        return findall(x -> x != 0, locs)
    end
end


"""
    testParamCausality(model::VARmodel, param_idx::Int64, T::Int64; siglvl::Float64 = 0.01)

Returns a tuple `(is_causal, p_value)` where `is_causal` is a boolean indicating whether or not the parameter with 
    index `param_idx` is causal for the observables in the VAR model `model` and `p_value` is the p-value of the 
    test. In order to compute the test statistic it is necessary to give the length `T` of the timeseries, that was
    used to estimate the VAR model. The test is performed at significance level `siglvl`.
"""
function testParamCausality(model::VARmodel, param_idx::Int64, T::Int64; siglvl::Float64 = 0.01)
    
    # compute the test statistic λ_F  
    param_loc = getParamLocations(model, param_idx)
    C_arr = C(param_loc, model.d_x, model.d_y, model.p)
    Γ = model.Γ_hat
    Σ_u = model.Σ_hat_u
    β = vec(model.B_hat)
    n = length(param_loc)
    λ_F = 1/n * (C_arr * β)' * (C_arr * ((Γ .* T)^(-1) ⊗ Σ_u) * C_arr')^(-1) * (C_arr * β)

    # now test the statistic against an F(n, T - model.d_y * model.p - 1) distribution
    if λ_F > cquantile(FDist(n, T - model.d_y * model.p - 1), siglvl)
        return true, ccdf(FDist(n, T - model.d_y * model.p - 1), λ_F)
    else
        return false, ccdf(FDist(n, T - model.d_y * model.p - 1), λ_F)
    end
end