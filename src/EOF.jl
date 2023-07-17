# This file introduces functions that allow us to compute the EOFs of a given system

using .VARmodel_module
using LinearAlgebra


"""
    EOF(progressor::Function, 
        v_init::Vector{Float64}, 
        p_init::Vector{Float64};
        random_vec_length::Union{Int64,Nothing}=nothing,
        T_max = 10000)

This function computes the EOFs of the dynamical system whos evolution is defined by 
    the progressor function `progressor` with 
    given initial state variables `v_init` and parameter variables `p_init`. The EOFs are computed by first 
    determining the slowest timescale of the system and then integrating the system for a time that is 
    demermined according to the timescale. The EOFs are then computed by performing a PCA on the 
    resulting trajectory. The function returns the EOFs and the corresponding eigenvalues.
"""
function EOF(progressor::Function, 
            v_init::Vector{Float64}, 
            p_init::Vector{Float64};
            random_vec_length::Union{Int64,Nothing}=nothing,
            T_max = 10000)

    DS = DynamicalSystem(progressor, (v,p) -> v, v_init, p_init, random_vec_length = random_vec_length)
    return EOF(DS, T_max = T_max)
end


"""
    EOF(DS::DynamicalSystem;
        T_max::Int64=10000)

This function computes the EOFs of the dynamical system `DS`. The EOFs are computed by first
    determining the slowest timescale of the system and then integrating the system for a time that is 
    demermined according to the timescale. The EOFs are then computed by performing a PCA on the 
    resulting trajectory. The function returns the EOFs and the corresponding eigenvalues.
"""
function EOF(DS::DynamicalSystem;
                T_max::Int64=10000)
    
    # Determine the timescale of the system
    ts, err = timeScale(DS, T_max = T_max)
    if ts - 2*err <= 1
        @warn "The timescale could not be determined with the desired accuracy. The EOFs will now be determined with a sample of length T_max = $T_max. This could potentially lead to inaccurate results."
    end
    # By subtracting two times the error on the timescale we ensure that the timescale is not overestimated
    ts -= 2*err

    # Integrate the system for a time scale dependent time (The integration time T is choses s.t. a system with
    # ts = 1/0.9 is integrated for a length of 200 which seems to ba a reasonable choise and this result is then
    # scaled by the actual ts)
    T = Int64(ceil(- 200 * log(0.9)/log(ts)))
    x_arr = integrateTraj(DS, T)[3]

    # Determine the EOFs
    mean = sum(x_arr, dims = 2)/size(x_arr, 2)
    mean_adjusted = x_arr .- mean
    cov = mean_adjusted * mean_adjusted' / size(x_arr, 2)
    Λ, U = eigen(cov)

    return Λ, U
end


"""
    EOF_to_obs(EOFs::Matrix{Float},
                eigenvals::Vector{Float},
                numb_of_obs::Int64)

This function returns a suggested function that computes reasonable observables from the EOFs of the system, by
    returning the `numb_of_obs` most relevant principal components.
"""
function EOF_to_obs(EOFs::Matrix{Float64},
                    eigenvals::Vector{Float64},
                    numb_of_obs::Int64)

    if numb_of_obs > length(eigenvals)
        throw(ArgumentError("The number of observables must be smaller than the number of EOFs."))
    end
    if numb_of_obs < 1
        throw(ArgumentError("The number of observables must be larger than zero."))
    end
    if size(EOFs, 2) != length(eigenvals)
        throw(ArgumentError("The number of EOFs must be equal to the number of eigenvalues."))
    end
    if size(EOFs, 1) != size(EOFs, 2)
        throw(ArgumentError("The EOFs must be a square matrix."))
    end

    order = sortperm(eigenvals, rev = true)
    obs(v::Vector{Float64},p::Vector{Float64}) = EOFs[:, order[1:numb_of_obs]]' * v
    return obs
end
