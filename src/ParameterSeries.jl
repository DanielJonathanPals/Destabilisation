# Here a function is presented that allows to generate the parameter variations given the time scale of the system,
# the initial parameter values as well as the noise amplitudes.

using Random


"""
    parameterSeriesGenerator(ts::Float64, p_init::Vector{Float64}, noise::Vector{Float64}; n::Union{Int64,Nothing} = nothing)

Generate a parameter series given the time scale of the system, the initial parameter values as well as the noise amplitudes.
    The series has the form of a VAR(1) process with the given time scale and the given noise amplitudes varying around the
    given initial parameter values. The length of the parameter series is determined by the time scale if not specified.

# Arguments
- `ts::Float64`: The time scale of the system which for instance can be determined using the `timeScale` function.
- `p_init::Vector{Float64}`: The initial parameter values.
- `noise::Vector{Float64}`: The noise amplitudes.
- `n::Union{Int64,Nothing} = nothing`: The length of the parameter series. If not specified the length is determined by the time scale.
"""
function parameterSeriesGenerator(ts::Float64, p_init::Vector{Float64}, noise::Vector{Float64}; n::Union{Int64,Nothing} = nothing)

    # check that the input is valid
    if length(p_init) != length(noise)
        throw(ArgumentError("The length of p_init and noise must be the same."))
    end
    if ts <= 1
        throw(ArgumentError("The time scale must be larger than 1."))
    end
    if ts > 1.5
        @warn "The time scale is larger than 1.5. Maybe the discritisation is too coarse."
    end

    # The length of the parameter series is determined by the time scale if not specified otherwise. To this end we assume that
    # for a time scale of 1/0.9 a parameter series of length 100 is sufficient for the analysis.
    if n === nothing
        n = Int(round(-log(0.9)/log(ts) * 100))
    end

    # generate the parameter series
    p_series = zeros(length(p_init),n)
    p_series[:,1] = p_init
    for i in 2:n
        p_series[:,i] = (ts - 1) / ts * p_init + 1 / ts * p_series[:,i-1] + noise .* randn(length(p_init))
    end

    return p_series
end