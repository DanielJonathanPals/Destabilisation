# Here we introduce the functions which can be used to fit a VAR model to a trajectory

using .FormatTests
using Kronecker

include("Matrices.jl")   
include("VARmodel.jl")


"""
    VARorder(traj::Matrix; p_max::Int64=10, criterion::String="AIC")

This function computes the optimal order of a VAR model for a given trajectory `x_traj` by
    minimizing the AIC, FPE, HQ or SC criterion. The optimal order is returned as an integer.
    The `p_max` argument denotes the maximal order which is considered. If the system is forced
    by varying parameter values the respective time series should be given in the `p_traj` argument
    which must have the same number of datapoints as the `x_traj` argument. The `h` argument allows
    for non-linearities to be included in the model. To this end `h` must be a vector valued function
    which takes two arguments x and p if `p_traj !== nothing` and one argument x if `p_traj === nothing`.
    So in total the optimal VAR order of a model of 
    the form `x_traj_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'` is computed.
"""
function VARorder(x_traj::Matrix;
                    p_traj::Union{Matrix,Nothing}=nothing,
                    h::Union{Function,Nothing}=nothing,
                    p_max::Int64=3, 
                    criterion::String="AIC")

    # check the input
    d_x,_,_ = check_xph(x_traj,p_traj=p_traj,h=h)
    if p_max < 1
        error("The maximal order `p_max` must be at least 1")
    end
    if criterion ∉ ["AIC", "FPE", "HQ", "SC"]
        error("The criterion `criterion` must be one of `AIC`, `FPE`, `HQ` or `SC`")
    end

    # compute the optimal order
    l = zeros(p_max)
    T = size(x_traj,2)
    for p in 1:p_max
        model = VARmodel(x_traj, p_traj=p_traj, h=h, p=p)
        d_y = model.d_y
        if criterion == "AIC"
            l[p] = log(det(model.Σ_tilde_u)) + 2*p*d_x*d_y/T
        end
        if criterion == "FPE"
            l[p] = det(model.Σ_tilde_u) * ((T+d_y*p+1)/(T-d_y*p-1))^d_x
        end
        if criterion == "HQ"
            l[p] = log(det(model.Σ_tilde_u)) + 2*p*d_x*d_y*log(log(T))/T
        end
        if criterion == "SC"
            l[p] = log(det(model.Σ_tilde_u)) + log(T)*p*d_x*d_y/T
        end
    end
    if argmin(l) == p_max
        @warn "The optimal order is the maximal order `p_max`"
    end
    return argmin(l)
end


"""
    VARmodel(x_traj::Matrix; p_traj::Union{Matrix,Nothing}=nothing, h::Union(Function,Nothing)=nothing, p::Int64=1)

Fits a VAR(p) model of the form `x_traj_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'` to the given input data. It returns
    a `VARmodel` object which contains the estimated parameters and the estimated covariance matrices.
"""
function VARmodel(x_traj::Matrix;
    p_traj::Union{Matrix,Nothing}=nothing,
    h::Union{Function,Nothing}=nothing,
    p::Int64=1)

    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj=p_traj,h=h)
    if p < 1
        error("The order of the VAR model must be at least 1")
    end
    
    x_traj = convert(Matrix{Float64}, x_traj) 
    (p_traj === nothing) || (p_traj = convert(Matrix{Float64}, p_traj))

    # Create y_traj from x_traj, p_traj and h
    y_traj = create_y_traj(x_traj,p_traj=p_traj,h=h)
    d_y = size(y_traj,1)

    # set up the matrices needed for estimating the VAR model
    T = size(x_traj,2) - p
    X_arr = X(x_traj,p=p,T=T)
    Y_arr = Y(y_traj,p=p,T=T)

    # estimate the VAR model
    B_hat = X_arr * Y_arr' * (Y_arr * Y_arr')^-1
    Σ_tilde_u = (X_arr - B_hat * Y_arr) * (X_arr - B_hat * Y_arr)' / T
    Σ_hat_u = Σ_tilde_u * T / (T - p * d_y - 1)
    Γ_hat = Y_arr * Y_arr' / T
    Σ_β_hat = Matrix(Γ_hat^(-1) ⊗ Σ_hat_u) / T
    Σ_x1_hat = (T + d_y * p + 1) / T * Σ_hat_u

    return VARmodel(p,d_x,d_p,d_h,d_y,B_hat,Σ_hat_u,Σ_tilde_u,Σ_β_hat,Γ_hat,Σ_x1_hat)
end


function testPredictions(model,x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing,h::Union{Function,Nothing}=nothing)
    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj=p_traj,h=h)
    if (d_x, d_p, d_h) != (model.d_x, model.d_p, model.d_h)
        error("The dimensions of the input data do not match the dimensions of the model")
    end
    if size(x_traj,2) <= model.p
        error("The trajectory must contain more than `p`=$(model.p) datapoints")
    end
    
    x_traj = convert(Matrix{Float64}, x_traj) 
    (p_traj === nothing) || (p_traj = convert(Matrix{Float64}, p_traj))

    # Create y_traj from x_traj, p_traj and h
    y_traj = create_y_traj(x_traj,p_traj=p_traj,h=h)
end