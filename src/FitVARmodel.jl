# Here we introduce the functions which can be used to fit a VAR model to a trajectory

using .FormatTests
using Kronecker
using .VARmodel_module


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
    d_x,_,_ = check_xph(x_traj,p_traj,h)
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
        model = fitVARmodel(x_traj, p_traj=p_traj, h=h, p=p)
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
    fitVARmodel(x_traj::Matrix; p_traj::Union{Matrix,Nothing}=nothing, h::Union(Function,Nothing)=nothing, p::Int64=1)

Fits a VAR(p) model of the form `x_traj_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'` to the given input data. It returns
    a `VARmodel` object which contains the estimated parameters and the estimated covariance matrices.
"""
function fitVARmodel(x_traj::Matrix;
    p_traj::Union{Matrix,Nothing}=nothing,
    h::Union{Function,Nothing} = nothing,
    p::Int64=1)

    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj,h)
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

    B_hat,Σ_hat_u,Σ_tilde_u,Σ_β_hat,Γ_hat,Σ_x1_hat = auxFitVARmodel(h,p,d_y,T,X_arr,Y_arr)
    return VARmodel(h,p,d_x,d_p,d_h,d_y,B_hat,Σ_hat_u,Σ_tilde_u,Σ_β_hat,Γ_hat,Σ_x1_hat)
end


# Dokumentation and testing missing
function fitVARmodel(x_traj::Matrix, 
                        x_ref_traj::Matrix, 
                        p_traj::Matrix, 
                        p_init::Vector, 
                        fixed_param_model::VARmodel)
    h(x::Vector,p::Vector) = reshape(Array(x ⊗ p), length(x)*length(p))

    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj,h)
    d_y = d_x + d_p + d_h
    check_xph(x_ref_traj,p_traj,h)
    
    x_traj = convert(Matrix{Float64}, x_traj) 
    x_ref_traj = convert(Matrix{Float64}, x_ref_traj)
    p_traj = convert(Matrix{Float64}, p_traj)

    T = size(x_traj,2) - 1
    X_arr = X(x_traj-x_ref_traj,p=1,T=T)
    Y_arr = Y(x_traj, p_traj, p_init, h, p = 1)

    B_hat_diff,Σ_hat_u_diff,Σ_tilde_u_diff,Σ_β_hat_diff,Γ_hat_diff,Σ_x1_hat_diff = auxFitVARmodel(h,1,d_y,T,X_arr,Y_arr)

    B_hat = zeros(Float64, d_x, 1 + d_y)
    B_hat[:,1] = fixed_param_model.B_hat[:,1] - B_hat_diff[:,1:d_p] * p_init
    B_hat[:,2:1+d_x] = fixed_param_model.B_hat[:,2:1+d_x] - B_hat_diff[:,d_p+1:end] * sparse(Id(d_x) ⊗ p_init)
    B_hat[:,2+d_x:end] = B_hat_diff

    Σ_hat_u = fixed_param_model.Σ_hat_u + Σ_hat_u_diff
    Σ_tilde_u = fixed_param_model.Σ_tilde_u + Σ_tilde_u_diff

    Σ_β_hat = zeros(Float64, length(B_hat), length(B_hat))
    Σ_β_hat[1:d_x,1:d_x] = fixed_param_model.Σ_β_hat[1:d_x,1:d_x]
    for (i,p_i) in enumerate(p_init), (j,p_j) in enumerate(p_init)
        Σ_β_hat[1:d_x,1:d_x] += p_i * p_j .* Σ_β_hat_diff[d_x*(i-1)+1:d_x*i,d_x*(j-1)+1:d_x*j]
    end
    Σ_β_hat[d_x+1:d_x+d_x^2,d_x+1:d_x+d_x^2] = fixed_param_model.Σ_β_hat[d_x+1:d_x+d_x^2,d_x+1:d_x+d_x^2]
    for i in 1:d_x^2, j in 1:d_x^2
        k_i = i % d_x
        k_i == 0 && (k_i = d_x)
        n_i = Int64((i-k_i)/d_x)
        k_j = j % d_x
        k_j == 0 && (k_j = d_x)
        n_j = Int64((j-k_j)/d_x)
        for (l_i,p_l_i) in enumerate(p_init), (l_j,p_l_j) in enumerate(p_init)
            Σ_β_hat[d_x+i,d_x+j] += p_l_i * p_l_j * Σ_β_hat_diff[d_x*d_p+d_x*(n_i*d_p+l_i-1)+k_i,d_x*d_p+d_x*(n_j*d_p+l_j-1)+k_j]
        end
    end
    Σ_β_hat[1:d_x,1+d_x:d_x + d_x^2] = fixed_param_model.Σ_β_hat[1:d_x,1+d_x:d_x+d_x^2]
    for n in 1:d_x
        for (i,p_i) in enumerate(p_init), (j,p_j) in enumerate(p_init)
            Σ_β_hat[1:d_x,d_x+(n-1)*d_x+1:d_x+n*d_x] += p_i * p_j * Σ_β_hat_diff[1+d_x*(i-1):d_x*i,d_x*d_p+d_x*d_p*(n-1)+d_x*(j-1)+1:d_x*d_p+d_x*d_p*(n-1)+d_x*j]
        end
    end
    Σ_β_hat[1+d_x:d_x+d_x^2,1:d_x] = copy(Σ_β_hat[1:d_x,1+d_x:d_x+d_x^2])'
    Σ_β_hat[1+d_x+d_x^2:end,1+d_x+d_x^2:end] = Σ_β_hat_diff
    Σ_β_hat[1:d_x+d_x^2,1+d_x+d_x^2:end] = zeros(Float64, d_x+d_x^2, d_x*d_p + d_p*d_x^2)
    for n in 1:d_x+1
        for (i,p_i) in enumerate(p_init)
            Σ_β_hat[d_x*(n-1)+1:d_x*n,d_x+d_x^2+1:end] += p_i * Σ_β_hat_diff[d_x*d_p*(n-1)+d_x*(i-1)+1:d_x*d_p*(n-1)+d_x*i,:]
        end
    end
    Σ_β_hat[1+d_x+d_x^2:end,1:d_x+d_x^2] = copy(Σ_β_hat[1:d_x+d_x^2,1+d_x+d_x^2:end])'

    return B_hat, Σ_hat_u, Σ_tilde_u, Σ_β_hat
end


# Auxiliary function for fitVARmodel
function auxFitVARmodel(h::Union{Function,Nothing},
                    p::Int64,
                    d_y::Int64,
                    T::Int64,
                    X_arr::Matrix{Float64},
                    Y_arr::Matrix{Float64})

    # estimate the VAR model
    B_hat = X_arr * Y_arr' * (Y_arr * Y_arr')^-1
    Σ_tilde_u = (X_arr - B_hat * Y_arr) * (X_arr - B_hat * Y_arr)' / T
    Σ_hat_u = Σ_tilde_u * T / (T - p * d_y - 1)
    Γ_hat = Y_arr * Y_arr' / T
    Σ_β_hat = Matrix(Γ_hat^(-1) ⊗ Σ_hat_u) / T
    Σ_x1_hat = (T + d_y * p + 1) / T * Σ_hat_u

    return B_hat,Σ_hat_u,Σ_tilde_u,Σ_β_hat,Γ_hat,Σ_x1_hat
end
