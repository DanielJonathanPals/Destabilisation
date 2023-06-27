# An element of the struct VARmodel contains all the information needed to define a VAR model

module VARmodel_module

export VARmodel
export oneStepPred

include("FormatTests.jl")
include("Matrices.jl")

"""
    VARmodel

An element of the struct VARmodel contains all the information needed to define a VAR model.

# Fields
- `h::Union{Function,Nothing}`: function that computes the regressor of the form `(x_t',p_t',h(x_t,p_t)')'`
- `p::Int64`: order of the VAR model
- `d_x::Int64`: dimension of the observable variable
- `d_p::Int64`: dimension of the parameter variable
- `d_h::Int64`: dimension of the hidden variable
- `d_y::Int64`: dimension of the augmented variable
- `B_hat::Matrix{Float64}`: estimated coefficient matrix
- `Σ_hat_u::Matrix{Float64}`: estimated covariance matrix of the residuals (unbiased estimator)
- `Σ_tilde_u::Matrix{Float64}`: estimated covariance matrix of the residuals (maximum likelihood)
- `Σ_β_hat::Matrix{Float64}`: estimated covariance matrix of the coefficient matrix (unbiased estimator)
- `Γ_hat::Matrix{Float64}`: estimator for YY'/T
- `Σ_x1_hat::Matrix{Float64}`: estimated MSE matrix of the one step forecast (unbiased estimator)
"""
struct VARmodel
    h::Union{Function,Nothing}      # is needed to compute the regressor of the form `(x_t',p_t',h(x_t,p_t)')'`
    p::Int64                        # order of the VAR model
    d_x::Int64                      # dimension of the observable variable
    d_p::Int64                      # dimension of the parameter variable
    d_h::Int64                      # dimension of the hidden variable
    d_y::Int64                      # dimension of the augmented variable
    B_hat::Matrix{Float64}          # estimated coefficient matrix
    Σ_hat_u::Matrix{Float64}        # estimated covariance matrix of the residuals (unbiased estimator)
    Σ_tilde_u::Matrix{Float64}      # estimated covariance matrix of the residuals (maximum likelihood)
    Σ_β_hat::Matrix{Float64}        # estimated covariance matrix of the coefficient matrix (unbiased estimator)
    Γ_hat::Matrix{Float64}          # estimator for YY'/T
    Σ_x1_hat::Matrix{Float64}       # estimated MSE matrix of the one step forecast (unbiased estimator)

    function VARmodel(h,p,dx,dp,dh,dy,B,Σ_hat_u,Σ_tilde_u,Σ_β,Γ,Σ_x1)
        if (h === nothing) 
            if dh != 0
                error("The dimension of the hidden variable `d_h` is not zero but no function h is given")
            end
        else
            if dp == 0
                check_h(h,ones(1,dx),nothing)
                if dh != length(h(ones(dx)))
                    error("The dimension of the hidden variable `d_h` does not match the dimension of the function h")
                end
            else
                check_h(h,ones(1,dx),ones(1,dp))
                if dh != length(h(ones(dx),ones(dp)))
                    error("The dimension of the hidden variable `d_h` does not match the dimension of the function h")
                end
            end
        end
        if p < 1
            error("The order of the VAR model must be at least 1")
        end
        if dy != dx + dp + dh
            error("The dimension of the augmented variable `d_y` must be equal to `d_x + d_p + d_h`")
        end
        if size(B) != (dx,p*dy+1)
            error("The dimension of the coefficient matrix `B` must be equal to `d_x x p*d_y+1`")
        end
        if size(Σ_hat_u) != (dx,dx)
            error("The dimension of the covariance matrix of the residuals `Σ_hat_u` must be equal to `d_x x d_x`")
        end
        if size(Σ_tilde_u) != (dx,dx)
            error("The dimension of the covariance matrix of the residuals `Σ_tilde_u` must be equal to `d_x x d_x`")
        end
        if size(Σ_β) != (dx*(p*dy+1),dx*(p*dy+1))
            error("The dimension of the covariance matrix of the coefficient matrix `Σ_β` must be equal to `(p*dy+1)*d_x x (p*dy+1)*d_x`")
        end
        if size(Γ) != (p*dy+1,p*dy+1)
            error("The dimension of the estimator for YY'/T `Γ` must be equal to `d_y+1 x d_y+1`")
        end
        if size(Σ_x1) != (dx,dx)
            error("The dimension of the MSE matrix of the one step forecast `Σ_x1` must be equal to `d_x x d_x`")
        end
        new(h,p,dx,dp,dh,dy,B,Σ_hat_u,Σ_tilde_u,Σ_β,Γ,Σ_x1)
    end
end


"""
    oneStepPred(model::VARmodel,x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing)

Compute the one step prediction of the VAR model given by `model` for the trajectory `x_traj` and the parameter
trajectory `p_traj`. If `p_traj` is not given, it is assumed that the VAR model is time invariant.
"""
function oneStepPred(model::VARmodel,x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing)

    # Check that the input trajectory has the correct dimensions
    d_x,d_p,_ = check_xph(x_traj,p_traj,model.h)
    if (d_x, d_p) != (model.d_x, model.d_p)
        error("The dimensions of the input data do not match the dimensions of the model")
    end
    if size(x_traj,2) < model.p
        error("The input trajectory is too short")
    end

    x_traj = convert(Matrix{Float64}, x_traj) 
    (p_traj === nothing) || (p_traj = convert(Matrix{Float64}, p_traj))

    # Create y_traj from x_traj, p_traj and h (here in order to later use the function Y we need to add 
    # a dummy data point at the end of the trajectory)
    y_traj = [create_y_traj(x_traj,p_traj=p_traj,h=model.h) ones(model.d_y)]

    # Compute the one step prediction
    return model.B_hat*Y(y_traj[:,end-model.p:end],p=model.p,T=1)
end

end