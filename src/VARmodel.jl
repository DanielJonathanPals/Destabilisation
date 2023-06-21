# An element of the struct VARmodel contains all the information needed to define a VAR model


"""
    VARmodel

An element of the struct VARmodel contains all the information needed to define a VAR model.

# Fields
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
end