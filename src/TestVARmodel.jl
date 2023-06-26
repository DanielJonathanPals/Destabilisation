# Here some methods are prsented that are used to check the validity of a given estimated VAR model

using Distributions
using Kronecker

include("FitVARmodel.jl")


"""
    testPredictions(model::VARmodel,
        x_traj::Matrix;
        p_traj::Union{Matrix,Nothing}=nothing,
        siglvl::Float64=0.05)

This function tests the validity of the given VAR model by checking the prediction error of the model. 
    The prediction error is computed as the sum of the squared nomalized, uncorrlated prediction errors 
    of the model. 
    The prediction error is then compared to the critical value of the χ² distribution with degrees 
    of freedom equal to the number of observations times the number of variables.
    If the percentil of the prediction error is smaller than the significance level `siglvl`, the model
    is rejected.
    The function returns a tuple consisting of a boolean value whether the model is accepted or not
    and the percentil of the prediction error. 
"""
function testPredictions(model::VARmodel,
    x_traj::Matrix;
    p_traj::Union{Matrix,Nothing}=nothing,
    siglvl::Float64=0.05)

    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj,model.h)
    if (d_x, d_p, d_h) != (model.d_x, model.d_p, model.d_h)
        error("The dimensions of the input data do not match the dimensions of the model")
    end
    if size(x_traj,2) <= model.p
        error("The trajectory must contain more than `p`=$(model.p) datapoints")
    end

    x_traj = convert(Matrix{Float64}, x_traj)
    (p_traj === nothing) || (p_traj = convert(Matrix{Float64}, p_traj)) 
    T = size(x_traj,2)

    # compute the statistic
    preds = zeros(d_x,T-model.p)
    if p_traj !== nothing
        for t in model.p+1:T
            preds[:,t-model.p] = oneStepPred(model,x_traj[:,t-model.p:t-1],p_traj=p_traj[:,t-model.p:t-1])
        end
    else
        for t in model.p+1:T
            preds[:,t-model.p] = oneStepPred(model,x_traj[:,t-model.p:t-1])
        end
    end
    errors = x_traj[:,model.p+1:end] - preds
    statistic = sum([errors[:,i]' * model.Σ_x1_hat^(-1) * errors[:,i] for i in 1:T-model.p])

    # check the statistic
    if statistic > cquantile(Chisq((T-model.p)*d_x), siglvl)
        return false, ccdf(Chisq((T-model.p)*d_x), statistic)
    else
        return true, ccdf(Chisq((T-model.p)*d_x), statistic)
    end
end


"""
    LMtest(model::VARmodel,h::Int64,x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing,siglvl::Float64=0.05)

This function is the implementation of the Lagrange multiplier test for testing the whiteness of the residuals.
    Here `x_traj` and `p_traj` are the trajectories used to estimate the model. `h` determines up to how many lags
    we test for uncorrelated residualts. `siglvl` is the significance level of the test.
"""
function LMtest(model::VARmodel,h::Int64,x_traj::Matrix;p_traj::Union{Matrix,Nothing}=nothing,siglvl::Float64=0.05)

    # check input data
    d_x,d_p,d_h = check_xph(x_traj,p_traj,model.h)
    if (d_x, d_p, d_h) != (model.d_x, model.d_p, model.d_h)
        error("The dimensions of the input data do not match the dimensions of the model")
    end
    if size(x_traj,2) <= model.p + h
        error("The trajectory must contain more than `p+h`=$(model.p+h) datapoints")
    end

    # compute the statistic
    T = size(x_traj,2) - model.p
    y_traj = create_y_traj(x_traj,p_traj=p_traj,h=model.h)
    Y_arr = Y(y_traj,p=model.p)
    U_hat = x_traj[:,model.p+1:end] - model.B_hat * Y_arr
    F_arr = F(h,T)
    mathcal_U_hat = (Id(h) ⊗ U_hat) * F_arr'
    λ_LM = vec(U_hat*mathcal_U_hat')' * ((mathcal_U_hat * mathcal_U_hat' - mathcal_U_hat * Y_arr' * (Y_arr * Y_arr')^(-1) * Y_arr * mathcal_U_hat')^(-1) ⊗ (model.Σ_hat_u)^(-1)) * vec(U_hat * mathcal_U_hat')

    # check the statistic
    if λ_LM > cquantile(Chisq(h*d_x^2), siglvl)
        return false, ccdf(Chisq(h*d_x^2), λ_LM)
    else
        return true, ccdf(Chisq(h*d_x^2), λ_LM)
    end
end