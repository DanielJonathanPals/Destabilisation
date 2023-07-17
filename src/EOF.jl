using .VARmodel_module
using LinearAlgebra


function EOF(progressor::Function; components::Int64=1)
    cov = model.Γ_hat[2:model.d_x+1,2:model.d_x+1]
    Λ, U = eigen(cov)
end

