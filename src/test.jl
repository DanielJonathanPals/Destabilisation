using .Destabilisation
using GLMakie

prog(v,p) = [p[1] + p[2]*v[1] + 5e-1*randn()]
obs(v,p) = v
v_init = [0.]
p_init = [0.,0.9]

DS = DynamicalSystem(prog,obs,v_init,p_init)
p_traj = parameterSeriesGenerator(1/0.9,p_init,[1e-3,1e-3],n=100)
p_tr = zeros(Float64,2,10000)
n = size(p_traj,2)
for i in 1:n
    p_tr[1,10*(i-1)+1:100*i] .= p_traj[1,i]
    p_tr[2,10*(i-1)+1:100*i] .= p_traj[2,i]
end

g(p) = p_init*0.1 + 0.9*p + 0.002 .*randn(2)
#g(p) = p

#v_tr, p_tr, x_tr = integrateTraj(DS,g,101,v_init,p_init)
#v_tr, p_tr, x_tr = integrateTraj(DS,g,1001,v_init,p_init)
_,_, x_tr = integrateTraj(DS,p_tr)

h(x,p) = x[1] .* p

model = fitVARmodel(x_tr,p_traj=p_tr,h=h,p=1)
#testParamCausality(model,1,101)
#model = fitVARmodel(v_tr,p=1)
timeScale(model)

#LMtest(model,5,v_tr,p_traj=p_tr)

.√([model.Σ_β_hat[i,i] for i in 1:6])
model.B_hat


"""
x = zeros(1,1003)
for i in 4:1003
    x[1,i] = 0.5*x[1,i-1] + 0.01*x[1,i-2] + 0.3*x[1,i-3] + 0.2*randn()
end

x1 = ones(2,1003) .* 100
for i in 4:1003
    x1[:,i] = 10 .+ 0.9.*x1[:,i-1] + 0.2*randn(2)
end

x2 = zeros(1,1003)
for i in 4:1003
    x2[1,i] = 0.9*x2[1,i-1] + 0.2*randn()
end
"""


"""
p = VARorder(x, criterion="AIC",p_max=5)
model = VARmodel(x[:,1:900],p=2)
testPredictions(model,x[:,900:end])
"""

.√([model.Σ_β_hat[i,i] for i in 1:6])
model.B_hat

testPredictions(model,x_tr,p_traj=p_tr)

lines(x_tr[1,:])
lines!(p_tr[1,:]./(1 .-p_tr[2,:]))


prog(v,p,r) = [p[1] + (p[1]+p[2])*v[1] + 5e-1*r[1]]
prog(v,p,r) = [-0.01*v[1]^3 + p[1]*v[1]^2 + 1.01*v[1] + p[2] + 1e-4*r[1]]
#prog(v,p,r) = [-0.02 + 0.98*v[1] - p[1] + p[2] - 2*v[1]*p[1] + 1e-4*r[1]]
obs(v,p) = v
v_init = [-1.]
p_init = [0.,0.]
DS = DynamicalSystem(prog,obs,v_init,p_init;random_vec_length=1)
p_traj = parameterSeriesGenerator(1/0.9,p_init,[1e-10,1e-10],n=1000, keep_const = 20)

_,_,x_tr = integrateTraj(DS,1000)
ref_model = fitVARmodel(x_tr,p=1)
v_arr, p_traj, x_arr, v_ref_arr, x_ref_arr, noise = integrateTraj(DS, p_traj; include_reference_traj=true)
B_hat, Σ_hat_u, Σ_tilde_u, Σ_β_hat = fitVARmodel(x_arr,x_ref_arr,p_traj,p_init,ref_model)

.√([Σ_β_hat[i,i] for i in 1:6])
B_hat
