using .Destabilisation
using GLMakie

prog(v,p) = [p[1] + p[2]*v[1] + 1e-4*randn()]
obs(v,p) = v
v_init = [0.]
p_init = [0.,0.9]

DS = DynamicalSystem(prog,obs,v_init,p_init)

#g(p) = p_init*0.1 + 0.9*p + 0.002*randn(2)
g(p) = p

#v_tr, p_tr, x_tr = integrateTraj(DS,g,1001,v_init,p_init)
v_tr, p_tr, x_tr = integrateTraj(DS,g,1001,v_init,p_init)

#h(x,p) = x[1] .* p

#model = fitVARmodel(v_tr,p_traj=p_tr,h=h,p=1)
model = fitVARmodel(v_tr,p=1)
timeScale(model)

#LMtest(model,5,v_tr,p_traj=p_tr)

.√([model.Σ_β_hat[i,i] for i in 1:2])
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