using .Destabilisation
using GLMakie

prog(v,p) = [p[1] + p[2]*v[1] + 1e-10*randn()]
obs(v,p) = v
v_init = [0.]
p_init = [0.,0.9]

DS = DynamicalSystem(prog,obs,v_init,p_init)

g(p) = p_init*0.1 + 0.9*p + 0.002*randn(2)

v_tr, p_tr, x_tr = integrateTraj(DS,g,101,v_init,p_init)

h(x,p) = x[1] .* p


x = zeros(1,23)
for i in 4:13
    x[1,i] = 0.5*x[1,i-1] + 0.01*x[1,i-2] + 0.3*x[1,i-3] + 0.002*randn()
end


model = VARmodel(v_tr,p_traj=p_tr,h=h,p=1)
p = VARorder(x, criterion="AIC",p_max=5)

.√([model.Σ_β_hat[i,i] for i in 1:6])
model.B_hat

lines(x_tr[1,:])