using .Destabilisation
using GLMakie

prog(v,p) = [v[1] + 0.01*(-p[1]*(v[1]-p[2]) + 0.1*randn())]
obs(v,p) = v
v_init = [2.]
p_init = [3.,2.]

DS = DynamicalSystem(prog,obs,v_init,p_init)

g(p) = p

v_tr, p_tr, x_tr = integrateTraj(DS,g,100,v_init,p_init)

lines(x_tr[1,:])