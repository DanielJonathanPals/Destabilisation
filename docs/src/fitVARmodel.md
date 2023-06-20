```@meta
CurrentModule = Destabilisation
```

# Fitting a VAR model to observed data

The VAR model we use in this case is of the following form:
$x_t = ν + A_1 (x_traj_{t-1}', p_traj_{t-1}', h(x_traj_{t-1}, p_traj_{t-1})')'
    + ... + A_p (x_traj_{t-p}', p_traj_{t-p}', h(x_traj_{t-p}, p_traj_{t-p})')'$

where $x_traj_t$ is the trajectory of the observables at time $t$, $p_traj_t$ is the trajectory of the parameters at time $t$ and $h$ is a function that allows for non-linearities in the model.

A suitable order $p$ for the VAR model can be determined using the function:
```@docs
VARorder
```

After selecting a suitable order the VAR model can be fitted to the data using the function:
```@docs
VARmodel
```