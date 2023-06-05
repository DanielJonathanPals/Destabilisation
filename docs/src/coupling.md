```@meta
CurrentModule = Destabilisation
```

# Coupling the destabilataions package to the dynamical system

The struct `DynamicalSystem` allows acces to the parameter dependent dynamics of the underlying system and also gives information about how to comupte the observables of interest.

```@docs
DynamicalSystem
DynamicalSystem(progressor,observable,x_init,p_init;p_bounds=nothing)
```