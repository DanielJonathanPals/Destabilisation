```@meta
CurrentModule = Destabilisation
```

# Test Parameter Dependence

The following function can be used to find the locations of the parameter dependent coefficients in the regressor of a given VAR model

```@docs
getParamLocations
```

In order to determine if the observables significantly depend on the parameters, the following function can be used

```@docs
testParamCausality
```