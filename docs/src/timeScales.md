```@meta
CurrentModule = Destabilisation
```

# Determine time scales of VAR models

The following method is used to determine the time scales of a VAR model. The time scales are determined by the eigenvalues of the VAR model.
```@docs
timeScale
```

This function uses the following auxiliary function which turns a given function into a polynomial
```@docs
toPolynomial
```