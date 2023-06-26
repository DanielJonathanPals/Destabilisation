```@meta
CurrentModule = Destabilisation
```

# Test estimated VAR models

Here we pesent some tests that we can perform on estimated VAR models in order to check their validity.

In this first test we present a function that quantifies the validity of a model by comparing the one step predictions of the model with the actual values of the observables by performing a hypothesis test.
```@docs
testPredictions
```

The whiteness i.e. the absence of autocorrelation in the residuals of the model is tested by the following function.
```@docs
LMtest
```