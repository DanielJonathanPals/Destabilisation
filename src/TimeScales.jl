# The methods presented here can be used to determine the timescale of a given VAR model. These also indicate
# the stability of the VAR model.

using .VARmodel_module
using LinearAlgebra
using Polynomials

"""
    toPolynomial(f::Function; deg::Int64 = 1, tolerance::Float64 = 1e-4)

Fit a polynomial of degree `deg` to the function `f`. The polynomial is fitted using the least squares method. The
    accuracy of the polynomial is checked by comparing the function `f` to the polynomial on 100 equidistant points
    in an appropriately chosen interval. If the maximum error is larger than `tolerance`, a warning is issued.
    Finally the polynomial is returned.
"""
function toPolynomial(f::Function; deg::Int64 = 1, tolerance::Float64 = 1e-4)
    # fit a polynomial of degree deg to the function f
    x = Array(LinRange(-1,1,deg+1))
    X = ones(deg+1,deg+1)
    for i in 1:deg+1
        X[:,i] = x.^(i-1)
    end
    y = X^-1 * f.(x)
    poly = Polynomial(y)

    # determine the interval relevant for checking the accuracy of the polynomial

    # At first we find the degrees of the coeficients which are larger than the coeficient of the highest degree
    idx = findall(x -> abs(x) > y[end], y[1:end-1])

    # if all coeficients are smaller than the coeficient of the highest degree, we use the interval [-3,3] since then the
    # highest order dominates the polynomial by a factor of at least 3.
    if idx == []
        test_interval = (-3,3) 

    # otherwise we determine an interval that includes [-3,3] for which the highest order dominates the polynomial at
    # the boundaries of the interval by a factor of at least 3.
    else
        x_bound = max(3, (abs(3*y[idx[end]]/y[end]))^(1/(deg-idx[end])))
        test_interval = (-x_bound,x_bound)
    end

    # check that the polynomial is accurate enough
    warn = false
    crit = zeros(2)
    for i in Array(LinRange(test_interval[1],test_interval[2],100))
        if abs(f(i) - poly(i)) > tolerance
            warn = true
            if abs(f(i) - poly(i)) > crit[2]
                crit = [i, abs(f(i) - poly(i))]
            end
        end
    end
    if warn
        @warn "Polynomial approximation is not accurate enough. Maximum error is $(crit[2]) at $(crit[1]) in [$(test_interval[1]),$(test_interval[2])]."
    end

    return poly
end


"""
    timeScale(model::VARmodel)

Determines a proxy of the slowest timescale of the VAR model `model` by returning the absolute value of the smallest
    root of the characteristic polynomial of the VAR model. This is done by fitting a polynomial of degree `deg` to
    the function `f(z) = det(Id - sum([z^i * B_i for i in 1:p]))` where `B_i` is the `i`-th block of the matrix `B_hat`
    of the VAR model `model`. The root of this polynomial with the smallest absolute value (which should be larger than
    one for the model to be stable) is returned.
    Note that this function can only be applied to VAR models without parameter variables and without hidden variables.
"""
function timeScale(model::VARmodel)

    # check that model has no parameter variables and no hidden variables
    if model.d_p != 0
        error("This function only works for VAR models without parameter variables")
    end
    if model.d_h != 0
        error("This function only works for VAR models without hidden variables")
    end

    function f(z)
        dx = model.d_x
        B = model.B_hat
        return det(Id(dx) - sum([z^i .* B[:,2+dx*(i-1):dx*i+1] for i in 1:model.p]))
    end

    poly = toPolynomial(f,deg=model.d_x*model.p)

    rts = roots(poly)

    return min(abs.(rts)...)
end


function timeScale(DS::DynamicalSystem; T_init::Int64=200)
    
end