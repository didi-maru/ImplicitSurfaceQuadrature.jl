# ImplicitSurfaceQuadrature.jl
A Julia package implementing the high-order quadrature method for implicitly defined domains described in the article
[R. I. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015)](http://dx.doi.org/10.1137/140966290), a C++ implementation of the same method exists, see the [Algoim GitHub page](https://algoim.github.io/).

The `generatequadrature` function returns quadrature nodes as `Array{SVector{N,T},1}` weights as `Array{T,1}` and 

**This Julia implementation only works for 2D domains.**

## Overview of the method
This Julia package can be used to produce a high-order quadrature scheme for a curved surface implicitly defined by a function ϕ: ℝᵈ → ℝ, and contained in a hyperrectangle U.
The integration domain is U ∩ Γ where Γ = {x: ϕ(x) = 0}.

The produced quadrature scheme is based on multiples `q` points Gauss-Legendre quadrature schemes.

## Example usage
Integration of f(x,y) = x² + y² on the circle of radius `0.5`
```julia
f(x) = sum(x.^2)

ϕ(x) = sum(x.^2) - 0.5 # Defining the isosurface (see Restrictions on the isosurface below)
a, b = (-1., -1.), (1., 1.) # Defining the hyperrectangle containing the integration domain

nodes, weights = generatequadrature(4, a, b, ϕ) # Compute the quadrature scheme, with 4 points Gauss-Legendre quadratures
int = dot( weights, f.(nodes) ) # Compute the integral
```

## Restrictions on the isosurface
The isosurface function `ϕ` should support the method `ϕ(x::NTuple{N,T}) where {N,T<:Real}`.

The quadrature method uses an automatic computation of first-order Taylor in order to compute estimations of bounds of functions.
Thus, the isosurface function should only be defined with operations that are also defined for `Linearization` objects in `src/linerization.jl`, which are:
- `+`, `-`, `*` between function arguments and reals (ex: `ϕ(x) = x[1] + 3.`)
- `+`, `-`, `*` between multiples function arguments (ex: `ϕ(x) = x[1] + x[2]`)
- `abs`, `^n`, `/c` on function arguments with `isa(n,Integer)` and `isa(c,Real)` (ex: `ϕ(x) = abs(x[1]^2 + x[2]/2)`)


