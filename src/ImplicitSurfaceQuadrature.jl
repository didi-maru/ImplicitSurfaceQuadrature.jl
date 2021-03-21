module ImplicitSurfaceQuadrature

using Roots
using ForwardDiff
using StaticArrays
using LinearAlgebra
using FastGaussQuadrature

export generatequadrature, HyperRectangle, rectanglecoords, Linearization

include("hyperrectangle.jl")
include("linearization.jl")
include("binarytree.jl")
include("quadrature.jl")

"""
    generatequadrature(q::Integer, lowercorner::NTuple{N,T}, uppercorner::NTuple{N,T},
                       isosurface::Function [, gradient=nothing; minrecsize::Real=Inf ,maxtreedepth=25]) where {N,T<:Real}

Compute the quadrature nodes and weights in the integration domain U ∩ {x; ϕ(x) = 0},
with `isosurface` = ϕ and U the hyperrectangle (x₁ᴸ,x₁ᵁ)×⋯×(xₙᴸ,xₙᵁ) where `lowercorner` = (x₁ᴸ,…,xₙᴸ), `uppercorner` = (x₁ᵁ,…,xₙᵁ). 

# Arguments
- `q::Int`: the order of the Gauss-Legendre quadrature used
- `lowercorner::NTuple{N,T}`: the lower corner of the hyperrectangle
- `uppercorner::NTuple{N,T}`: the upper corner of the hyperrectangle
- `iosurface::Function`: the isosurface function which defined the integration domain
- `gradient::Function`: the gradient of the isosurface, will be computed with `ForwardDiff.jl` is none is given

# Keyword Arguments
- `minrecsize::Real=Inf`: minimum size for hyperrectangles before the computation of the quadrature
- `maxtreedepth=25`: the maximum depth for the binary tree of hyperrectangles
"""
function generatequadrature(q::Integer, lowercorner::NTuple{N,T}, uppercorner::NTuple{N,T}, isosurface::Function, gradient=nothing; minrecsize::Real=Inf ,maxtreedepth=25) where {N,T<:Real}
    rec = HyperRectangle(lowercorner, uppercorner)
    gradient = gradient === nothing ? x -> ForwardDiff.gradient(isosurface, x) : gradient
    glnodes, glweights = gausslegendre(q)

    quadraturetree = maketree!(TreeNode(rec), isosurface, gradient, glnodes, glweights, minrecsize, maxtreedepth, 0)
    get_nodes_weights(quadraturetree)
end

function generatequadrature_and_tree(q::Integer, lowercorner::NTuple{N,T}, uppercorner::NTuple{N,T}, isosurface::Function, gradient=nothing; minrecsize::Real=Inf ,maxtreedepth=25) where {N,T<:Real}
    rec = HyperRectangle(lowercorner, uppercorner)
    gradient = gradient === nothing ? x -> ForwardDiff.gradient(isosurface, x) : gradient
    glnodes, glweights = gausslegendre(q)

    quadraturetree = maketree!(TreeNode(rec), isosurface, gradient, glnodes, glweights, minrecsize, maxtreedepth, 0)
    nodes, weigths = get_nodes_weights(quadraturetree)
    nodes, weigths, quadraturetree
end

end # module
