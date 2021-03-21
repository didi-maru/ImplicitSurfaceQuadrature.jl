"""
    Linearization{R::Hyperrectangle,N,T<:Real}(α::T, β::SVector{N,T}, δ::SVector{N,T}, ϵ::T) <: Real

A first-order Taylor series with bounded remainder term for some function f: ℝᵈ -> ℝ
centered at the midpoint point xc of the hyperrectangle `R`, with `δ = δ(R)`.

f (xc + y) = α + β · y + r(y) for all y ∈ ℝᵈ

    Linearization(rec::Rectangle{N,T}, axis::Integer)

Linearization for the hyperrectangle `rec` and the base function f:(x₁,...,xₙ) -> xᵢ where `i == axis`.
"""
struct Linearization{R,N,T<:Real} <: Real
    α::T
    β::SVector{N,T}
    δ::SVector{N,T}
    ϵ::T
end

function Linearization{R,N,T}(α::T, β::SVector{N,T}, ϵ::T) where {R,N,T}
    Linearization{R,N,T}(α, β, δ(R), ϵ)
end

function Linearization(R::HyperRectangle{N,T}, axis::Integer) where {N,T}
    0 < axis <= N || error("axis should be in {1,...,$N}")
    α = center(R)[axis]
    β = SVector{N}(i == axis ? one(T) : zero(T) for i = 1:N)
    Linearization{R,N,T}(α, β, δ(R), zero(T))
end


"""
    linearize(f::Function, rec::HyperRectangle{N,T}) where {N,T<:Real} -> Linearization{rec,N,T}

Compute a linearization for the fonction `f(x::NTuple{N,T}) <: T`
or `f(x::NTuple{N,T}) <: NTuple{M,T}` in a HyperRectangle `rec`.

The function `f` should be composed with the following operations :
`+`, `-`, `*`, `/c` and `^c` with `c` a real constant.
`Base.abs` is also supported

# Examples
For f: ℝᵈ -> ℝ
```julia
f(vec) = sum(vec.^2) - 1. # Defining the function f:(x,y) -> x² + y² - 1 
rec = HyperRectangle((0.,0.), (1.,1.)) # Defining a 2D hyperrectangle
lin = linearize(f, rec) # get the linearization of f on rec
```

For f: ℝᵈ -> ℝ²
```julia
f(x) = (vec[1]^2, 2vec[2]). # Defining the function f:(x,y) -> (x², 2y) 
rec = HyperRectangle((0.,0.), (1.,1.)) # Defining a 2D hyperrectangle
linx, liny = linearize(f, rec) # linx is the linearization of x -> x² on rec
                               # liny is the linearization of y -> 2y on rec
```
"""
linearize(f::Function, rec::HyperRectangle) = f(SVector{dim(rec)}(Linearization(rec, i) for i = 1:dim(rec)))


"""
    maxgap(lin::Linearization{rec,N,T})
    maxgap(f::Function, rec::HyperRectangle)

Compute an approximation of sup_{x ∈ rec}|f(x) - f(xc)|
with xc the center of the hyperrectangle `rec`
"""
maxgap(lin::Linearization) = dot(abs.(lin.β), lin.δ) + lin.ϵ
maxgap(f::Function, rec::HyperRectangle) = maxgap(linearize(f, rec))


# Operations on the type ::Type{Linearization}

Base.one(::Type{Linearization{R,N,T}}) where {R,N,T} = Linearization{R,N,T}(one(T), zero(SVector{N,T}), zero(T))
Base.zero(::Type{Linearization{R,N,T}}) where {R,N,T} = Linearization{R,N,T}(zero(T), zero(SVector{N,T}), zero(T))


# Operations on a single Linearization

Base.:(-)(lin::Linearization{R,N,T}) where {R,N,T} = Linearization{R,N,T}(-lin.α, -lin.β, lin.δ, lin.ϵ)
Base.abs(lin::Linearization{R,N,T}) where {R,N,T} = Linearization{R,N,T}(abs(lin.α), sign(lin.α)lin.β, lin.δ, lin.ϵ)


# Operations between Linearization and Real

Base.:(+)(lin::Linearization{R,N,T}, c::Real) where {R,N,T} = Linearization{R,N,T}(lin.α + c, lin.β, lin.δ, lin.ϵ)
Base.:(+)(c::Real, lin::Linearization{R,N,T}) where {R,N,T}  = lin + c

Base.:(-)(lin::Linearization{R,N,T}, c::Real) where {R,N,T} = lin + (-c)
Base.:(-)(c::Real, lin::Linearization{R,N,T}) where {R,N,T} = lin + (-c)

Base.:(*)(lin::Linearization{R,N,T}, c::Real) where {R,N,T} = Linearization{R,N,T}(c*lin.α, c*lin.β, lin.δ, abs(c)lin.ϵ)
Base.:(*)(c::Real, lin::Linearization{R,N,T}) where {R,N,T} = lin * c

Base.:(/)(lin::Linearization{R,N,T}, c::Real) where {R,N,T} = lin * (1/c)


# Operations between 2 Linearizations

Base.:(+)(lin1::Linearization{R,N,T}, lin2::Linearization{R,N,T}) where {R,N,T} = Linearization{R,N,T}(lin1.α+lin2.α, lin1.β+lin2.β, lin1.δ, lin1.ϵ+lin2.ϵ)
Base.:(-)(lin1::Linearization{R,N,T}, lin2::Linearization{R,N,T}) where {R,N,T} = lin1 + -lin2

function Base.:(*)(lin1::Linearization{R,N,T}, lin2::Linearization{R,N,T}) where {R,N,T}
    l1 = dot(abs.(lin1.β), lin1.δ)
    l2 = dot(abs.(lin2.β), lin2.δ)
    Linearization{R,N,T}(
        lin1.α*lin2.α,
        lin1.α*lin2.β + lin2.α*lin1.β,
        lin1.δ,
        l1*l2 + (abs(lin1.α)+l1)lin2.ϵ + (abs(lin2.α)+l2)lin1.ϵ + lin1.ϵ*lin2.ϵ
    )
end


# Others operations

function Base.:(^)(lin::Linearization, n::Integer)
    if n < 1
        error("Exposant should be positive")
    elseif n == 1
        return lin
    else
        return lin * lin^(n-1)
    end
end