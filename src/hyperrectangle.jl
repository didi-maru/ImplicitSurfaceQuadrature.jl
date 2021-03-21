"""
    HyperRectangle{N,T}(lower::SVector{N,T}, upper::SVector{N,T})
    HyperRectangle{N,T}(lower::NTuple{N,T}, upper::NTuple{N,T})

A simple hyperrectangle of dimension N
with coordinates of type T.

# Arguments
- `lower`: the lower corner of the hyperrectangle
- `upper`: the upper corner of the hyperrectangle

The assumption that `all(lower .< upper)` is true is made

# Example
A 2D hyperrectangle (i.e. a rectangle)

```julia
rec = HyperRectangle((-2., -1.), (2., 1.))
````
"""
struct HyperRectangle{N,T<:Real}
    lower::SVector{N,T} # Lower corner
    upper::SVector{N,T} # Upper corner
end

HyperRectangle(lower::NTuple{N,T}, upper::NTuple{N,T}) where {N,T<:Real} = HyperRectangle(SVector(lower), SVector(upper))

# Base.:(==)(rec1::HyperRectangle{N,T}, rec2::HyperRectangle{N,T}) where {N,T} = all(rec1.lower .== rec2.lower) && all(rec1.upper .== rec2.upper)

dim(rec::HyperRectangle) = length(rec.upper)
δ(rec::HyperRectangle) = (rec.upper .- rec.lower)/2
center(rec::HyperRectangle) = (rec.lower .+ rec.upper)/2


"Return the rectangle 4 corners coordinates in an Array for plotting"
function rectanglecoords(rec::HyperRectangle{2,<:Real})
    x1, y1 = rec.lower
    x2, y2 = rec.upper
    [x1 y1; x2 y1; x2 y2; x1 y2; x1 y1]
end


"""
    split(rec::HyperRectangle[, axis::Int]) -> (HyperRectangle, HyperRectangle)

Split a hyperrectangle in 2 halfs along a given axis,
split along the longer axis if none is given.
    
Construct and return the 2 half-hyperrectangles.
# """
function split(rec::HyperRectangle, axis::Int)
    middle = (rec.upper[axis] + rec.lower[axis])/2
    newlower = SVector{dim(rec)}(i == axis ? middle : rec.lower[i] for i = 1:dim(rec))
    newupper = SVector{dim(rec)}(i == axis ? middle : rec.upper[i] for i = 1:dim(rec))
    HyperRectangle(rec.lower, newupper), HyperRectangle(newlower, rec.upper)
end

split(rec::HyperRectangle) = split(rec, argmax(rec.upper .- rec.lower))
