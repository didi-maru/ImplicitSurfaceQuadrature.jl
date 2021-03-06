using ImplicitSurfaceQuadrature
using LinearAlgebra

a, b = (-1.5, -1.5), (2.8, 3.2)

radius = 1.3

ϕ(x,y) = x^2 + y^2 - radius^2
ϕ(x) = ϕ(x[1], x[2])
∇ϕ(x) = (2x[1], 2x[2])

f(x,y) = x^2 + y^2
f(x) = f(x[1], x[2])

solution = 2π*radius^3

println("Without specifying the gradient (uses ForwardDiff.jl)")
nodes, weights = generatequadrature(4, a, b, ϕ)
int = dot( weights, f.(nodes) )
println("approx: ", int, ", exact: ", solution, ", relative error: ", abs((int-solution)/solution), "\n")


println("Specifying the gradient")
nodes, weights = generatequadrature(4, a, b, ϕ, ∇ϕ)
int = dot( weights, f.(nodes) )
println("approx: ", int, ", exact: ", solution, ", relative error: ", abs((int-solution)/solution), "\n")
