using ImplicitSurfaceQuadrature

using LinearAlgebra
using Plots

rectangleshape(rec::HyperRectangle) = Shape(rectanglecoords(rec)[:,1], rectanglecoords(rec)[:,2])

function displaynodes(node::ImplicitSurfaceQuadrature.TreeNode)
    if node.pruning == ImplicitSurfaceQuadrature.inside
        plot!(rectangleshape(node.rec), c=:blue, opacity = 0.4, label=:none)
    elseif node.pruning == ImplicitSurfaceQuadrature.outside
        plot!(rectangleshape(node.rec), c=:orange, opacity = 0.4, label=:none)
    elseif node.levelsetaxis > 0
        plot!(rectangleshape(node.rec), c=:green, opacity = 0.4, label=:none)
        nodesx = [node.quadraturenodes[i][1] for i = 1:length(node.quadraturenodes)]
        nodesy = [node.quadraturenodes[i][2] for i = 1:length(node.quadraturenodes)]
        length(nodes) > 0 && scatter!(nodesx, nodesy, label=:none)
    else
        node.childleft  !== nothing && displaynodes(node.childleft)
        node.childright !== nothing && displaynodes(node.childright)
    end
end


a, b = (-1.5, -1.5), (2.8, 3.2)

radius = 1.3

ϕ(x,y) = x^2 + y^2 - radius^2
ϕ(x) = ϕ(x[1], x[2])
∇ϕ(x) = (2x[1], 2x[2])

f(x) = 1.

nodes, weights, treeroot = ImplicitSurfaceQuadrature.generatequadrature_and_tree(4, a, b, ϕ)
int = dot( weights, f.(nodes) )
println("approx: ", int, ", exact: ", 2π*radius, ", relative error: ", abs((int-2π*radius)/2π*radius))


plot(aspect_ratio = 1)
displaynodes(treeroot)

θ = range(0., 2π, length=100)
plot!(radius*cos.(θ), radius*sin.(θ), lw=2, label="ϕ=0")
