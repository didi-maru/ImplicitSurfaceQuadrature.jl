const C = 20 # see article 3.2.3 (ii) (and p.14)

@enum Pruning none inside outside

"""
    TreeNode(rec::HyperRectangle{N,T}) where {N,T} -> TreeNode{N,T}

A binary tree element used to store the succesive cutting of the integration domain and others related variable.

# Struct fields
- `rec::HyperRectangle{N,T}`: the hyperrectangle represented by this `TreeNode` 
- `quadratureweights::Vector{T}`: the quadrature weigths in the hyperrectangle `rec` (stay empty if there is none)
- `quadraturenodes::Vector{SVector{N,T}}`: the quadrature nodes in the hyperrectangle `rec` (stay empty if there is none)
- `pruning::Pruning`: the status of the pruning (see `ImplicitSurfaceQuadrature.Pruning`)
- `levelsetaxis::Int`: an admissible levelset axis in the hyperrectangle `rec` (`0` if there is none)
- `childleft::Union{TreeNode, Nothing}`: the first child node of this node in the binary tree
- `childright::Union{TreeNode, Nothing}`: the second child node of this node in the binary tree
"""
mutable struct TreeNode{N,T<:Real}
    rec::HyperRectangle{N,T}
    quadratureweights::Vector{T}
    quadraturenodes::Vector{SVector{N,T}}
    pruning::Pruning
    levelsetaxis::Int
    childleft::Union{TreeNode, Nothing}
    childright::Union{TreeNode, Nothing}
end

TreeNode(rec::HyperRectangle{N,T}) where {N,T} = TreeNode(rec, T[], SVector{N,T}[], none, 0, nothing, nothing)

Base.show(io::IO, node::TreeNode) = print(io, "node: $(node.rec)")


"""
    displaytree(node::TreeNode)

Pretty print function for a `TreeNode` object
"""
function displaytree(node::TreeNode, depth::Int=0)
    println("  "^depth, depth > 0 ? "$depth|-- " : "",
            node.pruning == inside ? "Inside Domain : " : "", node.pruning == outside ? "Outside Domain : " : "",
            node.levelsetaxis > 0 ? "Levelset axis k = $(node.levelsetaxis), " : "", node.rec
    )
    node.childleft  !== nothing && displaytree(node.childleft, depth+1)
    node.childright !== nothing && displaytree(node.childright, depth+1)
    return
end


"""
    ispruned(node::TreeNode, ϕ::Function)

Do the pruning, in the hypperrectangle stored in `node`, with the isosurface `ϕ`
return true if the hyperrectangle is pruned and update `node.pruning`
"""
function ispruned(node::TreeNode, ϕ)
    δ = maxgap(ϕ, node.rec)
    ϕc = ϕ(center(node.rec))
    if abs(ϕc) >= δ
        node.pruning = ϕc <= 0 ? inside : outside
        return true
    end
    false
end

"""
    getlevelsetaxis(rec::HyperRectangle, ∇ϕ::Function)

Return the axis index of an existing level set for the hyperrectangle `rec` and `∇ϕ`,
the gradient of the isosurface.
return `0` if none exists
"""
function getlevelsetaxis(rec::HyperRectangle, ∇ϕ::Function)
    xc = center(rec)
    g = ∇ϕ(xc)
    δ = maxgap.(linearize(∇ϕ, rec))
    norm∇ϕ = sum((g .+ δ).^2)
    for k = 1:dim(rec)
        if abs(g[k]) > δ[k] && norm∇ϕ / (g[k] - δ[k])^2 < C
            return k
        end
    end
    0
end


"""
    maketree!(node::TreeNode, ϕ::Function, ∇ϕ::Function, nodes, weights, minrecsize, maxdepth::Int, depth::Int)

Recursively construct the binary tree for the computation of the quadrature nodes and weights.
"""
function maketree!(node::TreeNode, ϕ::Function, ∇ϕ::Function, nodes, weights, minrecsize, maxdepth::Int, depth::Int)
    if ispruned(node, ϕ)
        return node
    elseif maximum(δ(node.rec)) <= minrecsize
        k = getlevelsetaxis(node.rec, ∇ϕ)
        if k > 0
            node.levelsetaxis = k
            node.quadratureweights, node.quadraturenodes = levelset_quadrature2D(node.rec, ϕ, ∇ϕ, k, nodes, weights)
            return node
        end
    end
    if depth <= maxdepth
        recA, recB = split(node.rec)
        node.childleft  = maketree!(TreeNode(recA), ϕ, ∇ϕ, nodes, weights, minrecsize, maxdepth, depth+1)
        node.childright = maketree!(TreeNode(recB), ϕ, ∇ϕ, nodes, weights, minrecsize, maxdepth, depth+1)
    else
        println("maximum detph reached (depth = ", depth, ").")
    end
    node
end

function get_nodes_weights(node::TreeNode{N,T}) where {N,T}
    nodes = SVector{N,T}[]
    weights = T[]
    update_nodes_weights!(nodes, weights, node)
    nodes, weights
end

function update_nodes_weights!(n, w, node::TreeNode)
    if node.levelsetaxis > 0
        append!(w, node.quadratureweights)
        append!(n, node.quadraturenodes)
    elseif node.childleft !== nothing && node.childright !== nothing
        update_nodes_weights!(n, w, node.childleft)
        update_nodes_weights!(n, w, node.childright)
    end
end
