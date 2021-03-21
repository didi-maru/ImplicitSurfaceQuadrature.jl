function levelset_quadrature2D(rec::HyperRectangle{N,T}, ϕ, ∇ϕ, k, glnodes, glweights) where {N,T}
    e1 = @SVector [one(T), zero(T)]; e2 = @SVector [zero(T), one(T)]
    k̃ = k == 1 ? 2 : 1
    ek = k == 1 ? e1 : e2
    ek̃ = k̃ == 1 ? e1 : e2

    ϕL(x̃) = ϕ(x̃*ek̃ + rec.lower[k]*ek)
    ϕU(x̃) = ϕ(x̃*ek̃ + rec.upper[k]*ek)

    # Integration segments (along the k\tilde direction)
    rootsL = find_zeros(ϕL, rec.lower[k̃], rec.upper[k̃])
    rootsU = find_zeros(ϕU, rec.lower[k̃], rec.upper[k̃])
    roots = sort(vcat(rec.lower[k̃], rec.upper[k̃], rootsL, rootsU))

    s = sign(∇ϕ(center(rec))[k]) # sign of ∂ϕ/∂xᵢ in the rectangle
    
    allweights = T[]
    allnodes = SVector{N,T}[]
    for j = 1:length(roots)-1 # Integration via quadrature on each segments

        xc = (roots[j] + roots[j+1])/2 # Current segment center

        if s*ϕL(xc) <= 0 && s*ϕU(xc) >= 0

            m = (roots[j+1] - roots[j])/2 # Current segment length/2
            x̃ = m*glnodes .+ xc # Quadrature nodes along the k\tilde direction

            for i = 1:length(glnodes)
                if sign(ϕL(x̃[i])) != sign(ϕU(x̃[i]))
                    # Get level set value x = h(x\tilde)
                    levelsetroots = find_zeros(y -> ϕ(x̃[i]*ek̃ + y*ek), rec.lower[k], rec.upper[k])
                    r = length(levelsetroots) == 1 ? levelsetroots[1] : error("more than 1 zero")
                    
                    x = x̃[i]*ek̃ + r*ek # i^th quadrature node
                    append!(allnodes, [x])

                    g = ∇ϕ(x) # ∇ϕ evaluation on the quadrature point
                    append!(allweights, glweights[i] * m * norm(g)/abs(g[k]))
                end
            end
        end
    end
    allweights, allnodes
end
