using ImplicitSurfaceQuadrature
using ForwardDiff
using Test

@testset "ForwardDiff and Linearization" begin
    f(x) = 2x[1]^2 + x[1]*x[2]^3 + 2.5 # f: ℝ² -> ℝ

    ∇f_fd(x) = ForwardDiff.gradient(f, x)    # ∇f with ForwardDiff
    ∇f_ex(x) = (4x[1] + x[2]^3, 3x[1]x[2]^2) # ∇f exact

    rec = HyperRectangle((0., 0.), (1., 1.))
    
    lin_fd = ImplicitSurfaceQuadrature.linearize(∇f_fd, rec)
    lin_ex = ImplicitSurfaceQuadrature.linearize(∇f_ex, rec)

    @test all(lin_fd .== lin_ex)
end
