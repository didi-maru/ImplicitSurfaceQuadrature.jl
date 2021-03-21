using ImplicitSurfaceQuadrature
using Test

@testset "2D and 3D hyperrectangle" begin
    rec = HyperRectangle((1.,1.), (2.,4.))
    recA, recB = ImplicitSurfaceQuadrature.split(rec)
    @test recA === HyperRectangle((1.,1.),(2.,2.5))
    @test recB === HyperRectangle((1.,2.5),(2.,4.))

    rec = HyperRectangle((1.,1.,1.), (2.,4.,2.))
    recA, recB = ImplicitSurfaceQuadrature.split(rec)
    @test recA === HyperRectangle((1.,1.,1.),(2.,2.5,2.))
    @test recB === HyperRectangle((1.,2.5,1.),(2.,4.,2.))
end
