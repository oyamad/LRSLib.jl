using GLPKMathProgInterface
const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))
lpsolver = tuple()
@testset "Polyhedra tests" begin
    polyhedratest(LRSLibrary(GLPKSolverLP()), ["empty", "cubedecompose", "largedecompose"])
end
