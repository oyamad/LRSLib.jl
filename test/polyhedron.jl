const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "alltests.jl"))
@testset "Polyhedra tests" begin
    simplextest(LRSLibrary())
    #runtests(LRSLibrary())
end
