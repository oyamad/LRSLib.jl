const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "alltests.jl"))
basicpolyhedrontests(LRSLibrary())
@testset "Polyhedra tests" begin
    runtests(LRSLibrary())
end
