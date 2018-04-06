const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))
@testset "Polyhedra tests" begin
    polyhedratest(LRSLibrary(), ["empty"])
end
