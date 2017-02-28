using LRSLib
using Polyhedra
using Base.Test
using Clp
lpsolver = Clp.ClpSolver()

include("cube.jl")
include("simplex.jl")

include("polyhedron.jl")
