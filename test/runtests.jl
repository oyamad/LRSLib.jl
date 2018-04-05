using LRSLib
using Polyhedra
using Base.Test
using Clp
lpsolver = tuple(Clp.ClpSolver())

include("cube.jl")
include("simplex.jl")

include("polyhedron.jl")
