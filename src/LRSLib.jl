module LRSLib

using BinDeps
using Polyhedra

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("LRSLib not properly installed. Please run Pkg.build(\"LRSLib\")")
end

macro lrs_ccall(f, args...)
  quote
    ret = ccall(($"lrs_$f", liblrs), $(map(esc,args)...))
    ret
  end
end

if 0 == (@lrs_ccall init Clong (Ptr{UInt8},) C_NULL)
  error("Initialization of LRS failed")
end

include("lrstypes.jl")
include("matrix.jl")
include("description.jl")

end # module
