export LRSLibrary

type LRSLibrary <: PolyhedraLibrary
end

type LRSPolyhedron{N} <: Polyhedron{N, Rational{BigInt}}
  ine::Nullable{HRepresentation{N, Rational{BigInt}}}
  inem::Nullable{LRSInequalityMatrix{N}}
  ext::Nullable{VRepresentation{N, Rational{BigInt}}}
  extm::Nullable{LRSGeneratorMatrix{N}}
  hlinearitydetected::Bool
  vlinearitydetected::Bool
  noredundantinequality::Bool
  noredundantgenerator::Bool

  function LRSPolyhedron(ine::HRepresentation{N, Rational{BigInt}}, ext::VRepresentation{N, Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool)
    new(ine, nothing, ext, nothing, hld, vld, nri, nrg)
  end
  function LRSPolyhedron(ine::HRepresentation{N, Rational{BigInt}})
    new(ine, nothing, nothing, nothing, false, false, false, false)
  end
  function LRSPolyhedron(ext::VRepresentation{N, Rational{BigInt}})
    new(nothing, nothing, ext, nothing, false, false, false, false)
  end
end

eltype{N}(::Type{LRSPolyhedron{N}}) = Rational{BigInt}
eltype(::LRSPolyhedron) = Rational{BigInt}

(::Type{LRSPolyhedron}){N}(ine::HRepresentation{N, Rational{BigInt}}) = LRSPolyhedron{N}(HRepresentation{N, Rational{BigInt}}(ine))
(::Type{LRSPolyhedron}){N}(ext::VRepresentation{N, Rational{BigInt}}) = LRSPolyhedron{N}(VRepresentation{N, Rational{BigInt}}(ext))
(::Type{LRSPolyhedron}){N, T}(rep::Representation{N, T}) = LRSPolyhedron{N}(Representation{N, Rational{BigInt}}(rep))

# Helpers
function getine(p::LRSPolyhedron)
  if isnull(p.ine)
    if !isnull(p.inem) && checkfreshness(get(p.inem), :Fresh)
      p.ine = p.inem
    else
      p.ine = getextm(p, :Fresh)
      p.hlinearitydetected = true
      p.noredundantinequality = true
    end
  end
  get(p.ine)
end
function getinem(p::LRSPolyhedron, fresh::Symbol=:AnyFreshness)
  if isnull(p.inem) || !checkfreshness(get(p.inem), fresh)
    p.inem = LRSMatrix(getine(p))
  end
  get(p.inem)
end
function getext(p::LRSPolyhedron)
  if isnull(p.ext)
    if !isnull(p.extm) && checkfreshness(get(p.extm), :Fresh)
      p.ext = p.extm
    else
      p.ext = getinem(p, :Fresh)
      p.vlinearitydetected = true
      p.noredundantgenerator = true
    end
  end
  get(p.ext)
end
function getextm(p::LRSPolyhedron, fresh::Symbol=:AnyFreshness)
  if isnull(p.extm) || !checkfreshness(get(p.extm), fresh)
    p.extm = LRSMatrix(getext(p))
  end
  get(p.extm)
end

function clearfield!(p::LRSPolyhedron)
  p.ine = nothing
  p.inem = nothing
  p.ext = nothing
  p.extm = nothing
  hlinearitydetected = false
  vlinearitydetected = false
  noredundantinequality = false
  noredundantgenerator = false
end
function updateine!{N}(p::LRSPolyhedron{N}, ine::HRepresentation{N, Rational{BigInt}})
  clearfield!(p)
  p.ine = ine
end
function updateext!{N}(p::LRSPolyhedron{N}, ext::VRepresentation{N, Rational{BigInt}})
  clearfield!(p)
  p.ext = ext
end


# Implementation of Polyhedron's mandatory interface
polyhedron{N}(repr::Representation{N}, ::LRSLibrary) = LRSPolyhedron{N}(repr)

getlibraryfor{T<:Union{Int,Rational}}(p::LRSPolyhedron, ::Type{T}) = LRSLibrary()

function Base.copy{N}(p::LRSPolyhedron{N})
  ine = nothing
  if !isnull(p.ine)
    ine = copy(get(p.ine))
  end
  ext = nothing
  if !isnull(p.ext)
    ext = copy(get(p.ext))
  end
  LRSPolyhedron{N}(ine, ext, p.hlinearitydetected, p.vlinearitydetected, p.noredundantinequality, p.noredundantgenerator)
end
function Base.push!(p::LRSPolyhedron, ine::HRepresentation)
  updateine!(p, intersect(getine(p), ine))
end
function Base.push!(p::LRSPolyhedron, ext::VRepresentation)
  updateext!(p, getext(p) + ext)
end
function inequalitiesarecomputed(p::LRSPolyhedron)
  !isnull(p.ine)
end
function getinequalities(p::LRSPolyhedron)
  copy(getine(p))
end
function generatorsarecomputed(p::LRSPolyhedron)
  !isnull(p.ext)
end
function getgenerators(p::LRSPolyhedron)
  copy(getext(p))
end
#eliminate(p::Polyhedron, delset::IntSet)                     = error("not implemented")
function detecthlinearities!(p::LRSPolyhedron)
  if !p.hlinearitydetected
    getext(p)
    p.inem = nothing
    p.ine = nothing
    getine(p)
    # getine sets hlinearities as detected and no redundant ineq.
  end
end
function detectvlinearities!(p::LRSPolyhedron)
  if !p.vlinearitydetected
    getine(p)
    p.extm = nothing
    p.ext = nothing
    getext(p)
    # getext sets vlinearities as detected and no redundant gen.
  end
end
function removeredundantinequalities!(p::LRSPolyhedron)
  #if !p.noredundantinequality
    ine = getine(p)
    inem = getinem(p, :AlmostFresh) # FIXME does it need to be fresh ?
    linset = getinputlinsubset(inem)
    redset = redund(inem)
    nonred = setdiff(IntSet(1:size(ine.A, 1)), redset)
    nonred = collect(setdiff(nonred, linset))
    lin = collect(linset)
    ine.A = [ine.A[lin,:]; ine.A[nonred,:]]
    ine.linset = IntSet(1:length(linset))
    p.noredundantinequality = true
  #end
end
function removeredundantgenerators!(p::LRSPolyhedron)
  if !p.noredundantgenerator
    detectvlinearities!(p)
    ext = getext(p)
    extm = getextm(p, :AlmostFresh) # FIXME does it need to be fresh ?
    redset = redund(extm)
    nonred = setdiff(IntSet(1:size(ext.R, 1)), redset)
    nonred = collect(setdiff(nonred, ext.linset))
    lin = collect(ext.linset)
    ext.R = [ext.R[lin,:]; ext.R[nonred,:]]
    ext.linset = IntSet(1:length(ext.linset))
    p.noredundantgenerator = true
  end
end
#function getredundantinequalities(p::LRSPolyhedron)
#  redund(getinem(p, :AlmostFresh))
#end
function isredundantinequality(p::LRSPolyhedron, i::Integer)
  redundi(getinem(p, :AlmostFresh), i) # FIXME does it need to be fresh ?
end
function isredundantgenerator(p::LRSPolyhedron, i::Integer)
  redundi(getextm(p, :AlmostFresh), i) # FIXME does it need to be fresh ?
end
# Optional interface
function Polyhedra.loadpolyhedron!(p::LRSPolyhedron, filename::AbstractString, ::Type{Val{:ext}})
  clearfield!(p)
  p.extm = LRSGeneratorMatrix(string(filename, ".ext"))
end
