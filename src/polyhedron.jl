export LRSLibrary
importall Polyhedra

type LRSLibrary <: PolyhedraLibrary
end

type LRSPolyhedron{N} <: Polyhedron{N, Rational{BigInt}}
  ine::Nullable{HRepresentation{Rational{BigInt}}}
  inem::Nullable{LRSInequalityMatrix{N}}
  ext::Nullable{VRepresentation{Rational{BigInt}}}
  extm::Nullable{LRSGeneratorMatrix{N}}
  hlinearitydetected::Bool
  vlinearitydetected::Bool
  noredundantinequality::Bool
  noredundantgenerator::Bool

  function LRSPolyhedron(ine::HRepresentation{Rational{BigInt}}, ext::VRepresentation{Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool)
    if fulldim(ine) != fulldim(ext)
      error("dimension does not match")
    end
    new{fulldim(ine)}(ine, nothing, ext, nothing, hld, vld, nri, nrg)
  end
  function LRSPolyhedron(ine::HRepresentation{Rational{BigInt}})
    new{fulldim(ine)}(ine, nothing, nothing, nothing, false, false, false, false)
  end
  function LRSPolyhedron(ext::VRepresentation{Rational{BigInt}})
    new{fulldim(ine)}(nothing, nothing, ext, nothing, false, false, false, false)
  end
end

call{N, T}(::Type{LRSPolyhedron{N}}, repr::Representation{T}) = LRSPolyhedron{N}(Representation{Rational{BigInt}}(repr))

# Helpers
function getine(p::LRSPolyhedron)
  if isnull(p.ine)
    p.ine = Representation(vertexenumend(getextm(p)))
    # Now p.extm is no more just after the getfirstbasis so it is no more valid (I think)
    p.extm = nothing
  end
  get(p.ine)
end
function getinem(p::LRSPolyhedron)
  if isnull(p.inem)
    ine = getine(p)
    m = LRSMatrix(ine)
    getfirstbasis(m)
    p.inem = m
  end
  get(p.inem)
end
function getext(p::LRSPolyhedron)
  if isnull(p.ext)
    p.ext = Representation(vertexenumend(getinem(p)))
    # Now p.inem is no more just after the getfirstbasis so it is no more valid (I think)
    p.inem = nothing
  end
  get(p.ext)
end
function getextm(p::LRSPolyhedron)
  if isnull(p.extm)
    ext = getext(p)
    m = LRSMatrix(ext)
    getfirstbasis(m)
    p.extm = m
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
# TODO check dim of HRepresentation
function updateine!{N}(p::LRSPolyhedron{N}, ine::HRepresentation)
  if N != fulldim(ine)
    error("dimension does not match")
  end
  clearfield!(p)
  p.ine = ine
end
function updateext!{N}(p::LRSPolyhedron{N}, ext::LRSGeneratorMatrix)
  if N != fulldim(ext)
    error("dimension does not match")
  end
  clearfield!(p)
  p.ext = ext
end


# Implementation of Polyhedron's mandatory interface
polyhedron(repr::Representation, ::LRSLibrary) = LRSPolyhedron{fulldim(repr)}(repr)
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
  updateine!(intersect(getine(p), ine))
end
function Base.push!(p::LRSPolyhedron, ext::VRepresentation)
  updateext!(getext(p) + ext)
end
function inequalitiesarecomputed(p::LRSPolyhedron)
  !isnull(p.ine)
end
function getinequalities(p::LRSPolyhedron)
  copy(getine(p.ine))
end
function generatorsarecomputed(p::LRSPolyhedron)
  !isnull(p.ext)
end
function getgenerators(p::LRSPolyhedron)
  copy(getext(p))
end
#eliminate(p::Polyhedron, delset::IntSet)                     = error("not implemented")
function detectlinearities!(p::LRSPolyhedron)
  ine = getine(p)
  inem = getinem(p)
  ine.linset = extractlinset(inem)
end
function removeredundantinequalities!(p::LRSPolyhedron)
  ine = getine(p)
  inem = getinem(p)
  redset = redundend(inem)
  redvec = collect(redset)
  ine.A = ine.A[redvec,:]
  ine.b = ine.b[redvec]
  # TODO Update linset (put them at beginning)
end
function removeredundantgenerators!(p::LRSPolyhedron)
  ext = getext(p)
  extm = getextm(p)
  redset = redundend(extm)
  redvec = collect(redset)
  ext.V = ext.V[redvec,:] # TODO if ext is split...
  # TODO Update linset (put them at beginning)
end
function isredundantinequality(p::LRSPolyhedron, i::Integer)
  redundendi(getinem(p), i)
end
function isredundantgenerator(p::LRSPolyhedron, i::Integer)
  redundendi(getextm(p), i)
end
