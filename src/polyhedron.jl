export LRSLibrary

mutable struct LRSLibrary <: PolyhedraLibrary
end

mutable struct LRSPolyhedron{N} <: Polyhedron{N, Rational{BigInt}}
    ine::Nullable{HRepresentation{N, Rational{BigInt}}}
    inem::Nullable{LRSInequalityMatrix{N}}
    ext::Nullable{VRepresentation{N, Rational{BigInt}}}
    extm::Nullable{LRSGeneratorMatrix{N}}
    hlinearitydetected::Bool
    vlinearitydetected::Bool
    noredundantinequality::Bool
    noredundantgenerator::Bool

    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}, ext::VRepresentation{N, Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool) where {N}
        new{N}(ine, nothing, ext, nothing, hld, vld, nri, nrg)
    end
    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}, ::Void, hld::Bool, vld::Bool, nri::Bool, nrg::Bool) where {N}
        new{N}(ine, nothing, nothing, nothing, hld, vld, nri, nrg)
    end
    function LRSPolyhedron{N}(::Void, ext::VRepresentation{N, Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool) where {N}
        new{N}(nothing, nothing, ext, nothing, hld, vld, nri, nrg)
    end
    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}) where {N}
        new{N}(ine, nothing, nothing, nothing, false, false, false, false)
    end
    function LRSPolyhedron{N}(ext::VRepresentation{N, Rational{BigInt}}) where {N}
        new{N}(nothing, nothing, ext, nothing, false, false, false, false)
    end
end

# ine may decompose fast but if ine is nothing I do not want to ask to compute it to see the type it is
# saying false normally do not give troubles
decomposedhfast{N}(::Type{LRSPolyhedron{N}}) = false
decomposedvfast{N}(::Type{LRSPolyhedron{N}}) = false
decomposedhfast(p::LRSPolyhedron{N}) where {N} = decomposedhfast(LRSPolyhedron{N})
decomposedvfast(p::LRSPolyhedron{N}) where {N} = decomposedvfast(LRSPolyhedron{N})

eltype{N}(::Type{LRSPolyhedron{N}}) = Rational{BigInt}
eltype(::LRSPolyhedron) = Rational{BigInt}

# Helpers
function getine(p::LRSPolyhedron)
    if isnull(p.ine)
        if !isnull(p.inem) && checkfreshness(get(p.inem), :Fresh)
            p.ine = p.inem
        else
            p.ine = LiftedHRepresentation(getextm(p, :Fresh))
            p.inem = nothing
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
            p.ext = LiftedVRepresentation(getinem(p, :Fresh))
            p.extm = nothing
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
function updateine!(p::LRSPolyhedron{N}, ine::HRepresentation{N, Rational{BigInt}}) where N
    clearfield!(p)
    p.ine = ine
end
function updateext!(p::LRSPolyhedron{N}, ext::VRepresentation{N, Rational{BigInt}}) where N
    clearfield!(p)
    p.ext = ext
end


# Implementation of Polyhedron's mandatory interface
polyhedron(repit::Union{Representation{N},HRepIterator{N},VRepIterator{N}}, ::LRSLibrary) where {N} = LRSPolyhedron{N}(repit)

getlibraryfor(p::LRSPolyhedron, n::Int, ::Type{T}) where {T<:Union{Integer,Rational}} = LRSLibrary()
Polyhedra.changefulldim{N}(::Type{LRSPolyhedron{N}}, n::Int)= LRSPolyhedron{n}

LRSPolyhedron{N}(it::HRepIterator{N,T}) where {N, T} = LRSPolyhedron{N}(LRSInequalityMatrix{N}(it))
LRSPolyhedron{N}(it::VRepIterator{N,T}) where {N, T} = LRSPolyhedron{N}(LRSGeneratorMatrix{N}(it))

function LRSPolyhedron{N}(eqs::EqIterator, ineqs::IneqIterator) where N
    LRSPolyhedron{N}(LRSInequalityMatrix{N}(eqs, ineqs))
end
function LRSPolyhedron{N}(points::PointIterator, rays::RayIterator) where N
    LRSPolyhedron{N}(LRSGeneratorMatrix{N}(points, rays))
end

function Base.copy(p::LRSPolyhedron{N}) where N
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
function Base.push!(p::LRSPolyhedron{N}, ine::HRepresentation{N}) where N
    updateine!(p, intersect(getine(p), changeeltype(ine, Rational{BigInt})))
end
function Base.push!(p::LRSPolyhedron{N}, ext::VRepresentation{N}) where N
    updateext!(p, convexhull(getext(p), changeeltype(ext, Rational{BigInt})))
end
function hrepiscomputed(p::LRSPolyhedron)
    !isnull(p.ine)
end
function hrep(p::LRSPolyhedron)
    copy(getine(p))
end
function vrepiscomputed(p::LRSPolyhedron)
    !isnull(p.ext)
end
function vrep(p::LRSPolyhedron)
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
function removehredundancy!(p::LRSPolyhedron)
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
function removevredundancy!(p::LRSPolyhedron)
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
function ishredundant(p::LRSPolyhedron, i::Integer; strongly=false, cert=false, solver=Polyhedra.defaultLPsolverfor(p))
    @assert !strongly && !cert
    redundi(getinem(p, :AlmostFresh), i) # FIXME does it need to be fresh ?
end
function isvredundant(p::LRSPolyhedron, i::Integer; strongly=false, cert=false, solver=Polyhedra.defaultLPsolverfor(p))
    @assert !strongly && !cert
    redundi(getextm(p, :AlmostFresh), i) # FIXME does it need to be fresh ?
end
# Optional interface
function Polyhedra.loadpolyhedron!(p::LRSPolyhedron, filename::AbstractString, ::Type{Val{:ext}})
    clearfield!(p)
    p.extm = LRSGeneratorMatrix(string(filename, ".ext"))
end

for f in [:hashreps, :nhreps, :starthrep, :hasineqs, :nineqs, :startineq, :haseqs, :neqs, :starteq]
    @eval $f(p::LRSPolyhedron) = $f(getine(p))
end
for f in [:donehrep, :nexthrep, :doneineq, :nextineq, :doneeq, :nexteq]
    @eval $f(p::LRSPolyhedron, state) = $f(getine(p), state)
end

for f in [:hasvreps, :nvreps, :startvrep, :haspoints, :npoints, :startpoint, :hasrays, :nrays, :startray]
    @eval $f(p::LRSPolyhedron) = $f(getext(p))
end
for f in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval $f(p::LRSPolyhedron, state) = $f(getext(p), state)
end
