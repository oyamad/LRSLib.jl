export LRSLibrary
import MathProgBase
const MPB = MathProgBase
import JuMP

struct LRSLibrary <: PolyhedraLibrary
    solver::MPB.AbstractMathProgSolver
    function LRSLibrary(solver=JuMP.UnsetSolver())
        new(solver)
    end
end
Polyhedra.similar_library(::LRSLibrary, ::FullDim, ::Type{T}) where T<:Union{Integer,Rational} = LRSLibrary()
Polyhedra.similar_library(::LRSLibrary, d::FullDim, ::Type{T}) where T = Polyhedra.default_library(d, T)

mutable struct LRSPolyhedron{N} <: Polyhedron{N, Rational{BigInt}}
    ine::Nullable{HRepresentation{N, Rational{BigInt}}}
    inem::Nullable{LRSInequalityMatrix{N}}
    ext::Nullable{VRepresentation{N, Rational{BigInt}}}
    extm::Nullable{LRSGeneratorMatrix{N}}
    hlinearitydetected::Bool
    vlinearitydetected::Bool
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver::MPB.AbstractMathProgSolver

    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}, ext::VRepresentation{N, Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver::MPB.AbstractMathProgSolver) where N
        new{N}(ine, nothing, ext, nothing, hld, vld, nri, nrg, solver)
    end
    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}, ::Void, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver::MPB.AbstractMathProgSolver) where N
        new{N}(ine, nothing, nothing, nothing, hld, vld, nri, nrg, solver)
    end
    function LRSPolyhedron{N}(::Void, ext::VRepresentation{N, Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver::MPB.AbstractMathProgSolver) where N
        new{N}(nothing, nothing, ext, nothing, hld, vld, nri, nrg, solver)
    end
    function LRSPolyhedron{N}(ine::HRepresentation{N, Rational{BigInt}}, solver::MPB.AbstractMathProgSolver) where N
        new{N}(ine, nothing, nothing, nothing, false, false, false, false, solver)
    end
    function LRSPolyhedron{N}(ext::VRepresentation{N, Rational{BigInt}}, solver::MPB.AbstractMathProgSolver) where N
        new{N}(nothing, nothing, ext, nothing, false, false, false, false, solver)
    end
end
LRSPolyhedron{N}(h::HRepresentation{N}, solver::MPB.AbstractMathProgSolver) where N = LRSPolyhedron{N}(HRepresentation{N, Rational{BigInt}}(h), solver)
LRSPolyhedron{N}(v::VRepresentation{N}, solver::MPB.AbstractMathProgSolver) where N = LRSPolyhedron{N}(VRepresentation{N, Rational{BigInt}}(v), solver)

Polyhedra.library(::LRSPolyhedron) = LRSLibrary()
Polyhedra.default_solver(p::LRSPolyhedron) = p.solver
Polyhedra.supportssolver(::Type{<:LRSPolyhedron}) = true

Polyhedra.arraytype(::Union{LRSPolyhedron, Type{<:LRSPolyhedron}}) = Vector{Rational{BigInt}}
Polyhedra.similar_type(::Type{<:LRSPolyhedron}, ::FullDim{N}, ::Type{Rational{BigInt}}) where N = LRSPolyhedron{N}
Polyhedra.similar_type(::Type{<:LRSPolyhedron}, d::FullDim, ::Type{T}) where T = Polyhedra.default_type(d, T)

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
polyhedron(rep::Representation{N}, lib::LRSLibrary) where N = LRSPolyhedron{N}(rep, lib.solver)

function LRSPolyhedron{N}(hits::Polyhedra.HIt{N}...; solver=JuMP.UnsetSolver()) where N
    LRSPolyhedron{N}(LRSInequalityMatrix{N}(hits...), solver)
end
function LRSPolyhedron{N}(vits::Polyhedra.VIt{N}...; solver=JuMP.UnsetSolver()) where N
    LRSPolyhedron{N}(LRSGeneratorMatrix{N}(vits...), solver)
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
    LRSPolyhedron{N}(ine, ext, p.hlinearitydetected, p.vlinearitydetected, p.noredundantinequality, p.noredundantgenerator, p.solver)
end
function Base.intersect!(p::LRSPolyhedron{N}, ine::HRepresentation{N}) where N
    updateine!(p, intersect(getine(p), HRepresentation{N, Rational{BigInt}}(ine)))
end
function Polyhedra.convexhull!(p::LRSPolyhedron{N}, ext::VRepresentation{N}) where N
    updateext!(p, convexhull(getext(p), VRepresentation{N, Rational{BigInt}}(ext)))
end
function hrepiscomputed(p::LRSPolyhedron)
    !isnull(p.ine)
end
function hrep(p::LRSPolyhedron)
    getine(p)
end
function vrepiscomputed(p::LRSPolyhedron)
    !isnull(p.ext)
end
function vrep(p::LRSPolyhedron)
    getext(p)
end
#eliminate(p::Polyhedron, delset::IntSet)                     = error("not implemented")
function detecthlinearity!(p::LRSPolyhedron)
    if !p.hlinearitydetected
        getext(p)
        p.inem = nothing
        p.ine = nothing
        getine(p)
        # getine sets hlinearity as detected and no redundant ineq.
    end
end
function detectvlinearity!(p::LRSPolyhedron)
    if !p.vlinearitydetected
        getine(p)
        p.extm = nothing
        p.ext = nothing
        getext(p)
        # getext sets vlinearity as detected and no redundant gen.
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
        detectvlinearity!(p)
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
_getrepfor(p::LRSPolyhedron, ::Polyhedra.HIndex, status::Symbol) = getinem(p, status)
_getrepfor(p::LRSPolyhedron, ::Polyhedra.VIndex, status::Symbol) = getextm(p, status)
function isredundant(p::LRSPolyhedron, idx::Polyhedra.Index; strongly=false, cert=false, solver=Polyhedra.solver(p))
    @assert !strongly && !cert
    redundi(_getrepfor(p, idx, :AlmostFresh), idx.value) # FIXME does it need to be fresh ?
end
# Optional interface
function Polyhedra.loadpolyhedron!(p::LRSPolyhedron, filename::AbstractString, ::Type{Val{:ext}})
    clearfield!(p)
    p.extm = LRSGeneratorMatrix(string(filename, ".ext"))
end
