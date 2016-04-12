import Base.convert, Polyhedra.HRepresentation, Polyhedra.VRepresentation

# LRSMatrix -> Representation
HRepresentation{N}(matrix::LRSInequalityMatrix{N}) = LiftedHRepresentation{N,Rational{BigInt}}(matrix)
VRepresentation{N}(matrix::LRSGeneratorMatrix{N}) = LiftedVRepresentation{N,Rational{BigInt}}(matrix)

# converters Representation -> LRSMatrix
function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::LiftedHRepresentation{N,Rational{BigInt}})
  (P, Q) = initmatrix(ine.A, ine.linset, true)
  m = LRSInequalityMatrix{N}(P, Q)
  #setdebug(m, true)
  #debugA(m)
  m
end
function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::HRepresentation{N})
  LRSInequalityMatrix{N}(Base.convert(LiftedHRepresentation{N,Rational{BigInt}}, ine))
end

# FIXME
Base.convert{N}(::Type{LRSMatrix{N}}, ine::HRepresentation{N, Rational{BigInt}}) = Base.convert(LRSInequalityMatrix{N}, ine)
#Base.convert{N}(::Type{LRSMatrix{N}}, ine::LiftedHRepresentation{N, Rational{BigInt}}) = LRSInequalityMatrix{N}(ine)
#Base.convert{N}(::Type{LRSMatrix{N}}, ine::SimpleHRepresentation{N, Rational{BigInt}}) = LRSInequalityMatrix{N}(ine)
#Base.convert{N}(::Type{LRSMatrix{N}}, ine::HRepresentation{N, Rational{BigInt}}) = LRSInequalityMatrix{N}(ine)
#call{N}(::Type{LRSMatrix{N}}, ine::HRepresentation{N, Rational{BigInt}}) = LRSInequalityMatrix{N}(ine)
#Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::HRepresentation{N, Rational{BigInt}}) = Base.convert(LRSInequalityMatrix{N}, ine)

function Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::LiftedVRepresentation{N, Rational{BigInt}})
  (P, Q) = initmatrix(ext.R, ext.linset, false)
  LRSGeneratorMatrix{N}(P, Q)
end
Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::VRepresentation{N}) = LRSGeneratorMatrix{N}(Base.convert(LiftedVRepresentation{N,Rational{BigInt}}, ext))

Base.convert{N}(::Type{LRSMatrix{N}}, ext::VRepresentation{N, Rational{BigInt}}) = Base.convert(LRSGeneratorMatrix{N}, ext)

# Specified N
Base.convert{N, S}(::Type{LRSMatrix{N}}, desc::Representation{N,S}) = Base.convert(LRSMatrix{N}, Base.convert(Representation{N,Rational{BigInt}}, desc))

Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::HRepresentation) = Base.convert(LRSMatrix{N}, ine)
Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::VRepresentation) = Base.convert(LRSMatrix{N}, ext)
# Unspecified N
Base.convert{N,S}(::Type{LRSMatrix}, desc::Representation{N,S}) = Base.convert(LRSMatrix{N}, desc)

Base.convert{N,T}(::Type{LRSInequalityMatrix}, ine::HRepresentation{N,T}) = Base.convert(LRSMatrix, ine)
Base.convert{N,T}(::Type{LRSGeneratorMatrix}, ext::VRepresentation{N,T}) = Base.convert(LRSMatrix, ext)


# converters LRSMatrix -> Representation

#FIXME could do BigInt actually
function Base.convert{N}(::Type{LiftedHRepresentation{N, Rational{BigInt}}}, matrix::LRSInequalityMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 0

  linset = extractinputlinset(Q)
  A = extractA(P, Q, 0)
  LiftedHRepresentation{N,Rational{BigInt}}(A, linset)
end

function Base.convert{N}(::Type{Representation{N, Rational{BigInt}}}, ine::LRSInequalityMatrix{N})
  Base.convert(LiftedHRepresentation{N, Rational{BigInt}}, ine)
end

function Base.convert{N}(::Type{LiftedVRepresentation{N, Rational{BigInt}}}, matrix::LRSGeneratorMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 1

  linset = extractinputlinset(Q)
  A = extractA(P, Q, 1)
  LiftedVRepresentation{N, Rational{BigInt}}(A, linset)
end

Base.convert{N}(::Type{Representation{N, Rational{BigInt}}}, ine::LRSGeneratorMatrix{N}) = Base.convert(LiftedVRepresentation{N, Rational{BigInt}}, ine)

Base.convert{N, S}(::Type{SimpleHRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(SimpleHRepresentation{N,S}, HRepresentation{N,S}(matrix))
Base.convert{N, S}(::Type{LiftedHRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(LiftedHRepresentation{N,S}, HRepresentation{N,S}(matrix))
Base.convert{N, S}(::Type{SimpleVRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(SimpleVRepresentation{N,S}, VRepresentation{N,S}(matrix))
Base.convert{N, S}(::Type{LiftedVRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(LiftedVRepresentation{N,S}, VRepresentation{N,S}(matrix))
Base.convert{N, S}(::Type{HRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(HRepresentation{N,S}, HRepresentation{N,Rational{BigInt}}(matrix))
Base.convert{N, S}(::Type{VRepresentation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(VRepresentation{N,S}, VRepresentation{N,Rational{BigInt}}(matrix))
Base.convert{N, S}(::Type{Representation{N,S}}, matrix::LRSMatrix{N}) = Base.convert(Representation{N,S}, Base.convert(Representation{N,Rational{BigInt}}, matrix))
Representation{N}(matrix::LRSMatrix{N}) = Base.convert(Representation{N, Rational{BigInt}}, matrix)
