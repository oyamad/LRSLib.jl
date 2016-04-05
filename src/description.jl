import Base.convert, Polyhedra.HRepresentation, Polyhedra.VRepresentation

# LRSMatrix -> Representation
HRepresentation(matrix::LRSInequalityMatrix) = HRepresentation{Rational{BigInt}}(matrix)
VRepresentation(matrix::LRSGeneratorMatrix) = VRepresentation{Rational{BigInt}}(matrix)

# converters Representation -> LRSMatrix
function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::HRepresentation{Rational{BigInt}})
  if N != fulldim(ine)
    error("N should be equal to the number of columns of A")
  end
  M = [ine.b -ine.A]
  (P, Q) = initmatrix(M, ine.linset, true)
  LRSInequalityMatrix{N}(P, Q)
end

Base.convert{N}(::Type{LRSMatrix{N}}, ine::HRepresentation{Rational{BigInt}}) = Base.convert(LRSInequalityMatrix{N}, ine)

function settoCarray(set::IntSet, m::Integer)
  s = zeros(Rational{BigInt}, m)
  for el in set
    s[el] = Base.convert(Rational{BigInt}, 1)
  end
  s
end

function Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::VRepresentation{Rational{BigInt}})
  if N != fulldim(ext)
    error("N should be equal to the number of columns of V and R")
  end
  mA = [ext.V; ext.R]
  b = settoCarray(ext.vertex, size(mA, 1))
  linset = ext.Rlinset
  if !isempty(ext.Vlinset)
    linset = copy(linset)
    for i in ext.Vlinset
      push!(linset, size(ext.V, 1)+i)
    end
  end
  (P, Q) = initmatrix([b mA], linset, false)
  LRSGeneratorMatrix{N}(P, Q)
end

function Base.convert{N}(::Type{LRSMatrix{N}}, ext::VRepresentation{Rational{BigInt}})
  Base.convert(LRSGeneratorMatrix{N}, ext)
end

# Specified N
Base.convert{N, S}(::Type{LRSMatrix{N}}, desc::Representation{S}) = Base.convert(LRSMatrix{N}, Base.convert(Representation{Rational{BigInt}}, desc))

Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::HRepresentation) = Base.convert(LRSMatrix{N}, ine)
Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::VRepresentation) = Base.convert(LRSMatrix{N}, ext)
# Unspecified N
Base.convert{S}(::Type{LRSMatrix}, desc::Representation{S}) = Base.convert(LRSMatrix{fulldim(desc)}, desc)

Base.convert{T}(::Type{LRSInequalityMatrix}, ine::HRepresentation{T}) = Base.convert(LRSMatrix, ine)
Base.convert{T}(::Type{LRSGeneratorMatrix}, ext::VRepresentation{T}) = Base.convert(LRSMatrix, ext)


# converters LRSMatrix -> Representation
function extractAb(P::Clrs_dic, Q::Clrs_dat, offset::Int)
  m = P.m
  d = P.d-offset
  #d = Q.n-offset-1 # FIXME when it is modified...
  #@show P.m
  #@show P.d
  #@show Q.n
  b = Vector{Rational{BigInt}}(m)
  A = Matrix{Rational{BigInt}}(m, d)
  for i in 1:m
    gcd = extractbigintat(Q.Gcd, 1+i) # first row is the objective
    lcm = extractbigintat(Q.Lcm, 1+i)
    row = unsafe_load(P.A, 1+i)
    extractthisrow(i::Int) = (extractbigintat(row, offset+i) * gcd) // lcm
    b[i] = extractthisrow(1)
    for j in 1:d
      A[i, j] = extractthisrow(1+j)
    end
  end
  (b, A)
end

function extractlinset(Q::Clrs_dat)
  linset = IntSet([])
  for i in 1:Q.nlinearity
    push!(linset, unsafe_load(Q.linearity, i))
  end
  linset
end

function Base.convert{N}(::Type{HRepresentation{Rational{BigInt}}}, matrix::LRSInequalityMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 0

  linset = extractlinset(Q)
  (b, A) = extractAb(P, Q, 0)
  HRepresentation(-A, b, linset)
end

Base.convert{N}(::Type{Representation{Rational{BigInt}}}, ine::LRSInequalityMatrix{N}) = Base.convert(HRepresentation{Rational{BigInt}}, ine)

# I don't want it to overwrite Base.convert behaviour
function myconvert(::Type{IntSet}, a::Vector{Rational{BigInt}})
  b = Array{Bool}(a)
  s = IntSet()
  for i = 1:length(a)
    if b[i]
      push!(s, i)
    end
  end
  s
end

function Base.convert{N}(::Type{VRepresentation{Rational{BigInt}}}, matrix::LRSGeneratorMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 1

  linset = extractlinset(Q)
  (b, A) = extractAb(P, Q, 1)
  VRepresentation(A, myconvert(IntSet, b), linset)
end

Base.convert{N}(::Type{Representation{Rational{BigInt}}}, ine::LRSGeneratorMatrix{N}) = Base.convert(VRepresentation{Rational{BigInt}}, ine)

Base.convert{N, S}(::Type{HRepresentation{S}}, matrix::LRSMatrix{N}) = Base.convert(Representation{S}, Base.convert(Representation{Rational{BigInt}}, matrix))
Base.convert{N, S}(::Type{VRepresentation{S}}, matrix::LRSMatrix{N}) = Base.convert(Representation{S}, Base.convert(Representation{Rational{BigInt}}, matrix))
Base.convert{N, S}(::Type{Representation{S}}, matrix::LRSMatrix{N}) = Base.convert(Representation{S}, Base.convert(Representation{Rational{BigInt}}, matrix))
Base.convert{N}(::Type{Representation}, matrix::LRSMatrix{N}) = Base.convert(Representation{Rational{BigInt}}, matrix)
