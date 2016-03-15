import Base.convert, Polyhedra.InequalityDescription, Polyhedra.GeneratorDescription

# LRSMatrix -> Description
InequalityDescription(matrix::LRSInequalityMatrix) = InequalityDescription{Rational{BigInt}}(matrix)
GeneratorDescription(matrix::LRSGeneratorMatrix) = GeneratorDescription{Rational{BigInt}}(matrix)

# converters Description -> LRSMatrix
function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::InequalityDescription{Rational{BigInt}})
  if N != fulldim(ine)
    error("N should be equal to the number of columns of A")
  end
  M = [ine.b -ine.A]
  (P, Q) = initmatrix(M, ine.linset, true)
  LRSInequalityMatrix{N}(P, Q)
end

Base.convert{N}(::Type{LRSMatrix{N}}, ine::InequalityDescription{Rational{BigInt}}) = Base.convert(LRSInequalityMatrix{N}, ine)

function settoCarray(set::IntSet, m::Integer)
  s = zeros(Rational{BigInt}, m)
  for el in set
    s[el] = Base.convert(Rational{BigInt}, 1)
  end
  s
end

function Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::GeneratorDescription{Rational{BigInt}})
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

function Base.convert{N}(::Type{LRSMatrix{N}}, ext::GeneratorDescription{Rational{BigInt}})
  Base.convert(LRSGeneratorMatrix{N}, ext)
end

# Specified N
Base.convert{N, S}(::Type{LRSMatrix{N}}, desc::Description{S}) = Base.convert(LRSMatrix{N}, Base.convert(Description{Rational{BigInt}}, desc))

Base.convert{N}(::Type{LRSInequalityMatrix{N}}, ine::InequalityDescription) = Base.convert(LRSMatrix{N}, ine)
Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, ext::GeneratorDescription) = Base.convert(LRSMatrix{N}, ext)
# Unspecified N
Base.convert{S}(::Type{LRSMatrix}, desc::Description{S}) = Base.convert(LRSMatrix{fulldim(desc)}, desc)

Base.convert{T}(::Type{LRSInequalityMatrix}, ine::InequalityDescription{T}) = Base.convert(LRSMatrix, ine)
Base.convert{T}(::Type{LRSGeneratorMatrix}, ext::GeneratorDescription{T}) = Base.convert(LRSMatrix, ext)


# converters LRSMatrix -> Description
function extractbigintat(array::Ptr{GMPInteger}, i::Int)
  tmp = BigInt(0)
  ccall((:__gmpz_set, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}), pointer_from_objref(tmp), array + (i-1) * sizeof(GMPInteger))
  tmp
end

function extractAb(P::Clrs_dic, Q::Clrs_dat, offset::Int)
  m = P.m
  d = P.d-offset
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

function Base.convert{N}(::Type{InequalityDescription{Rational{BigInt}}}, matrix::LRSInequalityMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 0

  linset = extractlinset(Q)
  (b, A) = extractAb(P, Q, 0)
  InequalityDescription(-A, b, linset)
end

Base.convert{N}(::Type{Description{Rational{BigInt}}}, ine::LRSInequalityMatrix{N}) = Base.convert(InequalityDescription{Rational{BigInt}}, ine)

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

function Base.convert{N}(::Type{GeneratorDescription{Rational{BigInt}}}, matrix::LRSGeneratorMatrix{N})
  P = unsafe_load(matrix.P)
  Q = unsafe_load(matrix.Q)
  @assert Q.hull == 1

  linset = extractlinset(Q)
  (b, A) = extractAb(P, Q, 1)
  GeneratorDescription(A, myconvert(IntSet, b), linset)
end

Base.convert{N}(::Type{Description{Rational{BigInt}}}, ine::LRSGeneratorMatrix{N}) = Base.convert(GeneratorDescription{Rational{BigInt}}, ine)

Base.convert{N, S}(::Type{InequalityDescription{S}}, matrix::LRSMatrix{N}) = Base.convert(Description{S}, Base.convert(Description{Rational{BigInt}}, matrix))
Base.convert{N, S}(::Type{GeneratorDescription{S}}, matrix::LRSMatrix{N}) = Base.convert(Description{S}, Base.convert(Description{Rational{BigInt}}, matrix))
Base.convert{N, S}(::Type{Description{S}}, matrix::LRSMatrix{N}) = Base.convert(Description{S}, Base.convert(Description{Rational{BigInt}}, matrix))
Base.convert{N}(::Type{Description}, matrix::LRSMatrix{N}) = Base.convert(Description{Rational{BigInt}}, matrix)
