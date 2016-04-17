export LRSMatrix, LRSInequalityMatrix, LRSGeneratorMatrix, setdebug, debugA
abstract LRSMatrix{N}

function checkfreshness(m::LRSMatrix, fresh::Symbol)
  fresh == :AnyFreshNess ||
  (fresh == :Fresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty]) ||
  (fresh == :AlmostFresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty, :RedundancyChecked])
end

type LRSLinearitySpace{N}
  Lin::Clrs_mp_matrix
  nlin::Int
  n::Int
  hull::Bool
  homogeneous::Bool

  function LRSLinearitySpace(Lin::Clrs_mp_matrix, nlin, n, hull, homogeneous)
    m = new(Lin, nlin, n, hull, homogeneous)
    finalizer(m, myfree)
    m
  end
end

function lrs_alloc_dat()
  @lrs_ccall alloc_dat Ptr{Clrs_dat} (Ptr{Cchar},) C_NULL
end

function lrs_alloc_dic(Q::Ptr{Clrs_dat})
  @lrs_ccall alloc_dic Ptr{Clrs_dic} (Ptr{Clrs_dat},) Q
end

function initmatrix(filename::AbstractString)
  Q = lrs_alloc_dat()
  # The first element does not matter
  argv = ["", filename]
  ok = Clrs_true == @lrs_ccall read_dat Clong (Ptr{Clrs_dat}, Cint, Ptr{Ptr{Cchar}}) Q length(argv) argv
  if !ok
    error("Invalid file $filename")
  end
  P = lrs_alloc_dic(Q)
  ok = Clrs_true == @lrs_ccall read_dic Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}) P Q
  if !ok
    error("Invalid file $filename")
  end
  (P,Q)
end

function initmatrix(M::Matrix{Rational{BigInt}}, linset, Hrep::Bool)
  m = Clong(size(M, 1))
  n = Clong(size(M, 2))
  Q = lrs_alloc_dat()
  @lrs_ccall init_dat Void (Ptr{Clrs_dat}, Clong, Clong, Clong) Q m n Clong(Hrep ? 0 : 1)
  P = lrs_alloc_dic(Q)
  #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
  for i in 1:m
#   num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), M[i,:])
#   den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), M[i,:])
    ineq = !(i in linset)
    setrow(P, Q, i, M[i,:], !(i in linset))
#   @lrs_ccall set_row_mp Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
  end
  # This is the objective. If I have no objective LRS might fail
  if Hrep
    setrow(P, Q, 0, ones(Rational{BigInt}, n), true)
  end
  (P, Q)
end


function myfree(l::LRSLinearitySpace)
  if l.nlin > 0
    @lrs_ccall clear_mp_matrix Void (Clrs_mp_matrix, Clong, Clong) l.Lin l.nlin l.n
  end
end

function myfree(m::LRSMatrix)
  @lrs_ccall free_dic_and_dat Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}) m.P m.Q
end

type LRSInequalityMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
  status::Symbol
  lin::Nullable{LRSLinearitySpace{N}}
  function LRSInequalityMatrix(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat})
    m = new(P, Q, :AtNoBasis, nothing)
    finalizer(m, myfree)
    m
  end
end

function LRSInequalityMatrix(filename::AbstractString)
  P, Q = initmatrix(filename)
  LRSInequalityMatrix{unsafe_load(P).d}(P, Q)
end

type LRSGeneratorMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
  status::Symbol
  lin::Nullable{LRSLinearitySpace{N}}
  function LRSGeneratorMatrix(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat})
    m = new(P, Q, :AtNoBasis, nothing)
    finalizer(m, myfree)
    m
  end
end

function LRSGeneratorMatrix(filename::AbstractString)
  P, Q = initmatrix(filename)
  LRSGeneratorMatrix{unsafe_load(m.P).d-1}(P, Q)
end

#I should also remove linearity (should I remove one if hull && homogeneous ?)
#getd{N}(m::LRSInequalityMatrix{N}) = N
#getd{N}(m::LRSGeneratorMatrix{N}) = N+1
#Let's do it the easy way
getd{N}(m::LRSInequalityMatrix{N}) = unsafe_load(m.P).d
getd{N}(m::LRSGeneratorMatrix{N}) = unsafe_load(m.P).d

function setrow(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, i::Int, row::Vector{Rational{BigInt}}, ineq::Bool)
  num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), row)
  den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), row)
  @lrs_ccall set_row_mp Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
end

function setdebug(m::LRSMatrix, debug::Bool)
  @lrs_ccall setdebug Void (Ptr{Clrs_dat}, Clong) m.Q (debug ? Clrs_true : Clrs_false)
end

function debugA(m::LRSMatrix)
  P = unsafe_load(m.P)
  Q = unsafe_load(m.Q)
  m = P.m
  d = P.d
  for i in 1:m+1
    row = unsafe_load(P.A, i)
    for j in 1:d+1
      x = extractbigintat(row, j)
      print(" $x")
    end
    println()
  end
end

function extractA(P::Clrs_dic, Q::Clrs_dat, offset::Int)
  m = P.m
  d = P.d-offset
  #d = Q.n-offset-1 # FIXME when it is modified...
  A = Matrix{Rational{BigInt}}(m, d+1)
  for i in 1:m
    gcd = extractbigintat(Q.Gcd, 1+i) # first row is the objective
    lcm = extractbigintat(Q.Lcm, 1+i)
    row = unsafe_load(P.A, 1+i)
    extractthisrow(i::Int) = (extractbigintat(row, offset+i) * gcd) // lcm
    for j in 1:d+1
      A[i, j] = extractthisrow(j)
    end
  end
  A
end

function extractinputlinset(Q::Clrs_dat)
  linset = IntSet([])
  for i in 1:Q.nlinearity
    push!(linset, unsafe_load(Q.linearity, i))
  end
  linset
end

function extractoutputlinset(Q::Clrs_dat)
  k = (Q.hull == Clrs_true && Q.homogeneous == Clrs_true) ? 1 : 0
  nredundcol = Q.nredundcol
  IntSet(1:(nredundcol-k))
end

# FIXME The only think that is done is that the linearities given that were redundant have been
# removed from the linearity set so that they can be marked as redundant inequalities.
# New linearities are detected but getinputlinsubset does not give them.
# I should check in redundcols
function getinputlinsubset(m::LRSMatrix)
  if m.status == :AtNoBasis
    getfirstbasis(m)
  end
  extractinputlinset(unsafe_load(m.Q))
end
function getoutputlinset(m::LRSMatrix)
  if m.status == :AtNoBasis
    getfirstbasis(m)
  end
  extractoutputlinset(unsafe_load(m.Q))
end


function convertoutput(x::Clrs_mp_vector, n, hull)
  first = extractbigintat(x, 1)
  rest = Vector{BigInt}(n-1)
  for i = 2:n
    rest[i-1] = extractbigintat(x, i)
  end
  if hull || first == 0
    Rational{BigInt}[first; rest]
  else
    [one(Rational{BigInt}); rest // first]
  end
end

function getmat{N}(lin::LRSLinearitySpace{N})
  startcol = lin.hull && lin.homogeneous ? 2 : 1 # col zero not treated as redundant
  A = Matrix{BigInt}(lin.nlin-startcol+1, lin.n)
  for col in startcol:lin.nlin # print linearity space */
    A[col-startcol+1,:] = convertoutput(unsafe_load(lin.Lin, col), lin.n, lin.hull)
  end
  A
end

function getfirstbasis{N}(m::LRSMatrix{N})
  Lin = Ref{Clrs_mp_matrix}(C_NULL)
  Pptr = Ref{Ptr{Clrs_dic}}(m.P)
  # The "Clrs_true" at the last argument since that it should not be verbose
  found = Clrs_true == (@lrs_ccall getfirstbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Ptr{Clrs_mp_matrix}, Clong) Pptr m.Q Lin Clrs_true)
  m.P = Pptr[]
  if !found
    # Note that I can have a basis found while the polyhedron is empty
    m.status = :Empty
    # FIXME in that case does redundancy checking with getindex still works ?
  else
    m.status = :AtFirstBasis
    # FIXME does this linearity also works if the first basis is not found ?
    #       I could say that there are linearities which are x_1 = 0 and x_1 = 1
    Q = unsafe_load(m.Q)
    if Q.nredundcol > 0
      # There may have been column redundancy
      # If so the linearity space is obtained and redundant
      # columns are removed. User can access linearity space
      # from lin dimensions nredundcol x d+1

      m.lin = LRSLinearitySpace{N}(Lin[], Q.nredundcol, Q.n, Q.hull == Clrs_true, Q.homogeneous == Clrs_true)
    end
  end
end

function getnextbasis(m::LRSMatrix)
  Pptr = Ref{Ptr{Clrs_dic}}(m.P)
  x = (@lrs_ccall getnextbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Clong) Pptr m.Q Clrs_false)
  found = Clrs_true == x
  m.P = Pptr[]
  m.status = :AtSomeBasis
  found
end

function getsolution(m::LRSMatrix, col::Int)
  Q = unsafe_load(m.Q)
  output = @lrs_ccall alloc_mp_vector Clrs_mp_vector (Clong,) Q.n
  found = Clrs_true == (@lrs_ccall getsolution Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clrs_mp_vector, Clong) m.P m.Q output col)
  if found
    out = convertoutput(output, Q.n, Q.hull == Clrs_true)
  else
    out = nothing
  end
  @lrs_ccall clear_mp_vector Void (Clrs_mp_vector, Clong) output Q.n
  out
end

function checkindex(m::LRSMatrix, index::Int)
  if m.status == :AtNoBasis
    getfirstbasis(m)
  end
  # FIXME if it is at some basis or last basis, does this still works ?
  ret = @lrs_ccall2 checkindex Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong) m.P m.Q index
  m.status = :RedundancyChecked
  [:nonredundant, :redundant, :linearity][ret+1]
end
