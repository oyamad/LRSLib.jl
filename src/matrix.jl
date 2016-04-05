export LRSMatrix, LRSInequalityMatrix, LRSGeneratorMatrix, setdebug
abstract LRSMatrix{N}

# TODO do this in this order
# lrs_free_dic (P,Q);           # deallocate lrs_dic
# lrs_free_dat (Q);             # deallocate lrs_dat


type LRSInequalityMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
end
type LRSGeneratorMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
end

getd{N}(m::LRSInequalityMatrix{N}) = N
getd{N}(m::LRSGeneratorMatrix{N}) = N+1

function initmatrix(M::Matrix{Rational{BigInt}}, linset, Hrep::Bool)
  m = Clong(size(M, 1))
  n = Clong(size(M, 2))
  # TODO not needed to pass linset like that
  linsetarray = Vector{Clong}(length(linset))
  i = 1
  for j in linset
    linsetarray[i] = j
    i += 1
  end
  Q = @lrs_ccall alloc_dat Ptr{Clrs_dat} (Ptr{Cchar},) C_NULL
  @lrs_ccall init_dat Void (Ptr{Clrs_dat}, Clong, Clong, Clong) Q m n Clong(Hrep ? 0 : 1)
  P = @lrs_ccall alloc_dic Ptr{Clrs_dic} (Ptr{Clrs_dat},) Q
  #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
  for i in 1:m
    num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), M[i,:])
    den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), M[i,:])
    ineq = !(i in linset)
    @lrs_ccall set_row_mp Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
  end
  (P, Q)
end

type LRSLinearitySpace{N}
  Lin::Clrs_mp_matrix
  nlin::Int
  n::Int
  hull::Bool

  function LRSLinearitySpace(Lin::Clrs_mp_matrix, nlin, n, hull)
    m = new(Lin, nlin, n, hull)
    finalizer(m, myfree)
  end
end

function setdebug(m::LRSMatrix, debug::Bool)
  @lrs_ccall setdebug Void (Ptr{Clrs_dat}, Clong) m.Q (debug ? Clrs_true : Clrs_false)
end

function myfree(l::LRSLinearitySpace)
  if l.nlin > 0
    @lrs_ccall clear_mp_matrix Void (Clrs_mp_matrix, Clong, Clong) l.lin l.nlin l.n
  end
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
  startcol = lin.hull ? 2 : 1 # col zero not treated as redundant
  A = Matrix{BigInt}(lin.nlin-startcol+1, lin.n)
  for col in startcol:lin.nlin # print linearity space */
    A[col-startcol+1,:] = convertoutput(unsafe_load(lin.Lin, col), lin.n, lin.hull)
  end
end

function getfirstbasis{N}(m::LRSMatrix{N})
  Lin = Ref{Clrs_mp_matrix}(C_NULL)
  Pptr = Ref{Ptr{Clrs_dic}}(m.P)
  # Clrs_true means no output
  found = Clrs_true == (@lrs_ccall getfirstbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Ptr{Clrs_mp_matrix}, Clong) Pptr m.Q Lin Clrs_true)
  if !found
    error("getfirstbasis failed")
  end
  m.P = Pptr[]
  Q = unsafe_load(m.Q)
  if Q.nredundcol > 0
    LRSLinearitySpace{N}(Lin[], Q.nredundcol, Q.n, Q.hull == Clrs_true)
  else
    nothing
  end
end

function getnextbasis(m::LRSMatrix)
  Pptr = Ref{Ptr{Clrs_dic}}(m.P)
  found = Clrs_true == (@lrs_ccall getnextbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Clong) Pptr m.Q Clrs_false)
  m.P = Pptr[]
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
