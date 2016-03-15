export LRSMatrix, LRSInequalityMatrix, LRSGeneratorMatrix
abstract LRSMatrix{N}

type LRSInequalityMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
end
type LRSGeneratorMatrix{N} <: LRSMatrix{N}
  P::Ptr{Clrs_dic}
  Q::Ptr{Clrs_dat}
end

function initmatrix(M::Matrix{Rational{BigInt}}, linset, inequality::Bool)
  m = Clong(size(M, 1))
  n = Clong(size(M, 2))
  linsetarray = Vector{Clong}(length(linset))
  i = 1
  for j in linset
    linsetarray[i] = j
    i += 1
  end
  Q = @lrs_ccall alloc_dat Ptr{Clrs_dat} (Cstring,) "coucou"
  @lrs_ccall init_dat Void (Ptr{Clrs_dat}, Clong, Clong, Clong, Clong, Ptr{Clong}) Q m n Clong(inequality ? 0 : 1) length(linset) linsetarray
  P = @lrs_ccall alloc_dic Ptr{Clrs_dic} (Ptr{Clrs_dat},) Q
  N = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), M)
  D = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), M)
  @lrs_ccall init_dic_vec Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clrs_mp_vector, Clrs_mp_vector) P Q N D
  (P, Q)
end
