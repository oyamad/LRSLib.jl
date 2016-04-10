export vertexenum
function vertexenum(m::LRSMatrix)
  # code from here is borrowed from lrs_main

  # Pivot to a starting dictionary
  lin = getfirstbasis(m)

  # There may have been column redundancy
  # If so the linearity space is obtained and redundant
  # columns are removed. User can access linearity space
  # from lin dimensions nredundcol x d+1

  if lin !== nothing
    println("LINEARITY")
    println(getmat(lin))
  end

  vertexenumend(m)
end

function vertexenumend{N}(m::LRSMatrix{N})
  # We initiate reverse search from this dictionary
  # getting new dictionaries until the search is complete
  # User can access each output line from output which is
  # vertex/ray/facet from the lrs_mp_vector output

  M = Matrix{Rational{BigInt}}(0, N+1)
  while true
    for col in 0:getd(m)
      output = getsolution(m, col)
      if output !== nothing
        M = [M; output']
      end
    end
    if !getnextbasis(m)
      break
    end
  end
  M
end
function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, m::LRSGeneratorMatrix{N})
  M = vertexenum(m)
  (P, Q) = initmatrix(M, IntSet([]), true)
  LRSInequalityMatrix{N}(P, Q)
end
LRSInequalityMatrix{N}(m::LRSGeneratorMatrix{N}) = LRSInequalityMatrix{N}(m)
function Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, m::LRSInequalityMatrix{N})
  M = vertexenum(m)
  (P, Q) = initmatrix(M, IntSet([]), false)
  LRSGeneratorMatrix{N}(P, Q)
end
LRSGeneratorMatrix{N}(m::LRSInequalityMatrix{N}) = LRSGeneratorMatrix{N}(m)
