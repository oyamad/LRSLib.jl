function redund(m::LRSMatrix)
  # if non-negative flag is set, non-negative constraints are not input
  # explicitly, and are not checked for redundancy

  # Step 2: Find a starting cobasis from default of specified order
  #         Lin is created if necessary to hold linearity space

  lin = getfirstbasis(m)

  # Pivot to a starting dictionary
  # There may have been column redundancy
  # If so the linearity space is obtained and redundant
  # columns are removed. User can access linearity space
  # from lrs_mp_matrix Lin dimensions nredundcol x d+1

  redundend(m)
end

function redundendi(m::LRSMatrix, ineq::Int)
  # TODO ineq -> index
  :redundend == checkindex(m, index)
end

function redundend(m::LRSMatrix)
  # Step 3: Test each row of the dictionary to see if it is redundant

  # note some of these may have been changed in getting initial dictionary
  m_A = unsafe_load(P).m_A
  d = unsafe_load(P).d
  lastdv = unsafe_load(Q).lastdv
  linset = extractlinset(m)

  # rows 0..lastdv are cost, decision variables, or linearities
  # other rows need to be tested

  redset = IntSet([])
  for index in (lastdv + 1):(m_A + d)
      ineq = unsafe_load(unsafe_load(Q).inequality, index - lastdv + 1) # the input inequality number corr. to this index

      status = checkindex(m, index)
      if status == :linearity && !(ineq in linset)
        error("please report this")
      elseif redundendi(m, index)
        push!(redset, ineq)
      end
  end
  redset
end
