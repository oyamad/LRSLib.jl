A = [1 0; 0 1; -1 0; 0 -1]
b = [1, 1, 0, 0]
ine = Polyhedra.HRepresentation(A, b)
inem1 = LRSMatrix(ine)
#setdebug(poly1, true)
ine1  = Polyhedra.Representation{Int}(inem1)
println(ine)
println(ine1)
