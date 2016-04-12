A = [1 0; 0 1; -1 0; 0 -1]
b = [1, 1, 0, 0]
ine = Polyhedra.SimpleHRepresentation(A, b)
inem1 = LRSMatrix(ine)
#setdebug(poly1, true)
ine1  = Polyhedra.Representation{2,Int}(inem1)
println(ine)
println(ine1)
