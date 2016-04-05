A = [1 0; 0 1; -1 0; 0 -1]
b = [1, 1, 0, 0]
ine = Polyhedra.HRepresentation(A, b)
poly1 = LRSMatrix(ine)
#setdebug(poly1, true)
ineout1  = Polyhedra.Representation{Int}(poly1)
println(ine)
println(ineout1)
vertexenum(poly1)
