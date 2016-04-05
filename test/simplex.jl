A = [1 1; -1 0; 0 -1]
b = [1, 0, 0]
linset = IntSet([1])
V = [0 1; 1 0]
vertex = IntSet([1,2])

ine = Polyhedra.HRepresentation(A, b, linset)
ext = Polyhedra.VRepresentation(V, vertex)
poly1 = LRSMatrix(ine)
poly2 = LRSGeneratorMatrix(poly1)
ineout1  = Polyhedra.Representation{Int}(poly1)
extout1  = Polyhedra.Representation{Int}(poly2)
println(ine)
println(ineout1)
println(ext)
println(extout1)

# poly2 = LRSMatrix(ext)
# extout2  = Polyhedra.Representation{Int}(poly2)
# println(ext)
# println(extout2)
# vertexenum(poly2)
