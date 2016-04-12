# A = [1 1; -1 0; 0 -1]
# b = [1, 0, 0]
# linset = IntSet([1])
# V = [0 1; 1 0]
# vertex = IntSet([1,2])

A = [1 1; -1 0; 0 -1]
b = [1, 0, 0]
linset = IntSet([])
V = [0 0; 0 1; 1 0]
vertex = IntSet([1,2,3])


ine = Polyhedra.SimpleHRepresentation(A, b, linset)
ext = Polyhedra.SimpleVRepresentation(V)

inem1 = LRSMatrix(ine)
setdebug(inem1, true)
extm1 = LRSGeneratorMatrix(inem1)
ine1  = Polyhedra.Representation{2,Int}(inem1)
ext1  = Polyhedra.Representation{2,Int}(extm1)
println(ine)
println(ine1)
println(ext)
println(ext1)

extm2 = LRSMatrix(ext)
setdebug(extm2, true)
println(LRSGeneratorMatrix(extm2))
inem2 = LRSInequalityMatrix(extm2)
ine2  = Polyhedra.Representation{2,Int}(inem2)
ext2  = Polyhedra.Representation{2,Int}(extm2)
println(ine)
println(ine2)
println(ext)
println(ext2)
