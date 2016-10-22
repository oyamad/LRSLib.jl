# A = [1 1; -1 0; 0 -1]
# b = [1, 0, 0]
# linset = IntSet([1])
# V = [0 1; 1 0]

A = [1 1; -1 0; 0 -1]
b = [1, 0, 0]
linset = IntSet([])
V = [0 0; 0 1; 1 0]

function minitest(ine::LRSInequalityMatrix)
    ine  = SimpleHRepresentation{2,Int}(ine)
    @test sortrows([ine.b -ine.A]) == sortrows([b -A])
    @test ine.linset == linset
end
function minitest(ext::LRSGeneratorMatrix)
    ext  = SimpleVRepresentation{2,Int}(ext)
    @test sortrows(ext.V) == V
    @test length(ext.R) == 0
    @test ext.Vlinset == IntSet()
    @test ext.Rlinset == IntSet()
end

ine = Polyhedra.SimpleHRepresentation(A, b, linset)
ext = Polyhedra.SimpleVRepresentation(V)

inem1 = LRSMatrix(ine)
minitest(inem1)
extm1 = LRSGeneratorMatrix(inem1)
minitest(extm1)

extm2 = LRSMatrix(ext)
minitest(extm2)
inem2 = LRSInequalityMatrix(extm2)
minitest(inem2)
