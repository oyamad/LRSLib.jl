facts("Test representation conversion with the simplex") do
    # A = [1 1; -1 0; 0 -1]
    # b = [1, 0, 0]
    # linset = IntSet([1])
    # V = [0 1; 1 0]

    A = [1 1; -1 0; 0 -1]
    b = [1, 0, 0]
    linset = IntSet()
    V = [0 0; 0 1; 1 0]

    function minitest(ine::LRSInequalityMatrix)
        ine  = SimpleHRepresentation{2,Int}(ine)
        @fact sortrows([ine.b -ine.A]) --> sortrows([b -A])
        @fact ine.linset --> linset
    end
    function minitest(ext::LRSGeneratorMatrix)
        ext  = SimpleVRepresentation{2,Int}(ext)
        @fact sortrows(ext.V) --> V
        @fact length(ext.R) --> 0
        @fact ext.Vlinset --> IntSet()
        @fact ext.Rlinset --> IntSet()
    end

    ine = Polyhedra.SimpleHRepresentation(A, b, linset)
    ext = Polyhedra.SimpleVRepresentation(V)

    inem1 = LRSInequalityMatrix("simplex.ine")
    minitest(inem1)
    extm1 = LRSGeneratorMatrix(inem1)
    minitest(extm1)

    inem2 = LRSMatrix(ine)
    minitest(inem2)
    extm2 = LRSGeneratorMatrix(inem2)
    minitest(extm2)

    extm3 = LRSMatrix(ext)
    minitest(extm3)
    inem3 = LRSInequalityMatrix(extm3)
    minitest(inem3)
end
