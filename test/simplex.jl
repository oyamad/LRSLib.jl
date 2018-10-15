@testset "Test representation conversion with the simplex" begin
    # A = [1 1; -1 0; 0 -1]
    # b = [1, 0, 0]
    # linset = IntSet([1])
    # V = [0 1; 1 0]

    A = [1 1; -1 0; 0 -1]
    b = [1, 0, 0]
    linset = IntSet()
    V = [0 0; 0 1; 1 0]

    function minitest(ine::LRSLib.HMatrix)
        ine  = convert(MixedMatHRep{Int}, ine)
        @test sortrows([ine.b -ine.A]) == sortslices([b -A], dims=1)
        @test ine.linset == linset
    end
    function minitest(ext::LRSLib.VMatrix)
        ext  = convert(MixedMatVRep{Int}, ext)
        @test sortrows(ext.V) == V
        @test length(ext.R) == 0
        @test ext.Rlinset == IntSet()
    end

    ine = hrep(A, b, linset)
    ext = vrep(V)

    inem1 = LRSLib.HMatrix("simplex.ine")
    minitest(inem1)
    extm1 = convert(LRSLib.VMatrix, inem1)
    minitest(extm1)

    inem2 = LRSLib.RepMatrix(ine)
    minitest(inem2)
    extm2 = convert(LRSLib.VMatrix, inem2)
    minitest(extm2)

    extm3 = LRSLib.RepMatrix(ext)
    minitest(extm3)
    inem3 = convert(LRSLib.HMatrix, extm3)
    minitest(inem3)
end
