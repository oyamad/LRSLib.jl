using LRSLib: liblrs

@testset "Tests for initialization" begin

    @testset "Test for lrs_init" begin
        @test unsafe_load(cglobal((:lrs_ofp, liblrs), Ptr{Cvoid})) != C_NULL
    end

end
