using LRSLib: RepMatrix, setobj, lpsolve
using LinearAlgebra


@testset "Tests for lp.jl" begin

    @testset "Simple test" begin
        a = [7//3, 1//3]
        b = 1//2
        c = [-1//1, 0]
        hr = hrep([HalfSpace(a, b),
                   HalfSpace([-1//1, 0], 0),
                   HalfSpace([0, -1//1], 0)])
        maximize = false

        status = :Optimal
        val = -3//14
        x = [3//14, 0//1]

        m = RepMatrix(hr)
        setobj(m, c, maximize)
        sol = lpsolve(m)

        @test sol.objval == val
        @test sol.sol == x
    end

    @testset "Tests from scipy/optimize/tests/test_linprog.py" begin

        @testset "test_linprog_unbounded" begin
            c = [1, 1]
            A = [-1 1; -1 -1; -1 0; 0 -1]
            b = [-1, 2, 0, 0]
            maximize = true

            status = :Unbounded

            hr = hrep(A, b)
            m = RepMatrix(hr)
            setobj(m, c, maximize)
            sol = lpsolve(m)

            @test sol.status == status
        end

        @testset "test_linprog_infeasible" begin
            c = [1, 1]
            A = [1 0; 0 1; -1 -1]
            b = [2, 2, -5]
            maximize = true

            status = :Infeasible

            hr = hrep(A, b)
            m = RepMatrix(hr)
            setobj(m, c, maximize)
            sol = lpsolve(m)

            @test sol.status == status
        end

        @testset "test_nontrivial_problem" begin
            c = [-1, 8, 4, -6]
            A_ub = [-7 -7 6 9;
                    1 -1 -3 0;
                    10 -10 -7 7;
                    6 -1 3 4]
            b_ub = [-3, 6, -6, 6]
            A_eq = [-10, 1, 1, -8]
            b_eq = -4
            n = length(c)
            nonneg = hrep(-Matrix{Int}(I, (n, n)), zeros(Int, n))
            maximize = false

            status = :Optimal
            desired_fun = 7083//1391
            desired_x = [101//1391, 1462//1391, 0, 752//1391]

            hr = intersect(hrep(A_ub, b_ub), HyperPlane(A_eq, b_eq), nonneg)
            m = RepMatrix(hr)
            setobj(m, c, maximize)
            sol = lpsolve(m)

            @test sol.status == status
            @test sol.objval == desired_fun
            @test sol.sol == desired_x
        end

    end  # Tests from scipy/optimize/tests/test_linprog.py
end
