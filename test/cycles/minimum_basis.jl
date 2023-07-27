@testset "Minimum Cycle Basis" begin
    function evaluate(x, y)
        x_sorted = sort(sort.(x))
        y_sorted = sort(sort.(y))
        @test x_sorted == y_sorted
    end

    # Graph with one cycle
    elist = [(1, 2), (2, 3), (3, 4), (4, 1), (1, 5)]
    ex = Graph(SimpleEdge.(elist))
    expected_cyclebasis = Array{Int64,1}[[1, 2, 3, 4]]
    @testset "one cycle" for g in test_generic_graphs(ex)
        ex_cyclebasis = minimum_cycle_basis(g)
        evaluate(ex_cyclebasis, expected_cyclebasis)
    end

    # test_dimensionality(self):
    # checks |MCB|=|E|-|V|+|NC|
    @testset "dimensionality" begin
        ntrial = 10
        for _ in 1:ntrial
            rg = erdos_renyi(10, 0.3)
            for g in test_generic_graphs(rg) 
                nnodes = nv(g)
                nedges = ne(g)
                ncomp = length(connected_components(g))
                dim_mcb = sum([length(basis) for basis in minimum_cycle_basis(g)])
                @test dim_mcb == nedges - nnodes + ncomp
            end
        end
    end

    # Tests that in a complete graph all cycles are of length 3
    cg = complete_graph(5)
    @testset "clique" for g in test_generic_graphs(cg)
        mcb = minimum_cycle_basis(g)
        @test sum([length(cycle) for comp in mcb for cycle in comp]) == sum(length(basis) for basis in (mcb))*3
    end

    # Tests that a tree has no cycles and thus no cycle basis
    tg = binary_tree(6)
    @testset "tree" for g in test_generic_graphs(tg)
        @test length(minimum_cycle_basis(g)) == 0
    end
end
