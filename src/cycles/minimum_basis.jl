# Code in this file inspired by NetworkX.

function generate_test()
    G = Graph()
    add_vertices!(G, 5)
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 3)
    add_edge!(G, 3, 1)
    add_edge!(G, 3, 4)
    add_edge!(G, 4, 5)
    add_edge!(G, 5, 3)
    return G
end

"""Returns a minimum weight cycle basis for G

    Minimum weight means a cycle basis for which the total weight
    (length for unweighted graphs) of all the cycles is minimum.

    Parameters
    ----------
    G : NetworkX Graph

    Returns
    -------
    A list of cycle lists.  Each cycle list is a list of nodes
    which forms a cycle (loop) in G. Note that the nodes are not
    necessarily returned in a order by which they appear in the cycle

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_cycle(G, [0, 1, 2, 3])
    >>> nx.add_cycle(G, [0, 3, 4, 5])
    >>> print([sorted(c) for c in nx.minimum_cycle_basis(G)])
    [[0, 1, 2, 3], [0, 3, 4, 5]]

    References:
        [1] Kavitha, Telikepalli, et al. "An O(m^2n) Algorithm for
        Minimum Cycle Basis of Graphs."
        http://link.springer.com/article/10.1007/s00453-007-9064-z
        [2] de Pina, J. 1995. Applications of shortest path methods.
        Ph.D. thesis, University of Amsterdam, Netherlands

    See Also
    --------
    simple_cycles, cycle_basis
"""
function minimum_cycle_basis(G)
    # We first split the graph in connected subgraphs
    basis = []
    for c in connected_components(G)
        comp, vmap = induced_subgraph(G,c)
        cmin = _min_cycle_basis(comp, vmap)
        cmin_cycle = edges_to_cycles(cmin)
        if length(cmin_cycle) > 0
            push!(basis, cmin_cycle)
        end
    end
    return basis
end

function edges_to_cycles(cycle_vect)
    paths = []
    # For each basis
    for vect in cycle_vect
        # for each cycle

        for eset in vect
            nodes = Set()
            # Get the nodes
            for e in eset
                for u in e
                    push!(nodes, u)
                end
            end
            nodes = sort(collect(nodes))
            # Construct a graph
            g = Graph()
            add_vertices!(g, maximum(nodes))
            for e in eset
                t = Tuple(sort(collect(e)))
                add_edge!(g, t[1], t[2])
            end
            path = non_backtracking_randomwalk(g, nodes[1], length(eset); seed=0)
            push!(paths, path)
        end
    end

    return paths
end


function _min_cycle_basis(comp, vmap)
    cb = []
    # Get a spanning tree by running BFS
    component_edges = Set([Set([src(e), dst(e)]) for e in collect(edges(comp))])
    spanning_tree_edges = Set([Set([src(e), dst(e)]) for e in edges(bfs_tree(comp, 1))])
    edges_excl = [e for e in component_edges if !(e in spanning_tree_edges)]
    N = length(edges_excl)

    # We maintain a set of vectors orthogonal to sofar found cycles
    set_orth = [Set([edge]) for edge in edges_excl]
    for k in 1:N
        # kth cycle is "parallel" to kth vector in set_orth
        new_cycle = _min_cycle(comp, set_orth[k])
        push!(cb, Vector([Set(new_cycle)]))
        # now update set_orth so that k+1,k+2... th elements are
        # orthogonal to the newly found cycle, as per [p. 336, 1]
        base = set_orth[k]
        for (idx,orth) in enumerate(set_orth[k+1:length(set_orth)])
            if length(intersect(orth, new_cycle)) % 2 != 0
                set_orth[k+idx] = union(setdiff(orth,base), setdiff(base,orth))
            else
                set_orth[k+idx] = orth
            end
        end
        #set_orth[k + 1 :] = Vector([if (length(orth & new_cycle) % 2) (orth ⊻ base) else orth end for orth in set_orth[k + 1:]])
    end
    return cb
end

"""
Computes the minimum length cycle in G,
orthogonal to the vector orth as per [p. 338, 1]
"""
function _min_cycle(G, orth)
    T = Graph()

    nodes_idx = Dict(node=>idx for (idx, node) in enumerate(vertices(G)))
    idx_nodes = Dict(idx=>node for (node, idx) in nodes_idx)

    nnodes = length(nodes_idx)
    add_vertices!(T, nnodes*2)
    # Add 2 copies of each edge in G to T. If edge is in orth, add cross edge;
    # otherwise in-plane edge
    for e in collect(edges(G))
        u = src(e)
        v = dst(e)
        uidx, vidx = nodes_idx[u], nodes_idx[v]
        if Set([u, v]) in orth
            add_edge!(T, uidx, nnodes + vidx)
            add_edge!(T, nnodes + uidx, vidx)
        else
            add_edge!(T, uidx, vidx)
            add_edge!(T, nnodes + uidx, nnodes + vidx)
        end
    end

    all_shortest_pathlens = Dict()
    for node in vertices(T)
        sp = dijkstra_shortest_paths(T, node)
        all_shortest_pathlens[node] = sp.dists
    end
    #all_shortest_pathlens = dijkstra_shortest_paths(T, collect(vertices(T)))
    cross_paths_w_lens = Dict(n=>all_shortest_pathlens[n][nnodes + n] for n in 1:nnodes)
    for n in sort(collect(keys(cross_paths_w_lens)))
    end

    # Now compute shortest paths in T, which translates to cyles in G
    start = argmin([cross_paths_w_lens[n] for n in 1:nnodes])
    last = nnodes + start
    sp = dijkstra_shortest_paths(T, start)
    min_path = reconstruct_path(sp, last)
    # Now we obtain the actual path, re-map nodes in T to those in G
    min_path_nodes = [if node <= nnodes node else node - nnodes end for node in min_path]

    # Now remove the edges that occur two times
    mcycle_pruned = _path_to_cycle(min_path_nodes)
    return Set(Set((idx_nodes[u], idx_nodes[v])) for (u, v) in mcycle_pruned)
end

function reconstruct_path(sp, target)
    path = Vector{Int64}()
    curr_node = target
    push!(path, curr_node)
    while true
        curr_node = sp.parents[curr_node]
        push!(path, curr_node)
        if sp.parents[curr_node] == 0
            break
        end
    end
    return reverse!(path)
end

function _path_to_cycle(path)
    """
    Removes the edges from path that occur even number of times.
    Returns a set of edges
    """
    edge_counts = Dict(Tuple([path[idx], path[idx+1]])=>0 for idx in 1:length(path)-1)
    for idx in 1:length(path)-1
        edge = Tuple([path[idx], path[idx+1]])
        edge_counts[edge] += 1
    end
    edges = Set()
    for (edge, count) in edge_counts
        if count % 2 != 0
            push!(edges, edge)
        end
    end
    return edges
end
