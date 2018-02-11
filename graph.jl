struct Graph
    n :: Int # |V|
    m :: Int # |E|
    nbr :: Array{Array{Int, 1}, 1} # neighbros of each vertex
end

# Graph -> Graph
function max_cc(G :: Graph)
    vis = fill(false, G.n)

    dfs(u :: Int) = begin
        vis[u] = true
        ret = [u]
        for v in G.nbr[u]
            if ! vis[v]
                append!(ret, dfs(v))
            end
        end
        return ret
    end

    max_cc = Vector{Int}(0)
    for u in 1 : G.n
        if ! vis[u]
            cc = dfs(u)
            if length(cc) > length(max_cc)
                max_cc = cc
            end
        end
    end

    id = Dict{Int,Int}( tuple.(max_cc, indices(max_cc,1)) )

    n = length(id)

    nbr = Array{Array{Int, 1}}(n)
    for u in filter( x -> haskey(id, x), indices(G.nbr, 1) )
        nbr[ id[u] ] = getindex.(id, filter( x -> haskey(id, x), G.nbr[u] ) )
    end

    m = div( sum( length.(nbr) ), 2 )

    return Graph(n, m, nbr)
end

# Graph -> Graph
function bfs_max_cc(G :: Graph)
    vis = fill(false, G.n)

    bfs(u :: Int) = begin
        vis[u] = true
        ret = [u]
        hinge = [u]
        while !isempty(hinge)
            v = pop!(hinge)
            for t in G.nbr[v]
                if !vis[t]
                    push!(ret, t)
                    push!(hinge,t)
                    vis[t] = true
                    #println(t)
                end
            end
        end
        return ret
    end
    max_cc = Vector{Int}(0)
    for u in 1 : G.n
        if ! vis[u]
            cc = bfs(u)
            if length(cc) > length(max_cc)
                max_cc = cc
            end
        end
    end

    id = Dict{Int,Int}( tuple.(max_cc, indices(max_cc,1)) )

    n = length(id)

    nbr = Array{Array{Int, 1}}(n)
    for u in filter( x -> haskey(id, x), indices(G.nbr, 1) )
        nbr[ id[u] ] = getindex.(id, filter( x -> haskey(id, x), G.nbr[u] ) )
    end

    m = div( sum( length.(nbr) ), 2 )

    return Graph(n, m, nbr)
end

# String -> Graph
function read_graph(str :: String)
    ints = parse.(Int, split(str))

    n = 0
    d = Dict{Int,Int}()
    edges = Set{Tuple{Int,Int}}()
    getid(x :: Int) = haskey(d, x) ? d[x] : d[x] = n += 1

    for i in 1 : 2 : length(ints)
        u = getid(ints[i])
        v = getid(ints[i + 1])
        if u == v continue end
        if u > v  u,v = v,u end
        push!(edges, (u,v))
    end
    m = length(edges)

    nbr = Array{Array{Int, 1}}(n)
    for i in indices(nbr,1) nbr[i] = [] end
    for (u,v) in edges
        push!(nbr[u], v)
        push!(nbr[v], u)
    end

    return Graph(n,m,nbr)
end

# String -> Graph
function read_file(filename :: String)
    return open(filename) do f
        #read_graph( readstring(f) )
        bfs_max_cc( read_graph( readstring(f) ) )
    end
end
