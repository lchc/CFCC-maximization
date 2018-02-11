include("graph.jl")

function adja(G :: Graph)
    A:: Matrix{Float64} = zeros(G.n,G.n) 
    for i in indices(G.nbr,1)
        A[CartesianIndex.( tuple.(i, G.nbr[i]) )] .=  1.0
    end
    #A[CartesianIndex.( tuple.(1 : G.n, 1 : G.n) )] .-= vec(sum(L,2))
    return A
end

function sparseAdja(G :: Graph)
    A :: SparseMatrixCSC{Float64} = spzeros(G.n,G.n)
    
    I, J, V :: Array{Float64, 1}, D :: Array{Float64, 1} = [], [], [], []
    #edges :: Array{Tuple{Int,Int}, 1} = []
    for i in indices(G.nbr,1)
        push!(D,size(G.nbr[i],1))
        append!(I,fill(i, size(G.nbr[i],1)))
        append!(J,G.nbr[i])
        append!(V,fill(1, size(G.nbr[i],1)))
        #append!(edges, tuple.(i, G.nbr[i]))
        #A[CartesianIndex.( tuple.(i, G.nbr[i]) )] .=  1.0
    end
    return sparse(I,J,V), sparse([I; 1:G.n], [J; 1:G.n], [-V; D])
    #A[CartesianIndex.(edges)] .= 1.0
    #return A
end


function sparseAdja2(G :: Graph)
    A :: SparseMatrixCSC{Float64} = spzeros(G.n,G.n)
    for i in indices(G.nbr,1)
        t = time()
        for j = 1: size(G.nbr[i],1)
            A[i,G.nbr[i][j]] .= 1.0
        end
        println("Adja2: ",i, " ", time()-t)
        #A[CartesianIndex.( tuple.(i, G.nbr[i]) )] .=  1.0
    end
    return A
end


function lapl(G :: Graph)
    L :: Matrix{Float64}= zeros(G.n,G.n) 
    for i in indices(G.nbr,1)
        L[CartesianIndex.( tuple.(i, G.nbr[i]) )] .=  -1.0
    end
    L[CartesianIndex.( tuple.(1 : G.n, 1 : G.n) )] .-= vec(sum(L,2))
    return L
end

function sparseLapl(G :: Graph)
    L :: SparseMatrixCSC{Float64} = spzeros(G.n,G.n)
    tot :: Int = 0
    for i in indices(G.nbr,1)
        tot += size(G.nbr[i],1)
        L[CartesianIndex.( tuple.(i, G.nbr[i]) )] .=  -1.0
    end
    println("tot = ",tot)
    L[CartesianIndex.( tuple.(1 : G.n, 1 : G.n) )] .=- vec(sum(L,2))
    return L
end

function mppinv(L :: Matrix{Float64})
    n = size(L,2)
    let L = copy(L)
        L .-= 1.0 / n
        L .= inv(L)
        L .+= 1.0 / n
        return L
    end
end

# Graph -> Dict{ (Int,Int), Float64 }
function er_mppinv(G :: Graph)
    H = mppinv( lapl(G) )
    er = Dict{ Tuple{Int,Int}, Float64 }()
    for u in 1 : G.n
        for v in G.nbr[u]
            if u < v
                er[ (u,v) ] = H[u,u] + H[v,v] - 2 * H[u,v]
            end
        end
    end
    return er
end
