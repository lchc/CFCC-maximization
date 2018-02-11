include("logw.jl")
include("graph.jl")
include("lapl.jl")
include("Linvdiag.jl")
include("Lpartinv.jl")
include("Lpartinv2.jl")
include("Lpartinv3.jl")
include("appxInvTrace.jl")

using Laplacians

function delnode2(L:: SparseMatrixCSC{Float64}, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function approx(G, k, w = open("log.txt", "w")) #v- the vertex chosed; G- graph; k- the number of edge added

    logw(w,"")
    logw(w,"running approx greedy...")
    n = G.n

    ans = zeros(3,1)
    t = time()
    A, L = sparseAdja(G)

    logw(w,"k=1 starts")
    start_time=time()

    t=time()
    Ldele = LinvdiagSS(A;JLfac=20)
    logw(w,"LaplSolve total time: ",time()-t, " (s)")

    u = indmin(Ldele)
    L = delnode2(L,u,n)
    A = delnode2(A,u,n)

    elapsed_time = time()-start_time

    start_time = time()
    for dep = 1:k-1
        logw(w,"k=",dep+1," starts")

        t=time()
        f = approxCholSddm(L,tol=1e-5, verbose=false);
        logw(w,"SDDM preconditioning time: ",time()-t, " (s)")
        t=time()
        Ldele = LpartinvSS(n - dep, f;JLfac=20)

        logw(w,"SDDM1 solve time: ",time()-t, " (s)")
        t=time()
        Ldele2 = Lpartinv2SS(n - dep, f, A;JLfac=20)
        logw(w,"SDDM2 solve time: ",time()-t, " (s)")
        t=time()
        Ldele3 = Lpartinv3SS(L, f;JLfac=20)
        logw(w,"SDDM3 solve time: ",time()-t, " (s)")

        Ldele ./= (Ldele2 .+ Ldele3)

        u = indmax(Ldele)

        L = delnode2(L,u,n-dep)
        A = delnode2(A,u,n-dep)
    end

    elapsed_time += time()-start_time
    logw(w,"approx_elapsed_time = ",elapsed_time, "s")
    ans[1] = elapsed_time
    logw(w,"calculating CFCC of the solution returned by approx greedy..")
    ans[3] = (n > 30000) ? (n/appxInvTrace(L;JLfac=200)) : ( n / trace( inv( full(L) ) ) )
    logw(w,"CFCC achieved by approx greedy: ", ans[3])
    logw(w,"")
    return ans
end
