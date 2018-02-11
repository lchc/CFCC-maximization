include("graph.jl")
include("lapl.jl")
include("logw.jl")

function delnode2(L, v, t)
    return L[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function delnode(Linv, v, t)
    Linv2 = (Linv - Linv[:,v]*Linv[:,v]'/(Linv[v,v]))
    return Linv2[[1:v-1;v+1:t], [1:v-1;v+1:t]]
end

function exact(G, k, w :: IOStream) #G- graph; k- the number of vertices in set S

    logw(w,"")
    logw(w,"running exact greedy...")

    n = G.n;
    ans = zeros(3,1);

    L = lapl(G)

    logw(w,"k=1 starts")
    start_time = time()
    Linv = mppinv(L)

    min = Linv[1,1]
    u = 1
    for i = 2:n
        if Linv[i,i] < min
            u = i
            min = Linv[i,i]
        end
    end

    L2 = delnode2(L,u,n)
    Linv = inv(L2)

    elapsed_time = time()-start_time

    logw(w,"k = 1 total time: ",time() - elapsed_time, " (s)")

    start_time = time()
    for dep = 1:k-1
        logw(w,"k=", dep+1, " starts")
        t = time()
        dec = [(sum(Linv[:,i].^2)/Linv[i,i]) for i = 1:n-dep]
        u = indmax(dec)
        
        Linv = delnode(Linv,u,n-dep)
        logw(w,"k = ", dep + 1 ," total time: ", time() - t, " (s)")
    end

    elapsed_time+=time()-start_time
    logw(w,"exact_elapsed_time = ", elapsed_time, " (s)")
    ans[1] = elapsed_time
    logw(w,"calculating CFCC of the solution returned by exact greedy..")
    ans[3] = n/trace(Linv)
    logw(w,"CFCC achieved by exact greedy: ", ans[3])
    logw(w,"")
    return ans
end
