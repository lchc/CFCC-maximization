include("logw.jl")
include("lapl.jl")
include("Approx.jl")
include("Exact.jl")

function run_exact()
    return length(ARGS) <= 1 || ARGS[2] == "both" || ARGS[2] == "exact"
end

function run_approx()
    return length(ARGS) <= 1 || ARGS[2] == "both" || ARGS[2] == "approx"
end

k = length(ARGS) >= 3 ? parse(Int,ARGS[3]) : 10

datadir = string(ARGS[1],"/")
outFName=string(ARGS[1],".txt")
w = open(outFName, "w")
srand(Int(round(time())))
for rFile in filter( x->!startswith(x, "."), readdir(string(datadir)))
    logw(w, "reading graph from edges list file ", rFile)
    G=read_file(string(datadir,rFile))
    logw(w,"finished reading graph")
    logw(w, "\t LCC.n: ", G.n, "\t LCC.m: ", G.m)
    logw(w, "LCC.n: ", G.n, "\t LCC.m: ", G.m)

    if run_exact() exact_rst = exact(G, k, w) end

    if run_approx() approx_rst = approx(G, k, w) end

    logw(w,"")
    if run_exact()
        logw(w,"exact_elapsed_time: ",  exact_rst[1],  " (s)\t cent[k]: ", exact_rst[3])
    end

    if run_approx()
        logw(w,"approx_elapsed_time: ", approx_rst[1], " (s)\t cent[k]: ", approx_rst[3])
    end

    logw(w,"")
    logw(w,String(fill('*',80)))
    logw(w,"")
end
close(w)
