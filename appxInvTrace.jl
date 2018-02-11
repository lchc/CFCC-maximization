using Laplacians

function appxInvTrace(parL::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=500.0)

  f = approxCholSddm(parL,tol=1e-5, verbose=false)

  n = size(parL,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  #U = wtedEdgeVertexMat(a)
  #m = size(U,1)
  #R = randn(m,k)
  #UR = U'*R;

  #V = zeros(n,k)
  er = zeros(n,1)
  
  rst = 0

  for i = 1:k
  #  V[:,i] .= f(UR[:,i])
    r = rand([1,-1],n,1)
    #ur = U'*r
    v = zeros(n,1)
    v = f(r[:])
    #println(typeof(v))
    #println(typeof(r))
    rst += sum(v.*r)
  end

  return rst/k;

end
