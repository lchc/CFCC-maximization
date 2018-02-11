using Laplacians

function LinvdiagSS(a::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

  f = approxCholLap(a,tol=1e-5);

  n = size(a,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  er = zeros(n,1)

  for i = 1:k
    r = randn(m,1)
    ur = U'*r
    v = zeros(n,1)
    v = f(ur[:])
    er.+= v.^2/k
  end

  return er;

end
