using Laplacians

function Lpartinv2SS(n :: Int, f, a::SparseMatrixCSC{Float64}; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

  k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  er = zeros(n,1)

  for i = 1:k
    r = randn(m,1)
    ur = U'*r
    v = f(ur[:], verbose=false);
    er.+= v.^2/k
  end

  return er;

end
