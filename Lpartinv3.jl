using Laplacians


function Lpartinv3SS(L::SparseMatrixCSC{Float64}, f; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

  n = size(L,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  et = ones(n,1)
  x1 = sqrt.(L * et)
  er= zeros(n,1)

  for i = 1:k
    xr = randn(n,1) .* x1
    v = f(xr[:], verbose=false)
    er.+= v.^2/k
  end

  return er;

end
