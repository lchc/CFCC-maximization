using Laplacians

function LpartinvSS(n :: Int, f; ep=0.3, matrixConcConst=4.0, JLfac=200.0)

  k = round(Int, JLfac*log(n)) # number of dims for JL

  er = zeros(n,1)

  for i in 1:k
    r = randn(n,1)
    v = f(r[:], verbose=false);
    er.+= v.^2/k
  end

  return er;

end
