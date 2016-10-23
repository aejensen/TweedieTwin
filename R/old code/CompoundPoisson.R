r.compound.poisson=function(n,p,mu,sigma2){
# From BJ: The Theory of Dispersion Models pp 140
# See my notes of 19 Sep 08
  if ((p>=2)|(p<=1))
    stop("Invalid index; Use 1<p<2");
  alpha=(p-2)/(p-1)
  lambda=sigma2^(1/(1-p))
  theta=(alpha-1)*((mu/lambda)^(1/(alpha-1)))
  pois.rate=(lambda*(alpha-1)/alpha)*(theta/(alpha-1))^alpha
  pois.N=rpois(n,pois.rate)
  rgamma(n,-pois.N*alpha,-theta)
}

