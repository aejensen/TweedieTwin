r.inv.gaussian=function(n,mu,sigma2){
# J.R. Michael, W.R. Schucany and R.W.Haas. 1976.
# Generating Random Variates Using Transformations with Multiple Roots
# The Am. Statistician, Vol. 30. No. 2 pp 88-90
# mu: mean, sigma2: dispersion (Var(X)=mu^3*sigma2)
  v=rchisq(n,1)
  y=runif(n)
  w=mu*v
  cx=0.5*mu*sigma2
  x1=mu+cx*(w-sqrt(w*(4/sigma2+w)))
  ifelse(y>=mu/(mu+x1),mu*mu/x1,x1)
}

d.inv.gaussian=function(x,mu,sigma2){
exp(-(x-mu)^2/(2*mu^2*sigma2*x))/sqrt(2*pi*sigma2*x^3)
}

p.inv.gaussian=function(q,mu,sigma2){
t1=1/sqrt(q*sigma2)
t2=q/mu
pnorm(t1*(t2-1))+pnorm(-t1*(t2+1))*exp(2/(mu*sigma2))
}

q.inv.gaussian=function(p,mu,sigma2){
f=function(x){p.inv.gaussian(x,mu,sigma2)-p}
nlm(f,1)$estimate
}
