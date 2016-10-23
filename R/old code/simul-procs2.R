my.r.tweedie=function(r.idx,mu,sigma2,n){
    if (r.idx == 0) {
            y = rnorm(n, mu, sqrt(sigma2))
        
    } # (r.idx == 0)
    else if (r.idx == 1) {
        if (sigma2 == 1) 
            y = rpois(n, mu)
        else if (sigma2 > 1) 
            y = rnegbin(n, mu , mu/(sigma2 - 1))
        else {
            prob = 1 - sigma2
            size = round(mu/prob)
            y = rbinom(n, size, prob)
        }
    } # 
    else if ((r.idx > 1) & (r.idx < 2)) 
            y = r.compound.poisson(n,r.idx,mu, sigma2)
    else if (r.idx == 2) 
            y = rgamma(n, 1/sigma2, scale = mu*sigma2)
    else if (((r.idx > 2) & (r.idx < 3)) | 
            (r.idx > 3)) 
            y = rtweedie(n, r.idx, mu, sigma2)
    else if (r.idx == 3) 
            y = r.inv.gaussian(n, mu, sigma2)
    else stop(paste("r=", r.idx, " not implemented"))
    return(y)  
}

sim.geep.twin.ACE=function(nMZ,nDZ,ACE,rACE=c(2,2,2),mu=1,seed){
#	browser()
if (missing(seed)) seed=sample(1:1000000,1)
set.seed(seed)
  if (any((rACE!=3)&(rACE>2)))
      require("tweedie")
  if ((ACE[3] > 1) & (rACE[3]==1)) 
      require("MASS")
alpha.MZ=1
alpha.A.DZ=sqrt(c(1,7+4*sqrt(3))/(4+2*sqrt(3)))
alpha.plus.A=sqrt(3)
if (rACE[1]>0){
Z.A.DZ=matrix(my.r.tweedie(rACE[1],1/sqrt(3),sqrt(3)^rACE[1]*ACE[1],2*nDZ),nc=2)
Z.A.MZ=matrix(my.r.tweedie(rACE[1],1/2,2^rACE[1]*ACE[1],2*nMZ),nc=2)
Q.A.DZ=cbind(alpha.A.DZ[1]*Z.A.DZ[,1]+alpha.A.DZ[2]*Z.A.DZ[,2],alpha.A.DZ[2]*Z.A.DZ[,1]+alpha.A.DZ[1]*Z.A.DZ[,2])
Q.A.MZ=cbind(Z.A.MZ[,1]+Z.A.MZ[,2],Z.A.MZ[,1]+Z.A.MZ[,2])
}
else {
Q.A.DZ=matrix(rep(1,2*nDZ),nc=2)
Q.A.MZ=matrix(rep(1,2*nMZ),nc=2)
}
if (rACE[2]>0){
Z.C.DZ=matrix(my.r.tweedie(rACE[2],1/2,2^rACE[2]*ACE[2],2*nDZ),nc=2)
Z.C.MZ=matrix(my.r.tweedie(rACE[2],1/2,2^rACE[2]*ACE[2],2*nMZ),nc=2)
Q.C.DZ=cbind(Z.C.DZ[,1]+Z.C.DZ[,2],Z.C.DZ[,1]+Z.C.DZ[,2])
Q.C.MZ=cbind(Z.C.MZ[,1]+Z.C.MZ[,2],Z.C.MZ[,1]+Z.C.MZ[,2])
}
else {
Q.C.DZ=matrix(rep(1,2*nDZ),nc=2)
Q.C.MZ=matrix(rep(1,2*nMZ),nc=2)
}

Q.prod.DZ=Q.A.DZ*Q.C.DZ
Q.prod.MZ=Q.A.MZ*Q.C.MZ

y.DZ=sapply(t(Q.prod.DZ),function(q)my.r.tweedie(rACE[3],q*mu,ACE[3]*q^(1-rACE[3]),1))
y.MZ=sapply(t(Q.prod.MZ),function(q)my.r.tweedie(rACE[3],q*mu,ACE[3]*q^(1-rACE[3]),1))
dat=data.frame(pair=rep(1:(nMZ+nDZ),each=2),twin=rep(1:2,nMZ+nDZ),zyg=factor(rep(c("MZ","DZ"),c(2*nMZ,2*nDZ))),y=c(y.MZ,y.DZ))
  return(list(data=dat,call=match.call(),seed=seed))
}

#sim.and.plot.ACE=function(nMZ,nDZ,disp,Tw.idx=c(2,2,0,2),mu,seed,return.dat=F){	
sim.and.plot.ACE=function(nMZ,nDZ,disp,Tw.idx=c(2,2,2),mu,seed,return.dat=F){	
	if (missing(seed)) 
	  seed=sample(1:1000000,1)
#	a=sim.geep.twin(nMZ,nDZ,ACDE=disp,rACDE=Tw.idx,mu,seed)
	a=sim.geep.twin.ACE(nMZ,nDZ,ACE=disp,rACE=Tw.idx,mu,seed)
	a2=matrix(a[[1]][a[[1]]$zyg=="DZ","y"],nr=2)
	a1=matrix(a[[1]][a[[1]]$zyg=="MZ","y"],nr=2)
	plot(a2[1,],a2[2,],col=2,pch=16,cex=0.8)
	points(a1[1,],a1[2,],col=4,pch=16,cex=0.8)
	if (return.dat)return(data=a)
}	
	

