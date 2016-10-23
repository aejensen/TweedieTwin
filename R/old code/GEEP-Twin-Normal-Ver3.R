    ##########################
    #                        #
    #  G E E P  -   T W I N  #
    #  N O R M A L  C A S E  #      
    #                        #
    #  A Univariate version  #
    #                        #
    ##########################

  ###################################
  ##                               ##
  ##  M A I N  -  F U N C T I O N  ##
  ##                               ##
  ###################################

    fit.geep.twin.norm=function(fixed, cluster, personid,zygo, dat, offset, 
                           init.beta,init.sigma,fit.sigma,tol=1e-6,max.iter=100,
                           bias.correction=F,emp.sandwich=F,bias.adj.emp.var=T){
    # Fits a one-dimensional geep's                                                                                                          
    # fixed:                    A formula or mean-value functions (for linear and non-linear models respectively)
    #                           (a formula implies a linear model)
    # cluster:                  A one-sided fomula with the rhs being a variable (factor) defining the clusters
    # zygo:                     A one-sided fomula with the rhs being a variable (factor) defining the zygosity with levels "DZ" and "MZ"
    #
    # data:                     A data frame
    # offset:                   A name used for offset
    # r3:                       Order of the fundamental, the middle and the top level Tweedies; Integers
    # bias.correction           Use bias corrected est func for estimation of nuisance paramter sigma
    # emp.sandwich              Sandwich estimator; S^-1 * V * S^-T. Here V is the empirical variability. 
    #                           If emp.sandwich=F: Replace the upper left corner of the V by -S_beta (which is equal to V_beta)
    # bias.adj.emp.var          If TRUE the empirical variance is based on the bias corrected psi_sigma; Has no effect if bias.correction=FALSE
    
    
    # Details:
       
       init.lm=function(formula,dat){
         lm(formula,data=dat)$coef
       }
       mean.and.derivs=function(beta){
        # mm a model matrix                  
          mu.mat<<-xtabs((mm%*%beta)~person.var+cluster.var)
          for (i in (1:nClusters))
            derivs[,(i-1)*dim.beta+(1:dim.beta)]<<-mm.by.cluster[,(i-1)*dim.beta+(1:dim.beta)]
       } # mean.and.derivs
       
       calc.varY.and.inv=function(){
         var.y.ce <- corr.C*sigma[2]+corr.E*sigma[3]
         var.Y.MZ<<-var.y.ce+corr.A.MZ*sigma[1]
         var.Y.DZ<<-var.y.ce+corr.A.DZ*sigma[1]
         inv.var.Y.MZ<<-solve(var.Y.MZ)
         inv.var.Y.DZ<<-solve(var.Y.DZ)
       } # calc.varY.and.inv
  
       calc.psi.beta=function(){
         psi.beta[]<<-0
         for (i in (1:nClusters)){
            posi<-pos[,i]
            D.i=derivs[posi,(i-1)*dim.beta+(1:dim.beta)]
#            browser()
            if (isMZ[i])
              psi.beta<<-psi.beta+t(D.i)%*%inv.var.Y.MZ[posi,posi]%*%(response.mat[posi,i]-mu.mat[posi,i])
            else  
              psi.beta<<-psi.beta+t(D.i)%*%inv.var.Y.DZ[posi,posi]%*%(response.mat[posi,i]-mu.mat[posi,i])
         }
       } # calc.psi.reg
       
       sensitivity.reg=function(){
         S.beta[,]<<-0
         for (i in 1:nClusters){
            posi<-pos[,i]
            D.i=derivs[posi,(i-1)*dim.beta+(1:dim.beta)]
 #           browser()
            if (isMZ[i])
              S.beta<<-S.beta-t(D.i)%*%inv.var.Y.MZ[posi,posi]%*%D.i
            else  
              S.beta<<-S.beta-t(D.i)%*%inv.var.Y.DZ[posi,posi]%*%D.i
         }
       } # sensitivity.reg         
    
      calc.psi.and.S.assoc=function(){
         # browser()
         psi.sigma[fit.sigma]<<-0
         S.sigma[fit.sigma,fit.sigma]<<-0
         resid.mat<-response.mat-mu.mat
         dCi.dsigma=array(0,c(dim.sigma,2,2))
         for (i in (1:nClusters)){
            # First prepare matrices for cluster i
            posi<-pos[,i]
            if (isMZ[i]){
              C.i<-var.Y.MZ[posi,posi]
              C.i.inv<-inv.var.Y.MZ[posi,posi]
              corr.A<-corr.A.MZ}
            else {
              C.i<-var.Y.DZ[posi,posi]
              C.i.inv<-inv.var.Y.DZ[posi,posi]
              corr.A<-corr.A.DZ}            
            resid.i<-resid.mat[posi,i]
            sigma.resid=(resid.i%*%t(resid.i))-C.i
            if (fit.sigma[1]==T)
              dCi.dsigma[1,posi,posi]=corr.A[posi,posi]
            if (fit.sigma[2]==T)                                                       
              dCi.dsigma[2,posi,posi]=corr.C[posi,posi]
            if (fit.sigma[3]==T)
              dCi.dsigma[3,posi,posi]=corr.E[posi,posi]
            if (adjust.bias==T){
               stop("Bias correction not implemented yet")
               Godambe=Godambe+t(D.i)%*%C.i.inv%*%D.i}
            # Now calculate psi and S
            for (j in (1:dim.sigma)[fit.sigma]){                                
                 W.ij=C.i.inv%*%dCi.dsigma[j,posi,posi]%*%C.i.inv
                 psi.sigma[j] <<- psi.sigma[j]+sum(diag(W.ij%*%sigma.resid))
                 for (k in (j:dim.sigma)[fit.sigma[j:dim.sigma]]) # Calc S.sigma and bias adjustment
                    S.sigma[j,k]<<-S.sigma[j,k]-sum(diag(W.ij%*%(dCi.dsigma[k,posi,posi])))
            } # for j
         } # for  i
         # browser()
         idx=unlist(sapply(1:(dim.sigma-1),function(i,maxi){((i-1)*maxi+i+1):(i*maxi)},dim.sigma))
         S.sigma[idx]<<-t(S.sigma)[idx]
#         browser()
    }
  
       calc.empirical.variability=function(use.emp.V.beta){
         stop("NOT CHEKKED")
         for (i in (1:nClusters)){
            emp.var<<-emp.var+emp.var.aux[,i]%*%t(emp.var.aux[,i])
         }
         if (use.emp.V.beta==F) 
            emp.var[1:dim.beta,1:dim.beta]<<- (-S.beta)
       }
         
       calc.work.godambe.inv=function(){
          stop("NOT CHEKKED")
#         Calculates a block diagonal inverse Godambe-proxy: 
#         Upper left block: -S_beta^-1, lower rigth S_sigma^-1*V_sigma*S_sigma^-1, where V_sigma is empirical
          calc.empirical.variability(T)
          idx=c(rep(T,dim.beta),fit.sigma)
          V.sigma=(emp.var[idx,idx])[-(1:dim.beta),-(1:dim.beta)]
          work.godambe.inv=matrix(0,ncol=sum(idx),nrow=sum(idx))
          work.godambe.inv[(1:dim.beta),(1:dim.beta)]=-S.beta.inv
          work.godambe.inv[-(1:dim.beta),-(1:dim.beta)]=S.sigma.inv%*%V.sigma%*%S.sigma.inv
          return(work.godambe.inv)
       } # calc.work.godambe
       
    ### MAIN BODY ###
    
    # INITIALIZATION 
    # INIT - Arrange Data
#browser()
        m<-m.rep<-match.call()
        nm <- names(m)[-1]
        keep <- is.element(nm, c("dat", "na.action")) # "weights", "data",  "na.action"
        for (i in nm[!keep]) m[[i]] <- NULL
        if ((class(fixed)!="formula")||(length(fixed)!=3))
          stop("'fixed' must be a double sided formula") 
        response.name<-as.name(as.character(fixed)[[2]])
        if ((class(cluster)!="formula")||(length(cluster)!=2))
          stop("'cluster' must be a one sided formula") 
        cluster.name<-as.name(as.character(cluster)[[2]])
        allvars <- c(all.vars(fixed),all.vars(cluster))
        if ((class(zygo)!="formula")||(length(zygo)!=2))
          stop("'zygo' must be a one sided formula") 
        zygo.name<-as.name(as.character(zygo)[[2]])
        if (missing(zygo))
          stop("'zygo' formula must be specified") 
        allvars <- c(allvars,all.vars(zygo))

        if ((class(personid)!="formula")||(length(personid)!=2))
          stop("'personid' must be a one sided formula") 
        person.name<-as.name(as.character(personid)[[2]])
        if (missing(personid))
          stop("'personid' formula must be specified") 
        allvars <- c(allvars,all.vars(personid))
        Terms <- if(missing(dat))
            terms(fixed)
        else terms(fixed, data = dat)
        if (length(off <- attr(Terms, "offset")))
            allvars <- c(allvars, as.character(attr(Terms, "variables"))[off+1])
        allvars <- unique(allvars)        
        m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
        environment(m$formula) <- environment(fixed)
        m$drop.unused.levels <- TRUE
        m[[1]] <- as.name("model.frame")
        mf <- eval.parent(m)
        response.var=eval(response.name,mf)
        zygo.var=eval(zygo.name,mf)
        if (!is.factor(zygo.var)||(levels(zygo.var)!=c("DZ","MZ")) )
          stop("zygosity variable must be a factor with levels 'MZ' and 'DZ'") 
        person.var=eval(person.name,mf)
        cluster.var=eval(cluster.name,mf)
        isMZ=xtabs(~zygo.var + cluster.var)[2,]>0
        response.mat<-xtabs(response.var ~ person.var + cluster.var)
#        browser()
        pos <- xtabs( ~ person.var + cluster.var) > 0
        nClusters <- ncol(pos)
        if (nrow(pos)!=2) stop("'personid' mis-specified") 
#browser()    
        mm <- model.matrix(fixed,mf)
        dim.beta=ncol(mm)
        mm.by.cluster <- derivs <- matrix(0,ncol=ncol(mm)*nClusters,nrow=2)
        for (i in (1:nClusters)){
          if (i==1) idx<<-(1:sum(pos[,1])) else idx <-(sum(pos[,1:(i-1)])+1):sum(pos[,1:i])
          mm.by.cluster[pos[,i],(i-1)*ncol(mm)+(1:ncol(mm))]<-mm[idx,]
        }  
        adjust.bias<-bias.correction
        psi.beta<-rep(0,dim.beta)
        S.beta<-matrix(0,ncol=dim.beta,nrow=dim.beta)
        var.Y<-inv.var.Y<-matrix(0,nrow=2,ncol=nClusters*2)
    
        corr.A.MZ<-corr.C<-matrix(rep(1,4),ncol=2)
        corr.A.DZ<-matrix(c(1,0.5,0.5,1),ncol=2)
        corr.E<-diag(c(1,1))
        dim.sigma=3 # sigma_A, sigma_C, sigma_E
 
        if (missing(fit.sigma))
          fit.sigma=rep(T,dim.sigma)
        if (missing(init.sigma))
           init.sigma<-c(0.5,0.3,0.2) # sigma2.A, sigma2.C, sigma2.E
        sigma=rep(0,length(fit.sigma))
        sigma[1]=init.sigma[1]
        sigma[2]=init.sigma[2]
        sigma[3]=init.sigma[3]
           
        S.sigma<-matrix(0,ncol=dim.sigma,nrow=dim.sigma)
        psi.sigma<-rep(0,dim.sigma)
        blup<-list(Z.i=rep(0,nClusters),Z.iA=matrix(0,ncol=2*nClusters,nrow=2),Q.iA=matrix(0,ncol=2*nClusters,nrow=2))
#        emp.var.aux=matrix(0,nrow=dim.beta+dim.sigma,ncol=nClusters) # i'th column: the value of (psi_beta,psi_sigma)^T
#        emp.var=matrix(0,ncol=dim.beta+dim.sigma,nrow=dim.beta+dim.sigma)               
        
    # INIT - Other Settings
      act.tol=10
      it=0
      if (missing(init.beta)){
           beta=init.lm(fixed,dat)
         }
      else   
        beta=init.beta
      phase<-1 # phase==1 : Update only beta  # # phase==1 : Update both beta and sigma
        
    # ITERATIONS    
       act.tol=1
    repeat{
#browser()
      cat("it: ",it," phase: ",phase," -  beta: ",beta,";     A, C, E: ",sigma[1:3],"    act.tol: ",act.tol,"\n")
      it=it+1
#      browser()
      mean.and.derivs(beta)
      calc.varY.and.inv()
      sensitivity.reg()
      S.beta.inv=solve(S.beta)
      calc.psi.beta()
      old.beta<-beta
      beta<-as.vector(old.beta-S.beta.inv%*%psi.beta)
      if (phase==1){
         act.tol=sqrt(sum((beta-old.beta)^2))
         phase=ifelse(act.tol<10e-2,2,1)
         act.tol=10 # Make sure an iteration in phase 2 is take
 #        phase=1 # DEBUG
       }
      else{ # if (phase==1)
#        browser()
        calc.psi.and.S.assoc()
       # browser()
        S.sigma.inv=solve(S.sigma[fit.sigma,fit.sigma])
        old.sigma=sigma
        sigma[fit.sigma]=sigma[fit.sigma]-as.vector(S.sigma.inv%*%(psi.sigma[fit.sigma]))
#        w.godambe.inv=calc.work.godambe.inv()
        x=c(beta,sigma[fit.sigma])-c(old.beta,old.sigma[fit.sigma])
#        act.tol=t(x)%*%w.godambe.inv%*%x    # 2nd criteria 
#         act.tol=max(abs(x)*sqrt(diag(w.godambe.inv))) # 1st criteria
          act.tol=sqrt(sum(x*x))
        if (it > max.iter){
          cat("\n Convergence not achieved in ",max.iter, " iterations\n")
          break
        }    
        
        if (act.tol<tol){    # 1st criteria 
#         if (act.tol<tol^2){ # 2nd criteria 
          se.beta=sqrt(diag(-S.beta.inv))
          se.sigma=sqrt(diag(-S.sigma.inv))   # ??? psi.sigma is not a quasi-score!!!
          calc.psi.and.S.assoc()
          disp=c(sigma2.A=sigma[1],sigma2.C=sigma[2],sigma2.E=sigma[3])
          return(c(beta=beta,disp=disp,psi.beta=as.double(psi.beta),psi.sigma=psi.sigma,it=it,tol=tol))
        }
      } # if (phase==2)
      } # repeat
      return(NULL)
    }
