    #############################
    #                                                      #
    #  G E E P  -   T W I N                     #
    #  T W E E D I E   C A S E              #      
    #                                                      #
    #  A Univariate version                    #
    #                                                      #
    #############################

  ###################################
  ##                                                              ##
  ##  M A I N  -  F U N C T I O N                  ##
  ##                                                              ##
  ###################################
  # Ver4b: log-link for r3  
  # Modified to enforce same marginal variances for MZ's and DZ's singletons 
  
    fit.geep.twin=function(fixed, cluster,personid,zygo, dat, offset, ACDE.model="ACE", fit.r3, init.beta,init.gamma,tol=1e-6,max.iter=100,
     bias.correction=F,emp.sandwich=F,bias.adj.emp.var=T){
    # Fits a one-dimensional geep's                                                                                                          
    # fixed:                    A formula or mean-value functions (for linear and non-linear models respectively)
    #                           (a formula implies a linear model)
    # cluster:                  A one-sided fomula with the rhs being a variable (factor) defining the clusters
    # zygo:                     A one-sided fomula with the rhs being a variable (factor) defining the zygosity with levels "DZ" and "MZ"
    #
    # data:                     A data frame
    # ACDE.model                A string consisting of some of the characters "A", "C", "D" and "E"
    # offset:                   A name used for offset
    # fit.r3:                       Order of the fundamental, the middle and the top level Tweedies; If missing r3 will be estimated
    #                           If emp.sandwich=F: Replace the upper left corner of the V by -S_beta (which is equal to V_beta)
    # bias.adj.emp.var          If TRUE the empirical variance is based on the bias corrected psi_gamma; Has no effect if bias.correction=FALSE
    
    
    # Details:
       
       init.lm=function(formula,dat){
         lm(formula,data=dat)$coef
       }
       
       decode.ACDE.model=function(){
       	 if (class(ACDE.model)!="character")
       	   stop("ACDE.model parameter must be a string value")
       	 check.str=c("A","C","D","E")
       	 fit.acde=(is.element(check.str,(model.str=strsplit(ACDE.model,"")[[1]])))
       	 fit.acde[4]=TRUE # Note: Estimation of E (rho^2) is imperative 
       	 return(fit.acde)
       }
       
       mean.and.derivs=function(beta){
        # mm a model matrix                  
#        browser()
          mu.mat<<-xtabs((as.double(exp(mm%*%beta)))~person.var+cluster.var)
          for (i in (1:nClusters))
             derivs[,(i-1)*dim.beta+(1:dim.beta)]<<-(mu.mat[,i]*mm.by.cluster[,(i-1)*dim.beta+(1:dim.beta)])
       } # mean.and.derivs
       
       calc.varY.and.inv=function(){
       	 meat.MZ=as.double((fit.gamma[1:3]%*%gamma[1:3])+ prod(gamma[fit.gamma[1:3]])*2)*V.MZ
         component.idx=(1:3)[fit.gamma[1:3]]
       	 meat.DZ=V.DZ[component.idx[1],,]*gamma[component.idx[1]]
       	 if (length(component.idx)>1)
       	 meat.DZ=meat.DZ+V.DZ[component.idx[2],,]*gamma[component.idx[2]]*(1+meat.DZ)
       	 mu.to.r3=mu.mat^r3
          for (i in (1:nClusters)){
              posi=pos[,i]
              diag.mu.i=diag(mu.mat[posi,i])
              if (isMZ[i]) temp.zygo=meat.MZ 
              else  temp.zygo=meat.DZ
              temp<-rho*(diag(mu.to.r3[posi,i]))+diag.mu.i%*%temp.zygo[posi,posi]%*%diag.mu.i
              var.Y[,(i-1)*2+(1:2)][posi,posi]<<-temp
              inv.var.Y[,(i-1)*2+(1:2)][posi,posi]<<-solve(temp)
         }
       } # calc.varY.and.inv
  
       calc.psi.beta=function(){
         psi.beta[]<<-0
         for (i in (1:nClusters)){
            posi<-pos[,i]
            D.i=derivs[posi,(i-1)*dim.beta+(1:dim.beta)]
            psi.beta<<-psi.beta+t(D.i)%*%(inv.var.Y[,(i-1)*2+(1:2)][posi,posi])%*%(response.mat[posi,i]-mu.mat[posi,i])
         }
       } # calc.psi.reg
       
       sensitivity.reg=function(){
         S.beta[,]<<-0
         for (i in 1:nClusters){
            posi<-pos[,i]
            D.i=derivs[posi,(i-1)*dim.beta+(1:dim.beta)]
            S.beta<<-S.beta-t(D.i)%*%(inv.var.Y[,(i-1)*2+(1:2)][posi,posi])%*%D.i
         }
       } # sensitivity.reg         
    
      calc.psi.and.S.assoc=function(){
#         browser(fit.gamma)
         psi.gamma[fit.gamma]<<-0
         S.gamma[fit.gamma,fit.gamma]<<-0
         resid.mat<-response.mat-mu.mat
         dCi.dgamma=array(0,c(5,2,2))           # Derivatives of C_i wrt omega_A, omega_C, omega_D, rho, r3
         dCi.dgamma.aux=array(0,c(3,2,2,2)) # 2 x 2 matrices for MZ and DZ for each of A, C and E
         ACD.idx=1:3
      #   browser()
         for (i in 1:3){
           if (fit.gamma[i]==T){     
         	   if (any(fit.gamma[ACD.idx[-i]])) {                                                   # Two components 
             	   	j=(ACD.idx[-i])[fit.gamma[ACD.idx[-i]]]                                     # Index of the other component
                	dCi.dgamma.aux[i,1,,]=V.MZ+2*V.MZ*gamma[j]                      # MZ two components
                	dCi.dgamma.aux[i,2,,]=V.DZ[i,,]+V.DZ[i,,]*V.DZ[j,,]*gamma[j]   # DZ two components
         	   	}
         	   else {
         	   	    dCi.dgamma.aux[i,1,,]=V.MZ                                                    # MZ single component
                	dCi.dgamma.aux[i,2,,]=V.DZ[i,,]                                               # DZ  single component
         	   }}}
         zygo.idx=2-isMZ     # MZ=1 DZ=2
         for (i in (1:nClusters)){
            # First prepare matrices for cluster i
            posi<-pos[,i]
            C.i<-var.Y[,(i-1)*2+(1:2)][posi,posi]
            C.i.inv<-inv.var.Y[,(i-1)*2+(1:2)][posi,posi]
            resid.i<-resid.mat[posi,i]
            gamma.resid=(resid.i%*%t(resid.i))-C.i
            diag.mu.i=diag(mu.mat[posi,i])
            log.diag.mu.i=diag(log(mu.mat[posi,i]))
            zi=zygo.idx[i]

            if (fit.gamma[1]==T)  # A
                dCi.dgamma[1,posi,posi]=diag.mu.i%*% dCi.dgamma.aux[1,zi,,]%*%diag.mu.i
            if (fit.gamma[2]==T)  # C
                dCi.dgamma[2,posi,posi]=diag.mu.i%*% dCi.dgamma.aux[2,zi,,]%*%diag.mu.i
            if (fit.gamma[3]==T)  # D
                dCi.dgamma[3,posi,posi]=diag.mu.i%*% dCi.dgamma.aux[3,zi,,]%*%diag.mu.i
            if (fit.gamma[4]==T)  # E
                dCi.dgamma[4,posi,posi]=diag.mu.i^r3
            if (fit.gamma[5]==T)  # r3
#                dCi.dgamma[5,posi,posi]=r3*gamma[4]*diag.mu.i^(r3-1)
                dCi.dgamma[5,posi,posi]=r3*gamma[4]*diag.mu.i^r3*log.diag.mu.i
            if (adjust.bias==T){
               stop("Bias correction not implemented yet")
               Godambe=Godambe+t(D.i)%*%C.i.inv%*%D.i}
            # Now calculate psi and S
            for (j in (1:dim.gamma)[fit.gamma]){                                
                 W.ij=C.i.inv%*%dCi.dgamma[j,posi,posi]%*%C.i.inv
                 psi.gamma[j] <<- psi.gamma[j]+sum(diag(W.ij%*%gamma.resid))
                 for (k in (j:dim.gamma)[fit.gamma[j:dim.gamma]]) # Calc S.gamma and bias adjustment
                    S.gamma[j,k]<<-S.gamma[j,k]-sum(diag(W.ij%*%(dCi.dgamma[k,posi,posi])))
            } # for j
         } # for  i
         # browser()
         idx=unlist(sapply(1:(dim.gamma-1),function(i,maxi){((i-1)*maxi+i+1):(i*maxi)},dim.gamma))
         S.gamma[idx]<<-t(S.gamma)[idx]
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
#         Upper left block: -S_beta^-1, lower rigth S_gamma^-1*V_gamma*S_gamma^-1, where V_gamma is empirical
          calc.empirical.variability(T)
          idx=c(rep(T,dim.beta),fit.gamma)
          V.gamma=(emp.var[idx,idx])[-(1:dim.beta),-(1:dim.beta)]
          work.godambe.inv=matrix(0,ncol=sum(idx),nrow=sum(idx))
          work.godambe.inv[(1:dim.beta),(1:dim.beta)]=-S.beta.inv
          work.godambe.inv[-(1:dim.beta),-(1:dim.beta)]=S.gamma.inv%*%V.gamma%*%S.gamma.inv
          return(work.godambe.inv)
       } # calc.work.godambe
       
    ### MAIN BODY ###
    
    # INITIALIZATION 
    # INIT - Arrange Data
        m<-m.rep<-match.call()
        nm <- names(m)[-1]
        keep <- is.element(nm, c("dat", "na.action")) # "weights", "data", "subset", "na.action"
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
    
#        corr.A.MZ<-corr.C<-matrix(rep(1,4),ncol=2)
#        corr.A.DZ<-matrix(c(1,0.5,0.5,1),ncol=2)
#        corr.E<-diag(c(1,1))
#        dim.sigma=3 # sigma_A, sigma_C, sigma_E
 
#        if (missing(fit.sigma))
#          fit.sigma=rep(T,dim.sigma)
#        if (missing(init.sigma))
#           init.sigma<-c(0.5,0.3,0.2) # sigma2.A, sigma2.C, sigma2.E
#        sigma=rep(0,length(fit.sigma))
#        sigma[1]=init.sigma[1]
#        sigma[2]=init.sigma[2]
#        sigma[3]=init.sigma[3]

        a2.to.a1.ratio.DZ=c(A=2+sqrt(3),C=1,D=4+sqrt(15))
        a1.coeffs.DZ=c(a1.A=1/sqrt(4+2*sqrt(3)),a1.C=1,a1.D=1/sqrt(16+4*sqrt(15))) # Row 1: a1, Row 2: a2, Col1-3: A, C, D
        V.MZ=matrix(c(2,2,2,2),nrow=2)
        t1=a1.coeffs.DZ^2*(1+a2.to.a1.ratio.DZ^2)
        t2=2*a2.to.a1.ratio.DZ*a1.coeffs.DZ^2
        V.DZ=array(0,c(3,2,2))
        V.DZ[1,,]=matrix(c(t1[1],t2[1],t2[1],t1[1]),nr=2)    # A
        V.DZ[2,,]=matrix(c(t1[2],t2[2],t2[2],t1[2]),nr=2)    # C
        V.DZ[3,,]=matrix(c(t1[3],t2[3],t2[3],t1[3]),nr=2)    # D
        
        dim.gamma=5 # gamma_A, gamma_C, gamma_D, gamma_E, r3
        fit.gamma=c(decode.ACDE.model(),F)
        if (missing(init.gamma))
           init.gamma<-c(0.5,0.3,0.4,0.2,2) # gamma2.A, gamma2.C, gamma2.E, r3
        if ((missing(fit.r3))||(fit.r3==T))
           fit.gamma[5]=T        
        else  
          init.gamma[5]=as.double(fit.r3) 
          
        init.gamma[5]=log(init.gamma[5])  
        gamma=init.gamma
        sigma_A=gamma[1]
        sigma_C=gamma[2]
        sigma_D=gamma[3]
        rho=gamma[4]
        r3=exp(gamma[5])

        
           
        S.gamma<-matrix(0,ncol=dim.gamma,nrow=dim.gamma)
        psi.gamma<-rep(0,dim.gamma)
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
      phase<-1 # phase==1 : Update only beta  # # phase==1 : Update both beta and gamma
        
    # ITERATIONS    
       act.tol=1
    repeat{
#browser()
      cat("it: ",it," phase: ",phase," -  beta: ",beta,";     A, C, E: ",gamma[fit.gamma],"    act.tol: ",act.tol,"\n")
      it=it+1
      mean.and.derivs(beta)
      calc.varY.and.inv()
      sensitivity.reg()
      S.beta.inv=solve(S.beta)
      calc.psi.beta()
      old.beta<-beta
 #     browser()
      beta<-as.vector(old.beta-S.beta.inv%*%psi.beta)
      if (phase==1){
         act.tol=sqrt(sum((beta-old.beta)^2))
         phase=ifelse(act.tol<10e-2,2,1)
         act.tol=10 # Make sure an iteration in phase 2 is take
 #        phase=1 # DEBUG
       }
      else{ # if (phase==1)
#       browser()
        calc.psi.and.S.assoc()
#       browser()
        S.gamma.inv=solve(S.gamma[fit.gamma,fit.gamma])
        old.gamma=gamma
        gamma[fit.gamma]=gamma[fit.gamma]-as.vector(S.gamma.inv%*%(psi.gamma[fit.gamma]))
        sigma_A=gamma[1]
        sigma_C=gamma[2]
        sigma_D=gamma[3]
        rho=gamma[4]
        r3=exp(gamma[5])
#        w.godambe.inv=calc.work.godambe.inv()
        x=c(beta,gamma[fit.gamma])-c(old.beta,old.gamma[fit.gamma])
#        act.tol=t(x)%*%w.godambe.inv%*%x    # 2nd criteria 
#         act.tol=max(abs(x)*sqrt(diag(w.godambe.inv))) # 1st criteria
          act.tol=sqrt(sum(x*x))
        if (it > max.iter){
          cat("\n Convergence not achieved in ",max.iter, " iterations\n")
          break
        }    
#        browser()
        if (act.tol<tol){    # 1st criteria 
#         if (act.tol<tol^2){ # 2nd criteria 
          se.beta=sqrt(diag(-S.beta.inv))
          se.gamma=sqrt(diag(-S.gamma.inv))   # ??? psi.gamma is not a quasi-score!!!
          calc.psi.and.S.assoc()
          disp=c(gamma2.A=gamma[1],gamma2.C=gamma[2],gamma2.D=gamma[3],gamma2.E=gamma[4])[fit.gamma[1:4]]
          return(c(beta=beta,disp=disp,r3=r3,psi.beta=as.double(psi.beta),psi.gamma=psi.gamma,it=it,tol=tol))
        }
      } # if (phase==2)
      } # repeat
      return(NULL)
    }
