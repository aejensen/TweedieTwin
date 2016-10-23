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
# ver 4.0 22. Januar 2011: Extended to allow for multivariate responses
#                          Removed optional offsets terms
# UTILS 
# Theser are taken from MCMCPack:
lower.tri=function (x, diag = FALSE){
    x <- as.matrix(x)
    if (diag) 
        row(x) >= col(x)
    else row(x) > col(x)
}

upper.tri=function(x, diag = FALSE){
    x <- as.matrix(x)
    if (diag) 
        row(x) <= col(x)
    else row(x) < col(x)
}

xpnd=function(x, nrow = NULL){
    dim(x) <- NULL
    if (is.null(nrow)) 
        nrow <- (-1 + sqrt(1 + 8 * length(x)))/2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag = TRUE)] <- 0
    output <- output + t(hold)
    return(output)}

vech=function(x) {
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
        stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)}

##############

    fit.geep.twin.norm.multi=function(fixed, cluster, personid,zygo, data, subset, 
                           init.beta,ACDE.model,tol=1e-6,max.iter=100, emp.sandwich=F){
    # Fits a one-dimensional geep's                                                                                                          
    # fixed:                    A list of formulae for the mean-values y_j~X beta_j
    # cluster:                  A one-sided fomula with the rhs being a variable (factor) defining the clusters
    # personid                  A one-sided fomula with the rhs being a variable (factor) defining the perons within clusters
    # zygo:                     A one-sided fomula with the rhs being a variable (factor) defining the zygosity with levels "DZ" and "MZ"
    # data:                     A data frame
    # subset:                   Not implemented; Consider to omit                                         
    # init.beta                 Start values for fixed part
    # ACDE.model                A list of one-sided formulae for linear functions for the A, C and D components,
    #                           with item names in c("A","C","D","E"); e.g. A=~age
    # tol                       Convergence criterion
    # max.iter                  Convergence criterion
    # emp.sandwich              Sandwich estimator; S^-1 * V * S^-T. Here V is the empirical variability. 
    #                           If emp.sandwich=F: Replace the upper left corner of the V by -S_beta (which is equal to V_beta)
    
    
    # Details:
       
       init.lm=function(fix,dat){
         coef.list=list()
         for (i in (1:length(fix))) coef.list[[i]]=lm(fix[[i]],data=dat)$coef
       }
       
       mean.and.derivs=function(beta){
        # mm a model matrix                  
browser()        
          mu.mat<<-xtabs((mm%*%beta)~person.var+cluster.var)
          for (i in (1:nClusters))
            derivs[,(i-1)*dim.beta+(1:dim.beta)]<<-mm.by.cluster[,(i-1)*dim.beta+(1:dim.beta)]
       } # mean.and.derivs
       
       calc.varY.and.inv=function(){
         sigma.A=xpnd(Z.A%*%alpha.A)
         sigma.C=xpnd(Z.C%*%alpha.C)
         sigma.D=xpnd(Z.D%*%alpha.D)
         var.Y.ce<-kronecker(corr.C,sigma.C)+kronecker(corr.E,sigma.E)
         var.Y.MZ<<-var.y.ce+kronecker(corr.A.MZ,sigma.A)
         var.Y.DZ<<-var.y.ce+kronecker(corr.A.DZ,sigma.A)
         inv.var.Y.MZ<<-solve(var.Y.MZ)
         inv.var.Y.DZ<<-solve(var.Y.DZ)
       } # calc.varY.and.inv
  
       calc.psi.beta=function(){
         psi.beta[]<<-0
         for (i in (1:nClusters)){
            posi<-pos[,i]
            D.i=derivs[posi,(i-1)*dim.beta+(1:dim.beta)]
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
            if (isMZ[i])
              S.beta<<-S.beta-t(D.i)%*%inv.var.Y.MZ[posi,posi]%*%D.i
            else  
              S.beta<<-S.beta-t(D.i)%*%inv.var.Y.DZ[posi,posi]%*%D.i
         }
       } # sensitivity.reg         
    
      calc.psi.and.S.assoc=function(){
         # browser()
         psi.sigma[]<<-0
         S.sigma[,]<<-0
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
            # Now calculate psi and S
            for (j in (1:dim.sigma)[fit.sigma]){                                
                 W.ij=C.i.inv%*%dCi.dsigma[j,posi,posi]%*%C.i.inv
                 psi.sigma[j] <<- psi.sigma[j]+sum(diag(W.ij%*%sigma.resid))
                 for (k in (j:dim.sigma)[fit.sigma[j:dim.sigma]]) # Calc S.sigma 
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

    validate.formulae=function(formulae){
      rnames=vars=NULL  
      if ((class(formulae)=="formula")&&(length(formulae)==3)){
          formulae<-list(formulae)}
      if (class(formulae)=="list")
        if (all(unlist(lapply(formulae,function(y){((class(y)=="formula")&&(length(y)==3))})))){
          rnames=unlist(lapply(formulae,function(y){as.character(y)[[2]]}))
          vars=unlist(lapply(formulae,function(y){all.vars(y)}))
        }
      vars=setdiff(vars,rnames)  
      return(list(rnames=rnames,vars=vars))
    }
    
    # INITIALIZATION 
    # INIT - Arrange Data
        m<-m.rep<-match.call()
        nm <- names(m)[-1]
        keep <- is.element(nm, c("data", "na.action")) # "weights", "data", "subset", "na.action"
        for (i in nm[!keep]) m[[i]] <- NULL
        split.fixed=validate.formulae(fixed)
        if (is.null(split.fixed[[1]]))
          stop("'fixed' must be a double sided formula or a list of double sided formulae") 
        split.varcomps=validate.formulae(ACDE.model)
        if ((class(cluster)!="formula")||(length(cluster)!=2))
          stop("'cluster' must be a one sided formula") 
        cluster.name<-as.name(as.character(cluster)[[2]])
        allvars <- c(split.fixed$vars,split.varcomps$vars,all.vars(cluster)) #union(split.fixed$vars,all.vars(cluster))
        if ((class(zygo)!="formula")||(length(zygo)!=2))
          stop("'zygo' must be a one sided formula") 
        zygo.name<-as.name(as.character(zygo)[[2]])
        if (missing(zygo))
          stop("'zygo' formula must be specified") 
        allvars <- union(allvars,all.vars(zygo))
        if ((class(personid)!="formula")||(length(personid)!=2))
          stop("'personid' must be a one sided formula") 
        person.name<-as.name(as.character(personid)[[2]])
        if (missing(personid))
          stop("'personid' formula must be specified") 
        allvars <- c(allvars,all.vars(personid))
        
        allvars <- unique(allvars)        
        m$drop.unused.levels <- TRUE
        m[[1]] <- as.name("model.frame")
        mf<-NULL
        for (y in split.fixed$rnames){
          m$formula <- as.formula(paste("~", paste(c(y,allvars), collapse = "+")))
          environment(m$formula) <- environment(fixed)
          if (is.null(mf))
            mf <- eval.parent(m)
          else
            mf<-merge(mf,eval.parent(m),all=T)
        }      
        mf=mf[order(mf$pair,mf$twin),]
        zygo.var=eval(zygo.name,mf)
        if (!is.factor(zygo.var)||(levels(zygo.var)!=c("DZ","MZ")) )
          stop("zygosity variable must be a factor with levels 'MZ' and 'DZ'") 
        person.var=eval(person.name,mf)
        cluster.var=eval(cluster.name,mf)
        isMZ=xtabs(~zygo.var + cluster.var)[2,]>0
        response.vars=sapply(split.fixed$rnames,function(y)eval(as.name(y),mf))
        dim.response=ncol(response.vars)
        response.mat<-xtabs(response.vars ~ person.var + cluster.var,na.action=na.pass) # response.mat[,,i] response matrix for i'th response
        valid.response<-xtabs(!is.na(response.vars) ~ person.var + cluster.var)==T
        nClusters <- dim(response.mat)[2]
        if (dim(response.mat)[1]!=2) stop("'personid' mis-specified") 
        if (class(fixed)=="fomula")
          model.matrices<-list(model.matrix(fixed,mf))
        else  
          model.matrices <- lapply(fixed,function(fix)model.matrix(fix,mf))
        dim.beta=unlist(lapply(model.matrices,function(X)ncol(X)))
        total.dim.beta=sum(dim.beta)
        model.matrices.map=sapply(split.fixed$rnames,function(y)model.matrix(as.formula(paste(y,"-1+pair+twin",sep="~")),mf))
        if (is.null(names(ACDE.model))||!all(is.element(names(ACDE.model),c("A","C","D","E"))))
          stop("The 'ACDE.model' parameter is not correctly specified")
        if (class(ACDE.model)=="formulae")                             # Note that the same design matrix is used for all
          varComp.model.matrices<-list(model.matrix(ACDE.model,mf))    # parameters in a variance component
        else  
          varComp.model.matrices <- lapply(ACDE.model,function(ACDE)model.matrix(ACDE,mf))
        fit.ACDE=is.element(c("A","C","D","E"),names(ACDE.model))
        names(fit.ACDE)=c("A","C","D","E")
        dim.gamma=unlist(lapply(varComp.model.matrices,function(Z){ncol(Z)*dim.response*(dim.response+1)/2}))
        total.dim.gamma=sum(dim.gamma)
        psi.beta<-rep(0,total.dim.beta)
        S.beta<-matrix(0,ncol=total.dim.beta,nrow=total.dim.beta)
        var.Y<-inv.var.Y<-matrix(0,nrow=2*dim.response,ncol=nClusters*2*dim.response)
                                                                                    
        corr.A.MZ<-corr.C<-matrix(rep(1,4),ncol=2)
        corr.A.DZ<-matrix(c(1,0.5,0.5,1),ncol=2)
        corr.E<-diag(c(1,1))
        
#        init.sigma<-c(0.5,0.3,0.2) # sigma2.A, sigma2.C, sigma2.E
#        sigma=rep(0,length(fit.sigma))
#        sigma[1]=init.sigma[1]
#        sigma[2]=init.sigma[2]
#        sigma[3]=init.sigma[3]
           
        S.sigma<-matrix(0,ncol=total.dim.gamma,nrow=total.dim.gamma)
        psi.sigma<-rep(0,total.dim.gamma)
#        blup<-list(Z.i=rep(0,nClusters),Z.iA=matrix(0,ncol=2*nClusters,nrow=2),Q.iA=matrix(0,ncol=2*nClusters,nrow=2))
#        emp.var.aux=matrix(0,nrow=dim.beta+dim.sigma,ncol=nClusters) # i'th column: the value of (psi_beta,psi_sigma)^T
#        emp.var=matrix(0,ncol=dim.beta+dim.sigma,nrow=dim.beta+dim.sigma)               
        
    # INIT - Other Settings
      act.tol=10
      it=0
      if (missing(init.beta)){
          init.beta=list()
            beta=init.lm(fixed,data)
         }
      else   
        beta=init.beta
      phase<-1 # phase==1 : Update only beta  # # phase==1 : Update both beta and sigma
        
    # ITERATIONS    
       act.tol=1
    repeat{
#browser()
      cat("it: ",it," phase: ",phase," -  beta: ",beta,";     A, C, E: ,,    act.tol: ",act.tol,"\n")
      it=it+1
 #     browser()
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

      