#Estimation functions
MLE_MSCMGARCH<-function(r, type, start_val=NULL, BEKK1=NULL,nc){
  if(ncol(r)>2){
    return("Wrong dimensions, series must be of dimension 2 \n")
  }
  if(type=="Clayton Gumbel"){
    
    if(is.null(BEKK1)){
      #claculate BEKK
      spec=bekk_spec()
      BEKK1=bekk_fit(spec,r)
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
    }
      if(is.null(start_val)){
        set.seed(4566)
        start_val=numeric(3)
        start_val[1]=runif(1,0,1)
        start_val[2]=runif(1,0,1)
        start_val[3]=runif(1,0,17)
        start_val[4]=runif(1,1,17)
        start_val=list(start_val)
        for(i in 2:nc){
          start_val_temp=numeric(3)
          start_val_temp[1]=runif(1,0,1)
          start_val_temp[2]=runif(1,0,1)
          start_val_temp[3]=runif(1,0,17)
          start_val_temp[4]=runif(1,1,17)
          start_val[[i]]<-start_val_temp
        }
        
        f=function(par){return(random_grid_search_clayton_gumbel(BEKK1,par,r,1))}
        theta<-parLapplyLB(X=start_val,fun=f)
        max_index <- which.max(sapply(theta, '[[', 'best_val'))
        theta <- theta[[max_index]]
        start_val <- as.vector(theta[[1]])
      }
      cat("optimize")
      result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,control=list(fnscale=-1),bekk=BEKK1,r=r, lower=c(0,0,0,1) ,upper=c(1,1,17,17))
      
    
    return(list(BEKK1,result))
    
  }
  
if(type=="Normal Clayton"){
  if(is.null(BEKK1)){
    #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
    
    BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  if(is.null(start_val)){
    set.seed(1)
    start_val=numeric(3)
    start_val[1]=0.9
    start_val[2]=0.2
    start_val[3]=12
    start_val=list(start_val)
    for(i in 2:nc){
      start_val_temp=numeric(3)
      start_val_temp[1]=runif(1,0,1)
      start_val_temp[2]=runif(1,0.1,1)
      start_val_temp[3]=runif(1,0,17)
      start_val[[i]]<-start_val_temp
    }
    if(nc>1){
    f=function(par){return(random_grid_search_normal_clayton(BEKK1,par,r,1))}
    theta<-parLapplyLB(X=start_val,fun=f)
    max_index <- which.max(sapply(theta, '[[', 'best_val'))
    theta <- theta[[max_index]]
    start_val <- as.vector(theta[[1]])
    } else{
      start_val= as.vector(random_grid_search_normal_clayton(BEKK1,start_val[[1]],r,1)[[1]])
      
    }
  }
  
  
  f<-function(par){return(loglike_Normal_Clayton(BEKK1,par,r))}
  #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
  result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
  Filter_Prob=FilterProbs_Normal_Clayton(BEKK1,result$par,r)

return(list(BEKK1,result,Filter_Prob))

}

if(type=="Normal Gumbel"){
  if(is.null(BEKK1)){
  #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
  
  BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  if(is.null(start_val)){
    start_val=numeric(3)
    start_val[1]=0.9
    start_val[2]=0.2
    start_val[3]=12
    start_val=list(start_val)
    for(i in 2:nc){
      start_val_temp=numeric(3)
      start_val_temp[1]=runif(1,0.5,1)
      start_val_temp[2]=runif(1,0.1,1)
      start_val_temp[3]=runif(1,1,17)
      start_val[[i]]<-start_val_temp
    }
    if(nc>1){
    f=function(par){return(random_grid_search_normal_gumbel(BEKK1,par,r,1))}
    theta<-parLapplyLB(X=start_val,fun=f)
    max_index <- which.max(sapply(theta, '[[', 'best_val'))
    theta <- theta[[max_index]]
    start_val <- as.vector(theta[[1]])
    }else{
      start_val = random_grid_search_normal_gumbel(BEKK1,start_val[[1]],r,1)[[1]]
    }
  }
  
  
  #f<-function(par){return(loglike_Normal_Gumbel(BEKK1,par,r))}
  f<-function(par){return(loglike_Normal_Gumbel(BEKK1,par,r))}
  #result<-optim(par=as.vector(start_val),fn=f, control=list(fnscale=-1), method = "BFGS" )
  #result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-100))
   result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
  Filter_Prob=FilterProbs_Normal_Gumbel(BEKK1,result$par,r)
  
  return(list(BEKK1,result,Filter_Prob))

}
if(type=="Normal Frank"){
    if(is.null(BEKK1)){
      #claculate BEKK
      spec=bekk_spec()
      BEKK1=bekk_fit(spec,r)
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
    }
    if(is.null(start_val)){
      set.seed(45)
      start_val=numeric(3)
      start_val[1]=runif(1,0,1)
      start_val[2]=runif(1,0,1)
      start_val[3]=runif(1,-35,35)
      start_val=list(start_val)
      for(i in 2:nc){
        start_val_temp=numeric(3)
        start_val_temp[1]=runif(1,0,1)
        start_val_temp[2]=runif(1,0,1)
        start_val_temp[3]=runif(1,-35,35)
        start_val[[i]]<-start_val_temp
      }
      if(nc>1){
      f=function(par){return(random_grid_search_normal_frank(BEKK1,par,r,1))}
      theta<-parLapplyLB(X=start_val,fun=f)
      max_index <- which.max(sapply(theta, '[[', 'best_val'))
      theta <- theta[[max_index]]
      start_val <- as.vector(theta[[1]])
      }else{
        start_val=random_grid_search_normal_frank(BEKK1,start_val[[1]],r,1)[[1]]
      }
    }
    
    
    f<-function(par){return(loglike_Normal_Frank(BEKK1,par,r))}
    result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
    
    
    
    Filter_Prob=FilterProbs_Normal_Frank(BEKK1,result$par,r)
    
    return(list(BEKK1,result,Filter_Prob))
    
  }
if(type=="Normal Clayton Survival"){
  if(is.null(BEKK1)){
    #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
    
    BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  if(is.null(start_val)){
    set.seed(1)
    start_val=numeric(3)
    start_val[1]=0.9
    start_val[2]=0.2
    start_val[3]=12
    start_val=list(start_val)
    for(i in 2:nc){
      start_val_temp=numeric(3)
      start_val_temp[1]=runif(1,0.5,1)
      start_val_temp[2]=runif(1,0.1,1)
      start_val_temp[3]=runif(1,1,17)
      start_val[[i]]<-start_val_temp
    }
    if(nc>1){
    f=function(par){return(random_grid_search_normal_claytonS(BEKK1,par,r,1))}
    theta<-parLapplyLB(X=start_val,fun=f)
    max_index <- which.max(sapply(theta, '[[', 'best_val'))
    theta <- theta[[max_index]]
    start_val <- as.vector(theta[[1]])
    }else{
      start_val=random_grid_search_normal_claytonS(BEKK1,start_val[[1]],r,1)[[1]]
    }
  }
  
  
   f<-function(par){return(loglike_Normal_ClaytonS(BEKK1,par,r))}
  #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-100), method = "BFGS", hessian = T)
   result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
  

  Filter_Prob=FilterProbs_Normal_ClaytonS(BEKK1,result$par,r)
  
  return(list(BEKK1,result,Filter_Prob))

}
if(type=="Normal Gumbel Survival"){
  if(is.null(BEKK1)){
    #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
    
    BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  if(is.null(start_val)){
    set.seed(4630)
    start_val=numeric(3)
    start_val[1]=0.9
    start_val[2]=0.2
    start_val[3]=12
    start_val=list(start_val)
    for(i in 2:nc){
      start_val_temp=numeric(3)
      start_val_temp[1]=runif(1,0.5,1)
      start_val_temp[2]=runif(1,0.1,1)
      start_val_temp[3]=runif(1,1,12)
      start_val[[i]]<-start_val_temp
    }
    if(nc>1){
    f=function(par){return(random_grid_search_normal_gumbelS(BEKK1,par,r,1))}
    theta<-parLapplyLB(X=start_val,fun=f)
    max_index <- which.max(sapply(theta, '[[', 'best_val'))
    theta <- theta[[max_index]]
    start_val <- as.vector(theta[[1]])
  }else{
    start_val=random_grid_search_normal_gumbelS(BEKK1,start_val[[1]],r,1)[[1]]
  }
  }
  
  
  
  f<-function(par){return(loglike_Normal_GumbelS(bekk=BEKK1,theta=par,r=r))}
  
  #result<-optim(par=as.vector(start_val),fn=f,method="BFGS",hessian = T)
   result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))

  Filter_Prob=FilterProbs_Normal_GumbelS(BEKK1,result$par,r)
  
  return(list(BEKK1,result,Filter_Prob))

}
if(type=="Normal Gumbel GumbelSurvival"){
  if(is.null(BEKK1)){
  #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
  
  BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  cat("Generate starting values for MS \n")
  if(is.null(start_val)){
    set.seed(4602)
    start_val=numeric(8)
    start_val[1]=0.8
    start_val[2]=0.1
    start_val[3]=0.2
    start_val[4]=0.6
    start_val[5]=0.1
    start_val[6]=0.4
    start_val[7]=12
    start_val[8]=12
    if(nc>1){
      start_val=list(start_val)
    for(i in 2:nc){
    start_val_temp=numeric(8)
    start_val_temp[1]=runif(1,0,1)
    start_val_temp[2]=runif(1,0,1-start_val_temp[1])
    start_val_temp[3]=runif(1,0,1)
    start_val_temp[4]=runif(1,0,1-start_val_temp[3])
    start_val_temp[5]=runif(1,0,1)
    start_val_temp[6]=runif(1,0,1-start_val_temp[5])
    start_val_temp[7]=runif(1,10,12)
    start_val_temp[8]=runif(1,10,12)
    
    start_val[[i]]<-start_val_temp
    }
      f=function(par){return(random_grid_search_normal_gumbel_gumbelsurvival(BEKK1,par,r))}
      theta<-parLapply(X=start_val,fun=f)
      max_index <- which.max(sapply(theta, '[[', 'best_val'))
      theta <- theta[[max_index]]
      start_val <- as.vector(theta[[1]])
    }    else{
    start_val=random_grid_search_normal_gumbel_gumbelsurvival(BEKK1,start_val,r)$theta_optim
    }
  }
 
  cat("Optimize")
  f=function(par){return(loglike_Normal_Gumbel_GumbelSurvival(bekk=BEKK1,theta=as.matrix(par),r=r))}
  
  result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
  
  Filter_Prob=FilterProbs_normal_gumbel_gumbelS(BEKK1,result$par,r)

return(list(BEKK1,result,Filter_Prob))

}
if(type=="Normal GumbelSurvival ClaytonSurvival")  {
  
  if(is.null(BEKK1)){
    #claculate BEKK
    spec=bekk_spec()
    BEKK1=bekk_fit(spec,r)
    
    BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
  }
  cat("Generate starting values for MS \n")
  if(is.null(start_val)){
    set.seed(4602)
    start_val=numeric(8)
    start_val[1]=0.8
    start_val[2]=0.1
    start_val[3]=0.2
    start_val[4]=0.6
    start_val[5]=0.1
    start_val[6]=0.4
    start_val[7]=12
    start_val[8]=12
    if(nc>1){
      start_val=list(start_val)
      for(i in 2:nc){
        start_val_temp=numeric(8)
        start_val_temp[1]=runif(1,0,1)
        start_val_temp[2]=runif(1,0,1-start_val_temp[1])
        start_val_temp[3]=runif(1,0,1)
        start_val_temp[4]=runif(1,0,1-start_val_temp[3])
        start_val_temp[5]=runif(1,0,1)
        start_val_temp[6]=runif(1,0,1-start_val_temp[5])
        start_val_temp[7]=runif(1,10,17)
        start_val_temp[8]=runif(1,10,17)
        
        start_val[[i]]<-start_val_temp
      }
      f=function(par){return(random_grid_search_normal_gumbelS_claytonS(BEKK1,par,r))}
      theta<-parLapply(X=start_val,fun=f)
      max_index <- which.max(sapply(theta, '[[', 'best_val'))
      theta <- theta[[max_index]]
      start_val <- as.vector(theta[[1]])
    }    else{
      start_val=random_grid_search_normal_gumbelS_claytonS(BEKK1,start_val,r)$theta_optim
    }
  }
  
  cat("Optimize")
  f=function(par){return(loglike_Normal_GumbelS_ClaytonS(bekk=BEKK1,theta=as.matrix(par),r=r))}
  
  result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
  
  Filter_Prob=FilterProbs_normal_gumbelS_claytonS(BEKK1,result$par,r)
  
  return(list(BEKK1,result,Filter_Prob))
  
  
}
  if(type=="Normal Gumbel Clayton")  {
    
    if(is.null(BEKK1)){
      #claculate BEKK
      spec=bekk_spec()
      BEKK1=bekk_fit(spec,r)
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
    }
    cat("Generate starting values for MS \n")
    if(is.null(start_val)){
      set.seed(4602)
      start_val=numeric(8)
      start_val[1]=0.8
      start_val[2]=0.1
      start_val[3]=0.2
      start_val[4]=0.6
      start_val[5]=0.1
      start_val[6]=0.4
      start_val[7]=12
      start_val[8]=12
      if(nc>1){
        start_val=list(start_val)
        for(i in 2:nc){
          start_val_temp=numeric(8)
          start_val_temp[1]=runif(1,0,1)
          start_val_temp[2]=runif(1,0,1-start_val_temp[1])
          start_val_temp[3]=runif(1,0,1)
          start_val_temp[4]=runif(1,0,1-start_val_temp[3])
          start_val_temp[5]=runif(1,0,1)
          start_val_temp[6]=runif(1,0,1-start_val_temp[5])
          start_val_temp[7]=runif(1,10,17)
          start_val_temp[8]=runif(1,10,17)
          
          start_val[[i]]<-start_val_temp
        }
        f=function(par){return(random_grid_search_normal_gumbel_clayton(BEKK1,par,r))}
        theta<-parLapply(X=start_val,fun=f)
        max_index <- which.max(sapply(theta, '[[', 'best_val'))
        theta <- theta[[max_index]]
        start_val <- as.vector(theta[[1]])
      }    else{
        start_val=random_grid_search_normal_gumbel_clayton(BEKK1,start_val,r)$theta_optim
      }
    }
    
    cat("Optimize")
    f=function(par){return(loglike_Normal_Gumbel_Clayton(bekk=BEKK1,theta=as.matrix(par),r=r))}
    
    result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
    
    Filter_Prob=FilterProbs_normal_gumbel_clayton(BEKK1,result$par,r)
    
    return(list(BEKK1,result,Filter_Prob))
    
    
  }
  if(type=="Normal Clayton ClaytonSurvival"){
    if(is.null(BEKK1)){
      #claculate BEKK
      spec=bekk_spec()
      BEKK1=bekk_fit(spec,r)
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
    }
    cat("Generate starting values for MS")
    if(is.null(start_val)){
      set.seed(42)
      start_val=numeric(8)
      start_val[1]=0.8
      start_val[2]=0.1
      start_val[3]=0.2
      start_val[4]=0.6
      start_val[5]=0.1
      start_val[6]=0.4
      start_val[7]=12
      start_val[8]=12
      if(nc>1){
        start_val=list(start_val)
        for(i in 2:nc){
          start_val_temp=numeric(8)
          start_val_temp[1]=runif(1,0.5,1)
          start_val_temp[2]=runif(1,0,1-start_val_temp[1])
          start_val_temp[3]=runif(1,0.1,1)
          start_val_temp[4]=runif(1,0,1-start_val_temp[3])
          start_val_temp[5]=runif(1,0.1,1)
          start_val_temp[6]=runif(1,0,1-start_val_temp[5])
          start_val_temp[7]=runif(1,1,17)
          start_val_temp[8]=runif(1,1,17)
          
          start_val[[i]]<-start_val_temp
        }
      
      
      f=function(par){return(random_grid_search_normal_clayton_claytonsurvival(BEKK1,par,r))}
      theta<-parLapplyLB(X=start_val,fun=f)
      max_index <- which.max(sapply(theta, '[[', 'best_val'))
      theta <- theta[[max_index]]
      start_val <- as.vector(theta[[1]])
    }else{
      start_val=random_grid_search_normal_clayton_claytonsurvival(BEKK1,start_val,r)$theta_optim
      
    }
    }
    
    cat("Optimize")
    f=function(par){return(loglike_Normal_Clayton_ClaytonSurvival(bekk=BEKK1,theta=as.matrix(par),r=r))}
    
    #result_MS<-optim(par=start_val,fn=g,control=list(fnscale=-1),method="BFGS",hessian = T)#,lower=c(0,0,0,0,0,0,0,1),upper=c(1,1,1,1,1,1,17,17),hessian = T)
    #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-100), method = "BFGS", hessian = T)
    result<-optim(par=start_val,fn=f , method = "BFGS",control=list(fnscale=-1))
    Filter_Prob=FilterProbs_normal_clayton_claytonS(BEKK1,result$par,r)
    
    return(list(BEKK1,result,Filter_Prob))
    
  }  

}

MLE_MSCMGARCH_3dim<-function(r, type, start_val=NULL, BEKK=NULL,asymmBEKK=FALSE,signs=NULL,nc){
  if(ncol(r)!=3){
    return("Wrong dimensions, series must be of dimension 3 \n")
  }
  if(nrow(type)==3){
  
    
  
   if(asymmBEKK==FALSE){
      #claculate BEKK
     if(is.null(BEKK)){
       spec=bekk_spec()
      BEKK1=bekk_fit(spec,r)
      
      BEKK=BEKK1$theta
     }
  
  
  
    if(is.null(start_val)){
      set.seed(456)
      start_val=numeric(5)
      start_val[1]=0.92
      start_val[2]=0.2
      start_val[3]=12
      start_val[4]=12
      start_val[5]=12
      start_val=list(start_val)
      if(nc>1){
      for(i in 2:nc){
        start_val_temp=numeric(5)
        start_val_temp[1]=runif(1,0.5,1)
        start_val_temp[2]=runif(1,0.1,1)
        start_val_temp[3]=runif(1,1,26)
        start_val_temp[4]=runif(1,1,26)
        start_val_temp[5]=runif(1,1,26)
        start_val[[i]]<-start_val_temp
      }
      
      f=function(par){return(random_grid_search_normal_copula_3(BEKK,par,type,r,1))}
      theta<-pblapply(X=start_val,FUN=f)
      max_index <- which.max(sapply(theta, '[[', 'best_val'))
      theta <- theta[[max_index]]
      start_val <- as.vector(theta[[1]])
    } else{
      f=function(par){return(random_grid_search_normal_copula_3(BEKK,par,type,r,1))}
      start_val=f(start_val[[1]])[[1]]
    }
    } 
   
    f<-function(par){return(loglike_Normal_Copula_3(BEKK,par, r,type))}
    #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
    
    result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
    Filter_Prob=FilterProbs_Normal_Copula_3(BEKK,result$par,r,type)
    
    return(list(BEKK,result,Filter_Prob))
    
  
  
  
   }
    if(asymmBEKK==TRUE){
      #claculate BEKK
      if(is.null(BEKK)){
      set.seed(1)
      spec=bekk_spec(mode=list(type="bekk",asymmetric=T),signs=signs,init_values="random")
      BEKK1=bekk_fit(spec,r)
      
      BEKK=BEKK1$theta
        }
      
      if(is.null(start_val)){
        set.seed(456)
        start_val=numeric(5)
        start_val[1]=0.92
        start_val[2]=0.2
        start_val[3]=12
        start_val[4]=12
        start_val[5]=12
        start_val=list(start_val)
        if(nc>1){
          for(i in 2:nc){
            start_val_temp=numeric(5)
            start_val_temp[1]=runif(1,0.5,1)
            start_val_temp[2]=runif(1,0.1,1)
            start_val_temp[3]=runif(1,1,26)
            start_val_temp[4]=runif(1,1,26)
            start_val_temp[5]=runif(1,1,26)
            start_val[[i]]<-start_val_temp
          }
          
          f=function(par){return(random_grid_search_normal_copula_3_asymm(BEKK,signs,par,type,r,1))}
          theta<-pblapply(X=start_val,FUN=f)
          max_index <- which.max(sapply(theta, '[[', 'best_val'))
          theta <- theta[[max_index]]
          start_val <- as.vector(theta[[1]])
        } else{
          f=function(par){return(random_grid_search_normal_copula_3_asymm(BEKK,signs,par,type,r,1))}
          start_val=f(start_val[[1]])[[1]]
        }
      } 
      
      
      f<-function(par){return(loglike_Normal_Copula_3_asymm(BEKK,signs,par, r,type))}
      #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
      
      result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
      Filter_Prob=FilterProbs_Normal_Copula_3_asymm(BEKK,signs,result$par,r,type)
      
      return(list(BEKK,result,Filter_Prob))
    }
  }
  else if(nrow(type)==6){
    
      
      
      if(asymmBEKK==FALSE){
        #claculate BEKK
        if(is.null(BEKK)){
        spec=bekk_spec()
        BEKK1=bekk_fit(spec,r)
        
        BEKK=BEKK1$theta
        }
        
        
        
        if(is.null(start_val)){
         
          set.seed(456)
          start_val=numeric(12)
          start_val[1]=0.9
          start_val[2]=0.05
          start_val[3]=0.9
          start_val[4]=0.05
          start_val[5]=0.05
          start_val[6]=0.05
          start_val[7]=4
          start_val[8]=3
          start_val[9]=3
          start_val[10]=4
          start_val[11]=3
          start_val[12]=3
          if(nc>1){
          start_val=list(start_val)
          for(i in 2:nc){
            start_val_temp=numeric(12)
            start_val_temp[1]=runif(1,0,1)
            start_val_temp[2]=runif(1,0,1-start_val_temp[1])
            start_val_temp[3]=runif(1,0,1)
            start_val_temp[4]=runif(1,0,1-start_val_temp[3])
            start_val_temp[5]=runif(1,0,1)
            start_val_temp[6]=runif(1,0,1-start_val_temp[5])
            start_val_temp[7]=runif(1,2,31)
            start_val_temp[8]=runif(1,2,31)
            start_val_temp[9]=runif(1,2,31)
            start_val_temp[10]=runif(1,2,31)
            start_val_temp[11]=runif(1,2,31)
            start_val_temp[12]=runif(1,2,31)
            start_val[[i]]<-start_val_temp
          }
          
          f=function(par){return(random_grid_search_normal_copula_copula_3(BEKK,par,type,r,1))}
          theta<-lapply(X=start_val,FUN=f)
          max_index <- which.max(sapply(theta, '[[', 'best_val'))
          theta <- theta[[max_index]]
          start_val <- as.vector(theta[[1]])
          }else{
            start_val=random_grid_search_normal_copula_copula_3(BEKK,start_val,type,r,1)[[1]]
        }
        }
        
        f<-function(par){return(loglike_Normal_Copula_Copula_3(BEKK,par, r,type))}
        #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
        
        result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
        Filter_Prob=FilterProbs_Normal_Copula_Copula_3(BEKK,result$par,r,type)
        
        return(list(BEKK,result,Filter_Prob))
        
        
        
        
      }
      if(asymmBEKK==TRUE){
        #claculate BEKK
        if(is.null(BEKK)){
          set.seed(1)
        spec=bekk_spec(model=list(type="bekk",asymmetric=T),init_values = "random",signs=signs)
        BEKK1=bekk_fit(spec,r)
        BEKK=BEKK1$theta
        }
        if(is.null(start_val)){
          set.seed(456)
          
          start_val=numeric(12)
          start_val[1]=0.9
          start_val[2]=0.05
          start_val[3]=0.9
          start_val[4]=0.05
          start_val[5]=0.05
          start_val[6]=0.05
          start_val[7]=4
          start_val[8]=3
          start_val[9]=3
          start_val[10]=4
          start_val[11]=3
          start_val[12]=3
          if(nc>1){
          start_val=list(start_val)
          for(i in 2:nc){
            start_val_temp=numeric(12)
            start_val_temp[1]=runif(1,0,1)
            start_val_temp[2]=runif(1,0,1-start_val_temp[1])
            start_val_temp[3]=runif(1,0,1)
            start_val_temp[4]=runif(1,0,1-start_val_temp[3])
            start_val_temp[5]=runif(1,0,1)
            start_val_temp[6]=runif(1,0,1-start_val_temp[5])
            start_val_temp[7]=runif(1,1,16)
            start_val_temp[8]=runif(1,1,16)
            start_val_temp[9]=runif(1,1,16)
            start_val_temp[10]=runif(1,1,16)
            start_val_temp[11]=runif(1,1,16)
            start_val_temp[12]=runif(1,1,16)
            start_val[[i]]<-start_val_temp
          }
          
          f=function(par){return(random_grid_search_normal_copula_copula_3_asymm(BEKK,signs,par,type,r,1))}
          theta<-parLapplyLB(X=start_val,fun=f)
          max_index <- which.max(sapply(theta, '[[', 'best_val'))
          theta <- theta[[max_index]]
          start_val <- as.vector(theta[[1]])
          } else{
            start_val=random_grid_search_normal_copula_copula_3_asymm(BEKK,signs,start_val,type,r,1)[[1]]
        }
          
        
        f<-function(par){return(loglike_Normal_Copula_Copula_3_asymm(BEKK,signs,par, r,type))}
        #result<-optim(par=as.vector(start_val),fn=f,control=list(fnscale=-1), method = "BFGS", hessian = T)
        print(f(start_val)) 
        result<-optim(par=start_val,fn=f , method = "BFGS",control = list(fnscale=-1))
        Filter_Prob=FilterProbs_Normal_Copula_Copula_3_asymm(BEKK,signs,result$par,r, type)
        
        return(list(BEKK,result,Filter_Prob))
      }
    }
}
}
MLE_CMGARCH<-function(r, type, start_val,signs=NULL, hessian=T){
  if(!is.numeric(type)){
  if(type=="Gumbel"){
    f=function(par){return(loglike_LL_Gumbel(par,r))}
    res=maxBFGS(f,start=start_val)
    return(res)
  }
  if(type=="Gumbel Survival"){
    f=function(par){return(loglike_LL_GumbelS(par,r))}
    res=maxBFGS(f,start=start_val)
    return(res)
  }
  if(type=="Clayton"){
    f=function(par){return(loglike_LL_Clayton(par,r))}
    res=maxBFGS(f,start=start_val)
    return(res)
  }
  if(type=="Clayton Survival"){
    f=function(par){return(loglike_LL_ClaytonS(par,r))}
    #res=optimParallel(par=start_val,fn=f, lower=c(rep(-1,11),0), upper=c(rep(1,11),17))
    res=maxBFGS(f,start=start_val)
    return(res)
  }
  if(type=="Frank"){
    f=function(par){return(loglike_LL_Frank(par,r))}
    res=maxBFGS(f,start=start_val)
    return(res)
  }
  }
   if(nrow(type)==3 && is.null(signs)){
    
    f=function(par){return(-LL_loglike_Copula_3(par,r,type))}
    start_val=random_grid_search_LL_copula_3(start_val,r,type,1)$theta_optim
    res=optim(fn=f,par=start_val, hessian = hessian)
    
    return(res)
  }
  if(nrow(type)==3 && !is.null(signs)){
    f=function(par){return(-LL_loglike_Copula_3_asymm(par,signs,r,type))}
    start_val=random_grid_search_LL_copula_3_asymm(start_val,signs,r,type,1)$theta_optim
    res=optim(fn=f,par=start_val, hessian = hessian)
    
    return(res)
  }
  
  
  
}  
