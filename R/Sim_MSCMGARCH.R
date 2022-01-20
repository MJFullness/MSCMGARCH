
Sim_MSCMGARCH<-function(type,n,amount,true_par,lower_bound,upper_bound,nc,cl=NULL,seed){
  
  
  amount_parallel=rep(1,amount)
 
  
  
 
  f=function(x){return(Test(type,n,x,true_par,lower_bound,upper_bound,seed,1))}
  #for replicability
  theta_list<-future_lapply(X=amount_parallel,FUN=f,future.seed=seed)
  
  list_final=unlist(theta_list[[1]])
   for(i in 2:length(theta_list)){
     list_final=rbind(list_final,unlist(theta_list[[i]]))
   }
  return(list_final)
}
