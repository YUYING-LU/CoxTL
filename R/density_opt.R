
density_opt<-function(Xs,Xt,lam){
  mm<-function(theta_0){
    theta_0<-as.matrix(theta_0,ncol=1)
    ff<-mean(exp(Xs%*%theta_0))-mean(Xt%*%theta_0)+lam*norm(theta_0,type = "2")^2
    return(ff)
  }
  return(mm)
}




