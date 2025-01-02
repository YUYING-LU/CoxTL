g<-function(u){
  return(1/(1+exp(-u)))
}

Generate_data<-function(N,Ns,Nt,setting,kappa_t,kappa_s,ratio_scale,p){

  mu<-rep(1,p-3)
  mu[seq(1,p-3,2)]<-0.5
  r<-0.1^seq(0,p-4)
  Sigma<-toeplitz(r)
  # generate Z based on latent W
  W<-mvrnorm(N,mu,Sigma)
  mu<-c(2,1)

  Sigma<-matrix(c(1,0,0,4),ncol=2)
  X<-cbind(mvrnorm(N,mu,Sigma),rbinom(N,1,0.4),W)

  tempr<-ratio_scale*(X[,2]-X[,6]/2)
  uu<-runif(N)
  R<-ifelse(uu<g(tempr),1,0)
  #mean(R)

  # separate the data to four subsets
  # first separate data from target (R=0) and source (R=1)
  Xs<-X[R==1,]
  Xs<-Xs[1:Ns,]
  Xt<-X[R==0,]
  Xt<-Xt[1:Nt,]

  XB_t<-0.5*Xt[,1]-2*Xt[,2]^2+Xt[,3]-Xt[,7]^2/2+3*Xt[,10]
  scale_t<-sqrt(2/kappa_t)*exp(-XB_t/2)
  shape_t<-2

  if(setting %in% c(1,3)){
    XB_s<-0.5*Xs[,1]-2*Xs[,2]^2+Xs[,3]-Xs[,7]^2/2+3*Xs[,10]
  } else if (setting %in% c(2,4)){
    XB_s<-2*Xs[,1]-2*Xs[,2]^2+Xs[,3]-Xs[,7]^2+2*Xs[,10]
  }

  if(setting %in% c(1,2)) {
    scale_s = sqrt(2/kappa_t)*exp(-XB_s/2)
    shape_s = 2
  } else if(setting %in% c(3,4)) {
    scale_s = sqrt(2/kappa_s)*exp(-XB_s/2)
    shape_s = 2
  }

  T_t <- rweibull(Nt, shape = shape_t, scale = scale_t)
  C_t <- runif(Nt, 0, 3.55)
  survTime_t <- ifelse(T_t < C_t, T_t, C_t)
  event_t <- ifelse(T_t < C_t, 1, 0)

  T_s <- rweibull(Ns, shape = shape_s, scale = scale_s)
  C_s <- runif(Ns, 0, 3.55)
  survTime_s <- ifelse(T_s < C_s, T_s, C_s)
  event_s <- ifelse(T_s < C_s, 1, 0)

  data_s=data.frame(Xs,time=survTime_s,status=event_s)
  data_t=data.frame(Xt,time=survTime_t,status=event_t)

  return(list(data_s=data_s,data_t=data_t))
}


