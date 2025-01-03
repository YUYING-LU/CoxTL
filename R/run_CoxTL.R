run_CoxTL<-function(data_t,data_s,weights=NULL,di,folds_num=5,lam_set=c(0.001,0.005,0.01,0.05,0.1,1:10)){
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]
  size_ratio<-num_t/num_s
  if (is.null(weights)){
    weights=rep(1,num_s)
  }
  cov<-colnames(data_t[,1:di])
  functext = paste0("cox_ini<- rms::cph(survival::Surv(time, status) ~ ", paste(cov, collapse = "+"), ",
                  data = data_s, weights = weights,x=TRUE,y=TRUE,surv=TRUE)")
  eval(parse(text = functext))
  split_set<-sample(1:folds_num,num_t,replace=TRUE)
  cv_record<-matrix(rep(0,folds_num*length(lam_set)),ncol=folds_num)
  colnames(cv_record)<-paste0('Fold:',1:folds_num)
  rownames(cv_record)<-lam_set
  for(lam_idx in 1:length(lam_set)){
    lambda<-size_ratio*lam_set[lam_idx]
    for(k in 1:folds_num){
      data_t1<-data_t[split_set==k,] # validation
      data_t2<-data_t[split_set!=k,] # train
      data_ts<-rbind(data_s,data_t2)
      functext = paste0("cox_train<- rms::cph(survival::Surv(time, status) ~ ", paste(c(cov,'strat(R)'), collapse = "+"), ",
                  data = data_ts,weights=c(lambda*weights,rep(1,dim(data_t2)[1])),x=TRUE,y=TRUE,surv=TRUE)")
      eval(parse(text = functext))
      cv_record[lam_idx,k]<-1-rcorr.cens(exp(as.matrix(data_t1[,1:di])%*%coef(cox_train)),Surv(data_t1$time,data_t1$status))['C Index']
    }
  }
  lam_best<-size_ratio*lam_set[which.max(rowMeans(cv_record))]
  functext = paste0("cox_final<- rms::cph(survival::Surv(time, status) ~ ", paste(c(cov,'strat(R)'), collapse = "+"), ",
                  data=rbind(data_s,data_t),
                 weights=c(lam_best*weights,rep(1,num_t)),x=TRUE,y=TRUE,surv=TRUE)")
  eval(parse(text = functext))
  return(cox_final)
}
