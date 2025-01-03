Tutorial for CoxTL
================
Yuying Lu
2025-1-1

This is a tutorial about how to use the CoxTL method.

If you haven’t install ‘CoxTL’ package, please use the following code:

``` r
library(devtools)
install_github("YUYING-LU/CoxTL")
```

## Import Necessary Packages

``` r
library(CoxTL)
library(tidyverse)
library(survival)
library(pec)
library(MASS)
library(Matrix)
library(rms)
library(survex)
```

## Data Generation

In our simulation part, we use function `Generate_data()` for data
generation. Here is an example of data generation.

``` r
set.seed(1)

N<-20000 # Total samples for generation
Ns<-4000 # sample size of source data
Nt<-200 # sample size of training target data
nt<-200 # sample size of test target data
N_orac<-5000 # sample size of target data that used in `Cox_Orac`

kt<-4 # coefficient of baseline hazard function for target sample 
      # also that for source sample with setting = 1, 2
ks<-1/2 # coefficient of baseline hazard function for source sample with setting = 3, 4

ratio_scale<-0.5 # 0: without covariate shift, 0.5: with covariate shift
                 # can be assigned with other values to mimic different degree of covariate shift 

p<-50 # dimension
setting <- 4 # setting can be 1, 2, 3 or 4


sim<-Generate_data(N,Ns,Nt+nt+N_orac,setting,kt,ks,ratio_scale,p)
data_test<-sim$data_t[(Nt+1):(Nt+nt),]
data_test$R<-as.factor(rep(0,nt))
data_orac<-sim$data_t[(Nt+nt+1):(Nt+nt+N_orac),]
data_orac$R<-as.factor(rep(0,N_orac))
data_both<-rbind(sim$data_s,sim$data_t[1:Nt,]) #combind source and target data
data_both$R<-as.factor(c(rep(1,Ns),rep(0,Nt)))
data_s<-data_both[1:Ns,]
data_t<-data_both[(Ns+1):(Nt+Ns),]

X_s<-as.matrix(data_s[1:Ns,1:p])
X_t<-as.matrix(data_t[1:Nt,1:p])
X_test<-as.matrix(data_test[,1:p])
cov<-colnames(data_t)[1:p]
```

## Model Fitting

Here we show how to apply `Cox_Orac`, `Cox_T`, `Cox_S`, `Cox_Str` and
our proposed `CoxTL` method to the simulated data. And we also calculate
the C-index and IBS for each method and compare their performance.

``` r
# Cox_Orac: use a large number of target samples 
functext = paste0("cox_orac<- rms::cph(survival::Surv(time, status) ~ ", paste(cov, collapse = "+"), ", 
                  data = data_orac,x=TRUE,y=TRUE,surv=TRUE)")
eval(parse(text = functext))


# Cox_T: use only target training data
functext = paste0("cox_t<- rms::cph(survival::Surv(time, status) ~ ", paste(cov, collapse = "+"), ", 
                  data = data_t,x=TRUE,y=TRUE,surv=TRUE)")
eval(parse(text = functext))


# Cox_S: use only source data
functext = paste0("cox_s<- rms::cph(survival::Surv(time, status) ~ ", paste(cov, collapse = "+"), ", 
                  data = data_s,x=TRUE,y=TRUE,surv=TRUE)")
eval(parse(text = functext))


# Cox_Str: stratified method use the target training data and the source data
functext = paste0("cox_str<- rms::cph(survival::Surv(time, status) ~ ", paste(c(cov,'strat(R)'), collapse = "+"), ", 
                  data = data_both,x=TRUE,y=TRUE,surv=TRUE)")
eval(parse(text = functext))

# CoxTL: our proposed method

## Step 1: estimate density ratio

theta_ini<-rep(0,p+1) # set initial value
Xs_new<-as.matrix(cbind(rep(1,Ns),X_s)) # add intercept term
Xt_new<-as.matrix(cbind(rep(1,Nt),X_t)) # add intercept term

w_opt <- nlm(density_opt(Xs_new,Xt_new,1), theta_ini, iterlim = 500)
theta<-w_opt$estimate
wei_ts<-exp(Xs_new%*%matrix(theta,ncol=1))

## Step 2: fit CoxTL
cox_tl<-run_CoxTL(data_t,data_s,weights=wei_ts,p)


time_test<-seq(0,1.5,0.01)
surv_test<-Surv(data_test$time,data_test$status)
surv_orac<-predictSurvProb(cox_orac,newdata=data_test,times=time_test)
surv_t<-predictSurvProb(cox_t,newdata=data_test,times=time_test)
surv_s<-predictSurvProb(cox_s,newdata=data_test,times=time_test)
surv_str<-predictSurvProb(cox_str,newdata=data_test,times=time_test)
surv_tl<-predictSurvProb(cox_tl,newdata=data_test,times=time_test)


Method<-c('Cox_Orac', 'Cox_T', 'Cox_S', 'Cox_Str','CoxTL')
beta_est<-rbind(coef(cox_orac),coef(cox_t),coef(cox_s),coef(cox_str),coef(cox_tl))

rownames(beta_est)<-Method

C_index<-apply(beta_est,1,function(x){1-rcorr.cens(exp(X_test%*%x),Surv(data_test$time,data_test$status))['C Index']})

IBS<-c(integrated_brier_score(surv_test, surv = surv_orac, times = time_test),
       integrated_brier_score(surv_test, surv = surv_t, times = time_test),
       integrated_brier_score(surv_test, surv = surv_s, times = time_test),
       integrated_brier_score(surv_test, surv = surv_str, times = time_test),
       integrated_brier_score(surv_test, surv = surv_tl, times = time_test))

cbind(Method=Method,data.frame(C_index=C_index,IBS=IBS,p=paste0('p=',p),Setting=setting),eta=ratio_scale)
```

    ##            Method   C_index       IBS    p Setting eta
    ## Cox_Orac Cox_Orac 0.7735005 0.1584045 p=50       4 0.5
    ## Cox_T       Cox_T 0.6933568 0.2069622 p=50       4 0.5
    ## Cox_S       Cox_S 0.7012025 0.2015833 p=50       4 0.5
    ## Cox_Str   Cox_Str 0.7089749 0.2172876 p=50       4 0.5
    ## CoxTL       CoxTL 0.7359584 0.1716043 p=50       4 0.5
