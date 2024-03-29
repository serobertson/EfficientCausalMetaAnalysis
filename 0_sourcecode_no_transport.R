#####Estimators for adjusting in S=0 only (equivalent to (1-IN_TRIAL))
no_transport_est<-function(DF, trialnum=0,str_covariates,keep){
  
  #redefine variable
  DF$IN_TRIAL<-ifelse(DF$S!=trialnum, 1,0)
  
  K<-calculate_estimators_S0only(DF,trialnum,str_covariates,keep)
  
  #---------------------
  #unadjusted analysis in trial, S=s
  
  DFS0<-subset(DF, DF$S==trialnum)
  unadjusted<- glm(Y~A, data=DFS0, family="gaussian")
  summary(unadjusted)
  
  (param_start_unadjusted<-c(unlist(coef(unadjusted))))
  param_start_unadjusted<-c(param_start_unadjusted[1]+param_start_unadjusted[2],
                            param_start_unadjusted[1],
                            param_start_unadjusted[2])
  UNADJUSTED_mest<-m_estimate(
    estFUN = Unadjusted_EE_S0only,
    data  = DF,
    root_control = setup_root_control(start = param_start_unadjusted),
    compute_roots = T,
    compute_vcov = T
  )
  
  #save variance + SE for ate
  param_num_mu<-length(param_start_unadjusted)
  UNADJUSTED<-UNADJUSTED_mest@estimates[param_num_mu]
  UNADJUSTED_sandwich_se <- diag(UNADJUSTED_mest@vcov)^0.5 
  UNADJUSTED_SE<-UNADJUSTED_sandwich_se[param_num_mu]
  
  #save variance + SE for EY1
  param_num_mu<-1
  UNADJUSTEDm1<-UNADJUSTED_mest@estimates[param_num_mu]
  UNADJUSTED_sandwich_se <- diag(UNADJUSTED_mest@vcov)^0.5 
  UNADJUSTED_SEm1<-UNADJUSTED_sandwich_se[param_num_mu]
  
  #save variance + SE for EY0
  param_num_mu<-2
  UNADJUSTEDm0<-UNADJUSTED_mest@estimates[param_num_mu]
  UNADJUSTED_sandwich_se <- diag(UNADJUSTED_mest@vcov)^0.5 
  UNADJUSTED_SEm0<-UNADJUSTED_sandwich_se[param_num_mu]
  
 
  print("Returns unadjusted estimates in one dataset (one trial or OBS study)")
  m <- matrix(c(UNADJUSTED,UNADJUSTED_SE), ncol = 2, nrow = 1)
  m <- data.frame(m)
  names(m) <- c("EST", "SE")
  m$trialnum<-trialnum
  m$ESTNAME<-c("UNADJUSTED")
  m$type<-"ATE"
  rownames(m) <- NULL #print without rownames
  print(m,row.names = FALSE)
  #return(m)
  
  m1 <- matrix(c(UNADJUSTEDm1,UNADJUSTED_SEm1), ncol = 2, nrow = 1)
  m1 <- data.frame(m1)
  names(m1) <- c("EST", "SE")
  m1$trialnum<-trialnum
  m1$ESTNAME<-c("UNADJUSTED")
  m1$type<-"EY1"
  rownames(m1) <- NULL #print without rownames
  print(m1,row.names = FALSE)
  
  m0 <- matrix(c(UNADJUSTEDm0,UNADJUSTED_SEm0), ncol = 2, nrow = 1)
  m0 <- data.frame(m0)
  names(m0) <- c("EST", "SE")
  m0$trialnum<-trialnum
  m0$ESTNAME<-c("UNADJUSTED")
  m0$type<-"EY0"
  rownames(m0) <- NULL #print without rownames
  print(m0,row.names = FALSE)
  
  mfinal<-rbind(m, m1, m0)
  return(mfinal)
  
}

calculate_estimators_S0only<-function(DF,trialnum,str_covariates,keep){
  
  DFtrial<-subset(DF, DF$S==trialnum)
  #1) Calculate working models
  
  #a) E[Y|X, IN_TRIAL=1, A=a]
  
  #model specification
  Amod<-as.formula(paste("A ~ ",str_covariates,sep = ""))
  Ymod<-as.formula(paste("Y ~ ",str_covariates,sep = ""))
  
  #OUTCOME MODEL FOR TREATED
  DFtrialA1 <-subset(DF, DF$S==trialnum & DF$A==1)
  Y_mod_treated <- glm(Ymod, data=DFtrialA1, family="gaussian")
  g_A1 <-predict(Y_mod_treated, newdata=DF)  #/*g_a(X)*/
  #OUTCOME MODEL FOR UNTREATED
  DFtrialA0 <-subset(DF, DF$S==trialnum & DF$A==0)
  Y_mod_untreated <- glm(Ymod, data=DFtrialA0, family="gaussian")
  g_A0 <-predict(Y_mod_untreated, newdata=DF)  #/*g_a(X)*/
  
  # b) Pr[A=a|X, IN_TRIAL=1]
  A_mod <- glm(Amod, data=DFtrial, family="binomial")
  e_a1 <- predict(A_mod, newdata=DF, type="response")
  e_a0<-1-e_a1
  
  w_S0<-(DF$S==trialnum)*((DF$A) /(e_a1) + ( (1-DF$A)) /((1-e_a1)))
  DF$w_S0<-w_S0
  DFtrial<-subset(DF, DF$S==trialnum)
  IPW2_mod <- glm(Y ~A, data=DFtrial,weights=w_S0, family="gaussian") #-43.68, normalized IPW
  
  IPW2<- coef(IPW2_mod)["A"]
  
  #2) Calculations for estimators
  
  S0_size <- sum(DF$S==trialnum)
  
  #-----------------------------
  #estimators for treated
  #OM1
  term1 = (DF$S == trialnum) * (g_A1)
  ESTIMATE_OM1=sum(term1)/S0_size
  
  #IPW1
  w_1= (1/(e_a1))
  term2_mod=(DF$S == trialnum) *(DF$A==1) * w_1 * (DF$Y)
  ESTIMATE_IPW1=sum(term2_mod)/S0_size
  #-----------------------------
  #estimators for untreated
  #OM0
  term1 = (DF$S == trialnum) * (g_A0)
  ESTIMATE_OM0=sum(term1)/S0_size
  
  #IPW0
  w_0= (1)/(e_a0)
  term2_mod=(DF$S == trialnum) *(DF$A==0) * w_0 * (DF$Y)
  ESTIMATE_IPW0=sum(term2_mod)/S0_size
  #-----------------------------
  #ATE
  IPW=ESTIMATE_IPW1-ESTIMATE_IPW0
  OM=ESTIMATE_OM1-ESTIMATE_OM0
  
  ESTIMATES<-data.frame(OM, IPW,IPW2)
  
  #return working models (useful for starting values for geex) and estimates
  newList <- list("Y_mod_treated"=Y_mod_treated,"Y_mod_untreated"=Y_mod_untreated,
                  "A_mod"=A_mod, "IPW2_mod"=IPW2_mod,"ESTIMATES"=ESTIMATES  )
  
  return(newList)
  
}


Unadjusted_EE_S0only <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  
  
  function(theta){
    
    mu_mean1<-theta[1]
    mu_mean0<-theta[2]
    mu_mean<-theta[3]
    
    #means
    summand1<-(A*(1-IN_TRIAL)*(Y - mu_mean1))
    summand0<-((1-A)*(1-IN_TRIAL)*(Y - mu_mean0)) 
    mean<-(1-IN_TRIAL)*(( mu_mean1-mu_mean0)- mu_mean)
    
    c(summand1,summand0,mean)
    
  }
}

OM_EE_S0only <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
 
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  Xdat$IN_TRIAL<-NULL
  
  X<-data.matrix(Xdat)
  
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  
  function(theta){
    
    num_cov<- ncol(X)
    #OUTCOME MODEL PIECE 
    beta1<-theta[1: num_cov]
    beta0<-theta[(num_cov+1):(2* num_cov)]
    mu<-theta[(2* num_cov+1)]
    
    m_A1 <-X %*% beta1
    
    #E[Y|X, A=1]
    ols_A1 <-crossprod(X, ((1-IN_TRIAL)*A)*(Y - m_A1))
    
    #E[Y|X, A=0]
    m_A0 <-X %*% beta0
    ols_A0 <-crossprod(X, ((1-IN_TRIAL)*(1-A))*(Y - m_A0))
    
    mean<-(1-IN_TRIAL)*(m_A1-m_A0-mu)
    
    c(ols_A1,ols_A0, mean)
  }
}

#unnormalized version
IPW_EE_S0only <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  Xdat$IN_TRIAL<-NULL
  
  
  X<-data.matrix(Xdat)
  
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  
  
  function(theta){
    
    num_cov<- ncol(X)
    
    mu_mean<-theta[num_cov+1]
    
    #TREATMENT MODEL PIECE
    lp2  <- X %*% theta[1: num_cov]
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,(1-IN_TRIAL)*(A - pa) )
    
    w = (A * (1-IN_TRIAL)/(pa)) + ((1 - A)*(1-IN_TRIAL)/((1-pa))) 
    
    #means
    summand1<-(w*A*(1-IN_TRIAL)*Y)  
    summand0<-(w*(1-A)*(1-IN_TRIAL)*Y)  
    mean<-(summand1-summand0)- (1-IN_TRIAL)*mu_mean
    c(score_eqns2, mean)
    
    
  }
}

#normalized version IPW
IPWnorm_EE_S0only <- function(data){
  
  A<-data$A
  S<- data$S
  Y <- data$Y
  IN_TRIAL<-data$IN_TRIAL
  
  Xdat<-data.frame(int=1,data) #matrix of only Xs
  Xdat$A<-NULL
  Xdat$S<-NULL
  Xdat$Y<-NULL
  Xdat$IN_TRIAL<-NULL
  
  X<-data.matrix(Xdat)
  
  matA <- cbind(1, data$A) 
  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  
  
  function(theta){
    
    num_cov<- ncol(X)
    
    #TREATMENT MODEL PIECE
    lp2  <- X %*% theta[1: num_cov]
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,(1-IN_TRIAL)*(A - pa) )
    
    w = (A * (1-IN_TRIAL)/(pa)) + ((1 - A)*(1-IN_TRIAL)/((1-pa))) 
    
    lp3  <- matA %*% theta[(num_cov+1):(num_cov+2)] 
    score_eqns3<-crossprod(matA,(1-IN_TRIAL)*w*(Y - lp3) )
    
    c(score_eqns2, score_eqns3)
    
    
  }
}

