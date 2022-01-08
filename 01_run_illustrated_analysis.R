#complete-case dataset for illustration

DF=read.csv("data_illustrated_example.csv")

#********************************************************
#********************************************************

#-------sandwich variance------#
library("geex")
#--------------------------------------------------------------------
#No transportability analysis:
#Returns unadjusted estimate and adjusted estimates from confounding from OM and IPW in the target"
#--------------------------------------------------------------------
source('~/0_sourcecode_no_transport.R')
Sresults = list()
#trial 0
Sresults[[1]]<-no_transport_est(DF=DF,trialnum=0,
                                str_covariates="X1+X2+X3",
                                keep=c("X1", "X2", "X3"))

#write.csv(Sresults[[1]], "tech_results_no_transport_illustrated_example.csv")

#--------------------------------------------------------------------
#Transporting the entire collection of trials to S=0
#--------------------------------------------------------------------
source('~/0_sourcecode_transport_one_trial.R')
head(DF)
DF$S<-ifelse(DF$S==0, 0, 1) #R
meta<-transport_onetrial_est(DF=DF,trialnum=1,
                             str_covariates="X1+X2+X3",
                             keep=c("X1", "X2", "X3"))

#write.csv(meta, "tech_results_transport_collection_illustrated_example.csv")
