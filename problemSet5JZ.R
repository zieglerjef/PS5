##################################
### Calculating fit statistics ###
##################################

### Setup ###

# clear environment
rm(list=ls())
# set seed
set.seed(12435)
# load stings as factors FALSE
options(stringsAsFactors=F)

# load packages
library(foreign)
library(mice)
library(plyr)
library(Hmisc)

# read in data
anes <- read.dta("~/Google Drive/WashU/Spring2016/appliedStats/problemSets/PS5/anes_timeseries_2012_stata12.dta")

## model Obama's feeling thermometer score as function
## of Clinton's feeling thermometer score
model1 <- lm(ft_dpc ~ ft_hclinton, anes)

## make a prediction for a single observation with
## hypothetical clinton score of 77
predict(model1, data.frame(ft_hclinton=77))
## we would expect a Obama score of 71.7


## Question 1
## randomly subset the data into two partitions
## use "training set" to build at least three models 
## of Obama's feeling thermometer score
## document carefully how you deal with missingness

# create new dataframe for relevent variables
anesSub <- NULL

# recode: ft_dpc
anesSub$ft_dpc <- as.numeric(revalue(as.character(anes$ft_dpc),  replace=c("-2"=NA, "-8"=NA, "-9"=NA)))
# which(is.na(anes$ft_dpc))
# resulting missingness n=15 

# recode: dem_agegrp_iwdate_x
anesSub$age <- as.factor(revalue(anes$dem_agegrp_iwdate_x, replace=c("-2. Missing, birthdate fields left blank"=NA)))

# recode: dem_edugroup_x
anesSub$education <- as.factor(revalue(anes$dem_edugroup_x, replace=c("-9. Refused"=NA,
                                       "-8. Don't know"=NA, "-2. Missing, other not codeable to 1-5"=NA)))

# recode: pid_self
anesSub$partyID <- as.factor(revalue(anes$pid_self, replace=c("-9. Refused"=NA,
                                                               "-8. Don't know"=NA)))

# recode: dem_racecps_black
anesSub$black <- as.factor(anes$dem_racecps_black)

# recode: dem_hisp
anesSub$hispanic <- as.factor(revalue(anes$dem_hisp, replace=c("-9. Refused"=NA,
                                       "-8. Don't know"=NA)))

# recode: dem_hsworkg
anesSub$employed <- as.factor(anes$dem_emptype_work)

# as.data.fame
anesSub <- as.data.frame(anesSub)

#################################
### Create training/test sets ###
#################################

# divide anesSub: trainingSet = 70% (n=4140), testSet = 30% (n=1774)
anesSample <- sample(x = seq_len(nrow(anesSub)), size=round(dim(anesSub)[1]*.70), replace=F)

# subset training set
trainingSet <- anesSub[anesSample, ]

# subset test set
testSet <- anesSub[-anesSample, ]

###################
### Missingness ###
###################

### get summary of missing data ###
sumMissing <- function(df){
  # load libraries
  # return proportion missing
  cat("Proportion missing:",
      sum(is.na(df))/prod(dim(df)), "\n")
  # get missingness summary by variables
  cat("\n Missingness by variable: \n",
      colSums(is.na(df)) / nrow(df), "\n")
  cat("\n Patterns of missingness: \n")
  print(md.pattern(df))
}

# missingness: trainingSet
sumMissing(trainingSet)
# missingness: trainingSet
sumMissing(testSet)

### create training sets with:
# 1) multiple impututation and 
# 2) case-wise deletion

# multiple impututation 
trainingSet.imp <- mice(trainingSet, m=5)
testSet.imp <- mice(testSet, m=5)

# case-wise deletion
trainingSet.caseWise <- na.omit(trainingSet)
testSet.caseWise <- na.omit(testSet)

##################
### Run models ###
##################

# model 1: ethnicity
model1.imputed <- lm.mids(ft_dpc ~ black + hispanic, data = trainingSet.imp)
model1.caseWise <- lm(ft_dpc ~ black + hispanic, data = trainingSet.caseWise)

# model 2: partisan affiliation and age
model2.imputed <- lm.mids(ft_dpc ~ partyID + age, data = trainingSet.imp)
model2.caseWise <- lm(ft_dpc ~ partyID + age, data = trainingSet.caseWise)

# model 3: economic standing
model3.imputed <- lm.mids(ft_dpc ~ education + employed, data = trainingSet.imp)
model3.caseWise <- lm(ft_dpc ~ education + employed, data = trainingSet.caseWise)

################################
### retrieve model estimates ###

### function to retrieve coefficient estimates from:
# 1) lm.mids models
# 2) casewise deleted models

# initiate function
outTableFunc <- function(modelImputed){

  # if normal lm:
  if(class(modelImputed)=="lm"){
    summary(modelImputed)[4]
  }
  # otherwise indicates imputed model:
  ## NOTE: warning message is fine because the function should only
  ## read the first argument, which is FALSE for imputed models
  else{
  # param=1 gives coefficient estimates, param=2 gives SEs
  lm.mids.vals <- function(obj,param) { 
    out.mat <- NULL
    for (i in 1:obj$call1$m) 
      out.mat <- rbind(out.mat,summary.lm(obj$analyses[[i]])$coef[,param])
    out.mat
  }

  ### get three required vectors ###
  # mean coefficient estimate
  impute.coef.vec <- apply(lm.mids.vals(modelImputed,1),2,mean)
  # variation in coefficient estimates
  between.var <- apply(lm.mids.vals(modelImputed,1),2,var)
  # standard error of coefficients
  within.var <- apply(lm.mids.vals(modelImputed,2)^2,2,mean)

  ### compute standard errors ###
  m.sets <- 5
  impute.se.vec <- sqrt(within.var + ((m.sets+1)/m.sets)*between.var)

  ### adjust degrees of freedom for t-statistic ###
  # see Little & Rubin (1987), p. 257
  impute.df <- (m.sets-1)*(1 + (1/(m.sets+1)) * within.var/between.var)^2

  ### print table ###
  # construct table
  out.table <- round(cbind(impute.coef.vec,impute.se.vec,impute.coef.vec/impute.se.vec,
                           1-pt(abs(impute.coef.vec/impute.se.vec), impute.df)),6)
  # add labeling
  colnames(out.table) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
  out.table
  }
}

### retrieve model output for:
# model1
outTableFunc(model1.imputed)
outTableFunc(model1.caseWise)

# model2
outTableFunc(model2.imputed)
outTableFunc(model2.caseWise)

# model3
outTableFunc(model3.imputed)
outTableFunc(model3.caseWise)
