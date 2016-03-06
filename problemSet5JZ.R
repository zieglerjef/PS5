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
anesSub$ft_dpc <- recode(anes$ft_dpc, "-2=NA; -8=NA; -9=NA")
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

# divide anesSub: trainingSet = 70% (n=4140), testSet = 30% (n=1774)
anesSample <- sample(x = seq_len(nrow(anesSub)), size=round(dim(anesSub)[1]*.70), replace=F)

# subset training set
trainingSet <- anesSub[anesSample, ]

# subset test set
testSet <- anesSub[-anesSample, ]
