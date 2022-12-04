# SOME WRITTEN VALUES MAY BE WRONG WITHIN THE COMMENTS DUE TO REWRITING CODE

library(rstan)
library(readxl)
library(posterior)
library(loo)
#This folder included the code file along with the stan files
setwd("~/notebooks/BDA")
set.seed(123)

normDate <- function(x) {
  compValue <- x[1]-x[1]%%86400
  ret <- x[1]%%86400/3600
  add <- 0
  for (i in 2:length(x)) {
    if (x[i]-86400 > compValue) {
      add <- 1
    }
    ret <- c(ret, x[i]%%86400/3600+add*24)
  }
  ret
}

ciCalc <- function(sample,ci) {
  ciVal <- (1-ci)/2
  quantile(sample,probs = c(ciVal,1-ciVal))
}

calcAcc <- function(true,sample,ci) {
  compQuant <- unname(ciCalc(sample,ci))
  (sum(true > compQuant[1] & true < compQuant[2])/length(true))
}

dflist <- list(as.data.frame(read_excel("~/notebooks/BDA/excelfiles/BayesData.xlsx",sheet = 1)))

for (i in 2:12) {
  dflist <- c(dflist,list(as.data.frame(read_excel(
    "~/notebooks/BDA/excelfiles/BayesData.xlsx",sheet = i))))
}

tslist <- dflist

for (i in 1:12) {
  tslist[[i]] <- dflist[[i]][dflist[[i]]$`TAC over 0.08`>0, ]
  tslist[[i]] <- normDate(tslist[[i]]$`timestamp`)
}

meanArr <- unlist(lapply(tslist,mean))

data_list_sep <- list(
  J=length(tslist),
  mA=meanArr
)

data_list_hier <- list(
  J=length(tslist),
  mA=meanArr,
  gI=1:12
)

# INITIALLY USED PRIORS
# SEP: mu ~ normal(24, 8); sigma ~ inv_chi_square(0.01);
# HIER: mu0 ~ normal(24,8); sigma0 ~ inv_chi_square(0.01); sigma ~ inv_chi_square(0.01);

sModelSep <- stan(file="projSep.stan",data=data_list_sep,iter=4000)
sModelHier <- stan(file="projHier.stan", data=data_list_hier,iter = 4000)

dfSep <- as.data.frame(sModelSep)
dfHier <- as.data.frame(sModelHier)

postSep <- dfSep[,c(1:12)]
postHier <- dfHier[,c(3:14)]

predSep <- dfSep[,c(25:36)]
predHier <- dfHier[,c(16:28)]

# Computing r-hat values

apply(postSep, 2, rhat)
apply(postHier, 2, rhat)

# Computing HMC diagnostics

check_hmc_diagnostics(sModelSep)
check_hmc_diagnostics(sModelHier)

# Computing ESS values

ess_basic(sModelSep)
ess_basic(sModelHier)

# The separate model seems to form too wide of an distribution for each dataset
# as all the real values are always included within the 90 percent middle
# quantile. The hierarchical model seems to more consistently form smaller
# distributions as a notable value of total values seem to be left out from the
# 90 percent middle quantile.

# Computing LOO values for further model comparasions:

loo(sModelSep)
loo(sModelHier)

# Sep elpd_loo = -56
# Hier elpd_loo -26

# While neither of the models seem to portray the datasets particularly well
# considering the notable amount of bad pareto-k values.

# The hierarchical model seems be better than the separate model just considering
# the amount of pareto-k values that are either bad or very bad.

# Plotting the pareto-k values for a visual measure

hist(pareto_k_values(loo(sModelSep)),20)
hist(pareto_k_values(loo(sModelHier)),20, main = "Histogram of pareto-k -values for the hierarchical model", xlab = "pareto-k")

# Most parameters are somewhat far from the perfection and the models currently
# consider the real data to consistently be outliers of the "real model",
# implying somewhat low accuracy.


# The hierarchical model may be further used to predict the distribution of drunkeness
# for a new person e.g. considering  possible predictive performance assessment.
# This is not possible for the separate model by default and thus it is not
# performed for the separate model.

predNewHier <- predHier$`ypred[13]`

hist(predNewHier,50)
mean(predNewHier)
quantile(predNewHier,c(0.05,0.95))

# On average it is predicted that the most likely time a new person would be
# intoxicated is at the 20.8 hour mark which equates to 20:48 in standard
# 24-hour time notation.

# Additionally, the 90 percent middle quantile would be [17.3, 24.1] (equating to
# the time interval [17:18, 0:06]) so if a new person would be drinking during
# some day, they could be measured to be drunk between the aforementioned interval
# with 90% certainty, according to the model.

# There is no proper way to further assess the performance of this predictive model
# as the whole dataset was used to create the distribution in question.


# We may consider further prior sensitivity analysis:

# Let us choose alternative priors, e.g. a weakly informative for the mean instead
# of a informative prior and a different prior for the standard deviation.

# Let us use the weakly informative prior for the mean normal(0,100) and a
# gamma(1,1) prior (separately) for the standard deviation of the model.
# to consider the sensitivity of the model.

sModelSepWIMean <- stan(file="projSepWIMean.stan",data=data_list_sep,iter=4000)
sModelHierWIMean <- stan(file="projHierWIMean.stan", data=data_list_hier,iter = 4000)
sModelSepDiffSD <- stan(file="projSepDiffSD.stan",data=data_list_sep,iter=4000)
sModelHierDiffSD <- stan(file="projHierDiffSD.stan", data=data_list_hier,iter = 4000)

postSepWI <- as.data.frame(sModelSepWIMean)[,c(1:12)]
postHierWI <- as.data.frame(sModelHierWIMean)[,c(3:14)]
postSepDSD <- as.data.frame(sModelSepDiffSD)[,c(1:12)]
postHierDSD <- as.data.frame(sModelHierDiffSD)[,c(3:14)]

# Following steps are made assuming sModelSepDiffSD is now the chosen model
# Additionally after the following sensitivity analysis, the loo values (and
# other similarly appropriate things) are calculated for the said model.

SepWIDiff <- matrix(apply(postSepWI,2,mean)-apply(postSepDSD,2,mean))
HierWIDiff <- matrix(apply(postHierWI,2,mean)-apply(postHier,2,mean))

SepDSDDiff <- matrix(apply(postSep,2,mean)-apply(postSepDSD,2,mean))
HierDSDDiff <- matrix(apply(postHierDSD,2,mean)-apply(postHier,2,mean))

SepWIDiffCI <- array(c(1,2))
HierWIDiffCI <- array(c(1,2))

SepDSDDiffCI <- array(c(1,2))
HierDSDDiffCI <- array(c(1,2))
for (i in 1:12) {
  SepWIDiffCI <- rbind(SepDSDDiffCI,unname(ciCalc(sample = postSepWI[,i],ci = 0.9))
                       -unname(ciCalc(sample = postSepDSD[,i], ci = 0.9)))
  HierWIDiffCI <- rbind(HierDSDDiffCI,unname(ciCalc(sample = postHierWI[,i],ci = 0.9))
                        -unname(ciCalc(sample = dfHier[,(2+i)],ci = 0.9)))
  SepDSDDiffCI <- rbind(SepDSDDiffCI,unname(ciCalc(sample = dfSep[,i],ci = 0.9))
                        -unname(ciCalc(sample = postSepDSD[,i],ci = 0.9)))
  HierDSDDiffCI <- rbind(HierDSDDiffCI,unname(ciCalc(sample = postHierDSD[,i],ci = 0.9))
                         -unname(ciCalc(sample = dfHier[,(2+i)],ci = 0.9)))
}

SepWIDiff <- cbind(SepWIDiff,SepWIDiffCI[-1,]); HierWIDiff <- cbind(HierWIDiff,HierWIDiffCI[-1,]);
SepDSDDiff <- cbind(SepDSDDiff,SepDSDDiffCI[-1,]); HierDSDDiff <- cbind(HierDSDDiff,HierDSDDiffCI[-1,]);

# Let us then consider the differences for means then

colnames(SepWIDiff) <- c("Mean","5%","95%")
colnames(HierWIDiff) <- c("Mean","5%","95%")
colnames(SepDSDDiff) <- c("Mean","5%","95%")
colnames(HierDSDDiff) <- c("Mean","5%","95%")

SepWIDiff; HierWIDiff; SepDSDDiff; HierDSDDiff


# /////////

dfSepNew <- as.data.frame(sModelSepDiffSD)
postSepNew <- dfSepNew[,c(1:12)]
predSepNew <- dfSepNew[,c(25:36)]
apply(postSepNew, 2, rhat)
check_hmc_diagnostics(sModelSepDiffSD)
ess_basic(sModelSepDiffSD)
loo(sModelSepDiffSD)

# Let us perform posterior predictive checks by generating a sample for each
# of the parameters and comparing a 90% middle quantile to the original dataset.

# This will be performed as a generic ratio of real values within the given
# quantile. Essentially the chosen middle quantile (0.9) is the optimal value
# for the accuracy.

# Checking this ratio is done through the calcAcc function.

sepAccuracy <- array(numeric(),c(0,1)) 
for (i in 1:12) {
  sepAccuracy <- c(sepAccuracy, calcAcc(true = tslist[[i]],sample = dfSepNew[,c(24+i)],ci = 0.9))
}
hierAccuracy <- array(numeric(),c(0,1)) 
for (i in 1:12) {
  hierAccuracy <- c(hierAccuracy, calcAcc(true = tslist[[i]],sample = dfHier[,c(16+i)],ci = 0.9))
}
sepAccuracy
hierAccuracy
