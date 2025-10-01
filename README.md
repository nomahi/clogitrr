
# clogitrr package


## Conditional Logistic Regression with Augmented Pseudo-Observations for Risk Ratio Estimation

Conditional logistic regression has been widely applied as a multivariable method for analyzing stratified binary outcome data. However, the resulting odds ratio estimator can only be interpreted as an approximation to the risk ratio under the rare-event assumption and cannot generally serve as a direct effect measure. To address this limitation, Noma (2025) proposed a novel approach that yields a consistent risk ratio estimator by incorporating pseudo-observations into conditional logistic regression. A key advantage of this method is that it can be implemented directly in standard statistical software for conditional logistic regression simply by modifying the dataset. This package provides computational tools for applying the proposed risk ratio estimation method within the conditional logistic regression framework.



## Installation

Please download "clogitrr_1.1-1.tar.gz" and install it by R menu: "packages" -> "Install package(s) from local files...".

Download: [please click this link](https://github.com/nomahi/clogitrr/raw/main/clogitrr_1.1-1.tar.gz)

Manual: [please click this link](https://github.com/nomahi/clogitrr/raw/main/clogitrr_1.1-1.pdf)





## Example code

```r

# "Conditional logistic regression with augmented pseudo-observations for risk ratio estimation"
#  by Hisashi Noma

#  R example code for implementing the modified conditional logistic regression analysis

# The "clogitrr" package
# GitHub webpage: https://github.com/nomahi/clogitrr/

###

# Download the R package file from the following URL:
# https://github.com/nomahi/clogitrr/raw/main/clogitrr_1.1-1.tar.gz

# Then, install the package (tar.gz format) by R menu: "packages" -> "Install package(s) from local files...".

###

library("clogitrr")				# load the "glmmrr" package
library("survival")				# load the "glmmML" package
library("medicaldata")				# load the "glmmML" package

##

# Analyzing the example dataset "indo_rct" using the modified conditional logistic regression analysis

data(indo_rct)
indo_rct$y <- as.numeric(indo_rct$outcome=="1_yes")

gm1 <- clogit(y ~ rx + age + gender + risk + pep + psphinc + amp + prophystent + strata(site), data=indo_rct)

indo_e <- adpdt(y, data=indo_rct)		# Adding pseudo-observations to the original dataset
gm2 <- clogit(y ~ rx + age + gender + risk + pep + psphinc + amp + prophystent + strata(site), data=indo_e)

scoef(gm1, eform=TRUE)			# Odds-ratio estimates
scoef(gm2, eform=TRUE)			# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)

###

library("doSNOW")					# load the "doSNOW" package
library("doParallel")				# load the "doParallel" package

cl <- makeSOCKcluster(max(detectCores()-1,1))
registerDoSNOW(cl)

B <- 1000		# number of resampling

opts <- list(progress = function(x) print(paste0(x,"th bootstrap is completed.")))

R1 <- foreach(b = 1:B, .combine = rbind, .options.snow = opts) %dopar% {

	library("clogitrr")
	library("survival")
	
	irct.b <- stboot(site, indo_rct)
	irct.t <- adpdt(y, irct.b) 

	clgt.b <- clogit(y ~ rx + age + gender + risk + pep + psphinc + amp + prophystent + strata(site), data=irct.t)		# Modified conditional logistic regression analysis with pseudo-observations
	clgt.b$coefficients

}

stopCluster(cl)

R2 <- scoef(gm2, eform=TRUE)					# Risk-ratio estimates; 95%CIs and P-values are incorrect (based on the naive model variances)

RR <- R2[,1]
cbind(RR,sumboot(R1,eform=TRUE))			# Estimates, 95%CIs and P-values of risk-ratios by bootstrap

