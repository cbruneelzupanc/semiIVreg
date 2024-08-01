## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", echo = TRUE, fig.retina = 2#,#fig.height=4, fig.width=5, fig.align='center'
)
if (!requireNamespace("ivreg", quietly = TRUE)) {
  install.packages("ivreg")
}
library(ivreg)
options(scipen=8)

## ----model1-------------------------------------------------------------------
# Model 1
library(semiIVreg)
N = 10000; set.seed(1234)
model_type = "homogenous"
param_error = c(1, 1.5, -0.6) # var_u, var_v, cov_uv # if homogenous
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9) # meanW0 state0, meanW1 state0, meanW0 state1, meanW1 state1, varW0, varW1, covW0W1
param_p = c(0, -0.5, 0.5, 0, 0, 0) # constant, alphaW0, alphaW1, alphaW0W1, effect of state, effect of parent educ
param_y0 = c(3.2, 0.8, 0, 0) # intercept, effect of Wd, effect of state, effect of parent educ;
param_y1 = c(3.2+0.4, 0.5, 0, 0) # the +0.2 = Average treatment effect; effect of W1, effect of state, effect of parent educ;
param_genX = c(0.4, 0, 2) # probability state=1 (instead of 0), mean_parenteduc, sd_parenteduc (parenteduc drawn as continuous)

data = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)

## ----model1true---------------------------------------------------------------
data$true_TE = data$y1 - data$y0
summary(lm(true_TE ~ w0 + w1, data))
true_param = c(param_y0[1], param_y1[1] - param_y0[1], param_y0[2], param_y1[2]); true_param
# constant y0, constant y1 - constant y0, effects of w0 on y0, effects of w1 on y1.

## ----naiveols1----------------------------------------------------------------
naive_ols = lm(y ~ d + w0 + w1, data); summary(naive_ols)

## ----selectionbias1-----------------------------------------------------------
mean(data$U1[which(data$d == 1)]); mean(data$U0[which(data$d==0)])

## ----ivreg1-------------------------------------------------------------------
library(ivreg) # remark: ivreg is not required otherwise in semiivreg, only in this vignette. 
valid_iv = ivreg(y ~ d | w0 + w1 + w0:w1, data=data); summary(valid_iv)

## ----semiiv1------------------------------------------------------------------
semiiv = semiivreg(y~d|w0|w1, data, est_method="homogenous", plotting=FALSE)
semiiv_homogenous = semiiv$estimate$est # extract the coefficients from the homogenous TE specification
summary_coeff = summary(semiiv_homogenous)
summary_coeff$coefficients[1:4,] # only print the first 4 coefficients, the other correspond to the control function of P
true_param

## ----semiivboot1, cache=TRUE--------------------------------------------------
semiivboot = semiivreg_boot(y~d|w0|w1, data, Nboot=200, est_method="homogenous", plotting=FALSE) # reduce the number of bootstrap simulation for speed;  
boot_se = semiivboot$estimate$coeff$std_error[1:4]
res = as.data.frame(cbind(summary_coeff$coefficients[1:4,1:2], boot_se)); colnames(res) = c("Estimate", "wrong analytic SE", "Bootstrapped SE")
res

## ----semiivqr-----------------------------------------------------------------
semiivqr = ivreg(y~d+I(1-d):w0 + I(d):w1|w0+w1+w0:w1, data=data)
summary(semiivqr)

## ----model2-------------------------------------------------------------------
# Model 2
N = 10000; set.seed(1234)
model_type = "homogenous"
param_error = c(1, 1.5, -0.6) # var_u, var_v, cov_uv # if homogenous
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9) # meanW0 state0, meanW1 state0, meanW0 state1, meanW1 state1, varW0, varW1, covW0W1
param_p = c(0, -0.3, 0.4, 0.3, 0, 0) # constant, alphaW0, alphaW1, alphaW0W1, effect of state, effect of parent educ
param_y0 = c(3.2, 0.8, 0, 0) # intercept, effect of Wd, effect of state, effect of parent educ;
param_y1 = c(3.2+0.4, 0.5, 0, 0) # the +0.2 = Average treatment effect; effect of W1, effect of state, effect of parent educ;
param_genX = c(0.4, 0, 2) # probability state=1 (instead of 0), mean_parenteduc, sd_parenteduc (parenteduc drawn as continuous)

data2 = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)

## ----semiivqr2----------------------------------------------------------------
semiivqr2 = ivreg(y~d+I(1-d):w0 + I(d):w1|w0+w1+w0:w1, data=data2); summary(semiivqr2)

## ----semiivreg2---------------------------------------------------------------
semiiv2 = semiivreg(y~d|w0|w1, data=data2, est_method="homogenous", plotting=FALSE)
semiiv2_homogenous = semiiv2$estimate$est # extract the coefficients from the homogenous TE specification
summary_coeff2 = summary(semiiv2_homogenous)
summary_coeff2$coefficients[1:4,] # only print the first 4 coefficients, the other correspond to the control function of P
true_param

