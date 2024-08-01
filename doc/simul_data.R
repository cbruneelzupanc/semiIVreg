## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", echo = TRUE, fig.retina = 2#,#fig.height=4, fig.width=5, fig.align='center'
)
options(scipen=8)

## ----model1-------------------------------------------------------------------
# Example of general model with heterogenous treatment effects
N = 100000; set.seed(1234)
model_type = "heterogenous"
param_error = c(1, 1, 0.6, 0.5) # var_u0, var_u1, cov_u0u1, var_cost (the mean cost = constant in D*) # if heterogenous
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9) # meanW0 state0, meanW1 state0, meanW0 state1, meanW1 state1, varW0, varW1, covW0W1
param_p = c(0, -0.7, 0.7, 0, 0, 0) # constant, alphaW0, alphaW1, alphaW0W1, effect of state, effect of parent educ
param_y0 = c(3.2, 0.8, 0, 0) # intercept, effect of Wd, effect of state, effect of parent educ;
param_y1 = c(3.2+0.4, 0.5, 0, 0) # the +0.2 = Average treatment effect; effect of W1, effect of state, effect of parent educ;
param_genX = c(0.4, 0, 2)

data = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)

## ----model2-------------------------------------------------------------------
# Model with homogenous treatment effects - not the same param_error to specify. 
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

## ----model3-------------------------------------------------------------------
# Example of generalized Roy Model where the semi-IVs are valid IVs
N = 50000; set.seed(1234)
model_type = "heterogenous"
param_error = c(1, 1, 0.6, 0.5) # var_u0, var_u1, cov_u0u1, var_cost (the mean cost = constant in D*) # if heterogenous
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9) # meanW0 state0, meanW1 state0, meanW0 state1, meanW1 state1, varW0, varW1, covW0W1
param_p = c(0, -0.7, 0.7, 0, 0, 0) # constant, alphaW0, alphaW1, alphaW0W1, effect of state, effect of parent educ
param_y0 = c(3.2, 0, 0, 0) # intercept, effect of Wd, effect of state, effect of parent educ;
param_y1 = c(3.2+0.4, 0, 0, 0) # the +0.2 = Average treatment effect; effect of W1, effect of state, effect of parent educ;
param_genX = c(0.4, 0, 2)

data = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)

param_y0[2]; # W0 is a valid IV because no direct effect on Y0
param_y1[2]; # W1 is a valid IV because no direct effect on Y1

