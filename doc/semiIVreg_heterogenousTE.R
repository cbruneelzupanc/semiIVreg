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
# Model
library(semiIVreg)
N = 50000; set.seed(1234)
# Specification
model_type = "heterogenous"
param_error = c(1, 1, 0.6, 0.5) # var_u0, var_u1, cov_u0u1, var_cost (the mean cost = constant in D*) # if heterogenous
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9) # meanW0 state0, meanW1 state0, meanW0 state1, meanW1 state1, varW0, varW1, covW0W1
param_p = c(0, -0.7, 0.7, 0, 0, 0) # constant, alphaW0, alphaW1, alphaW0W1, effect of state, effect of parent educ
param_y0 = c(3.2, 0.8, 0, 0) # intercept, effect of Wd, effect of state, effect of parent educ;
param_y1 = c(3.2+0.4, 0.5, 0, 0) # the +0.2 = Average treatment effect; effect of W1, effect of state, effect of parent educ;
param_genX = c(0.4, 0, 2)

data = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)

## ----raw, fig.height=3, fig.width=7, fig.align='center'-----------------------
semiiv = semiivreg(y~d|w0|w1, data, ref_indiv = data.frame(w0=0, w1=0))

## ----mtr, fig.height=4, fig.width=5, fig.align='center'-----------------------
mte_plot = semiiv$plot$mte; 
mtr_plot = semiiv$plot$mtr;
mtr_plot

## ----direct, fig.height=4, fig.width=5, fig.align='center'--------------------
summary(semiiv$estimate$est0)
summary(semiiv$estimate$est1);
# To be compared with:
param_y0[2]; param_y1[2]

## ----smooth, fig.height=3, fig.width=7, fig.align='center'--------------------
semiiv$bw
semiiv = semiivreg(y~d|w0|w1, data, ref_indiv = data.frame(w0=0, w1=0),
                   bw_method = "plug-in", bw_y0 = 0.20, bw_y1 = 0.20)

## ----propensity, fig.height=4, fig.width=5, fig.align='center'----------------
firststage = semiiv$estimate$propensity
# Cannot be compared with param_p directly if V gets rescaled -> but can compare the predicted P with the truth
Phat = predict(firststage, newdata=data, type="response")

summary(Phat - data$P) # almost perfect; 

## ----sieve, fig.height=3, fig.width=7, fig.align='center'---------------------
semiiv2 = semiivreg(y~d|w0|w1, data, ref_indiv = data.frame(w0=0, w1=0), 
                    est_method="sieve", pol_degree_sieve=5, 
                    plotting=FALSE)

mte_plot2 = semiiv2$plot$mte; mtr_plot2 = semiiv2$plot$mtr
grid.arrange(mte_plot2, mtr_plot2, ncol=2)


## ----sievecompare, fig.height=3, fig.width=7, fig.align='center'--------------
# If want to plot on the same plot, need some manipulation of the data

dat = semiiv$data$RES # take the original data
dat$V = dat$Phat

# for MTE:
mte_plot2 = mte_plot2 + geom_line(aes(x=V, y=mte), linetype="dashed", col="#db93a4", na.rm=TRUE, data=dat)


# for MTR need some manipulation
dat_plot = dat;
dat1 = dat_plot; dat1$mtr = dat_plot$mtr1; dat1$Treatment = 1
dat0 = dat_plot; dat0$mtr = dat_plot$mtr0; dat0$Treatment = 0;
dat2 = rbind(dat1, dat0)
dat2$Treatment = as.factor(dat2$Treatment)

mtr_plot2 = mtr_plot2 + geom_line(aes(x=V, y=mtr, col=Treatment, group=Treatment), linetype="dashed", na.rm=TRUE, data=dat2)

grid.arrange(mte_plot2, mtr_plot2, ncol=2)


## ----mte_true, fig.height=3, fig.width=7, fig.align='center'------------------
# True MTE and MTR estimations 
seq_p = seq(0, 1, by=0.001); 
w0 = 0; w1 = 0; 
sigma_V2 = param_error[1] + param_error[2] - 2*param_error[3] + param_error[4] # = var(V); var(data$V)
covU0V = param_error[1] - param_error[3] # = cov(data$U0, data$V)
covU1V = -param_error[2] + param_error[3] # = cov(data$U1, data$V)
ku0 = covU0V/sigma_V2*(qnorm(seq_p) - 0) # 0 = mean(V)
ku1 = covU1V/sigma_V2*(qnorm(seq_p) - 0) # 0 = mean(V)
true_mtr0 = param_y0[1] + param_y0[2]*w0 + ku0
true_mtr1 = param_y1[1] + param_y1[2]*w1 + ku1
true_mte = true_mtr1 - true_mtr0

newdata = data.frame(Ud=seq_p, w0=0, w1=0)
newdata$true_mte = true_mte; newdata$true_mtr1 = true_mtr1; newdata$true_mtr0 = true_mtr0

## # Remark: alternative estimation method if the truth does not have a simple closed form formula:
## data$diff = data$y1 - data$y0; pol_degree=5
## true_model_mte = lm(diff~w1 + w0 + poly(Ud, pol_degree, raw=TRUE), data); # MTE
## true_model_mtr1 = lm(y1 ~w1+poly(Ud, pol_degree, raw=TRUE), data) # MTR1
## true_model_mtr0 = lm(y0 ~w0+poly(Ud, pol_degree, raw=TRUE), data) # MTR0
## newdata$true_mte = predict(true_model_mte, newdata); 
## newdata$true_mtr1 = predict(true_model_mtr1, newdata); newdata$true_mtr0 = predict(true_model_mtr0, newdata)


# Comparison:
mte_plot2 = mte_plot2 + geom_line(data=newdata, aes(x=Ud, y=true_mte), linetype="dashed", col="red")
mtr_plot2 = mtr_plot2 + geom_line(data=newdata, aes(x=Ud, y=true_mtr1), linetype="dashed", col="blue") +
  geom_line(data=newdata, aes(x=Ud, y=true_mtr0), linetype="dashed", col="orange")

grid.arrange(mte_plot2, mtr_plot2, ncol=2)


## ----mte_trim1, fig.height=3, fig.width=7, fig.align='center'-----------------
semiiv1 = semiivreg(y~d|w0|w1, data, 
                    ref_indiv = data.frame(w0=0, w1=0), 
                    common_supp_trim = c(0.1, 0.90),
                    plotting=FALSE)
mte_plot = semiiv1$plot$mte; mtr_plot = semiiv1$plot$mtr;
grid.arrange(mtr_plot, mte_plot, ncol=2)

