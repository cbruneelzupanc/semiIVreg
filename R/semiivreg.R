#' @title Semi-IV Regression Function

#' @description Semi-IV regression function from \href{https://www.cbruneel.com/research}{Bruneel-Zupanc (2024)}.
#' Syntax inspired from ivreg. Returns MTE and MTR curves with confidence intervals.
#' The estimation is almost *instantaneous* (a few seconds at most).
#'
#' By default, return analytic standard errors not accounting for the fact that the propensity score is estimated in a first stage in `semiivreg`.
#' Use `semiivreg_boot` to obtain 'correct' bootstrapped confidence intervals (takes a bit longer).
#'
#' @param formula Formula of the regression, of the form outcome ~ treatment | semi-iv0 | semi-iv1 | commoncovariates. \cr
#' The treatment variable should be binary (0, 1). \cr
#' covariates with an effect that differs on D=1 and D=0 should be included in each semi-iv0 and semi-iv1. \cr
#' with est_method = "locpoly": cannot restrict covariates to have common effects (not implemented), so `commoncovariates` will just be estimated as having generally a different effect on Y0 and Y1.
#'
#' @param data Dataframe containing the data.
#' @param propensity_formula Formula for the 1st stage. If nothing specified, just runs a probit of d ~ semi-iv0 + semi-iv1 + covariates (removing the redundant variables).
#' @param ref_indiv Specify the reference individual (in terms of covariates) at which we will evaluate the function. \cr
#' by default takes the average value for the numerical covariates, and the reference level for factors. \cr
#  Can also specify a data.frame with several reference individuals. The MTR and MTE will be computed for all of them. But the plot only for the first one.
#' @param firststage_model By default, the first stage is a probit model. Can specify another model (e.g., "logit").
#' @param est_method Estimation method: default is "locpoly" for Robinson (1988) double residual regression for partially linear model. Other options include "sieve" to specify flexibly the control function as a polynomial with pol_degree_sieve, and "homogenous" which is a sieve where we also impose homogenous treatment effect.
#' @param bw0,bw1 Bandwidth of the first residual regressions of Wd and X on Phat. Need to be specified in the order of the covariates as specified in the model. Be very careful with factors. Default NULL and computed using the specified bw_method.
#' @param bw_y0,bw_y1 Bandwidth of the second regression of Y (net of the effects of the covariates) on Phat. Default NULL and computed using the specified bw_method.
#' @param bw_method Method to compute the bandwidth of the local polynomial regressions. Default is simple "plug-in" method.
#' @param pol_degree_locpoly1 Degree of the local polynomial regression of the covariates on Phat. Default is 1 as recommended by Fan and Gijbels (1996) because we want to estimate the regular function.
#' @param pol_degree_locpoly2 Degree of the local polynomial regression of Y (net of the effects of the covariates) on Phat. Default is 2 as recommended by Fan and Gijbels (1996) because we want to estimate the derivative function.
#' @param pol_degree_sieve Degree of the polynomial transformation for the control function.
#' @param conf_level Confidence level for the confidence intervals.
#' @param common_supp_trim Vector of two values indicating the set of propensity scores at which we will evaluate the function. \cr
#' Default is the full support \code{[0,1]}. But can be trimmed manually.
#' @param trimming_value Can either be a vector c(0.05, 0.95) indicating the quantile of the propensity score above which and below which we keep the observations for both D=0 and D=1. \cr
#' Can also be a single value, in which case symmetric trimming up and down. \cr
#' Inserting a trimming_value generates automatic_trim = TRUE automatically.
#' @param automatic_trim If TRUE, the estimation of the second stage is done on the common_support only.
#' @param plotting TRUE if wants to plot at the end of the function, FALSE otherwise.
#'
#'
#'
#'
#'


#' @return A list with the following elements:
#' \describe{
#'  \item{`$data`}{Returns data of output estimation used to plot the MTE and MTR. In details:
#'  \describe{
#'    \item{`$RES`}{Dataframe with the estimated MTE and MTR values (and their confidence intervals if `est_method="sieve"` or `"homogenous"`) for a sequence of unobservable resistance to treatment in the identifiable common support.}
#'    \item{`$data`}{Original data used for the estimation where we added the propensity score estimated, named `Phat`, and where we made the transformation of the eventual `factor` variables as dummies.}
#'    \item{`$ref_indiv`}{Reference individual(s) at which we evaluate the MTE and MTR.}
#'    \item{`$Xdat`}{Set of covariates (this output is used for the bootstrap).}
#'  }}
#'
#' \item{`$estimate`}{Returns the estimation of:
#' \describe{
#'    \item{`$est, or $est0 and $est1`}{If `est_method = "locpoly"`, `est0` and `est1` returns the second stage estimates of the effect of the covariates and semi-IVs on their respective potential outcomes.
#'    Coming out of the double residual regression à la Robinson (1988), running a no-intercept OLS of the residuals Y-E(Yd|P) on the residuals of every semi-IVs, Wd-E(Wd|P), and covariates, X-E(X|P). \cr
#'    If `est_method = "sieve"` or `"homogenous"`, `$est` returns the second stage estimates of E(Yd | D, Wd, X, P) where the propensity score is controlled for with a flexible control function (polynomial), Kd(P).}
#'    \item{`$propensity`}{First stage estimate of the propensity score.}
#'    \item{`$avg_MTE`}{Average of the MTE over the identified common support. If full common support, it is an estimate of the ATE(x, w0, w1). If `est_method="homogenous"`, the MTE is constant so it also gives the ATE(x, w0, w1).}
#' }}
#'
#' \item{`$bw`}{Returns the bandwidth used (or estimated via `bw_method`) in the Robinson double residual regression. \cr
#' `bw0` and `bw1` are the bandwidths of the first residual regressions of Yd, Wd and X on Phat. \cr
#' `bw_y0` and `bw_y1` are the bandwidths of the second regression of Y (net of the effects of the covariates) on Phat. These are the one that matters for the smoothness of the MTE and MTR estimates.
#' }
#'
#' \item{`$plot`}{Returns separately the following plot objects: `supp` (support), `mtr`, `mte`.}
#' \item{`$supp`}{Returns the common support of the propensity score `Phat` between the two treatment group.}
#' \item{`$call`}{Returns the call of the function and the covariates and semi-IVs used.}
#' }

#' @details
#' # The estimated model
#' `semiivreg` estimates the marginal treatment effect (MTE) and marginal treatment response (MTR) of a binary treatment variable using semi-IVs, W0 and W1.
#' As with standard IVs (see Andresen, 2018), we estimate a semi-parametric partially linear model, as described in Bruneel-Zupanc (2024).
#' For more details on the model and estimation procedure, see the vignette \code{vignette("semiIVreg", package = "semiIVreg")}, also available online \href{https://cbruneelzupanc.github.io/semiIVreg/articles/semiIVreg.html}{here}.
#' For more details on the use of the `semiivreg` function, see also the vignettes \code{vignette("semiIVreg_heterogenousTE", package = "semiIVreg")} and \code{vignette("semiIVreg_homogenousTE", package = "semiIVreg")}.
#' For more details about causal inference with semi-IVs in general, see Bruneel-Zupanc (2024).
#'
#' # Caution about the Estimated Standard errors
#' By default, `est_method="locpoly"` returns no standard errors. \cr
#' If `est_method="sieve"` or `est_method="homogenous"`, it returns \bold{analytic standard errors}: but these are wrong because they do not account for the fact that the propensity score is estimated. \cr
#' In any case, we recommend to use `semiivreg_boot` to obtain 'correct' bootstrapped confidence intervals. Implemented separately because the bootstrap takes more time, while the baseline `semiivreg` function is almost instantaneous.

#' @references
#' Bruneel-Zupanc, C. (2024). Don't (fully) exclude me, it's not necessary! Identification with semi-IVs. arXiv preprint arXiv:2303.12667. \cr
#'
#' For empirical applications of the estimation of Marginal Treatment Effects with standard IVs, see for example: \cr
#' Carneiro, P., Heckman, J. J., & Vytlacil, E. J. (2011). Estimating marginal returns to education. American Economic Review, 101(6), 2754-2781. \cr
#'
#' Brinch, C. N., Mogstad, M., & Wiswall, M. (2017). Beyond LATE with a discrete instrument. Journal of Political Economy, 125(4), 985-1039. \cr
#'
#' In particular, see Andresen, M. E. (2018). Exploring marginal treatment effects: Flexible estimation using Stata. The Stata Journal, 18(1), 118-158.
#'
#' For double residual estimation of partially Linear models, see
#' Robinson, P. M. (1988). Root-N-consistent semiparametric regression. Econometrica: Journal of the Econometric Society, 931-954. \cr
#'
#' For local polynomial regressions choice of degree:
#' Fan, J., & Gijbels, I. (1996). Local polynomial modelling and its applications.
#'

#'@author
#' Christophe Bruneel-Zupanc, \url{cbruneel.com}

#' @usage semiivreg(formula, data, propensity_formula=NULL,
#'                  ref_indiv =NULL, firststage_model = "probit",
#'                  est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
#'                  bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = "plug-in",
#'                  pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
#'                  pol_degree_sieve = 5, conf_level = 0.05,
#'                  common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
#'                  plotting=TRUE)


#' @examples
#' # Load data:
#' data(roydata)
#'
#' # Run the semi-IV regression
#' semiiv = semiivreg(y~d|w0|w1, data=roydata)
#' semiiv = semiivreg(y~d|w0|w1|Xbinary + Xcontinuous, data=roydata) # with covariates
#' semiiv = semiivreg(y~d|w0+Xbinary|w1+Xbinary|Xcontinuous, data=roydata)
#' # Xbinary has different effect on Y0 and Y1, Xcontinuous has the same.
#' semiiv = semiivreg(y~d|w0|w1, data=roydata, propensity_formula = d~w0+w1+w0:w1)
#' # if want to specify another first stage
#'
#' semiiv$plot$mtr # if want to plot mtr_plot



#' @rdname semiivreg
#' @export
semiivreg = function(formula, data, propensity_formula=NULL,
                     ref_indiv =NULL, firststage_model = "probit",
                     est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
                     bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = "plug-in",
                     pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
                     pol_degree_sieve = 5, conf_level = 0.05,
                     common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
                     plotting=TRUE) {


  # 0. Construction of the variables/objects to be used
  # ---------------------------------------------------
  # (i) Convert to data.frame + na.omit
  data = as.data.frame(data) # maybe problem with data loaded from Stata
  vars = all.vars(formula)
  data = subset(data, select=vars); data = na.omit(data)
  data = transform_factor(formula, data) # the variables which are specified as "factor(x)" are transformed into factors # done in order to ensure ordering with ref indiv
  data_orig = data;

  # (ii) Create reference individual if missing
  if(is.null(ref_indiv)) { ref_indiv_orig = create_ref_indiv(formula, data) } else { ref_indiv_orig = ref_indiv }

  # (iii) Construct transformed data
  # Transform factors into dummies for example
  res = construct_data(formula, data_orig)
  Xdat = res$data;
  var_outcome = res$var_outcome; var_treatment = res$var_treatment; var_w0 = res$var_w0; var_w1 = res$var_w1; var_covariates = res$var_covariates
  form_w0 = res$form_w0; form_w1 = res$form_w1; form_covariates = res$form_covariates
  formula_X_orig = res$formula_X_orig
  datayd = subset(data_orig, select=c(var_outcome, var_treatment))

  # Updated formula with dummies instead of factors
  new_formula = res$formula

  # Final data to be used:
  data = cbind(datayd, Xdat)


  # (iv) Transformed reference individual
  name_all_X_orig = all.vars(formula_X_orig)
  Xdat_orig = subset(data_orig, select=name_all_X_orig)
  ref_indiv_orig = subset(ref_indiv_orig, select=name_all_X_orig) # remove potential "id"
  # Transform into factors as in data:
  for(j in 1:ncol(Xdat_orig)) {
    if(is.factor(Xdat_orig[,j])) {
      ref_indiv_orig[,j] = factor(as.character(ref_indiv_orig[,j]), levels=levels(Xdat_orig[,j]))
    }
  }
  ref = rbind(ref_indiv_orig, Xdat_orig) # inflate the data just to have all the levels of the factors;
  res_ref = construct_data(formula, data=ref)
  ref_indiv = res_ref$data[1:nrow(ref_indiv_orig),]
  ref_indiv$id = 1:nrow(ref_indiv)

  Xdat$id = NA #


  # (v) Ensures that the treatment variable is binary 0, 1
  if(is.numeric(data[[var_treatment]])) {
    Nunique_values = length(unique(data[[var_treatment]]))
  } else {
    data[[var_treatment]] = as.factor(data[[var_treatment]])
    Nunique_values = length(levels(data[[var_treatment]]))
    # Transform to binary:
    data[[var_treatment]] = as.numeric(data[[var_treatment]]) - 1
  }
  if(Nunique_values > 2) { error = "Error: the treatment variable is not binary"; stop(error) }
  if(Nunique_values < 2) { error = "Error: the treatment variable takes less than 2 values"; stop(error)}





  # 1. First stage: Propensity score estimation
  # --------------
  if(is.null(propensity_formula)) { # by default, simple only include w0 and w1, additively
    # Formula:
    var_propensity = unique(c(var_w0, var_w1, var_covariates)) # unique to keep only one if some covariates in several potential outcomes with different effects
    var_propensity = var_propensity[which(var_propensity != "")]
    propensity_formula = as.formula(paste0(var_treatment, "~", paste0(var_propensity, collapse="+")))
  }

  # (i) Estimation:
  if(firststage_model == "probit") {
    propensity = glm(propensity_formula, family=binomial(link="probit"), data=data)
  }
  if(firststage_model == "logit") {
    propensity = glm(propensity_formula, family=binomial(link="logit"), data=data)
  }
  data$Phat = predict.glm(propensity, newdata=data, type="response")
  #Phat = propensity$fitted.values # No because possibly missing values will make everything wrong;
  Xdat$Phat = data$Phat


  # (ii) Common Support
  # Common Support + seq_u on which the MTE will be estimated;
  # Default common support using min and max:
  supp_P0 = c(min(data$Phat[which(data[[var_treatment]] == 0)]), max(data$Phat[which(data[[var_treatment]] == 0)]))
  supp_P1 = c(min(data$Phat[which(data[[var_treatment]] == 1)]), max(data$Phat[which(data[[var_treatment]] == 1)]))
  common_supp = c(max(supp_P0[1], supp_P1[1]), min(supp_P0[2], supp_P1[2]))

  # Could also do common support using quantiles of Phat, trim below a limit:
  if(!is.null(trimming_value)) {
    automatic_trim = TRUE # to apply the trimming after
    if(length(trimming_value) == 1) { trimming_value = c(trimming_value, 1-trimming_value) } # if one value, assumes same trimming on both sides;
    supp_P0 = c(quantile(data$Phat[which(data[[var_treatment]] == 0)], trimming_value[1]), quantile(data$Phat[which(data[[var_treatment]] == 0)], trimming_value[2]))
    supp_P1 = c(quantile(data$Phat[which(data[[var_treatment]] == 1)], trimming_value[1]), quantile(data$Phat[which(data[[var_treatment]] == 1)], trimming_value[2]))
    common_supp = c(max(supp_P0[1], supp_P1[1]), min(supp_P0[2], supp_P1[2]))
  }

  if(common_supp[1] > common_supp[2]) { print("No common support!"); return(common_supp) } # STOP if empty common support

  # (ii-a) Plot "true" common support:
  # Round down the lowest value and up the highest value (to make sense with the graphical bars of bin 0.01)
  common_supp_round = c( trunc(common_supp[1]*100)/100, trunc((common_supp[2]+0.01)*100)/100)
  data$Treatment = factor(data[[var_treatment]])
  supp_plot = supp_plot_fun(data, common_supp_round)

  data_not_trimmed = data;


  # (ii-b) Trim according to the common support
  # Trim the observations outside of the pre-specified common_support, for a better estimate
  data = data[which(data$Phat >= common_supp_trim[1] & data$Phat <= common_supp_trim[2]),]
  Xdat = Xdat[which(Xdat$Phat >= common_supp_trim[1] & Xdat$Phat <= common_supp_trim[2]),]


  # Sequence of u at which to evaluate the MTE:
  seq_u = seq(0, 1, by=0.001) # May trim a priori depending on the trimming arguments;
  seq_u = seq_u[which(seq_u >= common_supp_trim[1] & seq_u <= common_supp_trim[2])]
  if(automatic_trim == TRUE) { # can trim further based on estimated P
    data = data[which(data$Phat >= common_supp[1] & data$Phat <= common_supp[2]),] # only estimate on the common support
    Xdat = Xdat[which(Xdat$Phat >= common_supp[1] & Xdat$Phat <= common_supp[2]),]
    seq_u = seq_u[which(seq_u >= common_supp[1] & seq_u <= common_supp[2])]
  }
  # Update common support for later MTE and MTR plots
  common_supp[1] = max(common_supp[1], common_supp_trim[1])
  common_supp[2] = min(common_supp[2], common_supp_trim[2])
  common_supp_round_mte = c( trunc((common_supp[1]+0.001)*1000)/1000, trunc((common_supp[2])*1000)/1000)
  # for the first value take the first above (so that mimic the seq_u choice)





  # 2. Second Stage:
  # ----------------
  # Several estimation methods here. By default do Robinson (1988) double residual regression for partially linear model.

  # Method 1: Local Polynomial Regressions
  # ---------
  if(est_method == "locpoly") {
    #if(var_covariates != "") { print("-- Common covariates treated as having generally different effect on D=0 and D=1") }
    # with Robinson, we don't allow the covariates to have the same effect on both D=0 and D=1

    res0 = mtr_fun(d=0, data, ref_indiv, seq_u,
                   bwd = bw0, bw_y = bw_y0, bw_method = bw_method,
                   pol_degree1 = pol_degree_locpoly1, pol_degree2 = pol_degree_locpoly2,
                   var_outcome=var_outcome, var_treatment=var_treatment, var_w0=var_w0, var_w1=var_w1, var_covariates=var_covariates)
    res1 = mtr_fun(d=1, data, ref_indiv, seq_u,
                   bwd = bw1, bw_y = bw_y1, bw_method = bw_method,
                   pol_degree1 = pol_degree_locpoly1, pol_degree2 = pol_degree_locpoly2,
                   var_outcome=var_outcome, var_treatment=var_treatment, var_w0=var_w0, var_w1=var_w1, var_covariates=var_covariates)

    mtr0 = res0$RES; mtr1 = res1$RES
    mtr1 = subset(mtr1, select=c("id", "Phat", "mtr1"))
    RES = base::merge(mtr0, mtr1, by=c("id", "Phat")) # normal to have some NA

    RES$mte = RES$mtr1 - RES$mtr0;

    # Graphical Parameter:
    conf_band = FALSE # No confidence band with this method

    # Other outputs:
    est0 = res0$estd; est1 = res1$estd;
    bw0 = res0$bwd; bw_y0 = res0$bw_y;
    bw1 = res1$bwd; bw_y1 = res1$bw_y;

    var_cov_2nd = NULL
  }







  # Method 2: "Sieve" method -> Flexible polynomial specification for Kappa_d
  # ---------

  if(est_method %in% c("sieve", "homogenous")) {

    # Graphical parameter:
    conf_band = TRUE

    # Shape of the control functions - for now only polynomials of degree pol_degree
    pol_degree = pol_degree_sieve

    # Formula
    var_d0_base = paste0("I(1 - ", var_treatment, ")")
    var_d1_base = paste0("I(", var_treatment, ")")
    var_d0 = paste0("I(1 - ", var_treatment, "):", var_w0) #var_d0 = paste0("I( (1 - ", var_treatment, ")*", var_w0, ")")
    var_d1 = paste0("I(", var_treatment, "):", var_w1) #var_d1 = paste0("I(", var_treatment, "*", var_w1, ")")
    var_cov_2nd = c(var_d1_base, var_d0, var_d1, var_covariates) # var_d0_base is the intercept.
    var_cov_2nd = var_cov_2nd[which(var_cov_2nd != "")] # if no covariates
  }

  # Method 2.1. General Sieve, with heterogenous Treatment Effects
  # -----------
  if(est_method == "sieve") {

    # Sieve-1) Estimate E[Yd | D, Z, X]
    # ---------------------------------
    # -> gives estimate of Kappa_d(P) -> from kappa_d, obtains the MTR and MTE.
    # Estimate as a stacked regression (useful if want to impose same effect of covariates + practical afterwards)
    formula_control_function = paste0("I(1-", var_treatment, "):Kappa_fun(Phat, pol_degree) + I(", var_treatment, "):Kappa_fun(Phat, pol_degree)")
    formula_2nd_stage = as.formula(paste0(var_outcome, "~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function))

    # Estimation
    est = lm(formula_2nd_stage, data)


    # Sieve-2) Marginal Treatment Responses and MTE estimation
    # ---------------------------------------------------------
    # From the estimated coefficients, can construct any MTR and MTE now.

    # (i) Extract Coefficients from the main estimation;
    coeff = coefficients(est)
    vcov = vcov(est)
    t_value = qt(1-conf_level/2, df = df.residual(est))

    # (ii) Construct kd(u) functions -> see formula in Andresen for e.g.
    # Done outside of here; kdu need to correspond to Kappa transformation
    # Then: kdu is equal to: kdu_transform_fun(seq_u, d=rep(d, length(u)))%*%coeffkd

    # (iii) Marginal Treatment Responses: MTR(u, x) and MTE(u, x)
    RES = mtr_fun_sieve(coeff=coeff, vcov=vcov, var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, ref_indiv=ref_indiv,
                        seq_u=seq_u, homogenous=FALSE, pol_degree=pol_degree,
                        conf_level=conf_level, t_value=t_value, Xdat = Xdat)

  }


  # Method 2.2: "homogenous" method -> Homogenous treatment effects
  # -----------
  if(est_method == "homogenous") {

    # Homogenous-1) Estimation of the regression
    # Computation non trivial (see paper); note that it still allows for a flexible Kappa; only imposes that MTR1 - MTR0 is constant;
    # Formula
    formula_control_function_homogenous = paste0("I(-(1-", var_treatment, ")*Phat/(1-Phat) + ", var_treatment, "):Kappa_homogenous_fun(Phat, pol_degree)")
    formula_2nd_stage_homogenous = as.formula(paste0(var_outcome, "~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function_homogenous))

    # Estimation
    est = lm(formula_2nd_stage_homogenous, data)

    # Homogenous-2) MTR and MTE with homogenous TE
    # Remark: implies flat (constant) MTE, but MTR can still vary.

    # (i) Extract Coefficients from the main estimation;
    coeff = coefficients(est)
    vcov = vcov(est)
    t_value = qt(1-conf_level/2, df = df.residual(est))

    # (ii) MTR and MTE (MTE is constant here by construction)
    RES = mtr_fun_sieve(coeff=coeff, vcov=vcov, var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, ref_indiv=ref_indiv,
                        seq_u = seq_u, homogenous=TRUE, pol_degree=pol_degree, conf_level=conf_level, t_value=t_value, Xdat=Xdat)
  }






  # 3. Plotting reports
  # -------------------
  # MTR and MTE plots
  # For the MTE/MTR plots, can only plot for one reference individual
  if(nrow(ref_indiv) > 1) { RES_plot = RES[which(RES$id == 1),] } else { RES_plot = RES; }
  # for the main object, can still return the predictions for any individual;
  # Then can use again the function plot_mte for any of the reported individual data;
  dat_plot = RES_plot

  # MTR curves
  mtr_plot = mtr_plot_fun(dat_plot, common_supp_round_mte, conf_band)

  # MTE curves
  mte_plot = mte_plot_fun(dat_plot, common_supp_round_mte, conf_band)

  # By default: plot the supp_plot and the mte_plot:
  if(plotting == TRUE) { grid.arrange(supp_plot, mte_plot, ncol=2) }


  # Also estimate the average MTE for the ref individual (= homogenous TE if homogenous)
  # Not necessarily = ATE if not full support;
  avg_MTE = mean(RES_plot$mte, na.rm=TRUE)




  # 4. Return objects
  # -----------------

  output_data = list(RES, data_not_trimmed, ref_indiv, Xdat); names(output_data) = c("RES", "data", "ref_indiv", "Xdat")
  output_plot = list(supp_plot, mtr_plot, mte_plot); names(output_plot) = c("supp", "mtr", "mte")
  output_supp = common_supp
  output_call = list(new_formula, formula, var_treatment, var_outcome, var_w0, var_w1, var_covariates, var_cov_2nd); names(output_call) = c("formula", "formula_orig", "var_treatment", "var_outcome", "var_w0", "var_w1", "var_covariates", "var_cov_2nd")

  if(est_method %in% c("sieve", "homogenous")) {
    output_estimate = list(est, propensity, avg_MTE);
    names(output_estimate) = c("est", "propensity", "avg_MTE")

    output_bw = NA # not relevant for this method
  }

  if(est_method == "locpoly") {
    output_estimate = list(est0, est1, propensity, avg_MTE);
    names(output_estimate) = c("est0", "est1", "propensity", "avg_MTE")

    output_bw = list(bw0, bw1, bw_y0, bw_y1)
    names(output_bw) = c("bw0", "bw1", "bw_y0", "bw_y1")
  }

  output = list(output_data, output_estimate, output_bw, output_plot, output_supp, output_call);
  names(output) = c("data", "estimate", "bw", "plot", "supp", "call")

  return(output)
}
















# ------------------------------------------------------
# Part 2. Bootstrap standard errors/confidence intervals
# ------------------------------------------------------
# The main function does not take into account that the propensity score Phat is estimated in a first stage to compute the standard errors
# We try to account for this with corrected standard errors using bootstrap in this function.



#' @rdname semiivreg
#' @usage semiivreg_boot(formula, Nboot=500, data, propensity_formula=NULL, ref_indiv =NULL,
#'                firststage_model="probit",
#'                pol_degree_transform = 5, common_supp_trim=c(0,1), trimming_value = NULL,
#'                automatic_trim = FALSE, plotting=TRUE, conf_level = 0.05, CI_method = "delta")
#' @param Nboot Number of bootstrap samples.
#' @param CI_method "delta" for delta method, "curve" for bootstrap the MTE curves directly. With est_method = "locpoly", only "curve" method is possible.
#'
#' @export
semiivreg_boot = function(formula, Nboot=500, data, propensity_formula=NULL,
                          ref_indiv =NULL, firststage_model = "probit",
                          est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
                          bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = "plug-in",
                          pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
                          pol_degree_sieve = 5, conf_level = 0.05,
                          common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
                          plotting=TRUE, CI_method = "curve") {



  # 1. semiivreg on the full sample
  # -------------------------------
  # to extract main estimate and eventually the support:

  main_res = semiivreg(formula=formula, data=data, propensity_formula=propensity_formula,
                       ref_indiv = ref_indiv, firststage_model = firststage_model,
                       est_method = est_method, bw0 = bw0, bw1 = bw1, bw_y0 = bw_y0, bw_y1 = bw_y1, bw_method = bw_method,
                       pol_degree_locpoly1 = pol_degree_locpoly1, pol_degree_locpoly2 = pol_degree_locpoly2,
                       pol_degree_sieve = pol_degree_sieve, conf_level = conf_level,
                       common_supp_trim = common_supp_trim, trimming_value = trimming_value, automatic_trim = automatic_trim, plotting=FALSE)

  # Extract ref_indiv (if not specified)
  ref_indiv = main_res$data$ref_indiv

  # Updated formula
  transform_formula = main_res$call$formula

  # Updated data:
  orig_data = data;
  data = main_res$data$data

  # Extract the bandwidth (if est_method = "locpoly"):
  # same bandwidth for every bootstrap
  if(est_method == "locpoly") {
    bw0 = main_res$bw$bw0; bw1 = main_res$bw$bw1; bw_y0 = main_res$bw$bw_y0; bw_y1 = main_res$bw$bw_y1
  }

  # Extract Support:
  common_supp_trim_boot = common_supp_trim;
  if(automatic_trim == TRUE) { common_supp_trim_boot = main_res$supp; } # same supp for every simulation
  if(!is.null(trimming_value)) { common_supp_trim_boot = main_res$supp; } # same supp for every simulation
  # Remark: in every simulation, applies the same support for the estimates. So always put automatic_trim = FALSE and trimming_value = NULL
  #         The support will be determined by the common_supp_trim of the main regression.




  # 2. Bootstrap
  # ------------
  BOOT = list()
  for(k in 1:Nboot) {

    #set.seed(1234*k) # just set seed outside

    # Bootstrap sample
    bootstrap_indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    bootstrap_sample <- data[bootstrap_indices, ]

    # semi-iv estimation
    boot_res = invisible(semiivreg(formula=transform_formula, data=bootstrap_sample,
                                   propensity_formula=propensity_formula,
                                   ref_indiv = ref_indiv, firststage_model = firststage_model,
                                   est_method = est_method, bw0 = bw0, bw1 = bw1, bw_y0 = bw_y0, bw_y1 = bw_y1, bw_method = bw_method,
                                   pol_degree_locpoly1 = pol_degree_locpoly1, pol_degree_locpoly2 = pol_degree_locpoly2,
                                   pol_degree_sieve = pol_degree_sieve, conf_level = conf_level,
                                   common_supp_trim = common_supp_trim_boot, trimming_value=NULL, automatic_trim=FALSE, plotting=FALSE))
    BOOT[[k]] = boot_res
  }



  # 3. Standard Errors around the coefficients
  # ------------------------------------------
  if(est_method %in% c("sieve", "homogenous")) {

    BOOTest = list()
    for(k in 1:Nboot) { BOOTest[[k]] = coefficients(BOOT[[k]]$estimate$est) }
    EST = do.call('rbind', BOOTest)
    meanEST = apply(EST, 2, mean);
    main_res_coeff = coefficients(main_res$estimate$est)
    vcov = cov(EST)
    standarderror = sqrt(diag(vcov))
    COEFF = data.frame(estimate=main_res_coeff, std_error=standarderror, boot_mean_estimate=meanEST)

  }


  if(est_method == "locpoly") {

    # Main coefficients:
    main_coeff0 = coefficients(main_res$estimate$est0); names(main_coeff0) = paste0("Untreated_", names(main_coeff0))
    main_coeff1 = coefficients(main_res$estimate$est1); names(main_coeff1) = paste0("Treated_", names(main_coeff1))
    main_res_coeff = c(main_coeff0, main_coeff1)

    BOOTest = list()
    for(k in 1:Nboot) {
      # est0:
      coeff0 = coefficients(BOOT[[k]]$estimate$est0); names(coeff0) = paste0("Untreated_", names(coeff0))
      coeff1 = coefficients(BOOT[[k]]$estimate$est1); names(coeff1) = paste0("Treated_", names(coeff1))
      coeff = c(coeff0, coeff1)
      BOOTest[[k]] = coeff
    }

    EST = do.call('rbind', BOOTest)
    meanEST = apply(EST, 2, mean);
    vcov = cov(EST)
    standarderror = sqrt(diag(vcov))
    COEFF = data.frame(estimate=main_res_coeff, std_error=standarderror, boot_mean_estimate=meanEST)
  }







  # 4. Confidence intervals for heterogenous TE
  # -------------------------------------------
  # Two possibilities:
  # (i) Delta method: bootstrap the variance-covariance matrix of the estimates; then apply the usual estimation with this corrected (bootstrapped) vcov matrix.
  # (ii) Bootstrap the MTE curves directly.

  # Option 1: Delta Method
  # ---------
  if(CI_method == "delta") {
    if(! est_method %in% c("sieve", "homogenous")) { stop("Delta method only possible with est_method = 'sieve' or 'homogenous'")}

    t_value = qt(1-conf_level/2, df = df.residual(main_res$estimate$est))
    var_treatment = main_res$call$var_treatment; var_cov_2nd = main_res$call$var_cov_2nd
    seq_u = seq(0, 1, by=0.001); pol_degree = pol_degree_sieve
    Xdat = main_res$data$Xdat

    # Then simply return the main function results using this new vcov around the main estimates:
    if(est_method == "homogenous") { homogenous = TRUE } else { homogenous = FALSE }
    RES = mtr_fun_sieve(coeff=main_res_coeff,
                        vcov=vcov, var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, ref_indiv=ref_indiv,
                        seq_u=seq_u, homogenous=homogenous, pol_degree=pol_degree,
                        conf_level=conf_level, t_value=t_value, Xdat=Xdat)

  }

  # Remark: with est_method="locpoly"
  # -> could still report standard errors around the estimation of the coefficients in the propensity score stage
  #    and in the double residual regression stage.



  # Option 2: CI around the MTE curves directly
  # ---------
  if(CI_method == "curve") {

    BOOTDAT = list()
    for(k in 1:Nboot) { bootdat = BOOT[[k]]$data$RES; bootdat$boot_id = k; BOOTDAT[[k]] = bootdat }
    bootdat = do.call('rbind', BOOTDAT)

    bootdt = data.table(bootdat)
    if(est_method %in% c("sieve", "homogenous")) {
      bootdt = data.frame(bootdt)
      bootdt = bootdt[, which(!colnames(bootdt) %in% c("mtr0_lwr", "mtr0_upr", "mtr0_se", "mtr1_lwr", "mtr1_upr", "mtr1_se", "mte_lwr", "mte_upr", "mte_se"))]
      bootdt = data.table(bootdt)
    }

    CI = bootdt[, .(
      mte_avg = mean(mte, na.rm=TRUE),
      mte_lwr = quantile(mte, conf_level/2, na.rm=TRUE),
      mte_upr = quantile(mte, 1-conf_level/2, na.rm=TRUE),
      mtr0_avg = mean(mtr0, na.rm=TRUE),
      mtr0_lwr = quantile(mtr0, conf_level/2, na.rm=TRUE),
      mtr0_upr = quantile(mtr0, 1-conf_level/2, na.rm=TRUE),
      mtr1_avg = mean(mtr1, na.rm=TRUE),
      mtr1_lwr = quantile(mtr1, conf_level/2, na.rm=TRUE),
      mtr1_upr = quantile(mtr1, 1-conf_level/2, na.rm=TRUE)
    ), by = .(Phat, id)]
    # Missing values are normal -> because at the tails MTE and MTR will be missing for many obs


    # merge with main estimates for the middle point:
    mdat = main_res$data$RES;
    if(est_method %in% c("sieve", "homogenous")) {
      mdat = mdat[, which(!colnames(mdat) %in% c("mtr0_lwr", "mtr0_upr", "mtr0_se", "mtr1_lwr", "mtr1_upr", "mtr1_se", "mte_lwr", "mte_upr", "mte_se"))]
    }

    # # Keep analytic confidence intervals (to compare how far from bootstrap one they are) -> not possible with locpoly
    # mdat$mte_lwr_analytic = mdat$mte_lwr; mdat$mte_upr_analytic = mdat$mte_upr;
    # mdat$mtr0_lwr_analytic = mdat$mtr0_lwr; mdat$mtr0_upr_analytic = mdat$mtr0_upr;
    # mdat$mtr1_lwr_analytic = mdat$mtr1_lwr; mdat$mtr1_upr_analytic = mdat$mtr1_upr;
    # mdat = subset(mdat, select=c("id", "Phat", "mte", "mtr0", "mtr1",
    #                              "mte_lwr_analytic", "mte_upr_analytic",
    #                              "mtr0_lwr_analytic", "mtr0_upr_analytic", "mtr1_lwr_analytic", "mtr1_upr_analytic"))
    # mdat = merge(ref_indiv, mdat, by="id")

    RES = merge(mdat, CI, by=c("id", "Phat"))

  }




  # 5. Plotting
  # -----------
  common_supp_plot = common_supp_trim_boot #main_res$supp # purely a graphical option.

  # 5.1. Heterogenous TE
  if(nrow(ref_indiv) > 1) { RES_plot = RES[which(RES$id == 1),] } else { RES_plot = RES; } # if several ref indiv, only plot for the main one;

  # MTE and MTR plots
  mte_plot = mte_plot_fun(RES_plot, common_supp_plot, conf_band=TRUE)
  mtr_plot = mtr_plot_fun(RES_plot, common_supp_plot, conf_band=TRUE)

  # Plot results
  main_supp_plot = main_res$plot$supp
  if(plotting == TRUE) { grid.arrange(main_supp_plot, mte_plot, ncol=2) }



  # 6. Output
  # ---------
  output_main = main_res;
  output_plot = list(main_supp_plot, mtr_plot, mte_plot); names(output_plot) = c("supp", "mtr", "mte")
  output_estimate = list(RES, COEFF, vcov); names(output_estimate) = c("est", "coeff", "vcov")

  output = list(output_main, output_plot, output_estimate); names(output) = c("main", "plot", "estimate")
  return(output)
}


NULL
