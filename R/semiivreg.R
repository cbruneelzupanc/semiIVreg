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
#' @param propensity_data Data used to compute the 1st stage; ignore by default set to NULL and = data. Mainly useful for internal bootstrap function is the first stage formula is different from the default one.
#' @param ref_indiv Specify the reference individual (in terms of covariates) at which we will evaluate the function. \cr
#' By default takes the average value for all the covariates (on the trimmed dataset) to compute the average estimate. Remark: for factors, the average is computed on the dummy variables to get the proper average effect. \cr
#  Can also specify a data.frame with several reference individuals. The MTR and MTE will be computed for all of them. But the plot only for the first one.
#' @param firststage_model By default, the first stage is a probit model. Can specify another model (e.g., "logit").
#' @param est_method Estimation method: default is "locpoly" for Robinson (1988) double residual regression for partially linear model. Other options include "sieve" to specify flexibly the control function as a polynomial with pol_degree_sieve, and "homogenous" which is a sieve where we also impose homogenous treatment effect.
#' @param se_type Type of standard errors if sieve/homogenous estimation. By default = "HC1". Can otherwise me any of the possibilities of vcovHC in the sandwich package. Also possible to have se_type="nonrobust" for the non-robust (lm default).
#' @param bw0,bw1 Bandwidth of the first residual regressions of (Y, Wd and X) on Phat. \cr
#' Two possibilities: specify one value that is applied to all covariates (and Y), or specify a different bandwidth for the regression on each covariate. In the second case, need to be specified in the order of the covariates as specified in the model. Be very careful with factors. \cr
#' Default NULL and computed using the specified bw_method. Ideally, if one factor covariate, apply the same bandwidth to all of the dummies created from the factor.
#' @param bw_y0,bw_y1 Bandwidth of the second regression of Y (net of the effects of the covariates) on Phat. Default NULL and computed using the specified bw_method.
#' @param bw_method Method to compute the bandwidth of the local polynomial regressions (of the first-order derivative).
#' Default option is 1/5, which arbitrarily sets bw0, bw1, bw_y0 and bw_y1 to 1/5th of the support (rounded to the 3th digit). Can place any fraction < 1. \cr
#' Recommended alternatives include (global constant) bandwidth computations from nprobust package (Calonico, Cattaneo and Farrell, 2019)
#' (i) "mse-dpi": direct plug-in MSE optimal bandwidth from Fan and Gijbels (1996).
#' (ii) "mse-rot": rule-of-thumb implementation of the MSE-optimal bandwidth.
#' These two methods take long with large sample: use bw_subsamp_size to speed up the computation.
#' @param bw_subsamp_size Size of the subsample to use for the bandwidth selection. Default is 10,000. Use `bw_subsamp_size = NULL` to use the full sample (may take time). Otherwise, recommend to set a number around 20,000 at most for reasonable computation time (exponentially increasing with sample size). \cr
#' `bw_subsamp_size` introduces some randomness into the bandwidth selection procedure: recommended to set a seed before running semiivreg for reproducibility.
#' @param pol_degree_locpoly1 Degree of the local polynomial regression of the covariates on Phat. Default is 1 as recommended by Fan and Gijbels (1996) because we want to estimate the regular function.
#' @param pol_degree_locpoly2 Degree of the local polynomial regression of Y (net of the effects of the covariates) on Phat. Default is 2 as recommended by Fan and Gijbels (1996) because we want to estimate the derivative function.
#' @param kernel Kernel to use for the local polynomial regressions. Default is "gaussian" but can be "epanechnikov". Takes longer with Epanechnikov (cannot use fast locpoly implementation from KernSmooth).
#' @param `fast_robinson1` Default is TRUE to speed things up in a first stage (if many covariates in particular). If TRUE, will use the locpoly function from Kernsmooth library to speed up the computation of the Robinson double residual first stage. This is only possible if no external weights are used. Fast Locpoly will enforce a gaussian kernel.
#' @param `fast_robinson2` Default is FALSE. If TRUE, will use the locpoly function from Kernsmooth library to speed up the computation of the Robinson double residual second stage. This is only possible if no external weights are used. Fast Locpoly will enforce a gaussian kernel.
#' Default is FALSE for the second stage because fast_locpoly returns no standard errors and the gain in time is not so important for the second stage.
#' @param pol_degree_sieve Degree of the polynomial transformation for the control function.
#' @param conf_level Confidence level for the confidence intervals.
#' @param common_supp_trim Vector of two values indicating the set of propensity scores at which we will evaluate the function. \cr
#' Default is the full support \code{[0,1]}. But can be trimmed manually.
#' @param trimming_value Can either be a vector c(0.05, 0.95) indicating the quantile of the propensity score above which and below which we keep the observations for both D=0 and D=1. \cr
#' Can also be a single value, in which case symmetric trimming up and down. \cr
#' Inserting a trimming_value generates automatic_trim = TRUE automatically.
#' @param automatic_trim If TRUE, the estimation of the second stage is done on the common_support only.
#' @param weight_var A variable of weights to be applied to the observations. Default is NULL, apply equal weights to all observations. \cr
#' Implemented completely for `est_method = "sieve"` for now. For `locpoly`, the weights are not used when computing the "optimal bandwidth".
#' @param plotting TRUE if wants to plot at the end of the function, FALSE otherwise.
#' @param print_progress TRUE if wants to print the progress of the function, FALSE otherwise (default=FALSE).
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
#'    \item{`$deltaX`}{Returns the estimated effects of the covariates and semi-IVs (without intercept) for the specified reference individuals.}
#'  }}
#'
#' \item{`$estimate`}{Returns the estimation of:
#' \describe{
#'    \item{`$est, or $est0 and $est1`}{If `est_method = "locpoly"`, `est0` and `est1` returns the second stage estimates of the effect of the covariates and semi-IVs on their respective potential outcomes.
#'    Coming out of the double residual regression à la Robinson (1988), running a no-intercept OLS of the residuals Y-E(Yd|P) on the residuals of every semi-IVs, Wd-E(Wd|P), and covariates, X-E(X|P). \cr}
#'    \item{`$mtr0, $mtr1 and $mte`}{If `est_method = "sieve"` or `"homogenous"`, returns the functional form estimated for both MTR and MTE.}
#'    \item{`$kv`}{Returns the estimated k_d(v) (=E(Ud|V=v)). Includes the constant. If sums with the effect of covariates and semi-IVs (deltadX), gives the mtr_d.}
#'    \item{`$propensity`}{First stage estimate of the propensity score.}
#'    \item{`$est_kappa`}{If `est_method = "sieve"` or `"homogenous"`, this returns the estimated model for E(Y|D=d, X, Wd, P). From this, we extract Kappad(P) = E(Ud | D=d, P=p) from which we compute the kd(v) and mtrd(v, x, wd) functions. }
#'    \item{`$avg_MTE`}{Average of the MTE over the identified common support. If full common support, it is an estimate of the ATE(x, w0, w1). If `est_method="homogenous"`, the MTE is constant so it also gives the ATE(x, w0, w1).}
#' }}
#'
#' \item{`$bw`}{Returns the bandwidth used (or estimated via `bw_method`) in the Robinson double residual regression. \cr
#' `bw0` and `bw1` are the bandwidths of the first residual regressions of Yd, Wd and X on Phat. \cr
#' `bw_y0` and `bw_y1` are the bandwidths of the second regression of Y (net of the effects of the covariates) on Phat. These are the one that matters for the smoothness of the MTE and MTR estimates.
#' }
#'
#' \item{`$plot`}{Returns separately the following plot objects: `supp` (support), `mtr`, `mte` and `mte2`. `mte` reports the estimation from "local IV" approach, with standard errors from the Robinson 2nd stage. `mte2` reports the MTE estimated as the difference between the MTR (without standard errors).}
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
#' For local polynomial regressions choice of degree & Bandwidth computation
#' Fan, J., & Gijbels, I. (1996). Local polynomial modelling and its applications.
#' Calonico, S., Cattaneo, M. D., & Farrell, M. H. (2019). nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference. Journal of Statistical Software, 91(8), 1–33. https://doi.org/10.18637/jss.v091.i08
#'

#'@author
#' Christophe Bruneel-Zupanc, \url{cbruneel.com}

#' @usage semiivreg(formula, data, propensity_formula=NULL, propensity_data = NULL,
#'                  ref_indiv =NULL, firststage_model = "probit",
#'                  est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
#'                  se_type = "HC1",
#'                  bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = 1/5,
#'                  kernel="gaussian",
#'                  pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
#'                  pol_degree_sieve = 5, conf_level = 0.05,
#'                  common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
#'                  plotting=TRUE, print_progress=FALSE)


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
semiivreg = function(formula, data, propensity_formula=NULL, propensity_data = NULL,
                     ref_indiv =NULL, firststage_model = "probit",
                     est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
                     # Local polynomial parameters:
                     bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = 1/5,
                     kernel="gaussian", bw_subsamp_size = 10000, #NULL,
                     pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
                     fast_robinson1 = TRUE, fast_robinson2 = FALSE,
                     # Sieve parameters:
                     pol_degree_sieve = 5,
                     se_type = "HC1", conf_level = 0.05,
                     # Trim:
                     common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
                     weight_var = NULL,
                     plotting=TRUE, print_progress = FALSE) {


  # 0. Construction of the variables/objects to be used
  # ---------------------------------------------------
  # (i) Convert to data.frame + na.omit
  data = as.data.frame(data) # maybe problem with data loaded from Stata
  vars = all.vars(formula)
  if(!is.null(weight_var)) {
    data$WEIGHTS = data[[weight_var]]
  } else { data$WEIGHTS = rep(1, nrow(data)) }
  vars = c(vars, "WEIGHTS")
  data = subset(data, select=vars);
  data = na.omit(data)
  data = transform_factor(formula, data) # the variables which are specified as "factor(x)" are transformed into factors # done in order to ensure ordering with ref indiv
  data_orig = data;


  # Create reference individual if missing -> outdated
  # if(is.null(ref_indiv)) { ref_indiv_orig = create_ref_indiv(formula, data) } else { ref_indiv_orig = ref_indiv }

  # (ii) Construct transformed data
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


  # (iii) Transformed reference individual
  if(!is.null(ref_indiv)) { # If specific ref_indiv has been specified:

    ref_indiv_orig = ref_indiv # save the original version
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

  }
  # if no ref_indiv specified: compute the 'mean' effect in the sample.
  # Requires to compute the 'average' indiv.
  # But this is done AFTER eventual trimming of the first stage.

  # Remark: even for 'factors' just compute the mean, because it will give the mean effect in the sample.

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



  # (vi) re-add the weights
  data$WEIGHTS = data_orig$WEIGHTS
  Xdat$WEIGHTS = data_orig$WEIGHTS



  # 1. First stage: Propensity score estimation
  # --------------

  if(print_progress) { cat(sprintf("Estimating first stage... \r")) }

  if(is.null(propensity_formula)) { # by default, simple only include w0 and w1, additively
    # Formula:
    var_propensity = unique(c(var_w0, var_w1, var_covariates)) # unique to keep only one if some covariates in several potential outcomes with different effects
    var_propensity = var_propensity[which(var_propensity != "")]
    propensity_formula = as.formula(paste0(var_treatment, "~", paste0(var_propensity, collapse="+")))

    #if(is.null(propensity_data)) { propensity_data = data }

    # (i) Estimation:
    if(firststage_model == "probit") {
      propensity = glm(propensity_formula, family=binomial(link="probit"), data=data, weights = WEIGHTS)
    }
    if(firststage_model == "logit") {
      propensity = glm(propensity_formula, family=binomial(link="logit"), data=data, weights = WEIGHTS)
    }
    data$Phat = predict.glm(propensity, newdata=data, type="response")
    #Phat = propensity$fitted.values # No because possibly missing values will make everything wrong;
    Xdat$Phat = data$Phat

  } else {
    if(is.null(propensity_data)) { propensity_data = data_orig } # use the original data because maybe different transformation in 1st stage than in 2nd stage
    # /!\ crucial that propensity_data is ordered same way as data;

    if(firststage_model == "probit") { # if propensity formula is well provided, can run it directly without any modif to the original data
      propensity = glm(propensity_formula, family=binomial(link="probit"), data=propensity_data, weights=WEIGHTS)
    }
    if(firststage_model == "logit") {
      propensity = glm(propensity_formula, family=binomial(link="logit"), data=propensity_data, weights=WEIGHTS)
    }
    data$Phat = predict.glm(propensity, newdata=propensity_data, type="response") # add it to DATA -> should have the same ordering as data_orig
    Xdat$Phat = data$Phat
  }


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


  # (ii-c) Compute the 'average individual' (if not specified)
  if(is.null(ref_indiv)) { # If no specific ref_indiv -> Compute effects "at the mean"
    Xdat11 = Xdat[,which(colnames(Xdat) != "Phat")]
    ref_indiv = as.data.frame(t(colSums(Xdat11 * Xdat11$WEIGHTS) / sum(Xdat11$WEIGHTS))) # weighted.mean
    #ref_indiv = as.data.frame(t(colMeans(Xdat11)))  # no need for "na.rm=TRUE" because already na.omit on data at the beginning.
    ref_indiv$id = 1
    rm(Xdat11); #gc()
  }



  if(print_progress) {  cat(sprintf("First stage estimated.       \n")) }



  # 2. Second Stage:
  # ----------------
  # Several estimation methods here. By default do Robinson (1988) double residual regression for partially linear model.

  if(print_progress) {  cat(sprintf("2nd stage: Estimating MTR and MTE... \n")) }

  # Method 1: Local Polynomial Regressions
  # ---------
  if(est_method == "locpoly") {
    #if(var_covariates != "") { print("-- Common covariates treated as having generally different effect on D=0 and D=1") }
    # with Robinson, we don't allow the covariates to have the same effect on both D=0 and D=1

    # locpoly-1) Estimate covariates effects and kd(v)
    # ----------
    res = mtr_est_poly(data=data, seq_u=seq_u,
                       bw0 = bw0, bw1=bw1, bw_y0 = bw_y0, bw_y1=bw_y1, bw_method = bw_method, kernel=kernel,
                       bw_subsamp_size = bw_subsamp_size,
                       fast_robinson1 = fast_robinson1, fast_robinson2 = fast_robinson2,
                       se_type=se_type,
                       pol_degree1=pol_degree_locpoly1, pol_degree2=pol_degree_locpoly2,
                       var_outcome=var_outcome, var_treatment=var_treatment, var_w0=var_w0, var_w1=var_w1, var_covariates=var_covariates,
                       print_progress = print_progress)

    est0 = res$est0; est1 = res$est1;
    kv = res$kv;

    # Save bw as output as well:
    bw0 = res$bw0; bw1 = res$bw1;
    bw_y0 = res$bw_y0; bw_y1 = res$bw_y1; bw_mte = res$bw_mte;


    # locpoly-2) Compute MTR0, MTR1 and MTE from this for the reference individual
    # ----------
    eval_v = kv$v # same for k1.

    PREDICT = mtr_fun_poly(ref_indiv=ref_indiv, eval_v=eval_v, est0=est0, est1=est1, kv=kv, se_type=se_type, conf_level)
    RES = PREDICT$est
    deltaX = PREDICT$deltaX
    kv = PREDICT$kv

    # Graphical Parameter:
    if(fast_robinson2 == TRUE) { conf_band = FALSE } else { conf_band = TRUE }
    var_cov_2nd = NULL
  }







  # Method 2: "Sieve" method -> Flexible polynomial specification for Kappa_d
  # ---------
  # Either "sieve": general sieve with heterogenous treatment effects
  # Or "homogenous": sieve with homogenous treatment effects.

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



    # Sieve-1) Estimate E[Yd | D, Z, X, P]
    # ---------------------------------
    # -> gives estimate of Kappa_d(P) -> from kappa_d, obtains the MTR and MTE.
    # Estimate as a stacked regression (useful if want to impose same effect of covariates + practical afterwards)
    if(est_method == "sieve") {
      formula_control_function = paste0("I(1-", var_treatment, "):Kappa_fun(Phat, pol_degree) + I(", var_treatment, "):Kappa_fun(Phat, pol_degree)")
      formula_2nd_stage = as.formula(paste0(var_outcome, "~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function))
    }
    if(est_method == "homogenous") {
      formula_control_function = paste0("I(-(1-", var_treatment, ")*Phat/(1-Phat) + ", var_treatment, "):Kappa_homogenous_fun(Phat, pol_degree)")
      formula_2nd_stage= as.formula(paste0(var_outcome, "~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function))
    }

    # Estimation
    est = lm(formula_2nd_stage, data, weights = WEIGHTS)
    est_kappa = est;


    # Sieve-2) Marginal Treatment Responses and MTE estimation
    # ---------------------------------------------------------
    # From the estimated coefficients, can construct any MTR and MTE now.

    # (i) Extract Coefficients from the main estimation;
    coeff = coefficients(est)
    if(se_type == "nonrobust") {
      vcov = vcov(est)
    } else {
      vcov = vcovHC(est, type = se_type) # robust standard errors;
    }
    df = df.residual(est)
    t_value = qt(1-conf_level/2, df=df)
    if(length(which(is.na(coeff))) > 0) { warning("Some coefficients are NA. May be impossible to evaluate for some reference individual. ") }

    # Additional correction: if robust standard errors, no row/column for NA coeff. Add them back;
    ncoeff = length(coeff)
    index_nonNA = which(!is.na(coeff))
    if(se_type != "nonrobust") {
      vcov_base = vcov;
      vcov = matrix(NA, nrow=ncoeff, ncol=ncoeff)
      vcov[index_nonNA, index_nonNA] = vcov_base
    }

    coeff_kappa = coeff;
    vcov_kappa = vcov


    # (ii) Transform into the coefficients for the MTR0 and MTR1:
    # Known transformation from Kappa_d to kd(v) when Kappa is polynomial.
    mcoeffs = mtr_coeff(coeff=coeff_kappa, vcov=vcov_kappa,
                        var_treatment=var_treatment, est_method=est_method)
    coeff = mcoeffs$coeff; vcov = mcoeffs$vcov; names_var = mcoeffs$variables

    if(est_method == "homogenous") {
      exp_to_replace = paste0("I(-(1 - ", var_treatment, ") * Phat/(1 - Phat) + ", var_treatment, "):Kappa_homogenous_fun(Phat, pol_degree)")
      names_var = gsub(exp_to_replace, "kd(v): v^", names_var, fixed=TRUE)
    } # just for easy readability

    std_errors = sqrt(diag(vcov))
    t_values <- coeff / std_errors
    p_values <- 2 * pt(-abs(t_values), df = df)
    table_stacked = data.frame(Variable = names_var,
                               Estimate = unname(coeff),
                               Std_Error = unname(std_errors),
                               t_value = unname(t_values),
                               p_value = unname(p_values))


    # Also output a table of estimation results directly from here:
    res_est = mtr_est(coeff=coeff, vcov=vcov, names_var=names_var, var_treatment=var_treatment, df=df, est_method=est_method)
    est_mtr0 = res_est$table$mtr0
    est_mtr1 = res_est$table$mtr1
    est_mte = res_est$table$mte

    coeff_stacked = coeff
    coeff_mtr0 = res_est$coeff$mtr0; coeff_mtr1 = res_est$coeff$mtr1; coeff_mte = res_est$coeff$mte

    vcov_stacked = vcov
    vcov_mtr0 = res_est$vcov$mtr0; vcov_mtr1 = res_est$vcov$mtr1; vcov_mte = res_est$vcov$mte


    # (iii) Marginal Treatment Responses: MTRd(v, x) and MTE(v, x)
    # Using directly the coefficients from the correct functional form for the mtr.
    PREDICT = mtr_predict_sieve(coeff=coeff, vcov=vcov, ref_indiv=ref_indiv,
                                var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, pol_degree=pol_degree,
                                seq_u = seq_u, t_value = t_value, est_method=est_method, se_type=se_type)
    RES = PREDICT$est
    kv = PREDICT$kv
    deltaX = PREDICT$deltaX
  }




  if(print_progress) { cat(sprintf("Estimation complete.                                             \n")) }




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
  if(est_method == "locpoly") { conf_band1 = FALSE } else { conf_band1 = TRUE}
  mte_plot = mte_plot_fun(dat_plot, common_supp_round_mte, conf_band=conf_band1) # main curve, without standard errors - mte from diff between mtr1 - mtr0

  if(est_method == "locpoly" & fast_robinson2 == FALSE) { # Return a second MTE plot from the difference between MTR1 and MTR0.
    dat_plot2 = dat_plot;
    dat_plot2$mte = dat_plot2$mte2; conf_band2 = TRUE
    dat_plot2$mte_lwr = dat_plot2$mte2_lwr; dat_plot2$mte_upr = dat_plot2$mte2_upr
    mte_plot2 = mte_plot_fun(dat_plot2, common_supp_round_mte, conf_band2)
  } else { mte_plot2 = NULL }

  # By default: plot the supp_plot and the mte_plot:
  if(plotting == TRUE) {
    if(est_method == "locpoly") {
      grid.arrange(supp_plot, mte_plot2, ncol=2)
    } else {
      grid.arrange(supp_plot, mte_plot, ncol=2)
    }
  }


  # Also estimate the average MTE for the ref individual (= homogenous TE if homogenous)
  # Not necessarily = ATE if not full support;
  avg_MTE = mean(RES_plot$mte, na.rm=TRUE)




  # 4. Return objects
  # -----------------

  output_data = list(RES, deltaX, data_not_trimmed, ref_indiv, Xdat, data_orig, data); names(output_data) = c("RES", "deltaX", "data", "ref_indiv", "Xdat", "data_orig", "data_trimmed")
  output_plot = list(supp_plot, mtr_plot, mte_plot, mte_plot2); names(output_plot) = c("supp", "mtr", "mte", "mte2")
  output_supp = common_supp
  output_call = list(new_formula, formula, var_treatment, var_outcome, var_w0, var_w1, var_covariates,
                     var_cov_2nd, formula_X_orig, se_type, est_method, pol_degree_sieve,
                     pol_degree_locpoly1, pol_degree_locpoly2, kernel, fast_robinson1, fast_robinson2, bw_subsamp_size, conf_level);
  names(output_call) = c("formula", "formula_orig", "var_treatment", "var_outcome", "var_w0", "var_w1",
                         "var_covariates", "var_cov_2nd", "formula_X_orig", "se_type", "est_method",
                         "pol_degree_sieve", "pol_degree_locpoly1", "pol_degree_locpoly2",
                         "kernel", "fast_robinson1", "fast_robinson2", "bw_subsamp_size", "conf_level")

  if(est_method %in% c("sieve", "homogenous")) {
    output_estimate = list(est_mtr0, est_mtr1, est_mte, kv, propensity, avg_MTE, table_stacked, est_kappa);
    names(output_estimate) = c("mtr0", "mtr1", "mte", "kv", "propensity", "avg_MTE", "est", "est_kappa")

    output_coeff = list(coeff_stacked, coeff_mtr0, coeff_mtr1, coeff_mte, coeff_kappa)
    names(output_coeff) = c("coeff_stacked", "coeff_mtr0", "coeff_mtr1", "coeff_mte", "coeff_kappa")

    output_vcov = list(vcov_stacked, vcov_mtr0, vcov_mtr1, vcov_mte, vcov_kappa)
    names(output_vcov) = c("vcov_stacked", "vcov_mtr0", "vcov_mtr1", "vcov_mte", "vcov_kappa")

    output_bw = NA # not relevant for this method

  }

  if(est_method == "locpoly") {
    output_estimate = list(est0, est1, kv, propensity, avg_MTE);
    names(output_estimate) = c("est0", "est1", "kv", "propensity", "avg_MTE")
    # don't call it mtr0 and mtr1 because it's not the entire function. It's only the covariates effect

    output_bw = list(bw0, bw1, bw_y0, bw_y1, bw_mte, bw_method)
    names(output_bw) = c("bw0", "bw1", "bw_y0", "bw_y1", "bw_mte", "bw_method")

    output_coeff = NA; output_vcov = NA;
  }

  output = list(output_data, output_estimate, output_coeff, output_vcov,
                output_bw, output_plot, output_supp, output_call);
  names(output) = c("data", "estimate", "coeff", "vcov",
                    "bw", "plot", "supp", "call")

  if(est_method == "locpoly") {
    if (identical(parent.frame(), globalenv())) {
      message("Caution: the standard errors around the plot are not correct (too small). \nThese are standard errors around k1(v), k0(v) and k1(v) - k0(v). \nThey do not take the propensity score estimation, nor the Robinson 1st stage estimation of the effect of the covariates. \nFor proper standard errors, run the bootstrap in semiivreg_boot().")
    }
  }

  return(output)
}
















# ------------------------------------------------------
# Part 2. Bootstrap standard errors/confidence intervals
# ------------------------------------------------------
# The main function does not take into account that the propensity score Phat is estimated in a first stage to compute the standard errors
# We try to account for this with corrected standard errors using bootstrap in this function.



#' @rdname semiivreg
#' @usage semiivreg_boot(formula, Nboot=500, data, propensity_formula=NULL, ref_indiv =NULL,
#'                firststage_model="probit", est_method = "locpoly", se_type="HC1",
#'                bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = "rule-of-thumb",
#'                pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
#'                common_supp_trim=c(0,1), trimming_value = NULL,
#'                automatic_trim = FALSE, plotting=TRUE, conf_level = 0.05, CI_method = "curve", weight_var)
#' @param Nboot Number of bootstrap samples.
#' @param block_boot_var Variable on which to base the block bootstrap. By default, = NULL for standard bootstrap.
#' @param CI_method "delta" for delta method, "curve" for bootstrap the MTE curves directly. With est_method = "locpoly", only "curve" method is possible.
#' @param fast_robinson2 If TRUE, use locpoly from KernSmooth during the bootstrap. Default is TRUE to speed things up (because do not need to compute standard errors in bootstrap)
#' Set to FALSE if want to use epanechnikov kernel, or if want to use weights.
#' @param se_type Type of standard errors in main estimation and in each bootstrap replication. Can simplify by setting "nonrobust" which goes (slightly) faster.
#' @param print_progress_main Print progress of the main estimation or not.
#'
#' @export
semiivreg_boot = function(formula, Nboot=500, data, propensity_formula=NULL, propensity_data=NULL,
                          block_boot_var = NULL,
                          ref_indiv =NULL, firststage_model = "probit",
                          est_method = "locpoly", # "locpoly", "sieve", or "homogenous".
                          se_type="HC1", # "HC1",
                          bw0 = NULL, bw1 = NULL, bw_y0 = NULL, bw_y1 = NULL, bw_method = 1/5,
                          kernel="gaussian", bw_subsamp_size = 10000, #NULL,
                          pol_degree_locpoly1 = 1, pol_degree_locpoly2 = 2,
                          fast_robinson1 = TRUE, fast_robinson2 = TRUE,
                          pol_degree_sieve = 5, conf_level = 0.05,
                          common_supp_trim=c(0,1), trimming_value=NULL, automatic_trim=FALSE,
                          CI_method = "curve",
                          weight_var = NULL,
                          plotting=TRUE,
                          print_progress = TRUE, print_progress_main = TRUE) {


  # 1. semiivreg on the full sample
  # -------------------------------
  # to extract main estimate and eventually the support:

  cat(sprintf("Bandwidth and MTR/MTE estimation on main sample... \n"))

  captured_warnings <- character()  # To store captured warnings
  main_res <- withCallingHandlers(
    {semiivreg(formula=formula, data=data, propensity_formula=propensity_formula,
                          ref_indiv = ref_indiv, firststage_model = firststage_model,
                          est_method = est_method, se_type=se_type, bw0 = bw0, bw1 = bw1, bw_y0 = bw_y0, bw_y1 = bw_y1, bw_method = bw_method,
                          kernel=kernel, bw_subsamp_size = bw_subsamp_size,
                          pol_degree_locpoly1 = pol_degree_locpoly1, pol_degree_locpoly2 = pol_degree_locpoly2,
                          fast_robinson1 = fast_robinson1, fast_robinson2 = fast_robinson2,
                          pol_degree_sieve = pol_degree_sieve, conf_level = conf_level,
                          common_supp_trim = common_supp_trim, trimming_value = trimming_value, automatic_trim = automatic_trim,
                          weight_var = weight_var,
                          plotting=FALSE, print_progress = print_progress_main)
    },
    warning = function(w) {
      warning_message <- conditionMessage(w)
      captured_warnings <<- c(captured_warnings, warning_message)
    }
  )

  # If encounters error in the main regression, don't start with the bootstrap.
  if(length(grep("deficient fit", captured_warnings)) > 0) {
    # Stop the execution with a WARNING message, and still return an output (the main result), such that the user can still investigate what went wrong.
    message = "ERROR, EXECUTION STOPPED.\nMain model has a deficient fit.\nInvestigate with res$estimate$est0 and res$estimate$est1 and the corresponding res$data$ref_indiv."
    message(message)
    return(main_res)
  }



  # Extract ref_indiv (if not specified at the beginning, will always be the average indiv in the main (trimmed) sample)
  ref_indiv = main_res$data$ref_indiv

  # Updated formula
  transform_formula = main_res$call$formula # because use the transformed Xdata

  # Updated data:
  orig_data = main_res$data$data_orig
  propensity_data = orig_data
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

  seq_u = main_res$data$RES$Phat

  var_treatment = main_res$call$var_treatment; var_cov_2nd = main_res$call$var_cov_2nd


  cat(sprintf("\nBandwidth and MTR/MTE estimation on main sample: Done. \n"))



  # 2. Bootstrap
  # ------------
  #print("2/ Starting Bootstrap:")
  #progressbar <- txtProgressBar(min = 0, max = Nboot, style = 3)

  #BOOT = list()
  BOOTest = list()
  BOOTmtr0 = list(); BOOTmtr1 = list(); BOOTmte = list()
  BOOTDAT = list()

  k = 1; k1 = 1; k1_error = c()
  while(k <= Nboot) {

    #set.seed(1234*k) # just set seed outside


    if(is.null(block_boot_var)) {

      # (i) Standard Bootstrap
      bootstrap_indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE, prob=data$WEIGHTS)
      bootstrap_sample <- data[bootstrap_indices, ]
      bootstrap_propensity_sample = propensity_data[bootstrap_indices, ]
      # data and propensity_data are ordered exactly the same, does not cause problems

    } else {

      # (ii) Block Bootstrap
      # Based on block_boot_var:
      main_dat = propensity_data;
      unique_block_var = unique(main_dat[[block_boot_var]])
      N_block = length(unique_block_var)

      BOOTSTRAP_INDICES = list()
      for(b in 1:N_block) {
        lines_block = which(main_dat[[block_boot_var]] == unique_block_var[b])
        weights_block = main_dat$WEIGHTS[lines_block]
        size_block = length(lines_block)
        bootstrap_indices_b <- sample(1:size_block, size = size_block, replace = TRUE, prob=weights_block)
        bootstrap_indices_lines_b = lines_block[bootstrap_indices_b]
        BOOTSTRAP_INDICES[[b]] = bootstrap_indices_lines_b
      }
      bootstrap_indices = unlist(BOOTSTRAP_INDICES)

      bootstrap_sample <- data[bootstrap_indices, ]
      bootstrap_propensity_sample = propensity_data[bootstrap_indices, ]
    }


    bootstrap_sample$WEIGHTS = 1;
    bootstrap_propensity_sample$WEIGHTS = 1; # Reset the weights because already drawn with the weights for the bootstrap
    # So don't double weight!


    # semi-iv estimation
    captured_warnings_boot <- character()  # To store captured warnings
    boot_res <- withCallingHandlers(
      {invisible(semiivreg(formula=transform_formula, data=bootstrap_sample,
                           propensity_data = bootstrap_propensity_sample,
                           propensity_formula=propensity_formula,
                           ref_indiv = ref_indiv, firststage_model = firststage_model,
                           est_method = est_method, se_type=se_type, bw0 = bw0, bw1 = bw1, bw_y0 = bw_y0, bw_y1 = bw_y1, bw_method = bw_method,
                           kernel=kernel,
                           bw_subsamp_size = bw_subsamp_size, # irrelevant anyway because computed bw already;
                           fast_robinson1 = fast_robinson1, fast_robinson2 = fast_robinson2,
                           pol_degree_locpoly1 = pol_degree_locpoly1, pol_degree_locpoly2 = pol_degree_locpoly2,
                           pol_degree_sieve = pol_degree_sieve, conf_level = conf_level,
                           common_supp_trim = common_supp_trim_boot, trimming_value=NULL, automatic_trim=FALSE,
                           weight_var = NULL, # always NULL! Because the weights are placed INTO the bootstrap draw probabilities;
                           plotting=FALSE,
                           print_progress = FALSE))
      },
      warning = function(w) {
        warning_message <- conditionMessage(w)
        captured_warnings_boot <<- c(captured_warnings_boot, warning_message)
      }
    )
    #BOOT[[k]] = boot_res


    # If no error: store the results, update k. Otherwise, only update k1, check that not too many errors and go next.
    if(length(grep("deficient fit", captured_warnings_boot)) == 0) {

      # Storing the N_boot simulation would take too much memory if big datasets
      # Only extract the objects of interest:
      if(est_method %in% c("sieve", "homogenous")) {
        BOOTest[[k]] = boot_res$coeff$coeff_stacked
        BOOTmtr0[[k]] = boot_res$coeff$coeff_mtr0
        BOOTmtr1[[k]] = boot_res$coeff$coeff_mtr1
        BOOTmte[[k]] = boot_res$coeff$coeff_mte

      }
      if(est_method == "locpoly") {
        coeff0 = coefficients(boot_res$estimate$est0); names(coeff0) = paste0("Untreated_", names(coeff0))
        coeff1 = coefficients(boot_res$estimate$est1); names(coeff1) = paste0("Treated_", names(coeff1))
        coeff = c(coeff0, coeff1)
        BOOTest[[k]] = coeff
      }

      bootdat = boot_res$data$RES; bootdat$boot_id = k; BOOTDAT[[k]] = bootdat

      if(k %% 10 == 0) { gc() } # clean memory regularly

      cat(sprintf("Bootstrap Progress: %d/%d", k, Nboot), "\r")

      # Update k:
      k = k+1
      k1_error[k1] = 0;

    } else {
      # If error: only update the error list and check if not too many errors.
      k1_error[k1] = 1;

      # IF too many errors (more than 5% of the number of bootstrap required), stop the bootstrap:
      if(sum(k1_error) > Nboot/20) {
        message = paste0("ERROR, EXECUTION STOPPED.\nToo many deficiency fit in the Bootstrap samples (", sum(k1_error), " after ", k1, " attempts).\nLikely occurs because not both D are present with some covariates values.\nReturns the main model to investigate problems with res$main$estimate$est0 and res$main$estimate$est1 and the corresponding res$main$data$ref_indiv.\nCheck which combination of covariates have few observations with D=1 (or D=0) using res$main$data$data_trimmed.")
        message(message)
        return(list(main=main_res,
                    #BOOT=list(BOOTest=BOOTest, BOOTmte=BOOTmte, BOOTmtr0=BOOTmtr0, BOOTmtr1=BOOTmtr1),
                    error=k1_error))
      }

    }

    rm(boot_res)
    # General update of k1:
    k1 = k1+1

    #setTxtProgressBar(progressbar, k)
    cat(sprintf("Bootstrap Progress: %d/%d", k, Nboot), "\r")
  }



  # 3. Standard Errors around the coefficients
  # ------------------------------------------
  if(est_method %in% c("sieve", "homogenous")) {

    #BOOTest = list()
    #for(k in 1:Nboot) { BOOTest[[k]] = BOOT[[k]]$coeff$coeff_stacked }
    EST = do.call('rbind', BOOTest)
    meanEST = apply(EST, 2, mean);
    main_res_coeff = main_res$coeff$coeff_stacked #coefficients(main_res$estimate$est)
    vcov = cov(EST) # the main point is to obtain this vcov in order to do the delta method
    standarderror = sqrt(diag(vcov))
    COEFF = data.frame(estimate=main_res_coeff, std_error=standarderror, boot_mean_estimate=meanEST)
    # This will be used in the predict function after;


    # Table of coefficients for MTR and MTE directly:
    # if est_method = sieve or homogenous, can also DIRECTLY export table of results
    #BOOTmtr0 = list(); BOOTmtr1 = list(); BOOTmte = list()
    #for(k in 1:Nboot) {
    #  BOOTmtr0[[k]] = BOOT[[k]]$coeff$coeff_mtr0
    #  BOOTmtr1[[k]] = BOOT[[k]]$coeff$coeff_mtr1
    #  BOOTmte[[k]] = BOOT[[k]]$coeff$coeff_mte
    #}
    ESTmtr0 = do.call('rbind', BOOTmtr0)
    ESTmtr1 = do.call('rbind', BOOTmtr1)
    ESTmte = do.call('rbind', BOOTmte)

    meanESTmtr0 = apply(ESTmtr0, 2, mean)
    meanESTmtr1 = apply(ESTmtr1, 2, mean)
    meanESTmte = apply(ESTmte, 2, mean)

    main_res_mtr0 = main_res$coeff$coeff_mtr0
    main_res_mtr1 = main_res$coeff$coeff_mtr1
    main_res_mte = main_res$coeff$coeff_mte

    vcov_mtr0 = cov(ESTmtr0)
    vcov_mtr1 = cov(ESTmtr1)
    vcov_mte = cov(ESTmte)

    standarderror_mtr0 = sqrt(diag(vcov_mtr0))
    standarderror_mtr1 = sqrt(diag(vcov_mtr1))
    standarderror_mte = sqrt(diag(vcov_mte))

    MTR0 = data.frame(estimate=main_res_mtr0, std_error=standarderror_mtr0, boot_mean_estimate=meanESTmtr0)
    MTR1 = data.frame(estimate=main_res_mtr1, std_error=standarderror_mtr1, boot_mean_estimate=meanESTmtr1)
    MTE = data.frame(estimate=main_res_mte, std_error=standarderror_mte, boot_mean_estimate=meanESTmte)

  }



  if(est_method == "locpoly") {

    # Main coefficients:
    main_coeff0 = coefficients(main_res$estimate$est0); names(main_coeff0) = paste0("Untreated_", names(main_coeff0))
    main_coeff1 = coefficients(main_res$estimate$est1); names(main_coeff1) = paste0("Treated_", names(main_coeff1))
    main_res_coeff = c(main_coeff0, main_coeff1)

    #BOOTest = list()
    #for(k in 1:Nboot) {
    #  # est0:
    #  coeff0 = coefficients(BOOT[[k]]$estimate$est0); names(coeff0) = paste0("Untreated_", names(coeff0))
    #  coeff1 = coefficients(BOOT[[k]]$estimate$est1); names(coeff1) = paste0("Treated_", names(coeff1))
    #  coeff = c(coeff0, coeff1)
    #  BOOTest[[k]] = coeff
    #}

    EST = do.call('rbind', BOOTest)
    meanEST = apply(EST, 2, mean);
    vcov = cov(EST)
    standarderror = sqrt(diag(vcov))
    COEFF = data.frame(estimate=main_res_coeff, std_error=standarderror, boot_mean_estimate=meanEST)


    # Table of coefficients for MTR and MTE directly:
    # -> not possible with locpoly
    MTR0 = NULL; MTR1 = NULL; MTE = NULL
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

    # delta method uses the VCOV computed from the bootstrap above

    df = df.residual(main_res$estimate$est_kappa)
    t_value = qt(1-conf_level/2, df = df)
    #seq_u = seq(0, 1, by=0.001);
    pol_degree = pol_degree_sieve


    # Then simply return the main function results using this new vcov around the main estimates:
    PREDICT = mtr_predict_sieve(coeff=main_res_coeff,
                                vcov=vcov, ref_indiv=ref_indiv,
                                var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, pol_degree=pol_degree,
                                seq_u = seq_u, t_value = t_value, est_method=est_method, se_type=se_type)
    RES = PREDICT$est
  }

  # Remark: with est_method="locpoly"
  # -> could still report standard errors around the estimation of the coefficients in the propensity score stage
  #    and in the double residual regression stage.



  # Option 2: CI around the MTE curves directly
  # ---------
  if(CI_method == "curve") {

    #BOOTDAT = list()
    #for(k in 1:Nboot) { bootdat = BOOT[[k]]$data$RES; bootdat$boot_id = k; BOOTDAT[[k]] = bootdat }
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
    ## if(est_method %in% c("sieve", "homogenous")) {
    ##   mdat = mdat[, which(!colnames(mdat) %in% c("mtr0_lwr", "mtr0_upr", "mtr0_se", "mtr1_lwr", "mtr1_upr", "mtr1_se", "mte_lwr", "mte_upr", "mte_se"))]
    ## }
    mdat = mdat[, which(!colnames(mdat) %in% c("mtr0_lwr", "mtr0_upr", "mtr0_se", "mtr1_lwr", "mtr1_upr", "mtr1_se", "mte_lwr", "mte_upr", "mte_se"))]


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

  output_estimate = list(RES, COEFF, vcov, MTR0, MTR1, MTE); names(output_estimate) = c("est", "coeff", "vcov", "mtr0", "mtr1", "mte")
  output_boot_warnings = list(k1_error=k1_error) #

  output = list(output_main, output_plot, output_estimate, output_boot_warnings); names(output) = c("main", "plot", "estimate", "boot_warnings")
  return(output)
}










#' @rdname semiivreg
#' @usage semiiv_predict(semiiv, newdata, seq_v=NULL)
#' @param semiiv Object returns from a semiivreg estimation.
#' @param newdata New data for which to predict the MTE and MTR.
#' @param seq_v Sequence of v at which to predict the MTE and MTR. By default: NULL fits the default interval of the original semiivreg (equally space grid of proba, with step size of 0.001 on the common support).
#' @export
semiiv_predict = function(semiiv, newdata, seq_v=NULL) {
  # seq_v = NULL means that we use the original sequence of the sieve
  # can specify a single point for example if only want to estimate deltaX as an output!

  # 1. Format the newdata according to formula
  # This implies changing the factors, polynomials, log, etc. accordingly

  formula = semiiv$call$formula_orig
  formula_X_orig = semiiv$call$formula_X_orig
  data_orig = semiiv$data$data_orig

  name_all_X_orig = all.vars(formula_X_orig)
  Xdat_orig = subset(data_orig, select=name_all_X_orig)

  newdata1 = subset(newdata, select=name_all_X_orig) # remove potential "id"
  # Transform into factors as in data:
  for(j in 1:ncol(Xdat_orig)) {
    if(is.factor(Xdat_orig[,j])) {
      newdata1[,j] = factor(as.character(newdata1[,j]), levels=levels(Xdat_orig[,j]))
    }
  }
  ref = rbind(newdata1, Xdat_orig) # inflate the data just to have all the levels of the factors;
  res_ref = construct_data(formula, data=ref)

  ref_indiv = res_ref$data[1:nrow(newdata1),] # "ref_indiv" = newdata changed
  ref_indiv$id = 1:nrow(ref_indiv)
  rownames(ref_indiv) = 1:nrow(ref_indiv)




  # 2. Predict the MTR0, MTR1 and MTE for the new individuals
  est_method = semiiv$call$est_method;
  se_type = semiiv$call$se_type;
  if(is.null(seq_v)) {
    seq_u = sort(unique(semiiv$data$RES$Phat)) # the original sequence
  } else { seq_u = seq_v } # if pre-specified

  var_treatment = semiiv$call$var_treatment
  var_cov_2nd = semiiv$call$var_cov_2nd

  # 2.(i) Sieve or Homogenous

  if(est_method %in% c("sieve", "homogenous")) {
    coeff = semiiv$coeff$coeff_stacked
    vcov = semiiv$vcov$vcov_stacked
    pol_degree = semiiv$call$pol_degree_sieve
    conf_level = semiiv$call$conf_level
    df = df.residual(semiiv$estimate$est_kappa)
    t_value = qt(1-conf_level/2, df=df)

    # Support check:
    supp = semiiv$supp;
    if(length(which(seq_u < supp[1] | seq_u > supp[2])) > 0) {
      warning("Some values of seq_v are outside of the common support. kd(v) will be extrapolated for these values. ")
    } # Only returns a warning because can still compute with the functional form;

    PREDICT = mtr_predict_sieve(coeff=coeff, vcov=vcov, ref_indiv=ref_indiv,
                                var_treatment=var_treatment, var_cov_2nd=var_cov_2nd, pol_degree=pol_degree,
                                seq_u = seq_u, t_value = t_value, est_method=est_method, se_type=se_type)

  }


  # 2.(ii) Local Polynomial

  if(est_method == "locpoly") {
    conf_level = semiiv$call$conf_level

    est0 = semiiv$estimate$est0
    est1 = semiiv$estimate$est1
    kv = semiiv$estimate$kv

    # Support check:
    # Already done within the function (returns error if outside)

    PREDICT = mtr_fun_poly(ref_indiv=ref_indiv, eval_v=seq_u, est0=est0, est1=est1, kv=kv, se_type=se_type, conf_level=conf_level)
  }


  # 3. Output
  RES = PREDICT$est
  kv = PREDICT$kv
  deltaX = PREDICT$deltaX

  output = list(RES, deltaX, kv)
  names(output) = c("est", "deltaX", "kv")
  return(output)

}

NULL
