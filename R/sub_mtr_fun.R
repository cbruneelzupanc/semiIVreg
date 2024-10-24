#' @title MTE and MTR sub-estimation functions
#'
#' @description These functions allow to estimate the mte and mtr, and their confidence intervals, based on coefficients estimated from the main model in the main function.
#' More details can be found in Bruneel-Zupanc (2024). \cr
#' Different formulas must be applied depending on whether the treatment is homogenous or heterogenous.
#'
#'
#' @rdname mtr_fun
#' @usage mtr_fun = function(d, data, ref_indiv, seq_u,
#'                          bwd = NULL, bw_y = NULL, bw_method = "plug-in",
#'                          pol_degree1, pol_degree2, var_outcome, var_treatment, var_w0, var_w1, var_covariates)
#' @export
mtr_fun = function(d, data, ref_indiv, seq_u,
                   bwd = NULL, bw_y = NULL, bw_method = "plug-in",
                   pol_degree1, pol_degree2,
                   var_outcome, var_treatment, var_w0, var_w1, var_covariates) {

  # All regressions done on subsample D=d
  datad = data[which(data[[var_treatment]] == d),]
  if(d==0) { var_w = var_w0 } else { var_w = var_w1}
  var_covd = c(var_w, var_covariates); var_covd = unique(var_covd[which(var_covd != "")])
  formulad = as.formula(paste0("~ -1 + ", paste0(var_covd, collapse="+")))

  # Step 1. Double residual regression (Robinson, 1988)
  # -------
  Xd = model.matrix(formulad, data=datad)
  datd = cbind(subset(datad, select=var_outcome), Xd)
  Phatd = datad$Phat; yd = datd[[var_outcome]]
  residd = datd;

  BWD=c()
  for(j in 1:ncol(datd)) {
    #if(is.null(bwd)) { bwj = NULL } else { bwj = bwd[j] }
    bwj = bwd[j] # if is NULL will always put null
    resj = lpoly_regress(x=Phatd, y=datd[,j], bw=bwj, bw_method=bw_method, degree=pol_degree1, drv=0)
    # recommands degree = derivative order + 1; cf Fan and Gijbels (1996), footnote 19 of Carneiro et al 2011
    reg = resj$reg
    pred = approx(x=reg$x, y=reg$y, xout=Phatd)$y
    resid = datd[,j] - pred
    BWD[j] = resj$bw; residd[,j] = resid
  }

  formuladx = as.formula(paste0(var_outcome, "~ -1 + ."))
  estd = lm(formuladx, data=residd)

  # Step 2. Estimate Kappa_d(u) and then k_d(u)
  # -------
  # (i) Extract Yd net of the effect of the covariates
  ypredd = predict(estd, newdata=as.data.frame(Xd))
  ynetd = datd[[var_outcome]] - ypredd;

  # (ii) Estimate Kappa_d(u) with local polynomial regression
  est_Kappad = lpoly_regress(x=Phatd, y=ynetd, bw=bw_y, bw_method=bw_method, degree=pol_degree2, drv=0)
  # pol_degree2 because in the end I want to get the derivatives too.
  bw_y = est_Kappad$bw # if pre-specified will return the same value
  est_dKappad = lpoly_regress(x=Phatd, y=ynetd, bw = bw_y, degree=pol_degree2, drv=1) # same bw, no need to specify the method then
  # should be the same bandwidth in both, so don't recompute it

  # Estimation on the set of specified u:
  Kappad = approx(x=est_Kappad$reg$x, y=est_Kappad$reg$y, xout=seq_u)$y
  dKappad = approx(x=est_dKappad$reg$x, y=est_dKappad$reg$y, xout=seq_u)$y

  # (iii) Estimate kd(u):
  # See Andresen 2018, Table 2 for the formula; Can also derive manually;
  if(d==0) {
    ku = Kappad - (1-seq_u)*dKappad
  }
  if(d==1) {
    ku = Kappad + seq_u*dKappad
  }


  # Step 3. MTRd(u, wd, x)
  # ------
  # For each ref_indiv, predict the effect of the covariates - repeat it to have a data.frame as result
  REF = ref_indiv[rep(seq_len(nrow(ref_indiv)), each = length(seq_u)), ]
  REF$Phat = rep(seq_u, times=nrow(ref_indiv))

  # (i) Predict X beta_d + Wd alpha_d
  Walpha = predict(estd, newdata=REF)

  # (ii) Also repeat ku
  KU = rep(ku, nrow(ref_indiv))

  # (iii) MTRd
  mtrd = Walpha + KU


  # Results to return:
  RES = REF;
  RES[[paste0("mtr", d)]] = mtrd

  res_list = list(RES, estd, BWD, bw_y, Kappad, dKappad); names(res_list) = c("RES", "estd", "bwd", "bw_y", "Kappad", "dKappad")
  return(res_list)
}



#' @rdname mtr_fun
#' @param x Vector of x values
#' @param y Vector of y values
#' @param bw Pre-specified bandwidth
#' @param bw_method Method to compute the bandwidth (if bandwdith is NULL)
#' @param degree Degree of the polynomial: recommended to set to drv + 1
#' @param drv Derivative order of the function to be estimated.
#' @usage lpoly_regress(x, y, bw=NULL, bw_method="plug-in", degree, drv)
#' @export
lpoly_regress = function(x, y, bw=NULL, bw_method = "plug-in", degree, drv) {
  # To be adjusted if wants to change
  if(is.null(bw)) {
    if(bw_method == "plug-in") {
      xy <- cbind(x, y)
      xy <- xy[sort.list(xy[, 1L]), ]
      x1 <- xy[, 1L]
      y1 <- xy[, 2L]
      bw = dpill(x=x1, y=y1) # weirdly need to sort outside the function even though it's also done inside because can bug otherwise, no idea why!
    } # implement other methods later;
  }
  reg = locpoly(x=x, y=y, drv=drv, bandwidth=bw, degree=degree)
  res = list(reg, bw); names(res) = c("reg", "bw")
  return(res)
}





#'
#' @rdname mtr_fun
#' @usage mtr_fun_sieve(coeff, vcov, var_treatment, var_cov_2nd, ref_indiv, seq_u,
#'                homogenous=FALSE, pol_degree, conf_level, t_value, Xdat)
#' @export
mtr_fun_sieve = function(coeff, vcov, var_treatment, var_cov_2nd, ref_indiv, seq_u, homogenous=FALSE, pol_degree, conf_level, t_value, Xdat) {
  # Function computing the MTE and the MTR and reporting everything in a data.frame

  # Appropriate formula to construct the appropriate matrices
  if(homogenous==TRUE) {
    formula_control_function_mtr = paste0("I(1-", var_treatment, "+", var_treatment, "):ku_transform_homogenous_fun(Phat, pol_degree)")
    # The formula with homogenous effect derivation implies that k1(u) = k0(u) by construction;
  } else {
    formula_control_function_mtr = paste0("I(1-", var_treatment, "):kdu_transform_fun(Phat, 0, pol_degree) + I(", var_treatment, "):kdu_transform_fun(Phat, 1, pol_degree)") # general formula
  }
  #formula_MTR = as.formula(paste0("~ - 1 +", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function_mtr))
  formula_MTR = as.formula(paste0("~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function_mtr))


  # Corresponding data:
  REF = ref_indiv[rep(seq_len(nrow(ref_indiv)), each = length(seq_u)), ]
  REF$Phat = rep(seq_u, times=nrow(ref_indiv))

  # if includes polynomial transform or splines in the function, may lead to issues in the prediction here because ref_indiv usually contains only one level
  # small trick to avoid the issue: predict on an artificially inflated dataset, including all the initial covariates using Xdat
  Xdat = Xdat[,order(colnames(Xdat))]; REF = REF[, order(colnames(REF))] # ensures same order
  REFF = rbind(REF, Xdat)
  ntrue_rows = nrow(REF)
  REFF = droplevels(REFF) # very important to droplevels! because the lm dropped them as well;

  REF1 = REFF; REF1[[var_treatment]] = 1
  REF0 = REFF; REF0[[var_treatment]] = 0

  RES = REF # to save results

  # Possible that some coefficients (e.g., some fixed effects) are NA -  but not a problem if does not concern the reference individual
  coeff_base = coeff; vcov_base = vcov
  index_NA = which(is.na(coeff_base))
  if(length(index_NA > 0)) {
    coeff = coeff[-index_NA]
    vcov = vcov[-index_NA, -index_NA]
  }

  # (i) MTR computation
  # MTR1:
  mat1 = model.matrix(formula_MTR, data=REF1);
  mat1 = mat1[1:ntrue_rows, ] # reselect the proper rows once the bugs are avoided;
  if(length(index_NA) > 0) {
    # check that ref individual does not have a characteristics in the missing column:
    missing_col1 = mat1[, index_NA]
    if(length(which(missing_col1 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat1 = mat1[, -index_NA]
  }
  mtr1 = mat1%*%coeff;

  # MTR0:
  mat0 = model.matrix(formula_MTR, data=REF0);
  mat0 = mat0[1:ntrue_rows, ]
  if(length(index_NA) > 0) {
    # check that ref individual does not have a characteristics in the missing column:
    missing_col0 = mat0[, index_NA]
    if(length(which(missing_col0 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat0 = mat0[, -index_NA]
  }
  mtr0 = mat0%*%coeff;



  # (ii) Marginal Treatment Effect: MTE(u, x)
  mat_mte = mat1 - mat0 # element by element addition
  mte = mat_mte %*% coeff
  # equivalent to: mte = mtr1-mtr0

  # (iii) Confidence Intervals
  # MTR0:
  se0 = sqrt(rowSums((mat0 %*% vcov)*mat0)) # equivalent to doing, for each row i: r0 = matrix(mat0[i,], nrow=1); sqrt(r0 %*% vcov %*% t(r0))
  mtr0_lwr = mtr0 - t_value*se0;
  mtr0_upr = mtr0 + t_value*se0

  # MTR1:
  se1 = sqrt(rowSums((mat1 %*% vcov)*mat1))  # equivalent to doing, for each row i: r0 = matrix(mat0[i,], nrow=1); sqrt(r0 %*% vcov %*% t(r0))
  mtr1_lwr = mtr1 - t_value*se1;
  mtr1_upr = mtr1 + t_value*se1

  # MTE:
  se_mte = sqrt(rowSums((mat_mte %*% vcov)*mat_mte)) # equivalent to doing, for each row i: r0 = matrix(mat0[i,], nrow=1); sqrt(r0 %*% vcov %*% t(r0))
  mte_lwr = mte - t_value*se_mte;
  mte_upr = mte + t_value*se_mte

  # Save:
  RES$mtr1 = mtr1;
  RES$mtr0 = mtr0;
  RES$mte = mte
  RES$mtr0_lwr = mtr0_lwr; RES$mtr0_upr = mtr0_upr; RES$mtr0_se = se0
  RES$mtr1_lwr = mtr1_lwr; RES$mtr1_upr = mtr1_upr; RES$mtr1_se = se1;
  RES$mte_lwr = mte_lwr; RES$mte_upr = mte_upr; RES$mte_se = se_mte


  return(RES)
}

NULL
