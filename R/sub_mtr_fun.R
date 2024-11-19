#' @title MTE and MTR sub-estimation functions
#'
#' @description These functions allow to estimate the mte and mtr, and their confidence intervals, based on coefficients estimated from the main model in the main function.
#' More details can be found in Bruneel-Zupanc (2024). \cr
#' Different formulas must be applied depending on whether the treatment is homogenous or heterogenous.
#'
#'
#' @rdname mtr_fun
#' @usage mtr_est_poly(d, data, seq_u,
#'                     bwd = NULL, bw_y = NULL, bw_method = "plug-in",
#'                     pol_degree1, pol_degree2,
#'                     var_outcome, var_treatment, var_w0, var_w1, var_covariates)
#' @export
mtr_est_poly = function(d, data, seq_u,
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
    if(length(unique(datd[,j])) == 1) { stop(paste0("No variation in variable ", colnames(datd)[j], " in the subsample d = ", d, " of the trimmed subsample.")) }
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

  kv=data.frame(v=seq_u); # to return the function kd(v)
  kv[[paste0("k", d)]] = ku


  res_list = list(Kappad, dKappad, kv, estd, BWD, bw_y)
  names(res_list) = c("Kappad", "dKappad", "kv", "estd", "bwd", "bw_y")
  return(res_list)

}


#' @rdname mtr_fun
#' @usage mtr_fun_poly(ref_indiv, est0, est1, k0, k1)
#' @export
mtr_fun_poly = function(ref_indiv, eval_v, est0, est1, k0, k1) {

  # Main prediction: predict mtr, mte

  # For each ref_indiv, predict the effect of the covariates - repeat it to have a data.frame as result
  REF = ref_indiv[rep(seq_len(nrow(ref_indiv)), each = length(eval_v)), ]
  REF$Phat = rep(eval_v, times=nrow(ref_indiv))

  # Effect of covariates
  delta0x = predict(est0, newdata=REF)
  delta1x = predict(est1, newdata=REF)

  X0_pred = predict(est0, newdata=REF, se.fit=TRUE)
  X1_pred = predict(est1, newdata=REF, se.fit=TRUE)
  delta0x = X0_pred$fit
  delta1x = X1_pred$fit
  se_delta0x = X0_pred$se.fit
  se_delta1x = X1_pred$se.fit

  # Also repeat ku:
  range_v = c(min(k0$v), max(k0$v))
  if(length(which(eval_v > range_v[2] | eval_v < range_v[1])) > 0) {
    stop("Evaluation value for V is out of the estimated range.")
  }

  k0_est = approx(x=k0$v, y=k0$k0, xout=eval_v)$y
  k1_est = approx(x=k1$v, y=k1$k1, xout=eval_v)$y
  kv = data.frame(v=eval_v, k0=k0_est, k1=k1_est)

  # repeat it to fit the REF data;
  K0 = rep(k0_est, times=nrow(ref_indiv))
  K1 = rep(k1_est, times=nrow(ref_indiv))

  # MTRd:
  mtr0 = delta0x + K0
  mtr1 = delta1x + K1
  mte = mtr1 - mtr0

  # Save
  RES = REF;
  RES$mtr0 = mtr0; RES$mtr1 = mtr1
  RES$mte = mte
  RES$k0 = K0; RES$k1 = K1
  RES$delta0X = delta0x; RES$delta1X = delta1x
  RES$se_delta0X = se_delta0x; RES$se_delta1X = se_delta1x


  # Also return the effect of the covariates (but only one prediction by observation)
  deltaX = ref_indiv;
  delta0X_pred = predict(est0, newdata=ref_indiv, se.fit=TRUE)
  delta1X_pred = predict(est1, newdata=ref_indiv, se.fit=TRUE)
  deltaX$delta0X = delta0X_pred$fit; deltaX$delta1X = delta1X_pred$fit;
  deltaX$se_delta0X = delta0X_pred$se.fit; deltaX$se_delta1X = delta1X_pred$se.fit;


  res_list = list(RES, deltaX, kv); names(res_list) = c("est", "deltaX", "kv")
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
lpoly_regress = function(x, y, bw=NULL, bw_method = "rule-of-thumb", degree, drv) {
  # To be adjusted if wants to change
  if(is.null(bw)) {
    if(bw_method == "rule-of-thumb") {
      bw = thumbBw(x, y, deg=degree, kernel=gaussK, weig = rep(1, length(y))) # might yield weird results when degree = 2
      #if(is.nan(bw)) { stop("Error no variation in variable in this subsample.") }
      if(bw > max(x) - min(x)) {
        bw = thumbBw(x, y, deg=degree-1, kernel=gaussK, weig = rep(1, length(y))) # reduce the degree
        if(bw > max(x) - min(x)) { # if still greater:
          warning("Rule-of-Thumb bw too large... Setting to (max(x)-min(x))/100 arbitrarily")
          bw = (max(x) - min(x))/100
        }
      }
    }
    if(bw_method == "plug-in") {
      # Gets long if too many observations... By default compute on subsample
      sub_id = sample(1:length(x), min(4999, length(x))) # 5000 corresponds to .maxEvalPts in locpol
      suby = y[sub_id]; subx = x[sub_id]
      # requires degree to be odd number (not even)
      if(degree %% 2 == 0) {
        warning("Degree should be an odd number for plug-in bandwidth selection: reduced by 1 for bw computation.")
        degree = degree - 1
      }
      bw <- pluginBw(subx, suby, deg=degree, kernel=gaussK, weig = rep(1, length(suby)))
    }
    if(bw_method == "cv") {
      # Gets long if too many observations... By default compute on subsample
      sub_id = sample(1:length(x), min(4999, length(x))) # 5000 corresponds to .maxEvalPts in locpol
      suby = y[sub_id]; subx = x[sub_id]
      range = range(subx); upper =  (range[2] - range[1])/2
      lower = (range[2] - range[1])/1000
      interval = c(lower, (range[2] - range[1])/2)
      bw = regCVBwSelC(x=subx, y=suby, deg=degree, kernel=gaussK,weig=rep(1,length(suby)),
                       interval = interval)
      # bw very sensitive to interval choice;
    }
    ## # Remark: if want to allow for higher maxevalpts:
    ## unlockBinding(".maxEvalPts", asNamespace("locpol"))
    ## assign(".maxEvalPts", length(x) + 10, envir = asNamespace("locpol")) # Overwrite the value
    ## lockBinding(".maxEvalPts", asNamespace("locpol")) # Lock the binding again to avoid unintended modifications
  }
  reg = locpoly(x=x, y=y, drv=drv, bandwidth=bw, degree=degree)
  res = list(reg, bw); names(res) = c("reg", "bw")
  return(res)
}















#' @rdname mtr_fun
#' @param coeff Vector of coefficients from sieve/homogenous estimates of E(Yd|D=d, W, X, P)
#' @param vcov Covariance matrix of these coefficients
#' @param var_cov_2nd Names of the covariates and semi-IVs
#' @usage mtr_coeff(coeff, vcov, var_cov_2nd, est_method="sieve")
#' @return Returns the raw (stacked regression) coefficients and covariance matrix corresponding to mtr0, mtr1 and mte function.
#' @export
#'
#'
mtr_coeff = function(coeff, vcov, var_treatment, est_method="sieve") {
  # Only applicable for est_method "sieve" or "homogenous".
  # Returns the exact functional form coefficients estimated for mtr0, mtr1 and mte
  # with the corresponding covariance matrix around the estimation


  # 0. Identify the NA -> will be really important for the transformation and final output
  ncoeff = length(coeff); names_var = names(coeff)
  index_nonNA = which(!is.na(coeff)); index_NA = which(is.na(coeff))
  # pre-identify some lines as well:
  linekd = grep("Kappa_fun(Phat, pol_degree)", names_var, fixed=TRUE)
  if(length(which(linekd %in% index_NA)) > 0) {
    stop("The Kappa function is not defined for some coefficients. Please check the estimation.")
  }
  line_coeff1 = grep(paste0("I(", var_treatment, "):"), names_var, fixed=TRUE)
  line_coeff0 = grep(paste0("I(1 - ", var_treatment, "):"), names_var, fixed=TRUE)
  line_treated = which(names_var == paste0("I(", var_treatment, ")"))


  # 1. Create the transformation matrix A
  A = matrix(0, nrow=ncoeff, ncol=ncoeff) # transformation matrix
  diag(A) = 1 # by default, most coefficients are unchanged

  # 2. Update for the correct transformation depending on the specification
  # 2.A) Sieve polynomial:
  if(est_method == "sieve") {

    # (i) Change for MTR1:
    # If Kappa1(p) = alpha1*p + alpha2*p^2 + alpha3*p^3
    # Then k1(u) = 2alpha1*u + 3alpha2*u^2 + 4alpha3*u^3

    # What is the index of the coefficients of Kappa1?
    linek1 = linekd[which(linekd %in% line_coeff1)]

    multiplier1 = 2;
    for(j in 1:length(linek1)) {
      index_coeff_kappa1 = linek1[j]
      A[index_coeff_kappa1, index_coeff_kappa1] = multiplier1
      multiplier1 = multiplier1 + 1
    }

    # (ii) Change for MTR0 (and intercept):
    # E[Y0 | D=0, W0, X, P] = mu0 + W0 gamma0 + X beta0 + Kappa0(p).
    # Kappa0(p) = alpha1*p + alpha2*p^2 + alpha3*p^3 + ...
    # Then mtr0(u, x, w0) = mu0 - alpha1 +
    #                       w0 gamma0 + xbeta0 +
    #                       +2(alpha1-alpha2)*u + 3 (alpha2-alpha3)*u^2 + 4 alpha3*u^3.
    # The transformation is more complex!

    linek0 = linekd[which(linekd %in% line_coeff0)]
    multiplier0 = 2;

    # First coeff: alpha1. -> in the intercept and multiplying 2*p
    A[1, linek0[1]] = -1 # New intercept for mtr0
    # New coeff multiplying p (order 1):
    A[linek0[1], linek0[1]] = multiplier0

    multiplier0 = multiplier0 + 1
    if(length(linek0) >= 2) {
      for(j in 2:length(linek0)) {
        A[linek0[j-1], linek0[j]] = - (multiplier0-1) # coeff at the order p^{j-1}
        A[linek0[j], linek0[j]] = multiplier0 # coeff at the order p^j
        multiplier0 = multiplier0 + 1
      }
    }

    # (iv) Last adjustment: intercept (again)
    # Problem: the intercept does NOT change for mtr1...
    #          but it does here if you change the main coefficient!
    # To enforce the previous coefficient, need to readd what we subtracted from the intercept (to the I(d) coeff)
    A[line_treated, linek0[1]] = 1 # to "re-add" it to the intercept
  }



  # 2.B) Homogenous TE
  if(est_method == "homogenous") {
    # with homogenous, quite simple, k1(u) = k0(u) = d/dp Kappa1(p)*p, where Kappa1(p) includes a constant here

    linekd_homo = grep("Kappa_homogenous_fun(Phat, pol_degree)", names_var, fixed=TRUE)
    if(length(which(linekd_homo %in% index_NA)) > 0) {
      stop("The Kappa function is not defined for some coefficients. Please check the estimation.")
    }

    multiplier_homogenous = 1;
    for(j in 1:length(linekd_homo)) {
      A[linekd_homo[j], linekd_homo[j]] = multiplier_homogenous
      multiplier_homogenous = multiplier_homogenous + 1
    }

  }


  # 3. Updated coefficients and vcov to return
  # Need to work with non-NA lines
  A_nonNA = A[index_nonNA, index_nonNA]
  coeff_nonNA = coeff[index_nonNA]
  vcov_nonNA = vcov[index_nonNA, index_nonNA]

  new_coeff_nonNA = A_nonNA %*% coeff_nonNA
  new_vcov_nonNA = A_nonNA %*% vcov_nonNA %*% t(A_nonNA)

  # Merge in the global values including NA
  new_coeff = rep(NA, ncoeff); new_vcov = matrix(NA, nrow=ncoeff, ncol=ncoeff)
  new_coeff[index_nonNA] = new_coeff_nonNA;
  new_vcov[index_nonNA, index_nonNA] = new_vcov_nonNA
  names_coeff = names_var


  new_coeff = as.numeric(new_coeff); names(new_coeff) = names_coeff
  colnames(new_vcov) = names_coeff; rownames(new_vcov) = names_coeff

  res = list(new_coeff, new_vcov, names_coeff)
  names(res) = c("coeff", "vcov", "variables")
  return(res)

}





#' @rdname mtr_fun
#' @param coeff Vector of (stacked) coefficients for mtr0, mtr1.
#' @param vcov Covariance matrix of these coefficients.
#' @param names_var Names of the variables corresponding to the coefficients.
#' @param df Degrees of freedom for the p-values.
#' @usage mtr_est(coeff, vcov, names_var, df)
#' @return Returns the mtr0, mtr1 and mte estimates tables with their standard errors and p-values. Also exports the corresponding estimates and vcov matrices.
#' @export
#'
#'
mtr_est = function(coeff, vcov, names_var, var_treatment, df, est_method) {

  # I. Table formation for General heterogenous TE Sieve:
  if(est_method == "sieve") {

    names_var = gsub("Kappa_fun(Phat, pol_degree)", "kd(v): v^", names_var, fixed=TRUE)
    linekd = grep("kd(v): v^", names_var, fixed=TRUE)
    line_intercept = 1;


    # 1. MTR0 table

    line_coeff0 = grep(paste0("I(1 - ", var_treatment, "):"), names_var, fixed=TRUE)
    line_kd0 = linekd[which(linekd %in% line_coeff0)]
    line_coeff0 = line_coeff0[which(!line_coeff0 %in% line_kd0)]
    line0 = c(line_intercept, line_kd0, line_coeff0)

    names_var0 = names_var[line0]
    names_var0 = gsub(paste0("I(1 - ", var_treatment, "):"), "", names_var0, fixed=TRUE)
    coeff0 = unname(coeff[line0]); vcov0 = unname(vcov[line0, line0])


    std_errors0 = sqrt(diag(vcov0))
    t_values0 <- coeff0 / std_errors0
    p_values0 <- 2 * pt(-abs(t_values0), df = df)

    table0 = data.frame(Variable = names_var0, Estimate = coeff0, Std_Error = std_errors0, t_value = t_values0, p_value = p_values0)



    # 2. MTR1 table

    line_coeff1 = grep(paste0("I(", var_treatment, "):"), names_var, fixed=TRUE)
    line_coefftreated = which(names_var == paste0("I(", var_treatment, ")")) # Also the pure I(d)
    line_kd1 = linekd[which(linekd %in% line_coeff1)]
    line_coeff1 = line_coeff1[which(!line_coeff1 %in% line_kd1)]
    line1 = c(line_intercept, line_coefftreated, line_kd1, line_coeff1)

    names_var1 = names_var[line1]
    names_var1 = gsub(paste0("I(", var_treatment, "):"), "", names_var1, fixed=TRUE)
    coeff1 = coeff[line1]; vcov1 = vcov[line1, line1]



    # "Problem": need an additional transformation because two effects in the constant here
    ncoeff1 = length(coeff1)
    A1 = matrix(0, nrow=ncoeff1, ncol=ncoeff1)
    diag(A1) = 1
    line1_intercept = 1
    line1_treated = which(names_var1 == paste0("I(", var_treatment, ")"))
    A1[line1_intercept, line1_treated] = 1;

    #A1 = A1[-line1_treated,]
    #coeff1 = A1%*%coeff1
    #vcov1 = A1 %*% vcov1 %*% t(A1)


    # To apply the transformation, pay attention in case some missing coefficients
    index_nonNA1 = which(!is.na(coeff1));
    A1_nonNA = A1[index_nonNA1, index_nonNA1]
    coeff1_nonNA = coeff1[index_nonNA1]
    vcov1_nonNA = vcov1[index_nonNA1, index_nonNA1]
    names_var1_nonNA = names_var1[index_nonNA1]

    new_coeff1_nonNA = A1_nonNA %*% coeff1_nonNA
    new_vcov1_nonNA = A1_nonNA %*% vcov1_nonNA %*% t(A1_nonNA)
    new_coeff1 = rep(NA, ncoeff1); new_vcov1 = matrix(NA, nrow=ncoeff1, ncol=ncoeff1)
    new_coeff1[index_nonNA1] = new_coeff1_nonNA;
    new_vcov1[index_nonNA1, index_nonNA1] = new_vcov1_nonNA
    coeff1 = new_coeff1; vcov1 = new_vcov1

    # Remove the I(d) line/columns:
    coeff1 = coeff1[-line1_treated]
    vcov1 = vcov1[-line1_treated, -line1_treated]
    names_var1 = names_var1[-line1_treated]

    std_errors1 = sqrt(diag(vcov1))
    t_values1 <- coeff1 / std_errors1
    p_values1 <- 2 * pt(-abs(t_values1), df = df)

    table1 = data.frame(Variable = names_var1, Estimate = coeff1, Std_Error = std_errors1, t_value = t_values1, p_value = p_values1)



    # 3. MTE Table
    # Build the MTE as a transformation of the original coefficients
    # For each variable, need to determine if it is w0, w1, x or kd(u) related.

    ncoeff = length(coeff)

    line_coeff0 = grep(paste0("I(1 - ", var_treatment, "):"), names_var, fixed=TRUE)
    names0 = names_var[line_coeff0]
    names0 = gsub(paste0("I(1 - ", var_treatment, "):"), "", names0, fixed=TRUE)

    line_coeff1 = grep(paste0("I(", var_treatment, "):"), names_var, fixed=TRUE)
    names1 = names_var[line_coeff1]
    names1 = gsub(paste0("I(", var_treatment, "):"), "", names1, fixed=TRUE)

    line_w1 = line_coeff1[which(!names1 %in% names0)];
    names_w1 = names1[which(!names1 %in% names0)]
    names_w1 = paste0("+ W1: ", names_w1)
    line_w0 = line_coeff0[which(!names0 %in% names1)];
    names_w0 = names0[which(!names0 %in% names1)]
    names_w0 = paste0("- W0: ", names_w0) # to make it clear that it is multiplied by -1

    line_x1 = line_coeff1[which(names1 %in% names0)]
    line_x0 = line_coeff0[which(names0 %in% names1)]
    names_x_kd = names1[which(names1 %in% names0)]

    # First line: simply the difference of intercept, directly given by I(d) coeff
    line_treated = which(names_var == paste0("I(", var_treatment, ")"))
    mattreated = matrix(0, nrow=1, ncol=ncoeff)
    mattreated[, line_treated] = 1
    name_treated = "(Intercept)"

    # Lines for W1 and W0:
    matw1 = matrix(0, nrow=length(line_w1), ncol=ncoeff)
    for(j in 1:length(line_w1)) {matw1[j, line_w1[j]] = 1}
    matw0 = matrix(0, nrow=length(line_w0), ncol=ncoeff)
    for(j in 1:length(line_w0)) {matw0[j, line_w0[j]] = -1} # - 1 for w0.

    # Lines for covariates (and kdu): always find the matched line
    MATX = list();
    if(length(names_x_kd) > 0) {
      for(j in 1:length(names_x_kd)) {
        name_i = names_x_kd[j]
        line_x0_i = line_coeff0[which(names0 == name_i)]
        line_x1_i = line_coeff1[which(names1 == name_i)]

        mat = matrix(0, nrow=1, ncol=ncoeff)
        mat[, line_x1_i] = 1;
        mat[, line_x0_i] = -1;

        MATX[[j]] = mat
      }
    }
    matx = do.call('rbind', MATX)

    A_mte = rbind(mattreated, matw1, matw0, matx) # transformation matrix to apply to coeff
    names_var_mte = c(name_treated, names_w1, names_w0, names_x_kd)
    ncoeff_mte = length(names_var_mte)
    #coeff_mte = A_mte%*%coeff
    #vcov_mte = A_mte %*% vcov %*% t(A_mte)


    # Deal with NA coefficients:
    index_nonNA_mte = which(!is.na(coeff))
    index_NA_mte = which(is.na(coeff))
    index_row_ok = 1:nrow(A_mte); index_row_NA = list()
    if(length(index_NA_mte) > 0) {
      for(j in 1:length(index_NA_mte)) {
        index_row_NA[[j]] = which(A_mte[, index_NA_mte[j]] != 0)
      }
    }
    index_row_NA = unique(do.call('c', index_row_NA))
    index_row_ok = index_row_ok[which(!index_row_ok %in% index_row_NA)]

    A_mte2 = A_mte[index_row_ok, index_nonNA_mte]
    coeff_nonNA = coeff[index_nonNA_mte]
    vcov_nonNA = vcov[index_nonNA_mte, index_nonNA_mte]

    coeff_mte_nonNA = A_mte2%*%coeff_nonNA
    vcov_mte_nonNA = A_mte2 %*% vcov_nonNA %*% t(A_mte2)

    # Rebuild the full coefficients and vcovs
    new_coeff_mte = rep(NA, ncoeff_mte); new_vcov_mte = matrix(NA, nrow=ncoeff_mte, ncol=ncoeff_mte)
    new_coeff_mte[index_row_ok] = coeff_mte_nonNA;
    new_vcov_mte[index_row_ok, index_row_ok] = vcov_mte_nonNA
    coeff_mte = new_coeff_mte; vcov_mte = new_vcov_mte





    # Re-order to have kd(v) at the beginning:
    linekd_mte = grep("kd(v): v^", names_var_mte, fixed=TRUE)
    line_w1_mte = grep("W1:", names_var_mte, fixed=TRUE)
    line_w0_mte = grep("W0:", names_var_mte, fixed=TRUE)
    line_treated_mte = 1
    line_x_mte = which(!1:length(names_var_mte) %in% c(line_treated_mte, linekd_mte, line_w1_mte, line_w0_mte))

    order_mte = c(line_treated_mte, linekd_mte, line_w1_mte, line_w0_mte, line_x_mte)

    coeff_mte = coeff_mte[order_mte]
    vcov_mte = vcov_mte[order_mte, order_mte]
    names_var_mte = names_var_mte[order_mte]


    std_errors_mte = sqrt(diag(vcov_mte))
    t_values_mte <- coeff_mte / std_errors_mte
    p_values_mte <- 2 * pt(-abs(t_values_mte), df = df)

    table_mte = data.frame(Variable = names_var_mte, Estimate = coeff_mte, Std_Error = std_errors_mte, t_value = t_values_mte, p_value = p_values_mte)

  }








  # II. Table formation for Homogenous TE
  if(est_method == "homogenous") {

    exp_to_replace = paste0("I(-(1 - ", var_treatment, ") * Phat/(1 - Phat) + ", var_treatment, "):Kappa_homogenous_fun(Phat, pol_degree)")
    names_var = gsub(exp_to_replace, "kd(v): v^", names_var, fixed=TRUE)
    linekd = grep("kd(v): v^", names_var, fixed=TRUE)
    line_intercept = 1;


    # 1. MTR0 table

    line_coeff0 = grep(paste0("I(1 - ", var_treatment, "):"), names_var, fixed=TRUE) # don't include the kd in the homogenous case
    line_kd0 = linekd
    line0 = c(line_intercept, line_kd0, line_coeff0)

    names_var0 = names_var[line0]
    names_var0 = gsub(paste0("I(1 - ", var_treatment, "):"), "", names_var0, fixed=TRUE)
    coeff0 = coeff[line0]; vcov0 = vcov[line0, line0]

    # "Problem": need an additional transformation because 3 effects in the constant here
    # (Intercept) + I(d) + kd(v): v^const.
    ncoeff0 = length(coeff0)
    A0 = matrix(0, nrow=ncoeff0, ncol=ncoeff0)
    diag(A0) = 1
    line0_intercept = 1
    line0_const = which(names_var0 == "kd(v): v^const")
    A0[line0_intercept, line0_const] = 1;

    #A0 = A0[-line0_const,]
    #coeff0 = A0%*%coeff0
    #vcov0 = A0 %*% vcov0 %*% t(A0)

    # To apply the transformation, pay attention in case some missing coefficients
    index_nonNA0 = which(!is.na(coeff0));
    A0_nonNA = A0[index_nonNA0, index_nonNA0]
    coeff0_nonNA = coeff0[index_nonNA0]
    vcov0_nonNA = vcov0[index_nonNA0, index_nonNA0]
    names_var0_nonNA = names_var0[index_nonNA0]

    new_coeff0_nonNA = A0_nonNA %*% coeff0_nonNA
    new_vcov0_nonNA = A0_nonNA %*% vcov0_nonNA %*% t(A0_nonNA)
    new_coeff0 = rep(NA, ncoeff0); new_vcov0 = matrix(NA, nrow=ncoeff0, ncol=ncoeff0)
    new_coeff0[index_nonNA0] = new_coeff0_nonNA;
    new_vcov0[index_nonNA0, index_nonNA0] = new_vcov0_nonNA
    coeff0 = new_coeff0; vcov0 = new_vcov0

    # Remove the I(d) line/columns:
    coeff0 = coeff0[-line0_const]
    vcov0 = vcov0[-line0_const, -line0_const]
    names_var0 = names_var0[-line0_const]

    std_errors0 = sqrt(diag(vcov0))
    t_values0 <- coeff0 / std_errors0
    p_values0 <- 2 * pt(-abs(t_values0), df = df)

    table0 = data.frame(Variable = names_var0, Estimate = coeff0, Std_Error = std_errors0, t_value = t_values0, p_value = p_values0)





    # 2. MTR1 table

    line_coeff1 = grep(paste0("I(", var_treatment, "):"), names_var, fixed=TRUE) # don't include the kd in the homogenous case
    line_coefftreated = which(names_var == paste0("I(", var_treatment, ")")) # Also the pure I(d)
    line_kd1 = linekd
    line1 = c(line_intercept, line_coefftreated, line_kd1, line_coeff1)

    names_var1 = names_var[line1]
    names_var1 = gsub(paste0("I(", var_treatment, "):"), "", names_var1, fixed=TRUE)
    coeff1 = coeff[line1]; vcov1 = vcov[line1, line1]

    # "Problem": need an additional transformation because 3 effects in the constant here
    # (Intercept) + I(d) + kd(v): v^const.
    ncoeff1 = length(coeff1)
    A1 = matrix(0, nrow=ncoeff1, ncol=ncoeff1)
    diag(A1) = 1
    line1_intercept = 1
    line1_treated1 = which(names_var1 == paste0("I(", var_treatment, ")"))
    line1_constantkd = which(names_var1 == "kd(v): v^const")
    line1_treated = c(line1_treated1, line1_constantkd)
    A1[line1_intercept, line1_treated] = 1;

    #A1 = A1[-line1_treated,]
    #coeff1 = A1%*%coeff1
    #vcov1 = A1 %*% vcov1 %*% t(A1)

    # To apply the transformation, pay attention in case some missing coefficients
    index_nonNA1 = which(!is.na(coeff1));
    A1_nonNA = A1[index_nonNA1, index_nonNA1]
    coeff1_nonNA = coeff1[index_nonNA1]
    vcov1_nonNA = vcov1[index_nonNA1, index_nonNA1]
    names_var1_nonNA = names_var1[index_nonNA1]

    new_coeff1_nonNA = A1_nonNA %*% coeff1_nonNA
    new_vcov1_nonNA = A1_nonNA %*% vcov1_nonNA %*% t(A1_nonNA)
    new_coeff1 = rep(NA, ncoeff1); new_vcov1 = matrix(NA, nrow=ncoeff1, ncol=ncoeff1)
    new_coeff1[index_nonNA1] = new_coeff1_nonNA;
    new_vcov1[index_nonNA1, index_nonNA1] = new_vcov1_nonNA
    coeff1 = new_coeff1; vcov1 = new_vcov1

    # Remove the I(d) line/columns:
    coeff1 = coeff1[-line1_treated]
    vcov1 = vcov1[-line1_treated, -line1_treated]
    names_var1 = names_var1[-line1_treated]

    std_errors1 = sqrt(diag(vcov1))
    t_values1 <- coeff1 / std_errors1
    p_values1 <- 2 * pt(-abs(t_values1), df = df)

    table1 = data.frame(Variable = names_var1, Estimate = coeff1, Std_Error = std_errors1, t_value = t_values1, p_value = p_values1)






    # 3. MTE Table
    # Build the MTE as a transformation of the original coefficients
    # For each variable, need to determine if it is w0, w1, x or kd(u) related.
    # In the homogenous case: only a difference between w0, w1, x and a constant.
    # Since k0(v) = k1(v) (up to the constant)

    ncoeff = length(coeff)

    line_coeff0 = grep(paste0("I(1 - ", var_treatment, "):"), names_var, fixed=TRUE)
    names0 = names_var[line_coeff0]
    names0 = gsub(paste0("I(1 - ", var_treatment, "):"), "", names0, fixed=TRUE)

    line_coeff1 = grep(paste0("I(", var_treatment, "):"), names_var, fixed=TRUE)
    names1 = names_var[line_coeff1]
    names1 = gsub(paste0("I(", var_treatment, "):"), "", names1, fixed=TRUE)

    line_w1 = line_coeff1[which(!names1 %in% names0)];
    names_w1 = names1[which(!names1 %in% names0)]
    names_w1 = paste0("+ W1: ", names_w1)
    line_w0 = line_coeff0[which(!names0 %in% names1)];
    names_w0 = names0[which(!names0 %in% names1)]
    names_w0 = paste0("- W0: ", names_w0) # to make it clear that it is multiplied by -1

    line_x1 = line_coeff1[which(names1 %in% names0)]
    line_x0 = line_coeff0[which(names0 %in% names1)]
    names_x_kd = names1[which(names1 %in% names0)]

    # First line: simply the difference of intercept, directly given by I(d) coeff
    line_treated = which(names_var == paste0("I(", var_treatment, ")"))
    mattreated = matrix(0, nrow=1, ncol=ncoeff)
    mattreated[, line_treated] = 1
    name_treated = "(Intercept)"

    # Lines for W1 and W0:
    matw1 = matrix(0, nrow=length(line_w1), ncol=ncoeff)
    for(j in 1:length(line_w1)) {matw1[j, line_w1[j]] = 1}
    matw0 = matrix(0, nrow=length(line_w0), ncol=ncoeff)
    for(j in 1:length(line_w0)) {matw0[j, line_w0[j]] = -1} # - 1 for w0.

    # Lines for covariates (and kdu): always find the matched line
    MATX = list();
    if(length(names_x_kd) > 0) {
      for(j in 1:length(names_x_kd)) {
        name_i = names_x_kd[j]
        line_x0_i = line_coeff0[which(names0 == name_i)]
        line_x1_i = line_coeff1[which(names1 == name_i)]

        mat = matrix(0, nrow=1, ncol=ncoeff)
        mat[, line_x1_i] = 1;
        mat[, line_x0_i] = -1;

        MATX[[j]] = mat
      }
    }
    matx = do.call('rbind', MATX)

    A_mte = rbind(mattreated, matw1, matw0, matx) # transformation matrix to apply to coeff
    names_var_mte = c(name_treated, names_w1, names_w0, names_x_kd)
    ncoeff_mte = length(names_var_mte)
    #coeff_mte = A_mte%*%coeff
    #vcov_mte = A_mte %*% vcov %*% t(A_mte)


    # Deal with NA coefficients:
    index_nonNA_mte = which(!is.na(coeff))
    index_NA_mte = which(is.na(coeff))
    index_row_ok = 1:nrow(A_mte); index_row_NA = list()
    if(length(index_NA_mte) > 0) {
      for(j in 1:length(index_NA_mte)) {
        index_row_NA[[j]] = which(A_mte[, index_NA_mte[j]] != 0)
      }
    }
    index_row_NA = unique(do.call('c', index_row_NA))
    index_row_ok = index_row_ok[which(!index_row_ok %in% index_row_NA)]

    A_mte2 = A_mte[index_row_ok, index_nonNA_mte]
    coeff_nonNA = coeff[index_nonNA_mte]
    vcov_nonNA = vcov[index_nonNA_mte, index_nonNA_mte]

    coeff_mte_nonNA = A_mte2%*%coeff_nonNA
    vcov_mte_nonNA = A_mte2 %*% vcov_nonNA %*% t(A_mte2)

    # Rebuild the full coefficients and vcovs
    new_coeff_mte = rep(NA, ncoeff_mte); new_vcov_mte = matrix(NA, nrow=ncoeff_mte, ncol=ncoeff_mte)
    new_coeff_mte[index_row_ok] = coeff_mte_nonNA;
    new_vcov_mte[index_row_ok, index_row_ok] = vcov_mte_nonNA
    coeff_mte = new_coeff_mte; vcov_mte = new_vcov_mte


    std_errors_mte = sqrt(diag(vcov_mte))
    t_values_mte <- coeff_mte / std_errors_mte
    p_values_mte <- 2 * pt(-abs(t_values_mte), df = df)

    table_mte = data.frame(Variable = names_var_mte, Estimate = coeff_mte, Std_Error = std_errors_mte, t_value = t_values_mte, p_value = p_values_mte)
  }



  # 4. Form results
  coeff0 = as.numeric(coeff0); names(coeff0) = names_var0
  coeff1 = as.numeric(coeff1); names(coeff1) = names_var1
  coeff_mte = as.numeric(coeff_mte); names(coeff_mte) = names_var_mte

  colnames(vcov0) = names_var0; rownames(vcov0) = names_var0
  colnames(vcov1) = names_var1; rownames(vcov1) = names_var1
  colnames(vcov_mte) = names_var_mte; rownames(vcov_mte) = names_var_mte


  tables = list(table0, table1, table_mte); names(tables) = c("mtr0", "mtr1", "mte")
  coeffs = list(coeff0, coeff1, coeff_mte); names(coeffs) = c("mtr0", "mtr1", "mte")
  vcovs = list(vcov0, vcov1, vcov_mte); names(vcovs) = c("mtr0", "mtr1", "mte")
  #names_VAR = list(names_var0, names_var1, names_var_mte); names(names_VAR) = c("mtr0", "mtr1", "mte")

  output = list(tables, coeffs, vcovs)
  names(output) = list("table", "coeff", "vcovs")

  return(output)


}








#' @rdname mtr_fun
#' @param coeff Vector of (stacked) coefficients for mtr0, mtr1.
#' @param vcov Covariance matrix of these coefficients.
#' @param ref_indiv Newdata (reference individuals) at which to compute the predictions.
#' @param seq_u Sequence of v at which to compute the prediction for kd(v).
#' @param names_var Names of the variables corresponding to the coefficients.
#' @param df Degrees of freedom for the p-values.
#' @param est_method Either "sieve" or "homogenous".
#' @usage mtr_predict_sieve(coeff, vcov, ref_indiv, var_treatment, var_cov_2nd, pol_degree, seq_u, t_value, est_method, se_type)
#' @return Returns mtr0, mtr1, and mte estimates with confidence intervals for the specified ref_indiv.
#' @export
mtr_predict_sieve = function(coeff, vcov, ref_indiv, var_treatment, var_cov_2nd, pol_degree, seq_u,
                       t_value,
                       est_method="sieve", se_type="HC1") {

  # 0. Format fix
  # Possible that some coefficients (e.g., some fixed effects) are NA -  but not a problem if does not concern the reference individual
  coeff_base = coeff; vcov_base = vcov
  index_NA = which(is.na(coeff_base))
  if(length(index_NA > 0)) {
    coeff = coeff_base[-index_NA]
    vcov = vcov_base[-index_NA, -index_NA]
  }


  # 1. Predict effect of covariates+semi-IVs and kd(v) separately
  # -------------------------------------------------------------
  # 1.(i) Formula to compute the appropriate matrices to compute the coefficients
  if(est_method == "homogenous") {
    formula_control_function_mtr = paste0("I(1-", var_treatment, "+", var_treatment, "):Kappa_homogenous_fun(Phat, pol_degree)")
  }
  if(est_method == "sieve") {
    formula_control_function_mtr = paste0("I(1-", var_treatment, "):Kappa_fun(Phat, pol_degree) + I(", var_treatment, "):Kappa_fun(Phat, pol_degree)") # general formula
  }
  formula_MTR = as.formula(paste0("~ ", paste0(var_cov_2nd, collapse="+"), "+", formula_control_function_mtr))


  # 1.(ii) Predict effects of the covariates: X betad + Wd deltad
  # D=0: delta0(W0, X) = Xbeta0 + W0 delta0:
  REF_cov = ref_indiv; REF_cov$Phat = 0;
  REF0_cov = REF_cov; REF0_cov[[var_treatment]] = 0
  mat0_cov = model.matrix(formula_MTR, data=REF0_cov);
  mat0_cov[, "(Intercept)"] = 0 # removes the intercept, included in kd(u)
  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
    missing_col0 = mat0_cov[, index_NA]
    if(length(which(missing_col0 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat0_cov = mat0_cov[, -index_NA]
  }
  delta0X = mat0_cov%*%coeff
  se_delta0X = sqrt(rowSums((mat0_cov %*% vcov)*mat0_cov))



  # D=1: delta1(W1, X) = Xbeta1 + W1 delta1:
  REF_cov = ref_indiv; REF_cov$Phat = 0;
  REF1_cov = REF_cov; REF1_cov[[var_treatment]] = 1
  mat1_cov = model.matrix(formula_MTR, data=REF1_cov);
  mat1_cov[, "(Intercept)"] = 0 # removes the intercept, included in kd(u)
  line_treated = which(colnames(mat1_cov) == paste0("I(", var_treatment, ")"))
  mat1_cov[, line_treated] = 0 # also part of the intercept, put it at 0.

  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
    missing_col1 = mat1_cov[, index_NA]
    if(length(which(missing_col1 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat1_cov = mat1_cov[, -index_NA]
  }
  delta1X = mat1_cov%*%coeff;
  se_delta1X = sqrt(rowSums((mat1_cov %*% vcov)*mat1_cov))

  deltaX = data.frame(id=ref_indiv$id, delta0X = delta0X, delta1X = delta1X,
                      se_delta0X=se_delta0X, se_delta1X = se_delta1X)



  # 1.(iii) Predict kd(v) (including the constant)
  REF = ref_indiv[1,]; REF[1,] = 0 # set every values to 0.
  REF = REF[rep(1, times=length(seq_u)),]
  REF$Phat = seq_u;

  # k0(u):
  REF0 = REF; REF0[[var_treatment]] = 0
  mat0 = model.matrix(formula_MTR, data=REF0);
  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
    missing_col0 = mat0[, index_NA]
    if(length(which(missing_col0 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat0 = mat0[, -index_NA]
  }
  k0u = mat0%*%coeff; # numerical values output
  se_k0u = sqrt(rowSums((mat0 %*% vcov)*mat0))

  # k1(u):
  REF1 = REF; REF1[[var_treatment]] = 1
  mat1 = model.matrix(formula_MTR, data=REF1);
  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
    missing_col1 = mat1[, index_NA]
    if(length(which(missing_col1 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat1 = mat1[, -index_NA]
  }
  k1u = mat1%*%coeff; # numerical values output
  se_k1u = sqrt(rowSums((mat1 %*% vcov)*mat1))

  # return a table:
  kv = data.frame(v = seq_u, k0 = k0u, k1 = k1u, se_k0 = se_k0u, se_k1 = se_k1u)





  # 2. Predict MTR0, MTR1, MTE
  # Include proper confidence intervals (so recompute directly the entire thing to get them directly)

  REF = ref_indiv[rep(seq_len(nrow(ref_indiv)), each = length(seq_u)), ]
  REF$Phat = rep(seq_u, times=nrow(ref_indiv))
  REF1 = REF; REF1[[var_treatment]] = 1
  REF0 = REF; REF0[[var_treatment]] = 0
  RES = REF;

  # (i) MTR computation
  # MTR1:
  mat1 = model.matrix(formula_MTR, data=REF1);
  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
    missing_col1 = mat1[, index_NA]
    if(length(which(missing_col1 != 0)) > 0) { stop("Reference individual has a characteristic for which the coefficient could not be estimated")}
    mat1 = mat1[, -index_NA]
  }
  mtr1 = mat1%*%coeff;

  # MTR0:
  mat0 = model.matrix(formula_MTR, data=REF0);
  if(length(index_NA) > 0) { # check that ref individual does not have a characteristics in the missing column:
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



  # 3. Save
  # (i) Main result data
  RES$mtr1 = mtr1;
  RES$mtr0 = mtr0;
  RES$mte = mte
  RES$mtr0_lwr = mtr0_lwr; RES$mtr0_upr = mtr0_upr; RES$mtr0_se = se0
  RES$mtr1_lwr = mtr1_lwr; RES$mtr1_upr = mtr1_upr; RES$mtr1_se = se1;
  RES$mte_lwr = mte_lwr; RES$mte_upr = mte_upr; RES$mte_se = se_mte

  # Also add the "split" between k1u/k0u and "covariates" effects
  RES = merge(RES, deltaX, by="id")
  RES = merge(RES, kv, by.x="Phat", by.y="v")

  # re-order (the last merger messes up with the ordering)
  RES = RES[order(RES$Phat),]
  RES = RES[order(RES$id),]



  # (ii) Separately deltaX and KV
  output = list(RES, deltaX, kv)
  names(output) = c("est", "deltaX", "kv")

  return(output)
}





NULL
