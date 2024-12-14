#' @title MTE and MTR sub-estimation functions
#'
#' @description These functions allow to estimate the mte and mtr, and their confidence intervals, based on coefficients estimated from the main model in the main function.
#' More details can be found in Bruneel-Zupanc (2024). \cr
#' Different formulas must be applied depending on whether the treatment is homogenous or heterogenous.
#'
#'
#' @rdname mtr_fun
#' @param `fast_robinson1` Default is TRUE to speed things up in a first stage (if many covariates in particular). If TRUE, will use the locpoly function from Kernsmooth library to speed up the computation of the Robinson double residual first stage. This is only possible if no external weights are used. Fast Locpoly will enforce a gaussian kernel.
#' @param `fast_robinson2` Default is FALSE. If TRUE, will use the locpoly function from Kernsmooth library to speed up the computation of the Robinson double residual second stage. This is only possible if no external weights are used. Fast Locpoly will enforce a gaussian kernel.
#' Default is FALSE for the second stage because fast_locpoly returns no standard errors and the gain in time is not so important for the second stage.
#' @usage mtr_est_poly(data, seq_u,
#'                     bw0 = NULL, bw1=NULL, bw_y0 = NULL, bw_y1=NULL, bw_method = 1/5, kernel,
#'                     bw_subsamp_size = NULL, fast_locpoly = FALSE,
#'                     fast_robinson1 = TRUE, fast_robinson2 = FALSE,
#'                     pol_degree1, pol_degree2, var_outcome, var_treatment, var_w0, var_w1, var_covariates, print_progress)
#' @export
mtr_est_poly = function(data, seq_u,
                        bw0 = NULL, bw1=NULL, bw_y0 = NULL, bw_y1=NULL, bw_method = 1/5, kernel,
                        bw_subsamp_size = NULL,
                        fast_robinson1 = TRUE, fast_robinson2 = FALSE,
                        se_type,
                        pol_degree1, pol_degree2,
                        var_outcome, var_treatment, var_w0, var_w1, var_covariates,
                        print_progress=FALSE) {

  # -------------------------------------------
  # Double Residual Regression (Robinson, 1988)
  # -------------------------------------------

  supp = range(data$Phat)

  # 1st Stage. Local polynomial regressions of (Y, Wd, X) on P.
  # ----------
  for(d in c(0, 1)) {

    # All regressions done on subsample D=d
    datad = data[which(data[[var_treatment]] == d),]
    if(d==0) { var_w = var_w0; bwd = bw0; } else { var_w = var_w1; bwd = bw1; }
    var_covd = c(var_w, var_covariates); var_covd = unique(var_covd[which(var_covd != "")])
    formulad = as.formula(paste0("~ -1 + ", paste0(var_covd, collapse="+")))

    Xd = model.matrix(formulad, data=datad)
    datd = cbind(subset(datad, select=var_outcome), Xd)
    Phatd = datad$Phat; yd = datd[[var_outcome]]
    weightsd = datad$WEIGHTS
    residd = datd;

    BWD=c()

    # if specify only one value in bwd, apply it to ALL the variables (remark: if only one covariate, does not change anything)
    if(length(bwd) == 1) { bwd = rep(bwd, ncol(datd)) }

    if(print_progress) { cat(sprintf("D= %d, Robinson 1st residual regression of X on P...", d), "              \n") }
    for(j in 1:ncol(datd)) {

      # (i) Bandwidth Selection.
      bwj = bwd[j] # if is NULL will always put null
      if(is.null(bwj)) {
        bwj = lbw_select(x=Phatd, y=datd[,j], kernel=kernel, degree=pol_degree1, drv=0, bw_method = bw_method, bw_subsamp_size=bw_subsamp_size, supp=supp)$bw # drv = 0 because want to estimate the main function;
      }

      # (ii) Local Polynomial Regression
      #if(length(unique(datd[,j])) == 1) { warning(paste0("No variation in variable ", colnames(datd)[j], " in the subsample d = ", d, " of the trimmed subsample.")) }

      resj = lpoly(x=Phatd, y=datd[,j], bandwidth=bwj, degree = pol_degree1, drv = 0,
                   kernel = kernel, se_type=se_type, weights=weightsd, fast_locpoly=fast_robinson1)
      # recommands degree = derivative order + 1; cf Fan and Gijbels (1996), footnote 19 of Carneiro et al 2011

      pred = approx(x=resj$x, y=resj$y, xout=Phatd)$y
      resid = datd[,j] - pred
      BWD[j] = bwj; residd[,j] = resid

      if(print_progress & bw_method %in% c("mse-dpi", "mse_rot") & is.null(bwd)) { cat(sprintf("Progress: %d/%d", j, ncol(datd)), "                             \r") }
    }

    formuladx = as.formula(paste0(var_outcome, "~ -1 + ."))
    estd = lm(formuladx, data=residd, weights=weightsd)

    # Save:
    if(d==0) { est0 = estd; X0 = Xd; BW0 = BWD; dat0 = datd; Phat0 = Phatd; weights0 = weightsd; }
    if(d==1) { est1 = estd; X1 = Xd; BW1 = BWD; dat1 = datd; Phat1 = Phatd; weights1 = weightsd; }
  }


  # 2nd Stage. Estimate kd(v).
  # ----------

  # 1. Bandwidth Selection
  # Can either do two separate bandwidth selection for k0(v) and k1(v)
  # Or can do a bandwidth selection for the MTE(v) directly, i.e., k1(v) - k0(v).
  # Prefer the second option as it also has the advantage of applying the same bandwidth to both k0 and k1, that is optimal for the MTE.

  if(print_progress) { cat(sprintf("Robinson 2nd stage: Bandwidth Selection...                           \r")) }

  # To do the MTE approach, need to predict for every observation:
  yfit0 = predict(est0, newdata=data); yfit1 = predict(est1, newdata=data)
  weightsall = data$WEIGHTS
  xall = data$Phat
  yfit = data[[var_outcome]] - yfit0*(1-xall) - yfit1*xall


  bw_mte = bw_y0 # arbitrarily set it to bw_y0 -> gets updated if any missing bw
  if(is.null(bw_y0) | is.null(bw_y1)) { # Only if not pre-specified
    bw_y = lbw_select(x=xall, y=yfit, kernel=kernel, degree=pol_degree2, drv=1, # drv=1 because bw for the derivative of the fit.
                      bw_method = bw_method, bw_subsamp_size=bw_subsamp_size, supp=supp)$bw
    bw_mte = bw_y;
  }
  if(is.null(bw_y0)) { bw_y0 = bw_y }
  if(is.null(bw_y1)) { bw_y1 = bw_y }




  # 2. Robinson 2nd stage Local Polynomial Regression
  if(print_progress) { cat(sprintf("Robinson 2nd stage: Estimation of k0(v) and k1(v)... \n")) }

  # (i) Extract Yd net of the effect of the covariates
  # D=0:
  ypred0 = predict(est0, newdata=as.data.frame(X0))
  ynet0 = dat0[[var_outcome]] - ypred0
  x0 = Phat0;
  y0 = -ynet0*(1-x0) # this is the object we derive in order to obtain k0(v).

  # D=1:
  ypred1 = predict(est1, newdata=as.data.frame(X1))
  ynet1 = dat1[[var_outcome]] - ypred1
  x1 = Phat1;
  y1 = ynet1*x1 # this is the object we derive in order to obtain k1(v).


  # (ii) Local Polynomial regressions
  # D=0:
  k0_est = lpoly(x=x0, y=y0, bandwidth=bw_y0, degree = pol_degree2, drv = 1,
               kernel = kernel, se_type=se_type, weights=weights0, fast_locpoly=fast_robinson2)
  k0 = approx(x=k0_est$x, y=k0_est$y, xout=seq_u)$y
  if(length(which(!is.na(k0_est$se))) == 0) { # will be NA only if fast_robinson2==TRUE
    k0_se = rep(NA, length(seq_u))
  } else {
    k0_se = approx(x=k0_est$x, y=k0_est$se, xout=seq_u)$y
  }



  # D=1:
  k1_est = lpoly(x=x1, y=y1, bandwidth=bw_y1, degree = pol_degree2, drv = 1,
                 kernel = kernel, se_type=se_type, weights=weights1, fast_locpoly=fast_robinson2)
  k1 = approx(x=k1_est$x, y=k1_est$y, xout=seq_u)$y
  if(length(which(!is.na(k1_est$se))) == 0) { # will be NA only if fast_robinson2==TRUE
    k1_se = rep(NA, length(seq_u))
  } else {
    k1_se = approx(x=k1_est$x, y=k1_est$se, xout=seq_u)$y
  }


  # MTE:
  # MTE is = k1(v) - k0(v)
  # But if also wants standard errors around it, need to run a third lpoly regression - not exactly numerically equivalent though...
  deltak_est = lpoly(x=xall, y=yfit, bandwidth=bw_mte, degree = pol_degree2, drv = 1,
                 kernel = kernel, se_type=se_type, weights=weightsall, fast_locpoly=fast_robinson2)
  deltak = approx(x=deltak_est$x, y=deltak_est$y, xout=seq_u)$y

  if(length(which(!is.na(deltak_est$se))) == 0) {
    deltak_se = rep(NA, length(seq_u))
  } else {
    deltak_se = approx(x=deltak_est$x, y=deltak_est$se, xout=seq_u)$y
  }



  # (iii) Save:
  kv=data.frame(v=seq_u); # to return the function kd(v)
  kv$k0 = k0; kv$k1 = k1
  kv$k0_se = k0_se; kv$k1_se = k1_se
  kv$deltak = deltak; kv$deltak_se = deltak_se

  res_list = list(kv=kv, est0=est0, est1=est1,
                  bw0=BW0, bw1=BW1, bw_y0=bw_y0, bw_y1=bw_y1, bw_mte = bw_mte)
  return(res_list)

  #### Alternative computation (but longer + requires two locpolynomial per D + unknown optimal bw for these).
  ### (ii) Estimate Kappa_d(u) with local polynomial regression
  ##est_Kappad = lpoly_regress(x=Phatd, y=ynetd, bw=bw_y, bw_method=bw_method, degree=pol_degree2, drv=0, weights=weightsd)
  ### pol_degree2 because in the end I want to get the derivatives too.
  ### This regression is mostly ran to extract the bw in fact!
  ##bw_y = est_Kappad$bw # if pre-specified will return the same value
  ##est_dKappad = lpoly_regress(x=Phatd, y=ynetd, bw = bw_y, degree=pol_degree2, drv=1, weights=weightsd) # same bw, no need to specify the method then
  ### should be the same bandwidth in both, so don't recompute it
  ### Estimation on the set of specified u:
  ##Kappad = approx(x=est_Kappad$reg$x, y=est_Kappad$reg$y, xout=seq_u)$y
  ##dKappad = approx(x=est_dKappad$reg$x, y=est_dKappad$reg$y, xout=seq_u)$y
  ### (iii) Estimate kd(u):
  ### See Andresen 2018, Table 2 for the formula; Can also derive manually;
  ##if(d==0) {
  ##  ku = Kappad - (1-seq_u)*dKappad
  ##}
  ##if(d==1) {
  ##  ku = Kappad + seq_u*dKappad
  ##}

}






#' @rdname mtr_fun
#' @usage mtr_fun_poly(ref_indiv, eval_v, est0, est1, kv, se_type, conf_level)
#' @export
mtr_fun_poly = function(ref_indiv, eval_v, est0, est1, kv, se_type, conf_level) {

  # Main prediction: predict mtr, mte

  # For each ref_indiv, predict the effect of the covariates - repeat it to have a data.frame as result
  REF = ref_indiv[rep(seq_len(nrow(ref_indiv)), each = length(eval_v)), ]
  REF$Phat = rep(eval_v, times=nrow(ref_indiv))

  # (i) Effect of covariates

  if(se_type == "nonrobust") {
    X0_pred = predict(est0, newdata=REF, se.fit=TRUE)
    X1_pred = predict(est1, newdata=REF, se.fit=TRUE)
    delta0x = X0_pred$fit
    delta1x = X1_pred$fit
    se_delta0x = X0_pred$se.fit
    se_delta1x = X1_pred$se.fit
    # For the standard errors: cannot use predict default if do not use nonrobust.

    # If some problem in the reference individuals, will have warnings of the form:
    #1: In predict.lm(est0, newdata = REF, se.fit = TRUE) :
    #  prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases

  } else {
    # If specific standard errors (e.g., se_type = "HC1"), cannot use predict: need to do it manually

    coeff0 = coefficients(est0); coeff1 = coefficients(est1)
    vcov0 = vcovHC(est0, type = se_type) # robust standard errors;
    vcov1 = vcovHC(est1, type = se_type)

    # vcov does NOT include rows for NA coefficients
    index_NA0 = which(is.na(coeff0)); index_NA1 = which(is.na(coeff1))
    index_nonNA0 = which(!is.na(coeff0)); index_nonNA1 = which(!is.na(coeff1))
    # if(length(index_NA0) > 0 | length(index_NA1) > 0) { warning("Some coefficients are NA. May be impossible to evaluate for some reference individual. ") }
    # But the NA only causes problem if they are multiplied by something else than 0 in the ref_indiv.
    #   In which case, the prediction is NA -> stop the function.
    #   Otherwise: does not matter! It's missing because the trimming probably removed entirely one factor for example.

    # D=0:
    formula0y = paste0(deparse(formula(est0)), collapse="")
    formula0 = as.formula(gsub(".*~", "~", formula0y)) # remove Y
    mat0 = model.matrix(formula0, data=REF)

    # D=1:
    formula1y = paste0(deparse(formula(est1)), collapse="")
    formula1 = as.formula(gsub(".*~", "~", formula1y)) # remove Y
    mat1 = model.matrix(formula1, data=REF)

    # Check if no problem with NA: -> can occur if, for example, no D=1 for a given factor (even though some D=0 are present).
    if(length(index_NA0) > 0) {
      mat0_NA = mat0[, index_NA0]
      if(length(which(mat0_NA != 0)) > 0) { warning("Predicton from rank-deficient fit; Missing coefficients to evaluate the effect at the reference individuals when D=0.") }
    }
    if(length(index_NA1) > 0) {
      mat1_NA = mat1[, index_NA1]
      if(length(which(mat1_NA != 0)) > 0) { warning("Predicton from rank-deficient fit; Missing coefficients to evaluate the effect at the reference individuals when D=1.") }
    }

    # If passes the NA test:
    coeff0b = coeff0; coeff1b = coeff1; mat0b = mat0; mat1b = mat1; # save base
    coeff0 = coeff0b[index_nonNA0]; coeff1 = coeff1b[index_nonNA1]
    mat0 = matrix(mat0b[, index_nonNA0], nrow=nrow(mat0b)); # the "matrix" part is to ensure that it's still a matrix even if only one variable left.
    mat1 = matrix(mat1b[, index_nonNA1], nrow=nrow(mat1b))

    delta0x = as.vector(mat0 %*% coeff0)
    se_delta0x = sqrt(diag(mat0 %*% vcov0 %*% t(mat0)))

    delta1x = as.vector(mat1 %*% coeff1)
    se_delta1x = sqrt(diag(mat1 %*% vcov1 %*% t(mat1)))

  }


  # (ii) k(v):
  # Also repeat kv:
  range_v = range(kv$v)
  if(length(which(eval_v > range_v[2] | eval_v < range_v[1])) > 0) {
    stop("Evaluation value for V is out of the estimated range.")
  }

  kv_input = kv;

  # Reconstruct kv but on the new evaluation sample:
  k0 = approx(x=kv_input$v, y=kv_input$k0, xout=eval_v)$y
  k1 = approx(x=kv_input$v, y=kv_input$k1, xout=eval_v)$y
  deltak = approx(x=kv_input$v, y=kv_input$deltak, xout=eval_v)$y

  # S.E. (only if available, i.e., if fast_robinson2 == FALSE)
  if(length(which(!is.na(kv_input$k0_se))) == 0) {
    k0_se = rep(NA, length(eval_v));
    k1_se = rep(NA, length(eval_v));
    deltak_se = rep(NA, length(eval_v));
  } else {
    k0_se = approx(x=kv_input$v, y=kv_input$k0_se, xout=eval_v)$y
    k1_se = approx(x=kv_input$v, y=kv_input$k1_se, xout=eval_v)$y
    deltak_se = approx(x=kv_input$v, y=kv_input$deltak_se, xout=eval_v)$y
  }

  kv = data.frame(v=eval_v, k0, k1, k0_se, k1_se, deltak, deltak_se)


  # repeat it to fit the REF data;
  K0 = rep(k0, times=nrow(ref_indiv))
  K1 = rep(k1, times=nrow(ref_indiv))
  K0_SE = rep(k0_se, times=nrow(ref_indiv))
  K1_SE = rep(k1_se, times=nrow(ref_indiv))
  DELTAK = rep(deltak, times=nrow(ref_indiv))
  DELTAK_SE = rep(deltak_se, times=nrow(ref_indiv))

  # (iii) MTRd and MTE
  mtr0 = delta0x + K0
  mtr1 = delta1x + K1
  mte = mtr1 - mtr0 # measure of mte for direct comparison, not estimated "directly".

  # MTE from direct estimation: - to have standard errors around it...
  mte2 = delta1x - delta0x + DELTAK # Second measure of mte, not directly the same, but that's the one we use for s.e.


  # (iv) Confidence intervals - Not taking into account errors around delta_d X
  df0 = df.residual(est0); t_value0 = qt(1-conf_level/2, df=df0)
  df1 = df.residual(est1); t_value1 = qt(1-conf_level/2, df=df1)
  t_value_mte = qt(1-conf_level/2, df=df0+df1)

  mtr0_lwr = mtr0 - t_value0*K0_SE
  mtr0_upr = mtr0 + t_value0*K0_SE
  mtr1_lwr = mtr1 - t_value1*K1_SE
  mtr1_upr = mtr1 + t_value1*K1_SE
  mte2_lwr = mte2 - t_value_mte*DELTAK_SE
  mte2_upr = mte2 + t_value_mte*DELTAK_SE


  # Save
  RES = REF;
  RES$mtr0 = mtr0; RES$mtr1 = mtr1
  RES$mte = mte; RES$mte2 = mte2;
  RES$mtr0_lwr = mtr0_lwr; RES$mtr0_upr = mtr0_upr;
  RES$mtr1_lwr = mtr1_lwr; RES$mtr1_upr = mtr1_upr;
  RES$mte2_lwr = mte2_lwr; RES$mte2_upr = mte2_upr
  RES$k0 = K0; RES$k1 = K1; RES$deltak = DELTAK;
  RES$delta0X = delta0x; RES$delta1X = delta1x
  RES$se_k0 = K0_SE; RES$se_k1 = K1_SE; RES$se_deltak = DELTAK_SE
  RES$se_delta0X = se_delta0x; RES$se_delta1X = se_delta1x


  # (v) Also return the effect of the covariates (but only one prediction by observation)
  # remark: redundant to recompute, but quick..
  deltaX = ref_indiv;
  if(se_type == "nonrobust") {
    delta0X_pred = predict(est0, newdata=ref_indiv, se.fit=TRUE)
    delta1X_pred = predict(est1, newdata=ref_indiv, se.fit=TRUE)
    deltaX$delta0X = delta0X_pred$fit; deltaX$delta1X = delta1X_pred$fit;
    deltaX$se_delta0X = delta0X_pred$se.fit; deltaX$se_delta1X = delta1X_pred$se.fit;

  } else {
    mat0rb = model.matrix(formula0, data=deltaX)
    mat0r = matrix(mat0rb[, index_nonNA0], nrow=nrow(mat0rb))
    delta0X = mat0r %*% coeff0
    se_delta0X = sqrt(diag(mat0r %*% vcov0 %*% t(mat0r)))

    mat1rb = model.matrix(formula1, data=deltaX)
    mat1r = matrix(mat1rb[, index_nonNA1], nrow=nrow(mat1rb))
    delta1X = as.vector(mat1r %*% coeff1)
    se_delta1X = sqrt(diag(mat1r %*% vcov1 %*% t(mat1r)))

    deltaX$delta0X = delta0X; deltaX$delta1X = delta1X
    deltaX$se_delta0X = se_delta0X; deltaX$se_delta1X = se_delta1X
  }


  res_list = list(RES, deltaX, kv); names(res_list) = c("est", "deltaX", "kv")
  return(res_list)
}



#' @rdname mtr_fun
#' @description
#' Bandwidth selection for local polynomial regression. Different methods are available: "mse-dpi", "mse-rot" (from nprobust) or "arbitrary" (fixed bandwidth to a fraction of the support).
#' Can provide bandwidth for the main function or for any of its derivative order.
#' The bandwidth can be computed on a subsample of the data of size `bw_subsamp_size` to speed up the computation.
#' @param x Vector of x values
#' @param y Vector of y values
#' @param bw Pre-specified bandwidth
#' @param bw_method Method to compute the bandwidth (if bandwdith is NULL)
#' @param degree Degree of the polynomial: recommended to set to drv + 1
#' @param drv Derivative order of the function to be estimated.
#' @param supp Support of X in the complete data (not necessarily of the realized X in the current subsample), in order to compute the bandwidth if the rule is to take a fraction of the support. Ensure that same bandwidth on both samples D=0 and D=1. By default supp=NULL, in which case recompute the support of `x` directly.
#' @param bw_subsamp_size Size of the subsample to compute the bandwidth. Default is NULL (no subsample). If larger than sample size, it is ignored.
#' @usage lbw_select(x, y, kernel, degree, drv, bw_method = "arbitrary", bw_subsamp_size)
#' @export
lbw_select = function(x, y, kernel, degree, drv, bw_method = 1/5, bw_subsamp_size, supp=NULL) {

  # Rename kernel in nprobust formulation:
  kernel1 = ifelse(kernel == "gaussian", "gau",
                   ifelse(kernel == "epanechnikov", "epa", NA))
  if(is.na(kernel1)) { stop("Kernel not recognized.") }

  # 1. Arbitrary rule: fraction of the support
  if(is.numeric(bw_method)) {
    if(is.null(supp)) { supp = range(x) }
    bw = round(bw_method * (supp[2] - supp[1]), 3) # take a bandwidth equal to a fraction of the support
    if(bw_method < 0 | bw_method > 1) { stop("Bandwidth fraction of the support (bw_method) must be between 0 and 1.") }
  }


  # 2. Direct plug-in MSE:
  if(bw_method == "mse-dpi") {
    if(is.null(bw_subsamp_size)) {
      bw_est = lpbwselect(y=y, x=x, eval = NULL, neval = 30, p = degree, deriv = drv, kernel = kernel1, bwselect = "imse-dpi")
    } else {
      sub_id = sample(1:length(x), min(bw_subsamp_size, length(x))) # if subsampsize > length(x), just take x;
      suby = y[sub_id]; subx = x[sub_id]; #subweights = weights[sub_id]
      bw_est = lpbwselect(y=suby, x=subx, eval = NULL, neval = 30, p = degree, deriv = drv, kernel = kernel1, bwselect = "imse-dpi")
    }
    bw = bw_est$bws[,2]
  }

  # 3. Rule-of-thumb MSE:
  if(bw_method == "mse-rot") {
    if(is.null(bw_subsamp_size)) {
      bw_est = lpbwselect(y=y, x=x, eval = NULL, neval = 30, p = degree, deriv = drv, kernel = kernel1, bwselect = "imse-rot")
    } else {
      sub_id = sample(1:length(x), min(bw_subsamp_size, length(x))) # if subsampsize > length(x), just take x;
      suby = y[sub_id]; subx = x[sub_id]; #subweights = weights[sub_id]
      bw_est = lpbwselect(y=suby, x=subx, eval = NULL, neval = 30, p = degree, deriv = drv, kernel = kernel1, bwselect = "imse-rot")
    }
    bw = bw_est$bws[,2]
  }

  return(list(bw=as.numeric(bw)))
}


## ## OBSOLETE, DON T USE THE locpol PACKAGE ANYMORE -- too many errors + only computes BW for the main function, not derivative -> not what we want.
## if(bw_method == "rule-of-thumb") {
##   bw = thumbBw(x, y, deg=degree, kernel=gaussK, weig = weights) # might yield weird results when degree = 2
##   #if(is.nan(bw)) { stop("Error no variation in variable in this subsample.") }
##   if(bw > max(x) - min(x)) {
##     bw = thumbBw(x, y, deg=degree-1, kernel=gaussK, weig = weights) # reduce the degree
##     if(bw > max(x) - min(x)) { # if still greater:
##       warning("Rule-of-Thumb bw too large... Setting to (max(x)-min(x))/10 arbitrarily")
##       bw = (max(x) - min(x))/10
##     }
##   }
## }
## if(bw_method == "plug-in") {
##   # Gets long if too many observations... By default compute on subsample
##   sub_id = sample(1:length(x), min(4999, length(x))) # 5000 corresponds to .maxEvalPts in locpol
##   suby = y[sub_id]; subx = x[sub_id]; subweights = weights[sub_id]
##   # requires degree to be odd number (not even)
##   if(degree %% 2 == 0) {
##     warning("Degree should be an odd number for plug-in bandwidth selection: reduced by 1 for bw computation.")
##     degree = degree - 1
##   }
##   sub_bw <- pluginBw(subx, suby, deg=degree, kernel=gaussK, weig = subweights)
##   # Rescale bandwidth: Silverman (1986) rule: bw_full = bw_sub * (N_full/n_sub)^(-1/5)
##   # Quite arbitrary rule but saves a lot of time..
##   bw = sub_bw*(length(x)/length(subx))^(-1/5)
## }
## if(bw_method == "cv") {
##   # Gets long if too many observations... By default compute on subsample
##   sub_id = sample(1:length(x), min(4999, length(x))) # 5000 corresponds to .maxEvalPts in locpol
##   suby = y[sub_id]; subx = x[sub_id]; subweights = weights[sub_id]
##   range = range(subx); upper =  (range[2] - range[1])/2
##   lower = (range[2] - range[1])/1000
##   interval = c(lower, upper)
##   sub_bw = regCVBwSelC(x=subx, y=suby, deg=degree, kernel=gaussK,weig=subweights,
##                    interval = interval)
##   # bw very sensitive to interval choice;
##   # Rescale bandwidth: Silverman (1986) rule: bw_full = bw_sub * (N_full/n_sub)^(-1/5)
##   # Quite arbitrary rule but saves a lot of time..
##   bw = sub_bw*(length(x)/length(subx))^(-1/5)
## }
## ## # Remark: if want to allow for higher maxevalpts:
## ## unlockBinding(".maxEvalPts", asNamespace("locpol"))
## ## assign(".maxEvalPts", length(x) + 10, envir = asNamespace("locpol")) # Overwrite the value
## ## lockBinding(".maxEvalPts", asNamespace("locpol")) # Lock the binding again to avoid unintended modifications











#' @rdname mtr_fun
#' @description
#' The `lpoly` function estimates a (weighted) local polynomial regression of a specified degree at given evaluation points.
#' It supports derivative estimation and allows for heteroscedasticity-consistent standard errors. External weights (weights) are allowed.
#' @param x Vector of x values
#' @param y Vector of y values
#' @param bandwidth Pre-specified bandwidth
#' @param degree Degree of the polynomial: recommended to set to drv + 1
#' @param drv Derivative order of the function to be estimated.
#' @param se_type "HC1", "nonrobust" (for baseline homoscedastic), ...
#' @param weights Vector of external weights (in addition to the kernel weights)
#' @param x_eval Vector of evaluation points. If not pre-specified, use a grid (of gridsize). Default = NULL.
#' @param gridsize Size of the grid of evenly spaced points for x. Default is 201. Only relevant if x_eval = NULL.
#' @param fast_locpoly Default is FALSE. If `fast_locpoly` is TRUE, will use the locpoly function from Kernsmooth library to speed up the computation. This is only possible if no external weights are used. If the kernel is not set to Gaussian, the locpoly function will change it to Gaussian if fast_poly is TRUE.
#' @usage wlocpol(x, y, bandwidth, degree = 2, drv = 1, kernel = "gaussian", weights=NULL, x_eval=NULL, gridsize=201, fast_locpoly=FALSE)
#' @export
lpoly <- function(x, y, bandwidth, degree = 2, drv = 1,
                   kernel = "gaussian", se_type, weights=NULL, x_eval=NULL, gridsize=201, fast_locpoly=FALSE) {

  # Evaluation grid
  if(is.null(x_eval)) {
    x_eval <- seq(min(x), max(x), length.out = gridsize)
  }
  # Kernel function
  kernel_fn <- switch(kernel,
                      gaussian = function(u) exp(-0.5 * u^2) / sqrt(2 * pi),
                      epanechnikov = function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
                      uniform = function(u) ifelse(abs(u) <= 1, 0.5, 0),
                      stop("Unknown kernel"))
  # Weights:
  if(!is.null(weights)) { external_weights = weights } else { external_weights = rep(1, length(x)) }

  # 1. Fast version. Fast implementation in special case:
  # If fast_locpoly (& no weights & kernel = "gaussian"), use locpoly function from Kernsmooth library (speeds things up)
  weights_diff = external_weights - 1;

  if(max(abs(weights_diff)) < 0.0001 & fast_locpoly) {
    if(kernel != "gaussian") { warning("Fast locpolynomial regression specified -- Kernel Changed to Gaussian to estimate it.")} # no need to explicitly change it, no kernel specification in locpoly
    reg = locpoly(x=x, y=y, drv=drv, bandwidth=bandwidth, degree=degree)
    y1 = approx(x=reg$x, y=reg$y, xout=x_eval)$y
    return(list(x=x_eval, y=y1, se = rep(NA, length(y1))))
  }
  if(max(abs(weights_diff)) > 0.0001 & fast_locpoly) {
    warning("Fast locpolynomial regression impossible with external weights... Regular pace locpolynomial estimated instead.")
  }


  # 2. General Version. - takes longer, but applies in all cases and returns standard errors.
  # Initialize outputs
  n_eval <- length(x_eval)
  fit <- numeric(n_eval)  # Fitted values
  se <- numeric(n_eval)   # Standard errors

  # Loop over each evaluation point
  for (i in seq_along(x_eval)) {
    x0 <- x_eval[i]

    # Compute weights based on the kernel and bandwidth
    kernel_weights <- kernel_fn((x - x0) / bandwidth)
    weights <- kernel_weights * external_weights

    # Normalize weights
    weights <- weights / sum(weights)

    # Construct the design matrix
    X <- outer(x - x0, 0:degree, `^`) #make_design_matrix(x, x0, degree)

    # Weighted least squares
    XtW <- t(X * weights)         # Apply weights directly to rows of X
    XtWX <- XtW %*% X             # Weighted X'X
    XtWy <- XtW %*% y             # Weighted X'y
    beta <- solve(XtWX, XtWy)     # Coefficients

    # Variance estimation
    if(se_type == "nonrobust") {
      residuals <- y - X %*% beta
      sigma2 <- sum(weights * residuals^2) / (length(y) - (degree+1)) # length(y) = N # degree + 1 = number of coefficients to estimate with the intercept = degrees of freedom
      cov_matrix <- sigma2 * solve(XtWX)
    }
    if(se_type == "HC1") {
      residuals <- as.numeric(y - X %*% beta)  # Residuals
      XtWX_inv <- solve(XtWX)                  # (X'WX)^-1
      # Compute HC1 robust covariance matrix
      n <- nrow(X)                               # Number of observations
      k <- ncol(X)                               # Number of parameters
      diag_adjust <- (weights^2) * (residuals^2) # w_i^2 * r_i^2
      meat <- t(X) %*% (diag_adjust * X)         # Efficient computation of X' diag() X
      cov_matrix <- XtWX_inv %*% meat %*% XtWX_inv * (n / (n - k)) # HC1 correction
    }


    ## # Weighted Least Squares
    ## # via lm() directly - Quicker but takes a bit longer
    ## wls = lm(y~X-1, weights=weights)
    ## beta = coefficients(wls)
    ## if(se_type == "nonrobust") {
    ##   cov_matrix = vcov(wls)
    ## } else {
    ##   cov_matrix = vcovHC(wls, type = se_type)
    ## }


    # Extract coefficients of interest:
    if(drv == 0) { fit[i] = beta[1]; se[i] = cov_matrix[1,1] } # Raw function
    if(drv > 0) { # derivative:
      fit[i] = factorial(drv) * beta[drv + 1];
      se[i] = factorial(drv)*sqrt(cov_matrix[drv+1,drv+1])
    }

  }

  # Return results
  return(list(x=x_eval, y = fit, se = se))
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
