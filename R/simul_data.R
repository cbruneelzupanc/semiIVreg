#' Simulate data from the Generalized Roy Model with semi-IVs
#'
#'@description
#' This function simulates data from the Generalized Roy Model with semi-IVs, following the simulation specified in Bruneel-Zupanc (2024). \cr
#' For more details about the exact specification, see the vignettes \href{https://cbruneelzupanc.github.io/semiIVreg/reference/simul_data.html}{here} or by running \code{vignette("simul_data", package = "semiIVreg")}.

#' @details This function simulates data from the Generalized Roy Model with semi-IVs, following the simulation specified in Bruneel-Zupanc (2024). \cr
#' For more details about the exact specification, see the vignette \href{.https://cbruneelzupanc.github.io/semiIVreg/reference/simul_data.html}{here} or by running \code{vignette("simul_data", package = "semiIVreg")}.
#' One can use it to simulate general model with heterogenous treatment effects, but also restricted ones with homogenous treatment effects. \cr
#' `simul_data` was used to simulate the dataset available with this package, `data(roydata)` to obtain the simulated model with heterogenous treatment effect, and `data(roydata2)` to obtain the simulated model with homogenous treatment effect.



#' @param N Number of observations
#' @param model_type Type of model: "heterogenous" or "homogenous"
#' @param param_y0 Parameters for Y0 = (delta0, beta0, beta0X1, beta0X2) \cr
#' i.e., intercept, effects on w0, X_1, X_2 on Y0.
#' @param param_y1 Parameters for Y1: (delta1, beta1, beta1X1, beta1X2). \cr
#' i.e., intercept, effects w1, X1, X2 on Y1.
#' @param param_p Parameters for the selection: (alpha, alpha0, alpha1, alpha2, alphaX1, alphaX2)
#' i.e., intercept and effects of w0, w1, w0w1, Xbinary, Xcontinuous on the latent utility.
#' @param param_Z Parameters for the simulation of the semi-IVs: \cr
#' mean of W0 when X1=0, of W1 when X1=0, of W0 when X1=1, of W1 when X1=1; then variance of W0, W1, and covariance of W0 and W1.
#' @param param_genX Parameters for the covariates: p_X1, mu_X2, sigma_X2.
#' @param param_error Parameters for the error terms: depends on model_type: \cr
#' if heterogenous: variance of U0, U1, covariance of U0 and U1, variance of the cost (which has mean 0). \cr
#' if homogenous: variance of U, variance of V, covariance of U and V.

#' @return A data frame with the following columns:
#' \describe{
#'   \item{y}{The observed outcome.}
#'   \item{d}{The treatment.}
#'   \item{w0, w1}{The semi-IVs entering only D=0 and D=1.}
#'   \item{Xbinary, Xcontinuous}{Two covariates, one binary and one continuous.}
#'   \item{y0, y1}{The unobserved potential outcomes.}
#'   \item{P}{The unobserved true treatment probability.}
#'   \item{latent, V, Ud, U0, U1}{The unobserved shocks V. Ud is the normalized V ranks. U0 and U1 are the outcome shocks. `latent` gives the latent utility term in the selection equation.}
#' }

#' @examples
#' N = 10000; set.seed(12345)
#'
#' # Example 1: Heterogenous Treatment Effects.
#' model_type = "heterogenous"
#' param_error = c(1, 1, 0.6, 0.5) # var_u0, var_u1, cov_u0u1, var_cost
#' param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9)
#' # meanW0 Xbinary0, meanW1 Xbinary0, meanW0 Xbinary1, meanW1 Xbinary1, varW0, varW1, covW0W1
#' param_p = c(0, -0.7, 0.7, 0, 0, 0) # constant, W0, W1, W0xW1, Xbinary, Xcontinuous
#' param_y0 = c(3.2, 0.8, 0, 0) # intercept, W0, Xbinary, Xcontinuous;
#' param_y1 = c(3.2+0.4, 0.5, 0, 0) # the +0.4 = ATE; W1, Xbinary, Xcontinuous;
#' param_genX = c(0.4, 0, 2)
#'
#' data = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)
#'
#'
#' # Example 2: Homogenous Treatment Effects (constant MTE)
#' model_type = "homogenous"
#' param_error = c(1, 1.5, -0.6) # var_u, var_v, cov_uv
#' param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9)
#' param_p = c(0, -0.5, 0.5, 0, 0, 0) # the constant <=> mean_V
#' param_y0 = c(3.2, 0.8, 0, 0)
#' param_y1 = c(3.2+0.4, 0.5, 0, 0)
#' param_genX = c(0.4, 0, 2)
#'
#' data1 = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)
#'
#'
#' # Set the effects of w1 or w0 on its outcome to zero if want a valid IV, e.g.,
#' # param_y1 = c(3.2+0.4, 0, 0, 0) # w1 is a valid IV
#' # or: param_y0 = c(3.2, 0, 0, 0) # w0 is a valid IV


#' @references
#' Bruneel-Zupanc, C. (2023). Don't (fully) exclude me, it's not necessary! Identification with semi-IVs. arXiv preprint arXiv:2303.12667.
#'
#' Andresen, M. E. (2018). Exploring marginal treatment effects: Flexible estimation using Stata. The Stata Journal, 18(1), 118-158.
#'
#' Heckman, J. J., Urzua, S., & Vytlacil, E. (2006). Understanding instrumental variables in models with essential heterogeneity. The Review of Economics and Statistics, 88(3), 389-432.
#'
#' Heckman, J. J., & Vytlacil, E. J. (2007). Econometric evaluation of social programs, part II: Using the marginal treatment effect to organize alternative econometric estimators to evaluate social programs, and to forecast their effects in new environments. Handbook of econometrics, 6, 4875-5143.

#' @usage simul_data(N, model_type="heterogenous",
#'            param_y0, param_y1, param_p, param_Z, param_genX, param_error)

#' @export
simul_data = function(N, model_type="heterogenous", param_y0, param_y1, param_p, param_Z, param_genX, param_error) {

  # About the param error -> depends on model_type
  # param_error_heter = c(var_u0, var_u1, cov_u0u1, mean_cost var_cost) # if heterogenous
  # param_error_homo =  c(var_u, var_v, cov_uv, mean_v) # if homogenous

  # i) Shocks
  if(model_type == "heterogenous") {

    #param_error = param_error_heter

    # U0 and U1:
    mean_u0 = 0; mean_u1 = 0;
    var_u0 = param_error[1]; var_u1 = param_error[2]
    cov_u0u1 = param_error[3]
    Sigma = matrix(c(var_u0, cov_u0u1, cov_u0u1, var_u1), nrow=2, ncol=2, byrow=TRUE)
    mean_U = c(mean_u0, mean_u1) # to not get mixed with the intercept parameter
    U = mvrnorm(n=N, mu=mean_U, Sigma=Sigma)
    U0 = U[,1]; U1 = U[,2]

    # Cost shock
    mean_cost = 0 #param_error[4]; # mean cost is the same as the constant parameter in latent_p
    var_cost = param_error[4]
    C = rnorm(N, mean_cost, sqrt(var_cost))


    # Main proba shock
    V = -(U1 - U0 - C) # need to put the minus sign so that higher U1 increase proba that D=1

    # Normalized V ranks:
    # What is the distribution of V?
    # https://stats.stackexchange.com/questions/19948/what-is-the-distribution-of-the-sum-of-non-i-i-d-gaussian-variates
    #aX + bY ~ N(a mu_x + b mu_y, a^2 sigma^2_x + b^2 sigma^2_y + 2 ab sigma_(x,y))
    #V ~ N(mean_u0 - mean_u1 - mean_cost, var_u0 + var_u1 + 2*(-1)*cov_u0u1 + sd_cost^2)
    mean_V = mean_u0 - mean_u1 + mean_cost
    var_V = var_u0 + var_u1 + 2*(-1)*cov_u0u1 + var_cost
    Ud = pnorm(V, mean=mean_V, sd = sqrt(var_V))


  }

  # Homogenous TE if same U:
  if(model_type == "homogenous") {

    #param_error = param_error_homo

    # Simulate U and V directly:
    mean_u = 0; mean_v = 0 # param_error[4]; # Mean_V is the same parameter as the constant in latent P
    var_u = param_error[1]; var_v = param_error[2]
    cov_uv = param_error[3]
    Sigma = matrix(c(var_u, cov_uv, cov_uv, var_v), nrow=2, ncol=2, byrow=TRUE)
    mean_U = c(mean_u, mean_v) # to not get mixed with the intercept parameter
    ERROR = mvrnorm(n=N, mean_U, Sigma)
    U = ERROR[,1]; V = ERROR[,2] # the nonnormalized shocks

    U0 = U; U1 = U;

    # Normalized V ranks
    mean_V = mean_v; var_V = var_v
    Ud = pnorm(V, mean=mean_V, sd=sqrt(var_V))
  }


  # ii) Covariates
  # Generate two covariates -> so that can make "general estimation code" after.
  # Make one of the two correlated with semi-IVs,
  #   such that if do not include it in the regression, semi-IVs is not valid anymore? -> it depends?
  #   -> only if impacts the first stage and missing from it I think;

  # ii-a) Xbinary
  pXbinary = param_genX[1] # 2 Xbinarys with equal size
  Xbinary = rbinom(n=N, size=1, prob=pXbinary)

  # ii-b) Education of parents (continuous)
  mean_Xcontinuous = param_genX[2]; sd_Xcontinuous = param_genX[3]
  Xcontinuous = rnorm(N, mean=mean_Xcontinuous, sd=sd_Xcontinuous)


  # iii) semi-IVs
  mean_w0_s0 = param_Z[1]; mean_w1_s0 = param_Z[2]
  mean_w0_s1 = param_Z[3]; mean_w1_s1 = param_Z[4]
  var_w0 = param_Z[5]; var_w1 = param_Z[6]; cov_w0w1 = param_Z[7]
  SigmaZ = matrix(c(var_w0, cov_w0w1, cov_w0w1, var_w1), nrow=2, ncol=2, byrow=TRUE)

  # Binary X-specific semi-IVs draws (different mean)
  obs_Xbinary0 = which(Xbinary == 0); obs_Xbinary1 = which(Xbinary == 1)
  NXbinary0 = length(obs_Xbinary0); NXbinary1 = length(obs_Xbinary1)
  mean_Z0 = c(mean_w0_s0, mean_w1_s0)
  mean_Z1 = c(mean_w0_s1, mean_w1_s1)

  Z0 = mvrnorm(n=NXbinary0, mean_Z0, SigmaZ)
  Z1 = mvrnorm(n=NXbinary1, mean_Z1, SigmaZ)

  Z = matrix(rep(NA, N*2), ncol=2)
  Z[obs_Xbinary0,] = Z0; Z[obs_Xbinary1,] = Z1

  w0 = Z[,1]; w1 = Z[,2]



  # iv) Potential outcomes:
  Y0 = param_y0[1] + param_y0[2]*w0 + param_y0[3]*Xbinary + param_y0[4]*Xcontinuous + U0
  Y1 = param_y1[1] + param_y1[2]*w1 + param_y1[3]*Xbinary + param_y1[4]*Xcontinuous + U1


  # v) Discrete choice:
  # Logit formulation:
  latent = param_p[1] + w0*param_p[2] + w1*param_p[3] + w0*w1*param_p[4] + Xbinary*param_p[5] + Xcontinuous*param_p[6]

  Dstar = latent - V
  # Close to Heckman, Urzua and Vytlacil (2006) process. But change of sign for the latent
  # - latent to correspond to Heckman, Urzua and Vytlacil (2006) process.
  # Should be the same as Dstar = U1 - U0 - C - latent; since Ud = -(U1 - U0 - C)
  D = ifelse(Dstar > 0, 1, 0)

  P = pnorm(latent, mean=mean_V, sd=sqrt(var_V))


  # vi) Observable outcome:
  Y = rep(NA, N)
  Y[which(D==0)] = Y0[which(D==0)]
  Y[which(D==1)] = Y1[which(D==1)]

  ## netY0 = U0 + param_y0[1] #Y0 - param_y0[2]*w0;
  ## netY1 = U1 + param_y1[1] #Y1 - param_y1[2]*w1
  ## netY = rep(NA, N) # index sufficiency should hold for netY
  ## netY[which(D==0)] = netY0[which(D==0)]
  ## netY[which(D==1)] = netY1[which(D==1)]

  data = data.frame(id=1:N, y=Y, d=D, w0, w1, Xbinary, Xcontinuous,
                    P, y1=Y1, y0=Y0, latent, V, Ud, U0, U1) #, netY, netY0, netY1)
  return(data)
}

#' @import ggplot2
#' @import gridExtra
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm rbinom rnorm
#' @import data.table
#' @import KernSmooth
NULL
