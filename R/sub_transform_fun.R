#' @title Control functions transformations for selection probabilities
#'
#' @description These functions provides pre-specified transformations to control flexibly for the selection probabilities in the regression. \cr
#' These correspond to \eqn{\kappa_d(p)} and corresponding \eqn{k_d(u)} in Bruneel-Zupanc (2024). \cr
#' Special functions are used for homogenous treatment effect specifications because the code is different.
#' July 2024: for now, only polynomial transformations are encoded.
#'

#' @details See Andresen (2018) or Bruneel-Zupanc (2024) for computation details linking \eqn{\kappa_d(p)} and the corresponding corresponding \eqn{k_d(u)}. \cr
#' \eqn{\kappa_1(p) = E(U_1 | U_D \leq p) } and \eqn{\kappa_0(p) = E(U_0 | U_D>p) } while \eqn{k_d(u) = E(U_d | U_D=u)}. \cr
#'
#' In the case of homogenous treatment effects: \eqn{k_0(u) = k_1(u)}. This provides some restriction on \eqn{\kappa}, hence the special functions.

#' @param p Vector of propensity scores to transform into a flexible function.
#' @param d Which potential outcome to consider (only needed for \eqn{k_d(u)} with heterogenous treatment effects).
#' @param pol_degree Degree of the polynomial transformation.

#'
#' @examples
#' v = seq(0.1, 0.9, by=0.1)
#' # Transformations for general Heterogenous TE functions:
#' Kappa_fun(p=v, pol_degree=6)
#' k1u = kdu_transform_fun(v, d=1, pol_degree=6)
#'
#' # Transformations for Homogenous TE functions:
#' Kappa_homogenous_fun(p=v, pol_degree=6)
#' ku = ku_transform_homogenous_fun(v, pol_degree=6) # no d anymore, same for both d here;

#' @name Kappa_fun
#' @rdname Kappa_fun
#' @export
Kappa_fun = function(p, pol_degree=5) { return(poly(p, pol_degree, raw=TRUE)) } # Kappa_fun = E[U1 | V < p] if D = 1; E[U0 | V > p] if D=0


#' @rdname Kappa_fun
#' @export
kdu_transform_fun = function(u, d, pol_degree=5) {
  # E[U | V=v]: obtained from derivatives of Kappa (see Andresen 2018 for formula about what kdu corresponds to)
  # returns Kappa(u) + u*dKappa(u)/du if d=1, or Kappa(u)-(1-u)*dKappa(u)/du if d=0
  # If polynomial Kappa, trivial functional form

  pol = poly(u, pol_degree, raw=TRUE)
  multipliers = 2:(pol_degree+1)

  # D=1, return: 2u; 3u^2; 4u^3; ...
  res1 = sweep(pol, 2, multipliers, '*')

  # D=0, return: 2u - 1; 3u^2 - 2u; 4u^3 - 3u; ...
  if(pol_degree > 1) {
    pol2 = poly(u, pol_degree - 1, raw=TRUE)
    ppol2 = cbind(rep(1, length(u)), pol2)
  } else { ppol2 = cbind(rep(1, length(u))) }
  multipliers2 = 1:pol_degree

  res0_add = sweep(ppol2, 2, multipliers2, '*')
  res0 = res1 - res0_add #

  if(length(d) == 1 & length(u) > 1) { d=rep(d, length(u)) }
  res = res0;
  res[which(d==1),] = res1[which(d==1),]

  return(res)
}


#' @rdname Kappa_fun
#' @export
Kappa_homogenous_fun = function(p, pol_degree=5) {
  # Kappa functions for homogenous regression: need to include intercept
  # include intercept term for the homogenous regression
  pol = Kappa_fun(p, pol_degree); const = rep(1,length(p))
  return(cbind(const, pol))
}

#' @rdname Kappa_fun
#' @export
ku_transform_homogenous_fun = function(u, pol_degree=5) { # no d here, same for both, building on d=1 (given the specification I use)
  res1 = kdu_transform_fun(u, d=1, pol_degree) # builds off the d=1 case
  const = rep(1, length(u)) # just need to add a constant (cf derivation)
  return(cbind(const, res1))
}

NULL
