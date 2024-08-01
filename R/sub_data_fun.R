#' @title Data construction functions
#'
#' @description These functions are used to construct the data from a given formula. Handles change from factor into several dummies for example.
#' They also create reference individuals at which to evaluate the MTE and MTR (if no default is provided).
#' For the numerical variables, take the average of the variable; for the factors, take the first level.
#'
#' @param formula The formula of the model
#' @param data The original dataset
#'
#' @usage
#' construct_data(formula, data)
#'
#' @rdname construct_data
#' @export
construct_data = function(formula, data) {
  # Given a formula, in a given dataset, construct the appropriate dataset and subformulas
  # In particular, transform factors into subvariables;

  # 1. Extract sub-formulas
  form_string = deparse1(formula); form_string = gsub(" ", "", form_string) # remove useless spacing
  split_form = unlist(strsplit(form_string, split="|", fixed=TRUE))

  # Transform into formulas and variables:
  form_yd = split_form[1]
  form_yd2 = unlist(strsplit(form_yd, split="~", fixed=TRUE))
  form_w0 = split_form[2]
  form_w1 = split_form[3]
  if(length(split_form) > 3) { cov_included = TRUE; form_covariates = split_form[4] } else { cov_included = FALSE }
  if(cov_included==TRUE) {
    formula_X = as.formula(paste0("~", form_w0, " + ", form_w1, " + ", form_covariates))
  } else {formula_X = as.formula(paste0("~", form_w0, " + ", form_w1)) }

  # 2. Deal with possible "weird" covariates (e.g., factor(x))
  # W0:
  #dataw0 = as.data.frame(model.matrix(as.formula(paste0("~-1 +", form_w0)), data=data))
  dataw0 = as.data.frame(model.matrix(as.formula(paste0("~", form_w0)), data=data)) # include intercept, to drop ref factor levels;
  dataw0 = dataw0[,which(colnames(dataw0) != "(Intercept)"), drop=FALSE]
  colnames(dataw0) = make.names(colnames(dataw0)) # change automatic factor names if invalid (e.g., if a level includes a minus or a space).

  # W1:
  #dataw1 = as.data.frame(model.matrix(as.formula(paste0("~-1 +", form_w1)), data=data))
  dataw1 = as.data.frame(model.matrix(as.formula(paste0("~", form_w1)), data=data)) # include intercept, to drop ref factor levels;
  dataw1 = dataw1[,which(colnames(dataw1) != "(Intercept)"), drop=FALSE]
  colnames(dataw1) = make.names(colnames(dataw1)) # change automatic factor names if invalid (e.g., if a level includes a minus or a space).

  # Other covariates:
  if(cov_included==TRUE) {
    #datacov = as.data.frame(model.matrix(as.formula(paste0("~-1 +", form_covariates)), data=data))
    datacov = as.data.frame(model.matrix(as.formula(paste0("~", form_covariates)), data=data)) # include intercept, to drop ref factor levels;
    datacov = datacov[,which(colnames(datacov) != "(Intercept)"), drop=FALSE]
    colnames(datacov) = make.names(colnames(datacov)) # change automatic factor names if invalid (e.g., if a level includes a minus or a space).
  }

  # Recreate the formulas and dataset with only dummy/continuous variables:
  # Change the formula
  form_w0 = paste0(colnames(dataw0), collapse="+")
  form_w1 = paste0(colnames(dataw1), collapse="+")
  if(cov_included==TRUE) { form_covariates = paste0(colnames(datacov), collapse="+") } else { form_covariates = "" }


  # 3. Extract variables:
  var_outcome = form_yd2[1]
  var_treatment = form_yd2[2] #gsub(" ", "", form_yd2[2]) # gsub to remove useless spacing
  var_w0 = unlist(strsplit(form_w0, "+", fixed=TRUE))
  var_w1 = unlist(strsplit(form_w1, "+", fixed=TRUE))
  if(cov_included==TRUE) { var_covariates = unlist(strsplit(form_covariates, "+", fixed=TRUE)) } else { form_covariates = ""; var_covariates = "" }


  # 4. Update output dataset

  dataw1merge = subset(dataw1, select=colnames(dataw1)[which(!colnames(dataw1) %in% colnames(dataw0))]) # avoid repeating column names
  data_orig = data;
  Xdat = cbind(dataw0, dataw1merge)
  if(cov_included == TRUE) {
    datacovmerge = subset(datacov, select=colnames(datacov)[which(!colnames(datacov) %in% colnames(dataw0))])
    Xdat = cbind(Xdat, datacovmerge)
  }


  # Complete updated formula
  new_formula = paste0(var_outcome, "~", var_treatment, "|", form_w0, "|", form_w1)
  if(cov_included==TRUE) { new_formula = paste0(new_formula, "|", form_covariates) }
  new_formula=as.formula(new_formula)

  # Output:
  res = list(data=Xdat,
             formula = new_formula,
             form_w0=form_w0, form_w1=form_w1, form_covariates=form_covariates, formula_X_orig = formula_X,
             var_outcome=var_outcome, var_treatment=var_treatment, var_w0=var_w0, var_w1=var_w1, var_covariates=var_covariates)
  return(res)
}


#' @rdname construct_data
#' @usage create_ref_indiv(formula, data)
#' @export
create_ref_indiv = function(formula, data) {

  # Create the "reference individual" = mean value of the covariates + reference level if factors.

  # 0) Formulas
  form_string = deparse1(formula); form_string = gsub(" ", "", form_string) # remove useless spacing
  split_form = unlist(strsplit(form_string, split="|", fixed=TRUE))

  # Transform into formulas and variables:
  form_yd = split_form[1]
  form_yd2 = unlist(strsplit(form_yd, split="~", fixed=TRUE))
  form_w0 = split_form[2]
  form_w1 = split_form[3]
  if(length(split_form) > 3) { cov_included = TRUE; form_covariates = split_form[4] } else { cov_included = FALSE }

  var_outcome = form_yd2[1]
  var_treatment = form_yd2[2] #gsub(" ", "", form_yd2[2]) # gsub to remove useless spacing
  var_w0 = unlist(strsplit(form_w0, "+", fixed=TRUE))
  var_w1 = unlist(strsplit(form_w1, "+", fixed=TRUE))
  if(cov_included==TRUE) { var_covariates = unlist(strsplit(form_covariates, "+", fixed=TRUE)) } else { form_covariates = ""; var_covariates = "" }

  if(cov_included==TRUE) {
    formula_X = as.formula(paste0("~", form_w0, " + ", form_w1, " + ", form_covariates))
  } else {formula_X = as.formula(paste0("~", form_w0, " + ", form_w1)) }


  # 1) Create ref-individual
  # Special care on factor variables which are specified as "factor(x)" in the formula;
  name_all_X = all.vars(formula_X)
  Xdat = subset(data, select=name_all_X)

  # # -- subcode to deal with possible "factor(x) in the formulas --
  # # If factors are just informed as "factor(x)" in the formula, would recognize it as numerical, which we don't want;
  # # Use regular expression to identify factor() calls
  # formula_X_chars = paste0(as.character(formula_X), collapse="")
  # factor_pattern <- "factor\\(([^)]+)\\)"
  # factor_matches <- regmatches(formula_X_chars, gregexpr(factor_pattern, formula_X_chars))[[1]]
  # factor_vars = c()
  # for (match in factor_matches) {
  #   factor_var <- gsub("factor\\(|\\)", "", match)
  #   factor_vars <- c(factor_vars, factor_var)
  # }
  # # Ensure factor_vars are unique
  # factor_vars <- unique(factor_vars)
  # # Extract variable names inside factor()
  # for (match in factor_matches) {
  #   factor_var <- gsub("factor\\(|\\)", "", match)
  #   factor_vars <- c(factor_vars, factor_var)
  # }
  # # Ensure factor_vars are unique
  # factor_vars <- unique(factor_vars)
  # # -- end of subcode

  # If ref_indiv does not exists, create it:
  #if(is.null(ref_indiv)) { # only run if is null ref_indiv
  ref_indiv = Xdat[1,] # just to have the proper names
  for(j in 1:ncol(Xdat)) {
    var_selected = Xdat[,j]; name_selected = colnames(ref_indiv)[j]
    #if(name_selected %in% factor_vars) { var_selected = factor(var_selected) } # transform into factor for the reference individual.
    if(is.numeric(var_selected)) { ref_indiv[,j] = mean(var_selected) }
    if(is.factor(var_selected)) { ref_indiv[,j] = levels(var_selected)[1] }
    if(is.character(var_selected)) { var_selected = as.factor(var_selected); ref_indiv[,j] = levels(var_selected)[1] }
    if(!(is.numeric(var_selected) | is.factor(var_selected) | is.character(var_selected))) {
      stop(paste0("Error, the variable ", colnames(Xdat)[j], " is neither numerical nor a factor/character."))
    }
  }
  #}

  # Additionally, for factor data, ensures that these are saved as factor with the proper levels (very important for predictions later)
  for(j in 1:ncol(ref_indiv)) {
    name_refj = colnames(ref_indiv)[j]
    Xdatj = Xdat[, which(colnames(Xdat) == name_refj)]
    if(is.factor(Xdatj)) { ref_indiv[,j] = factor(ref_indiv[,j], levels=levels(Xdatj)) }
  }

  #ref_indiv$id = 1:nrow(ref_indiv)
  #Xdat$id = NA; # to use later just to avoid bugs in predictions if use splines or poly on covariates

  return(ref_indiv)
}



#' @rdname construct_data
#' @usage transform_factor(formula, data)
#' @export
transform_factor = function(formula, data) {
  # transform variables which are used as factor into a factor in the origin data;

  # 0) Formulas
  form_string = deparse1(formula); form_string = gsub(" ", "", form_string) # remove useless spacing
  split_form = unlist(strsplit(form_string, split="|", fixed=TRUE))

  # Transform into formulas and variables:
  form_yd = split_form[1]
  form_yd2 = unlist(strsplit(form_yd, split="~", fixed=TRUE))
  form_w0 = split_form[2]
  form_w1 = split_form[3]
  if(length(split_form) > 3) { cov_included = TRUE; form_covariates = split_form[4] } else { cov_included = FALSE }

  var_outcome = form_yd2[1]
  var_treatment = form_yd2[2] #gsub(" ", "", form_yd2[2]) # gsub to remove useless spacing
  var_w0 = unlist(strsplit(form_w0, "+", fixed=TRUE))
  var_w1 = unlist(strsplit(form_w1, "+", fixed=TRUE))
  if(cov_included==TRUE) { var_covariates = unlist(strsplit(form_covariates, "+", fixed=TRUE)) } else { form_covariates = ""; var_covariates = "" }
  if(cov_included==TRUE) {
    formula_X = as.formula(paste0("~", form_w0, " + ", form_w1, " + ", form_covariates))
  } else {formula_X = as.formula(paste0("~", form_w0, " + ", form_w1)) }

  # 1) Create ref-individual
  # Special care on factor variables which are specified as "factor(x)" in the formula;
  name_all_X = all.vars(formula_X)
  Xdat = subset(data, select=name_all_X)

  # -- subcode to deal with possible "factor(x) in the formulas --
  # If factors are just informed as "factor(x)" in the formula, would recognize it as numerical, which we don't want;
  # Use regular expression to identify factor() calls
  formula_X_chars = paste0(as.character(formula_X), collapse="")
  factor_pattern <- "factor\\(([^)]+)\\)"
  factor_matches <- regmatches(formula_X_chars, gregexpr(factor_pattern, formula_X_chars))[[1]]
  factor_vars = c()
  for (match in factor_matches) {
    factor_var <- gsub("factor\\(|\\)", "", match)
    factor_vars <- c(factor_vars, factor_var)
  }
  # Ensure factor_vars are unique
  factor_vars <- unique(factor_vars)
  # Extract variable names inside factor()
  for (match in factor_matches) {
    factor_var <- gsub("factor\\(|\\)", "", match)
    factor_vars <- c(factor_vars, factor_var)
  }
  # Ensure factor_vars are unique
  factor_vars <- unique(factor_vars)
  # -- end of subcode

  if(length(factor_vars) > 0) {
    for(j in 1:length(factor_vars)) {
      namej = factor_vars[j]
      column = which(colnames(data) == namej)
      data[, column] = as.factor(data[, column])
    }
  }

  return(data)
}



NULL
