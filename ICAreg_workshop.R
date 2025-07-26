# Load required packages
pacman::p_load(
  nlme,
  dplyr,
  ica,
  Matrix,
  pbapply
)


boot1 <- function(data, X, f1, f1X, has_intercept, formula, method,
                  dependent_var, independent_vars, 
                  independent_P_vars = NA, independent_X_vars = NA) {
  
  data <- data[, !colnames(data) %in% c("control_func")]
  
  full_formula <- as.formula(paste(all.vars(formula)[1], "~", 
                                   gsub("\\|", "+", as.character(formula)[3])))
  
  repeat {
    
    data_cleaned <- data %>%
      sample_n(size = nrow(data), replace = TRUE)
    design_matrix <- model.matrix(full_formula, data = data_cleaned)
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]
    
    if (rank_A == ncol(design_matrix)) { break }
    
  }
  
  if (length(f1) == 1) {
    
    # ICA
    d1 <- ica::ica(X = data_cleaned2, nc = ncol(data_cleaned2), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned$control_func <- control_func
    
    # Regressions with control function
    if (!has_intercept) {
      
      lm0 <- lm(Formula::as.Formula(formula), data_cleaned)
      lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . 
                                                          + control_func -1)), 
                data_cleaned)
      
    } else {
      
      lm0 <- lm(Formula::as.Formula(formula), data_cleaned)
      lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . 
                                                          + control_func)), 
                data_cleaned)
      
    }
    
    Estimates <- lm1$coefficients
    
  } else {
    
    # first-stage regression for depvar
    lm_Y <- lm(formula = as.formula(paste(dependent_var, as.character(f1X)[2], sep = " ~ ")), 
               data = data_cleaned)
    lm_Y_residuals <- residuals(lm_Y)
    
    # first-stage regression for endogs
    lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                             ncol = length(independent_P_vars))
    colnames(lm_P_residuals) <- independent_P_vars
    
    
    for (i in seq_along(independent_P_vars)) {
      
      p_var <- independent_P_vars[i]
      
      formula_P <- as.formula(paste(p_var, "~", paste(independent_X_vars, collapse = " + ")))
      lm_P <- lm(formula_P, data = data_cleaned)
      lm_P_residuals[, i] <- residuals(lm_P)
      
    }
    
    # Apply the ICA on the residuals
    data_residuals <- cbind(lm_Y_residuals, lm_P_residuals)
    
    # ICA
    d1 <- ica::ica(X = data_residuals, nc = ncol(data_residuals), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned$control_func <- control_func
    
    lm0 <- lm(Formula::as.Formula(formula), data_cleaned)
    lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . 
                                                        + control_func)), 
              data_cleaned)
    
    Estimates <- lm1$coefficients
    
  }
  
  return(Estimates)
  
}
boot2 <- function(data, X, f1, f1X, has_intercept, formula, method,
                  dependent_var, independent_vars, 
                  independent_P_vars = NA, independent_X_vars = NA) {
  
  data <- data[, !colnames(data) %in% c("control_func")]
  
  full_formula <- as.formula(paste(all.vars(formula)[1], "~", 
                                   gsub("\\|", "+", as.character(formula)[3])))
  
  repeat {
    
    data_cleaned <- data %>%
      sample_n(size = nrow(data), replace = TRUE)
    design_matrix <- model.matrix(full_formula, data = data_cleaned)
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]
    
    if (rank_A == ncol(design_matrix)) { break }
    
  }
  
  if (length(f1) == 1) {
    
    data_cleaned <- data_cleaned[, c(dependent_var, independent_vars)]
    
    # ICA
    d1 <- ica::ica(X = data_cleaned, nc = ncol(data_cleaned), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned$control_func <- control_func
    
    # first-stage regression for endogs
    lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                             ncol = length(independent_vars))
    colnames(lm_P_residuals) <- independent_vars
    
    
    for (i in seq_along(independent_vars)) {
      
      p_var <- independent_vars[i]
      
      formula_P <- as.formula(paste(p_var, "~", paste("control_func", collapse = " + ")))
      lm_P <- lm(formula_P, data = data_cleaned)
      lm_P_residuals[, i] <- residuals(lm_P)
      
    }
    
    data_cleaned[, independent_vars] <- lm_P_residuals
    
    # Regressions
    if (has_intercept) {
      
      lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -control_func")), data = data_cleaned)
      
    } else {
      
      lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -1 -control_func")), data = data_cleaned)
      
    }
    
    Estimates <- lm1$coefficients
    
  } else {
    
    # first-stage regression for depvar
    lm_Y <- lm(formula = as.formula(paste(dependent_var, as.character(f1X)[2], sep = " ~ ")), 
               data = data_cleaned)
    lm_Y_residuals <- residuals(lm_Y)
    
    # first-stage regression for endogs
    lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                             ncol = length(independent_P_vars))
    colnames(lm_P_residuals) <- independent_P_vars
    
    
    for (i in seq_along(independent_P_vars)) {
      
      p_var <- independent_P_vars[i]
      
      formula_P <- as.formula(paste(p_var, "~", paste(independent_X_vars, collapse = " + ")))
      lm_P <- lm(formula_P, data = data_cleaned)
      lm_P_residuals[, i] <- residuals(lm_P)
      
    }
    
    # Apply the ICA on the residuals
    data_residuals <- cbind(lm_Y_residuals, lm_P_residuals)
    
    # ICA
    d1 <- ica::ica(X = data_residuals, nc = ncol(data_residuals), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned$control_func <- control_func
    
    # first-stage regression for endogs
    lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                             ncol = length(independent_P_vars))
    colnames(lm_P_residuals) <- independent_P_vars
    
    for (i in seq_along(independent_P_vars)) {
      
      p_var <- independent_P_vars[i]
      
      formula_P <- as.formula(paste(p_var, "~", paste("control_func", collapse = " + ")))
      lm_P <- lm(formula_P, data = data_cleaned)
      lm_P_residuals[, i] <- residuals(lm_P)
      
    }
    
    data_cleaned[, independent_P_vars] <- lm_P_residuals
    
    # Regressions
    if (has_intercept) {
      
      lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -control_func")), data = data_cleaned)
      
    } else {
      
      lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -1 -control_func")), data = data_cleaned)
      
    }
    
    Estimates <- lm1$coefficients
    
  }
  
  return(Estimates)
  
}

ica_reg <- function(formula, data, method = "jade", CF = FALSE, nboots = 199) {
  
  # seperate endogenous and exogenous regressor(s)
  f1 <- nlme::splitFormula(formula, sep = "|")
  
  # check if intercept is removed
  has_intercept <- attr(terms(f1[[1]]), "intercept") == 1
  
  # no exogenous regressors
  if (length(f1) == 1) {
    
    f1P <- f1[[1]]
    
    # Check if all variables exist in the data
    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following variables are missing in the data:", 
                 paste(missing_vars, collapse=", ")))
      
    }
    
    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[variables], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)
    
    if (!all(numeric_vars)) {
      stop("The following variables are not numeric: ", paste(variables[!numeric_vars], collapse = ", "))
    }
    
    if (any(constant_vars)) {
      stop("The following variables are constant: ", paste(variables[constant_vars], collapse = ", "))
    }
    
    dependent_var <- all.vars(formula)[1]
    independent_vars <- all.vars(f1P)
    
    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_vars))) %>%
      na.omit() %>%
      as.data.frame()
    
    # Check if design matrix has full column rank
    design_matrix <- model.matrix(f1P, data = data_cleaned)
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]
    
    if(rank_A != ncol(design_matrix)) {
      
      stop(paste("Design matrix is rank deficient"))
      
    }
    
    # check if intercept is included
    if (has_intercept) { design_matrix1 <- design_matrix[, -1] } else { 
      design_matrix1 <- design_matrix }
    data_cleaned1 <- cbind(data_cleaned[, 1], design_matrix1)
    colnames(data_cleaned1)[1] <- colnames(data_cleaned)[1]
    
    # ICA
    d1 <- ica::ica(X = data_cleaned1, nc = ncol(data_cleaned1), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned1 <- as.data.frame(data_cleaned1)
    data_cleaned1$control_func <- control_func
    
    # Control function approach
    if (CF == TRUE) {
      
      # Regression
      if (has_intercept) {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned1)[1], "~ .")), data = data_cleaned1)
        
      } else {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned1)[1], "~ . -1")), data = data_cleaned1)
        
      }
      
      Estimates <- lm1$coefficients
      
      # Bootstrapping
      print("Estimation done. calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots, function(i) boot1(data = data_cleaned1, X = i, 
                                                      f1 = f1, has_intercept = has_intercept, 
                                                      formula = formula, method = method,
                                                      dependent_var = dependent_var,
                                                      independent_vars = independent_vars))
      
    } else if (CF == FALSE) {
      
      # first-stage regression for endogs
      lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                               ncol = length(independent_vars))
      colnames(lm_P_residuals) <- independent_vars
      
      
      for (i in seq_along(independent_vars)) {
        
        p_var <- independent_vars[i]
        
        formula_P <- as.formula(paste(p_var, "~", paste("control_func", collapse = " + ")))
        lm_P <- lm(formula_P, data = data_cleaned)
        lm_P_residuals[, i] <- residuals(lm_P)
        
      }
      
      data_cleaned[, independent_vars] <- lm_P_residuals
      
      # Regressions
      if (has_intercept) {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned1)[1], "~ .")), data = data_cleaned)
        
      } else {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned1)[1], "~ . -1")), data = data_cleaned)
        
      }
      
      Estimates <- lm1$coefficients
      
      # Bootstrapping
      print("Estimation done. calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots, function(i) boot2(data = data, X = i, 
                                                      f1 = f1, has_intercept = has_intercept, 
                                                      formula = formula, method = method,
                                                      dependent_var = dependent_var,
                                                      independent_vars = independent_vars))
      
    }
    
    if (is.numeric(trapped)) { ses <- sd(trapped) } else { ses <- apply(trapped, 1, sd) }
    
    Estimates1 <- cbind(Estimates, ses)
    colnames(Estimates1) <- c("Estimate", "Std. Error")
    
    ############################################################################
    
  } else {
    
    f1P <- f1[[1]]
    f1X <- f1[[2]]
    
    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following endogenous variables are missing in the data:", 
                 paste(missing_vars, collapse=", ")))
      
    }
    
    variables <- all.vars(f1X)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {
      
      stop(paste("The following exogenous variables are missing in the data:", 
                 paste(missing_vars, collapse=", ")))
      
    }
    
    dependent_var <- all.vars(formula)[1]
    independent_P_vars <- all.vars(f1P)
    independent_X_vars <- all.vars(f1X)
    
    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[independent_P_vars], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)
    
    if (!all(numeric_vars)) {
      stop("Only continuous variables can be endogenous. The following variables are not numeric: ", paste(variables[!numeric_vars], collapse = ", "))
    }
    
    if (any(constant_vars)) {
      stop("The following variables are constant: ", paste(variables[constant_vars], collapse = ", "))
    }
    
    # all variables
    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_P_vars, independent_X_vars))) %>%
      na.omit() %>%
      as.data.frame()
    
    # Check if design matrix has full column rank
    full_formula <- as.formula(paste(all.vars(formula)[1], "~", 
                                     gsub("\\|", "+", as.character(formula)[3])))
    
    design_matrix <- model.matrix(full_formula, data = data_cleaned)
    
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]
    
    if(rank_A != ncol(design_matrix)) {
      
      stop(paste("Design matrix is rank deficient"))
      
    }
    
    # first-stage regression for depvar
    lm_Y <- lm(formula = as.formula(paste(dependent_var, as.character(f1X)[2], sep = " ~ ")), 
               data = data_cleaned)
    lm_Y_residuals <- residuals(lm_Y)
    
    # first-stage regression for endogs
    lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                             ncol = length(independent_P_vars))
    colnames(lm_P_residuals) <- independent_P_vars
    
    for (i in seq_along(independent_P_vars)) {
      
      p_var <- independent_P_vars[i]
      
      formula_P <- as.formula(paste(p_var, "~", paste(independent_X_vars, collapse = " + ")))
      lm_P <- lm(formula_P, data = data_cleaned)
      lm_P_residuals[, i] <- residuals(lm_P)
      
    }
    
    # Apply the ICA on the residuals
    data_residuals <- cbind(lm_Y_residuals, lm_P_residuals)
    
    # ICA
    d1 <- ica::ica(X = data_residuals, nc = ncol(data_residuals), method = method)
    d2 <- d1$S
    
    # Extract the normal distribution
    ks_normality <- apply(d2, 2, function(x) {
      ks_test <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x)))
      return(ks_test$statistic)
    })
    
    control_func <- d2[, which.min(ks_normality)]
    data_cleaned$control_func <- control_func
    
    # Control function approach
    if (CF == TRUE) {
      
      # Regressions with control function
      lm0 <- lm(Formula::as.Formula(formula), data_cleaned)
      lm1 <- lm(Formula::as.Formula(stats::update.formula(stats::formula(lm0), ~ . 
                                                          + control_func)), 
                data_cleaned)
      Estimates <- lm1$coefficients
      
      # Bootstrapping
      print("Estimation done. calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots, function(i) boot1(data = data_cleaned, X = i, f1X = f1X,
                                                      f1 = f1, has_intercept = has_intercept, 
                                                      formula = formula, method = method,
                                                      dependent_var = dependent_var, 
                                                      independent_X_vars = independent_X_vars,
                                                      independent_P_vars = independent_P_vars))
      
    } else if (CF == FALSE) {
      
      # first-stage regression for endogs
      lm_P_residuals <- matrix(nrow = nrow(data_cleaned), 
                               ncol = length(independent_P_vars))
      colnames(lm_P_residuals) <- independent_P_vars
      
      for (i in seq_along(independent_P_vars)) {
        
        p_var <- independent_P_vars[i]
        
        formula_P <- as.formula(paste(p_var, "~", paste("control_func", collapse = " + ")))
        lm_P <- lm(formula_P, data = data_cleaned)
        lm_P_residuals[, i] <- residuals(lm_P)
        
      }
      
      data_cleaned[, independent_P_vars] <- lm_P_residuals
      
      # Regressions
      if (has_intercept) {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -control_func")), data = data_cleaned)
        
      } else {
        
        lm1 <- lm(formula = as.formula(paste(names(data_cleaned)[1], "~ . -1 -control_func")), data = data_cleaned)
        
      }
      
      Estimates <- lm1$coefficients
      
      # Bootstrapping
      print("Estimation done. calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots, function(i) boot2(data = data_cleaned, X = i, f1X = f1X,
                                                      f1 = f1, has_intercept = has_intercept, 
                                                      formula = formula, method = method,
                                                      dependent_var = dependent_var, 
                                                      independent_X_vars = independent_X_vars,
                                                      independent_P_vars = independent_P_vars))
      
    }
    
    ses <- apply(trapped, 1, sd)
    
    Estimates1 <- cbind(Estimates, ses)
    colnames(Estimates1) <- c("Estimate", "Std. Error")
    
  }
  
  ##############################################################################
  
  # Identification checks
  ks_test1 <- suppressWarnings(stats::ks.test(control_func, "pnorm"))
  warning_flag <- any(duplicated(control_func))
  
  if (ks_test1$p.value< .1) {warning("Joint component may not be normally distributed: Kolmogorov-Smirnov p = ", 
                                     paste(ks_test1$p.value, 
                                           collapse = ""), call. = FALSE)}
  
  if (warning_flag) {warning("Endogenous regressors contain ties (repeated values)", 
                             call. = FALSE)}
  
  return(list(Estimates1, control_func))
  
}
