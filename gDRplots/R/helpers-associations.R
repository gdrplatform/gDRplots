#' Calculate linear associations
#' 
#' Calculate the linear model associations between dependent variables and response variable(s) of interest.
#' 
#' @param X \code{matrix} Dependent variables data matrix (rows are samples, columns are features).
#' Must have the same length as \code{Y}.
#' @param Y \code{vector} or \code{matrix} Experimental response data (rows are samples).
#' Must have the same length as \code{X}.
#' 
#' @author James Hawley
#' @seealso \code{\link[cdsrmodels]{lin_associations}}
#' 
#' @keywords internal
#' 
#' @export 
calc_assoc <- function(X, Y) {

  checkmate::assert_matrix(X, mode = "numeric")
  checkmate::assert_names(rownames(X))
  checkmate::assert_multi_class(Y, c("matrix", "numeric"))
  if (is.matrix(Y)) checkmate::assert_names(rownames(Y))
  if (is.vector(Y)) checkmate::assert_names(names(Y))

  if (is.vector(Y)) {
    .calc_assoc_vector(X = X, Y = Y)
  } else if (is.matrix(Y)) {
    .calc_assoc_matrix(X = X, Y = Y)
  }
}

#' Calculate the linear model associations between dependent variables and some response variable of interest.
#' @param X \code{matrix} Dependent variables data matrix (rows are samples, columns are features).
#' Must have the same length as \code{Y}.
#' @param Y \code{vector} Experimental response data (rows are samples).
#' Must have the same length as \code{X}.
#' 
#' @return \code{data.table}
#' 
#' @author James Hawley
#' @seealso \code{\link[cdsrmodels]{lin_associations}}
#' 
#' @keywords internal 
.calc_assoc_vector <- function(X, Y) {
  
  checkmate::assert_matrix(X, mode = "numeric")
  checkmate::assert_names(rownames(X))
  checkmate::assert_class(Y, "numeric")
  checkmate::assert_names(names(Y))
  
  if (stats::sd(Y, na.rm = TRUE) == 0) {
    warning("Y has no variance, rendering all associations void. Please double check this is correct.")
    na_dt <- data.table::data.table(
      feature = names(Y)
      , est_beta = NA_real_
      , est_beta_se = NA_real_
      , posterior_mean = NA_real_
      , posterior_sd = NA_real_
      , prob_negative = NA_real_
      , prob_positive = NA_real_
      , rho = NA_real_
      , p_value = NA_real_
      , q_value = NA_real_
      , s_value = NA_real_
      , lfsr = NA_real_
      , lfdr = NA_real_
    )
    return(na_dt)
  }
  
  # use `cdsr_models` to calculate the linear model coefficients efficiently on large matrices
  res <- cdsrmodels::lin_associations(X = X, Y = Y)
  
  # convert results from a `matrix` to a `data.table`
  res_dt <- data.table::as.data.table(res$res.table)
  
  # store the correlation coefficient in the table
  available_genes <- rownames(res$res.table)
  res_dt[, rho := res$rho[available_genes, 1]]
  
  # return with a re-order columns for a more human-friendly output
  res_dt[
    , .(
      feature = ind.var
      , est_beta = betahat
      , est_beta_se = sebetahat
      , posterior_mean = PosteriorMean
      , posterior_sd = PosteriorSD
      , prob_negative = NegativeProb
      , prob_positive = PositiveProb
      , rho
      , p_value = p.val
      , q_value = qvalue
      , s_value = svalue
      , lfsr
      , lfdr
    )
  ]
}

#' Calculate the linear model associations between dependent variables and response variables of interest.
#' 
#' @param X \code{matrix} Dependent variables data matrix (rows are samples, columns are features).
#' Must have the same length as \code{Y}.
#' @param Y \code{matrix} Experimental response data (rows are samples).
#' Must have the same length as \code{X}.
#' 
#' @return \code{data.table}
#' 
#' @author James Hawley
#' @seealso \code{\link[cdsrmodels]{lin_associations}}
#' 
#' @keywords internal
.calc_assoc_matrix <- function(X, Y) {
  
  checkmate::assert_matrix(X, mode = "numeric")
  checkmate::assert_names(rownames(X))
  checkmate::assert_matrix(Y, mode = "numeric")
  checkmate::assert_names(rownames(Y))
  
  # calculate if any columns have zero variance
  zero_var_cols <- which(apply(Y, 2, stats::sd, na.rm = TRUE) == 0)
  zero_var_cols_len <- length(zero_var_cols)
  
  if (zero_var_cols_len > 0) {
    warning(
      paste0(
        "The following columns in Y have no variance, rendering associations void. 
        Please double check this is correct: "
        , paste(names(zero_var_cols), collapse = ", ")
      )
    )
    na_dt <- data.table::data.table(
      feature = rep(x = colnames(X), times = zero_var_cols_len)
      , response = rep(x = colnames(Y)[zero_var_cols], each = dim(X)[2])
      , est_beta = NA_real_
      , est_beta_se = NA_real_
      , posterior_mean = NA_real_
      , posterior_sd = NA_real_
      , prob_negative = NA_real_
      , prob_positive = NA_real_
      , rho = NA_real_
      , p_value = NA_real_
      , q_value = NA_real_
      , s_value = NA_real_
      , lfsr = NA_real_
      , lfdr = NA_real_
    )
    
    # remove the zero-variance columns from Y
    # The `as.matrix()` ensures that even if Y only has 1 column remaining,
    # it is still a matrix object, and doesn't get reduced to a vector.
    # If it was, it would cause later problems trying to name the `res_dt`
    # columns, correctly.
    y_colnames <- colnames(Y)
    Y <- as.matrix(Y[, -zero_var_cols])
    colnames(Y) <- y_colnames[-zero_var_cols]
  }
  
  # early return of the `na_dt` table if there are no columns of Y with non-zero variance
  if (dim(Y)[2] == 0) {
    return(na_dt)
  }
  
  # use `cdsr_models` to calculate the linear model coefficients efficiently on large matrices
  res <- cdsrmodels::lin_associations(X = X, Y = Y)
  
  # subset the data so that the correlations match the output of the table
  available_genes <- rownames(stats::na.omit(res$rho))
  
  # convert results from a `matrix` to a `data.table`
  res_dt <- data.table::as.data.table(res$res.table)
  
  # if Y only has 1 column (either originally, or because others were filtered out)
  # ensure that the `dep.var` column is listed, otherwise the renaming and 
  # re-ordering step below will fail
  if (!"dep.var" %in% names(res_dt)) {
    res_dt[, dep.var := rep(colnames(Y), each = NROW(available_genes))]
  }
  
  # combine each column of the correlation matrix into a single giant vector and
  # store it as the appropriate column in `res_dt`
  res_dt[, rho := vapply(res$rho[available_genes, ], as.vector, numeric(1))]
  
  # re-order columns for a more human-friendly output
  res_dt <- res_dt[
    , .(
      feature = ind.var
      , response = dep.var
      , est_beta = betahat
      , est_beta_se = sebetahat
      , posterior_mean = PosteriorMean
      , posterior_sd = PosteriorSD
      , prob_negative = NegativeProb
      , prob_positive = PositiveProb
      , rho
      , p_value = p.val
      , q_value = qvalue
      , s_value = svalue
      , lfsr
      , lfdr
    )
  ]
  
  # add in NA values from any columns that had 0 variance, if necessary
  if (zero_var_cols_len > 0) {
    res_dt <- data.table::rbindlist(list(res_dt, na_dt))
  }
  
  # return the combined table
  res_dt
}
