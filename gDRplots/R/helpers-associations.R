#' Calculate linear associations
#' 
#' Calculate the linear model associations between dependent variables and response variable(s) of interest.
#' 
#' @note inspired by the \code{calc_assoc} function written by James Hawley
#' 
#' @param X \code{matrix} dependent variables data matrix (rows are samples, columns are features).
#' Must have the same number of rows as matrix \code{Y} or equal to length of vector \code{Y}
#' @param Y \code{vector} or \code{matrix} experimental response data (rows are samples).
#' When \code{Y} is a matrix must have the same number of rows as matrix \code{X}; 
#' when \code{y} is a vector - its length has to be equal to number of rows in matrix \code{X}.
#' 
#' @return \code{data.table} with calculated linear associations
#' 
#' @examples
#' X <- matrix(rep(1:13, length.out = 42), nrow = 6, 
#'             dimnames = list(sprintf("row_%s", 1:6), sprintf("feat_%s", 1:7)))
#' Y <- matrix(c(10:15, 110:115, 210:215), ncol = 3,
#'             dimnames = list(sprintf("row_%s", 1:6), sprintf("met_%s", 1:3)))
#' tab_assoc <- calc_assoc(X, Y)
#' 
#' @seealso \code{\link[cdsrmodels:lin_associations]{cdsrmodels::lin_associations}}
#' 
#' @keywords internal
#' 
#' @export 
calc_assoc <- function(X, Y) {
  checkmate::assert_matrix(X, mode = "numeric")
  checkmate::assert_names(rownames(X))
  checkmate::assert_multi_class(Y, c("matrix", "numeric", "integer"))
  if (is.matrix(Y)) checkmate::assert_names(rownames(Y))
  if (is.vector(Y)) checkmate::assert_names(names(Y))
  
  stopifnot("The X and Y dimensions must match." = NROW(X) == NROW(Y))
  
  # prevent error with lack of "dep.var"
  if (is.matrix(Y) && is.null(colnames(Y))) colnames(Y) <- sprintf("var_%s", seq_len(NCOL(Y)))
  
  # when Y has no variance
  if (is.vector(Y) && stats::sd(Y, na.rm = TRUE) == 0) {
    warning("Y has no variance.
            Rendering all associations void. Please double check this is correct.")
    dt_na <- data.table::data.table(feature = colnames(X),
                                    est_beta = NA_real_,
                                    est_beta_se = NA_real_,
                                    posterior_mean = NA_real_,
                                    posterior_sd = NA_real_,
                                    prob_negative = NA_real_,
                                    prob_positive = NA_real_,
                                    rho = NA_real_,
                                    p_value = NA_real_,
                                    q_value = NA_real_,
                                    s_value = NA_real_,
                                    lfsr = NA_real_,
                                    lfdr = NA_real_)
    return(dt_na)
  } else if (is.matrix(Y) && any(apply(Y, 2, stats::sd, na.rm = TRUE) == 0, na.rm = TRUE)) {
    col_no_var <- which(apply(Y, 2, stats::sd, na.rm = TRUE) == 0)
    
    warning(sprintf(
        "The following columns in Y have no variance: %s.
        Rendering associations void. Please double check this is correct.",
        paste(names(col_no_var), collapse = ", ")
      ))
    
    dt_na <- data.table::data.table(expand.grid(feature = colnames(X), 
                                                response = names(col_no_var),
                                                stringsAsFactors = FALSE),
                                    est_beta = NA_real_,
                                    est_beta_se = NA_real_,
                                    posterior_mean = NA_real_,
                                    posterior_sd = NA_real_,
                                    prob_negative = NA_real_,
                                    prob_positive = NA_real_,
                                    rho = NA_real_,
                                    p_value = NA_real_,
                                    q_value = NA_real_,
                                    s_value = NA_real_,
                                    lfsr = NA_real_,
                                    lfdr = NA_real_)
    
    if (NROW(col_no_var) == NCOL(Y)) {
      return(dt_na)
    } else {
      Y <- Y[, -col_no_var, drop = FALSE]
    }
  } else {
    dt_na <- NULL
  }
  
  # use `cdsr_models` to calculate the linear model coefficients efficiently on large matrices
  res <- cdsrmodels::lin_associations(X = X, Y = Y)
  
  # convert results from a `matrix` to a `data.table`
  dt_res <- data.table::as.data.table(res$res.table)
  
  # fill lacking name in dt_res$ind.va
  if (!all(dt_res$ind.var %in% rownames(res$p.val))) {
    # finite values of res$p.val are used as the basis for the final result
    dt_pval <- data.table::as.data.table(stats::na.omit(res$p.val), keep.rownames = "ind.var")
    dt_res <-
      merge(dt_res[, -c("ind.var"), with = FALSE], dt_pval,
            by.x = "p.val", by.y = names(dt_pval)[names(dt_pval) != "ind.var"])
  }
  
  # add information about `rho`
  dt_rho <- data.table::melt(data.table::as.data.table(res$rho, keep.rownames = "ind.var"),
                             id.vars = "ind.var", variable.name = "dep.var", value.name = "rho")
  dt_res <- if (is.vector(Y)) {
    merge(dt_res, dt_rho[, -c("dep.var"), with = FALSE], by = "ind.var")
  } else {
    merge(dt_res, dt_rho, by = c("ind.var", "dep.var"))
  }
  
  # re-order columns for a more human-friendly output
  dt_res <- .order_assoc_result(dt_res) 
  data.table::setkey(dt_res, NULL)
  
  # add in NA values from any columns that had 0 variance, if necessary
  if (NROW(dt_na)) {
    dt_res <- rbind(dt_res, dt_na)
    data.table::setorder(dt_res, "feature", na.last = TRUE)
  }
  dt_res
}

#' Order and rename columns in associations results
#' 
#' @param res_dt \code{data.table} of associations results
#' 
#' @return \code{data.table} of associations results with reorder and rename columns
#'
#' @keywords internal
.order_assoc_result <- function(res_dt) {
  checkmate::assert_data_table(res_dt)
  
  ls_col <- c("ind.var", "dep.var", "betahat", "sebetahat", "PosteriorMean",
              "PosteriorSD",  "NegativeProb", "PositiveProb", "rho", 
              "p.val", "qvalue", "svalue", "lfsr", "lfdr")
  ls_col <- intersect(ls_col, names(res_dt)) # removre"dep.var" for vector
  res_dt <- res_dt[, (ls_col), with = FALSE]
  # rename
  data.table::setnames(res_dt, 
                       old = c("ind.var", "dep.var", "betahat", "sebetahat",
                               "PosteriorMean", "PosteriorSD",  "NegativeProb", "PositiveProb",
                               "p.val", "qvalue", "svalue"),
                       new = c("feature", "response", "est_beta", "est_beta_se",
                               "posterior_mean", "posterior_sd", "prob_negative", "prob_positive",
                               "p_value", "q_value", "s_value"),
                       skip_absent = TRUE)
  res_dt
}
