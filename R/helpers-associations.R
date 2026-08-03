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
        toString(names(col_no_var))
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

  res <- .lin_associations(X = X, Y = Y)

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


#' Compute pairwise linear associations between columns of X and Y
#'
#' Inline reimplementation of \code{lin_associations} from the
#' \href{https://github.com/cancerdatasci/cdsrmodels}{cdsrmodels} package
#' (MIT License, Broad Institute / DepMap Portal) to remove the GitHub-only
#' dependency.  The algorithm and output structure are unchanged; only the
#' package dependency is replaced.  Uses
#' \code{\link[WGCNA:cor]{WGCNA::cor}} for fast pairwise correlation on large
#' matrices and \code{\link[ashr:ash]{ashr::ash}} for empirical-Bayes
#' shrinkage of effect sizes.
#'
#' @param X \code{matrix} of independent variables (rows = samples, cols = features).
#' @param Y \code{matrix} or \code{vector} of response variables (rows = samples).
#' @param n.min integer; minimum number of finite paired observations required
#'   to compute a p-value (default 4).
#' @param shrinkage logical; apply \code{ashr} shrinkage (default \code{TRUE}).
#' @param alpha numeric; \code{ashr} alpha parameter (default 0).
#' @param MHC_direction character; \code{"x"} or \code{"y"} — direction of
#'   multiple-hypothesis correction.  Defaults to \code{"y"} when
#'   \code{ncol(Y) >= ncol(X)}, otherwise \code{"x"}.
#'
#' @return Named list with elements \code{N}, \code{rho}, \code{beta},
#'   \code{beta.se}, \code{p.val}, \code{q.val}, and \code{res.table}
#'   (a \code{data.frame} from \code{ashr}).
#'
#' @importFrom ashr ash
#' @importFrom dplyr bind_rows
#' @importFrom WGCNA cor
#' @importFrom stats pt p.adjust sd
#'
#' @keywords internal
.lin_associations <- function(X, Y, n.min = 4L, shrinkage = TRUE,
                              alpha = 0, MHC_direction = NULL) {
  if (is.null(MHC_direction)) {
    MHC_direction <- if (length(Y) >= length(X)) "x" else "y"
  }

  X.NA <- !is.finite(X)
  X[X.NA] <- NA
  Y.NA <- !is.finite(Y)
  Y[Y.NA] <- NA

  N <- (t(!X.NA) %*% (!Y.NA)) - 2L

  f_sd <- function(A) {
    if (is.matrix(A)) apply(A, 2L, stats::sd, na.rm = TRUE) else stats::sd(A, na.rm = TRUE)
  }
  sx <- f_sd(X)
  sy <- f_sd(Y)

  rho <- WGCNA::cor(X, Y, use = "pairwise.complete.obs")
  beta <- t(t(rho / sx) * sy)
  beta.se <- t(t(sqrt(1 - rho^2) / sx) * sy) / sqrt(N)
  p.val <- 2 * stats::pt(-abs(beta / beta.se), N)
  p.val[N < n.min] <- NA
  p.val[(sx == 0) | !is.finite(sx), ] <- NA
  p.val[, (sy == 0) | !is.finite(sy)] <- NA

  if (MHC_direction == "y") {
    q.val <- apply(p.val, 2L, function(x) stats::p.adjust(x, method = "BH"))
  } else {
    q.val <- t(apply(p.val, 1L, function(x) stats::p.adjust(x, method = "BH")))
  }

  if (shrinkage) {
    res.table <- vector("list", if (MHC_direction == "y") NCOL(p.val) else NROW(p.val))
    if (MHC_direction == "y") {
      for (ix in seq_len(NCOL(p.val))) {
        fin <- is.finite(p.val[, ix])
        feat_names <- rownames(p.val)[fin]
        res <- ashr::ash(beta[fin, ix], pmax(beta.se[fin, ix], 1e-10),
                         mixcompdist = "halfuniform", alpha = alpha)$result
        if (!is.null(res)) {
          res$dep.var <- colnames(Y)[ix]
          res$ind.var <- feat_names
          res$p.val <- p.val[fin, ix]
          res.table[[ix]] <- res
        }
      }
    } else {
      for (ix in seq_len(NROW(p.val))) {
        fin <- is.finite(p.val[ix, ])
        res <- tryCatch(
          ashr::ash(beta[ix, fin], pmax(beta.se[ix, fin], 1e-10),
                    mixcompdist = "halfuniform", alpha = alpha)$result,
          error = function(e) NULL
        )
        if (!is.null(res)) {
          res$ind.var <- rownames(p.val)[ix]
          y_col <- colnames(Y)[ix]
          if (!is.null(y_col) && !is.na(y_col)) res$dep.var <- y_col
          res$p.val <- p.val[ix, fin]
          res.table[[ix]] <- res
        }
      }
    }
    non_null <- Filter(Negate(is.null), res.table)
    if (length(non_null) == 0L) stop("Error: all input values are missing")
    res.table <- dplyr::bind_rows(non_null)
  } else {
    res.table <- NULL
  }

  list(N = N, rho = rho, beta = beta, beta.se = beta.se,
       p.val = p.val, q.val = q.val, res.table = res.table)
}
