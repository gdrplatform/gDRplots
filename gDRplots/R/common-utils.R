#' create a log-sequence
#'
#' Create a sequence of numbers growing in log-domain.
#'
#' The result is a numeric vector of length \code{length}.
#' Differences between items are constant in logarithmic domain
#' and therefore geometrically increase in linear domain.
#'
#' @param start,end numeric, lower and upper margins of the sequence
#' @param length integer, resulting sequence length
#' 
#' @keywords utils
#' @return A numeric vector, see \code{Details}.
#'
#' @export 
create_log_seq <- function(start, end, length) {

  checkmate::assert_number(start, lower = 0, finite = TRUE)
  checkmate::assert_number(end, lower = 0, finite = TRUE)
  checkmate::assert_number(length, lower = 1, finite = TRUE)
  
  limits <- c(start, end)
  limits_log <- log10(limits)
  sequence_log <- seq(from = limits_log[1], limits_log[2], length.out = length)
  sequence <- 10 ^ sequence_log
  return(sequence)
}                     

