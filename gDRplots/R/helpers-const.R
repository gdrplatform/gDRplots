#' Get names of metrics that should be presented on a log scale
#' 
#' @examples
#' get_metrics_to_transform()
#' 
#' @keywords utils
#' @return character vector of standard metric names
#' @export
get_metrics_to_transform <- function() {
  gDRutils::get_settings_from_json("METRICS_TRANSFORMED",
                                   system.file(package = "gDRplots", "settings.json"))
}
