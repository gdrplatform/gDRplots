.settings_cache <- new.env(parent = emptyenv())

#' @keywords internal
.get_setting <- function(key) {
  if (is.null(.settings_cache[[key]])) {
    .settings_cache[[key]] <- gDRutils::get_settings_from_json(
      key, system.file(package = "gDRplots", "settings.json")
    )
  }
  .settings_cache[[key]]
}
