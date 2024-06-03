#' Define tools in bar mode for plotly plot
#'
#' see: https://github.com/plotly/plotly.js/blob/master/src/plot_api/plot_config.js
#'
#' @param plt plotly object
#' @param editable logical determines whether the graph is editable or not
#' @param ... other parameters that may be passed to \code{\link[plotly]{config}}
#' 
#' @examples
#' plt <- plotly::plot_ly(mtcars, x = ~mpg, type = "histogram")
#' plt
#' gDR_plotly_config(plt) 
#'
#' @return plotly plot with gDR config
#' @keywords plotly_config
#' @export
gDR_plotly_config <- function(plt, 
                              editable = FALSE, 
                              ...) {
  checkmate::assert_class(plt, "plotly")
  checkmate::assert_flag(editable)
  
  plotly::config(
    plt,
    # options
    displaylogo = FALSE, # remove the plotly logo
    showSendToCloud = FALSE, # do not include the send data to cloud button (floppy icon) (dflt)
    showEditInChartStudio = FALSE, # do not display `Edit in Chart Studio` button (pencil icon) (dflt)
    showLink = FALSE, sendData = FALSE, # do not show "EditChart" link (dflt)
    # layout
    responsive = TRUE, # whether to change the layout size when window is resized (dflt)
    # modebar
    modeBarButtonsToRemove = list("toImage"),
    modeBarButtons = list(c("zoom2d", "resetViews", "autoScale2d",
                            "toggleSpikelines", "toggleHover", "hoverCompareCartesian")),
    # whether the graph is editable
    editable = editable,
    ...
  )
  
}


#' Get list of items that are editable in plotly plots
#' 
#' @examples
#' get_plotly_edits()
#' 
#' @keywords plotly_config
#' @return list of character vectors with editable items in plotly plots
#' @export
get_plotly_edits <- function() {
  gDRutils::get_settings_from_json("PLOTLY_EDITS",
                                   system.file(package = "gDRplots", "settings.json"))
}
