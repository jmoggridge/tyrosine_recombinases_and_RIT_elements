#' theme_bot()
#' ggplot(iris, aes(Species, Sepal.Width)) + geom_boxplot() + theme_roboto()

theme_bot <- function(base_size = 11,
                         strip_text_size = 12,
                         strip_text_margin = 5,
                         subtitle_size = 13,
                         subtitle_margin = 10,
                         plot_title_size = 16,
                         plot_title_margin = 10,
                         ...) {
  ret <- ggplot2::theme_minimal(
    base_family = "RobotoCondensed-Regular",
    base_size = base_size, 
    ...
    )
  ret$strip.text <- ggplot2::element_text(
    hjust = 0,
    size = strip_text_size,
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Roboto-Bold"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, 
    size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "RobotoCondensed-Regular"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, 
    size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Roboto-Bold"
  )
  ret
}


