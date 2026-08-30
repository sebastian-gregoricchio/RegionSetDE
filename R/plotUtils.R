#' @title .thinIndex
#'
#' @description Returns the positions of a regularly spaced subset of a vector, used to keep the point clouds drawable. The thinning is deterministic, so the same object always gives the same picture.
#'
#' @param n Number of available elements.
#' @param maxPoints Maximum number of elements to keep.
#'
#' @return An integer vector of positions.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.thinIndex <-
  function(n,
           maxPoints) {

    if (n <= maxPoints) {
      return(seq_len(n))
    }

    return(unique(round(seq(from = 1, to = n, length.out = maxPoints))))
  } # END function




#' @title .regionSetTheme
#'
#' @description Theme shared by the plotting functions of the package. It follows the look of a publication panel: white background, no grid, black axis lines, no box around the facet labels. Titles and axis titles are rendered as markdown, so that subscripts such as \code{log<sub>2</sub>} come out formatted.
#'
#' @param legendPosition String with the position of the legend, one among \code{"right"}, \code{"left"}, \code{"top"}, \code{"bottom"} or \code{"none"}. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot2} theme.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 theme_bw theme element_blank element_line element_text margin
#' @importFrom ggtext element_markdown
#'
#' @keywords internal

.regionSetTheme <-
  function(legendPosition = "right",
           baseSize = 12) {

    ggplot2::theme_bw(base_size = baseSize) +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
                     axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
                     axis.text.x = ggplot2::element_text(colour = "black"),
                     axis.text.y = ggplot2::element_text(colour = "black"),
                     axis.title.x = ggtext::element_markdown(colour = "black"),
                     axis.title.y = ggtext::element_markdown(colour = "black"),
                     plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5),
                     plot.subtitle = ggtext::element_markdown(face = "bold", hjust = 0.5),
                     plot.caption = ggplot2::element_text(colour = "gray30"),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold"),
                     legend.key = ggplot2::element_blank(),
                     legend.position = legendPosition,
                     plot.margin = ggplot2::margin(t = 6, r = 6, b = 6, l = 6))
  } # END function
