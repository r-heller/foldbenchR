#' Plot a crease pattern
#'
#' Mountain and valley creases are distinguished by linetype **and** colour, so
#' the figure survives greyscale printing and colour vision deficiency.
#'
#' @param pattern A `crease_pattern`.
#' @return A `ggplot` object.
#' @family plotting
#' @examples
#' plot_crease(miura_ori(3, 3))
#' @export
plot_crease <- function(pattern) {
  .check_pattern(pattern)
  v <- pattern$vertices
  seg <- data.frame(
    x    = v[pattern$creases$from, "x"],
    y    = v[pattern$creases$from, "y"],
    xend = v[pattern$creases$to, "x"],
    yend = v[pattern$creases$to, "y"],
    assignment = pattern$creases$assignment
  )
  ggplot2::ggplot(seg) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                   yend = .data$yend, colour = .data$assignment,
                   linetype = .data$assignment)) +
    ggplot2::scale_colour_manual(values = c(M = "#0d0887", V = "#cc4778")) +
    ggplot2::scale_linetype_manual(values = c(M = "solid", V = "dashed")) +
    ggplot2::coord_equal() +
    ggplot2::labs(colour = "Crease", linetype = "Crease", x = NULL, y = NULL) +
    ggplot2::theme_minimal()
}

#' Plot a folded manifold
#'
#' @param object A `fold_manifold` from [fold()].
#' @param ... Ignored.
#' @return A `ggplot` object showing the x-z profile of the fold.
#' @family plotting
#' @examples
#' autoplot.fold_manifold(fold(miura_ori(3, 3), pi / 6))
#' @export
autoplot.fold_manifold <- function(object, ...) {
  d <- as.data.frame(object$coords)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$y, y = .data$z)) +
    ggplot2::geom_point(size = 0.8, colour = "#0d0887") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "y", y = "z",
                  subtitle = sprintf("theta = %.3f", object$theta)) +
    ggplot2::theme_minimal()
}
