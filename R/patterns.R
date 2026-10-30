#' Construct a Miura-ori crease pattern
#'
#' A finite parallelogram sheet with free boundaries. Both `a` and `b` are
#' facet side lengths. Facets are wound counter-clockwise in the flat plane.
#' Mountain/valley labels describe the positive folding branch used by [fold()].
#' Boundary edges carry `"B"` and do not fold.
#'
#' @param n_rows,n_cols Positive integer numbers of cells.
#' @param a,b Positive facet side lengths.
#' @param alpha Acute facet angle, in radians, strictly between 0 and `pi/2`.
#' @return A `crease_pattern` with `vertices` (two columns), `facets` (ordered
#'   vertex indices), `creases` (`from`, `to`, `assignment`), `grid` and `params`.
#' @family pattern constructors
#' @examples
#' p <- miura_ori(3, 4)
#' table(p$creases$assignment)
#' @export
miura_ori <- function(n_rows = 4, n_cols = 4, a = 1, b = 1, alpha = pi / 3) {
  .check_grid(n_rows, n_cols, a, b)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= pi / 2) {
    stop("`alpha` must be a finite number strictly between 0 and pi/2.", call. = FALSE)
  }
  grid <- expand.grid(row = seq_len(n_rows + 1L), col = seq_len(n_cols + 1L))
  i <- grid$col - 1L; j <- grid$row - 1L
  vertices <- cbind(x = i * a + (j %% 2L) * b * cos(alpha), y = j * b * sin(alpha))
  idx <- function(row, col) (col - 1L) * (n_rows + 1L) + row
  facets <- list(); creases <- list()
  add <- function(from, to, assignment) {
    creases[[length(creases) + 1L]] <<- data.frame(from, to, assignment)
  }
  for (row in seq_len(n_rows)) for (col in seq_len(n_cols)) {
    facets[[length(facets) + 1L]] <- c(idx(row, col), idx(row, col + 1L),
                                     idx(row + 1L, col + 1L), idx(row + 1L, col))
  }
  for (row in seq_len(n_rows + 1L)) for (col in seq_len(n_cols)) {
    label <- if (row %in% c(1L, n_rows + 1L)) "B" else
      if ((col + row) %% 2L == 0L) "M" else "V"
    add(idx(row, col), idx(row, col + 1L), label)
  }
  for (col in seq_len(n_cols + 1L)) for (row in seq_len(n_rows)) {
    label <- if (col %in% c(1L, n_cols + 1L)) "B" else
      if (col %% 2L == 0L) "M" else "V"
    add(idx(row, col), idx(row + 1L, col), label)
  }
  structure(list(vertices = vertices, facets = facets, creases = do.call(rbind, creases),
                 grid = grid, family = "miura",
                 params = list(n_rows = n_rows, n_cols = n_cols, a = a, b = b, alpha = alpha)),
            class = "crease_pattern")
}

#' Construct a flat Yoshimura diamond pattern
#'
#' Returns a triangular tessellation for inspecting flat geometry. Interior
#' edges are `"U"` (unassigned): no folding branch is certified here. [fold()]
#' and [sample_manifold()] reject this family. This is not a second benchmark
#' mechanism. `b` is the row height for this constructor.
#'
#' @inheritParams miura_ori
#' @return A `crease_pattern`, with triangular facets.
#' @family pattern constructors
#' @examples
#' plot_crease(yoshimura(3, 4))
#' @export
yoshimura <- function(n_rows = 4, n_cols = 4, a = 1, b = 1) {
  .check_grid(n_rows, n_cols, a, b)
  grid <- expand.grid(row = seq_len(n_rows + 1L), col = seq_len(n_cols + 1L))
  vertices <- cbind(x = (grid$col - 1L) * a + ((grid$row - 1L) %% 2L) * a / 2,
                    y = (grid$row - 1L) * b)
  idx <- function(row, col) (col - 1L) * (n_rows + 1L) + row
  facets <- list()
  for (row in seq_len(n_rows)) for (col in seq_len(n_cols)) {
    q <- c(idx(row, col), idx(row, col + 1L), idx(row + 1L, col + 1L), idx(row + 1L, col))
    triangles <- if (row %% 2L == 1L) list(q[c(1, 2, 4)], q[c(2, 3, 4)]) else
      list(q[c(1, 2, 3)], q[c(1, 3, 4)])
    facets <- c(facets, triangles)
  }
  edges <- do.call(rbind, lapply(facets, function(v) {
    w <- c(v[-1], v[1]); cbind(from = pmin(v, w), to = pmax(v, w))
  }))
  key <- paste(edges[, 1], edges[, 2], sep = ":")
  unique_edges <- !duplicated(key)
  creases <- as.data.frame(edges[unique_edges, , drop = FALSE])
  creases$assignment <- ifelse(table(key)[key[unique_edges]] == 1L, "B", "U")
  structure(list(vertices = vertices, facets = facets, creases = creases,
                 grid = grid, family = "yoshimura",
                 params = list(n_rows = n_rows, n_cols = n_cols, a = a, b = b)),
            class = "crease_pattern")
}

#' Waterbomb research placeholder
#'
#' The prototype returned a relabelled staggered grid, not a waterbomb
#' tessellation. This entry point now fails explicitly until both its topology
#' and kinematics have been established. No waterbomb benchmark is available.
#' @inheritParams miura_ori
#' @return No value; raises an informative error.
#' @family pattern constructors
#' @export
waterbomb <- function(n_rows = 4, n_cols = 4, a = 1, b = 1) {
  stop("waterbomb geometry and kinematics are not implemented; use miura_ori() for the supported benchmark.", call. = FALSE)
}

#' @export
print.crease_pattern <- function(x, ...) {
  cat("<crease_pattern>", x$family, "\n")
  cat("  vertices:", nrow(x$vertices), " facets:", length(x$facets), "\n")
  print(table(x$creases$assignment))
  invisible(x)
}
