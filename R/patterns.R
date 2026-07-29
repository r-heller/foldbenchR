#' Construct a Miura-ori crease pattern
#'
#' Builds the vertex grid and crease assignment of a Miura-ori tessellation in
#' its flat, unfolded state. The returned flat coordinates are the exact
#' two-dimensional embedding the package exists to provide: folding is an
#' isometry, so these coordinates remain the ground truth for the folded
#' manifold at every fold angle.
#'
#' @param n_rows,n_cols Number of cells along each axis. Both must be >= 1.
#' @param a,b Cell side lengths. Both must be > 0.
#' @param alpha Acute angle of the parallelogram facet, in radians, strictly
#'   between 0 and `pi / 2`.
#'
#' @return An object of class `crease_pattern`: a list with `vertices` (an
#'   `n x 2` matrix of flat coordinates), `creases` (a data frame of vertex
#'   index pairs with an `assignment` column, `"M"` or `"V"`), `grid` (the
#'   row/column index of each vertex) and `params`.
#'
#' @family pattern constructors
#' @examples
#' p <- miura_ori(3, 4)
#' dim(p$vertices)
#' table(p$creases$assignment)
#' @export
miura_ori <- function(n_rows = 4, n_cols = 4, a = 1, b = 1, alpha = pi / 3) {
  .check_count(n_rows, "n_rows")
  .check_count(n_cols, "n_cols")
  .check_positive(a, "a")
  .check_positive(b, "b")
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= pi / 2) {
    stop("`alpha` must be a single number strictly between 0 and pi/2.",
         call. = FALSE)
  }

  nr <- n_rows + 1L
  nc <- n_cols + 1L

  # Miura-ori: every second row is sheared along x by the zig-zag offset.
  # The shear is what makes the pattern rigid-foldable with one degree of
  # freedom, and it is why the mountain/valley assignment alternates by row.
  offset <- b / tan(alpha)
  grid <- expand.grid(row = seq_len(nr), col = seq_len(nc))
  x <- (grid$col - 1) * a + ifelse(grid$row %% 2 == 0, offset, 0)
  y <- (grid$row - 1) * b
  vertices <- cbind(x = x, y = y)

  idx <- function(r, c) (c - 1L) * nr + r

  creases <- list()
  for (r in seq_len(nr)) {
    for (c in seq_len(nc - 1L)) {
      creases[[length(creases) + 1L]] <- data.frame(
        from = idx(r, c), to = idx(r, c + 1L),
        assignment = if (r %% 2 == 1) "M" else "V",
        stringsAsFactors = FALSE
      )
    }
  }
  for (c in seq_len(nc)) {
    for (r in seq_len(nr - 1L)) {
      creases[[length(creases) + 1L]] <- data.frame(
        from = idx(r, c), to = idx(r + 1L, c),
        assignment = if (c %% 2 == 1) "V" else "M",
        stringsAsFactors = FALSE
      )
    }
  }

  structure(
    list(
      vertices = vertices,
      creases  = do.call(rbind, creases),
      grid     = grid,
      params   = list(n_rows = n_rows, n_cols = n_cols,
                      a = a, b = b, alpha = alpha),
      family   = "miura"
    ),
    class = "crease_pattern"
  )
}


#' Construct a Yoshimura crease pattern
#'
#' A diamond (triangulated) tessellation. Facet shape and vertex degree differ
#' from [miura_ori()], so the same evaluation machinery meets a genuinely
#' different geometry.
#'
#' @inheritParams miura_ori
#' @return A `crease_pattern`, as described in [miura_ori()].
#' @family pattern constructors
#' @examples
#' p <- yoshimura(3, 4)
#' nrow(p$vertices)
#' @export
yoshimura <- function(n_rows = 4, n_cols = 4, a = 1, b = 1) {
  .check_count(n_rows, "n_rows")
  .check_count(n_cols, "n_cols")
  .check_positive(a, "a")
  .check_positive(b, "b")

  nr <- n_rows + 1L
  nc <- n_cols + 1L
  grid <- expand.grid(row = seq_len(nr), col = seq_len(nc))
  # Alternate rows offset by half a cell — the diamond stagger.
  x <- (grid$col - 1) * a + ifelse(grid$row %% 2 == 0, a / 2, 0)
  y <- (grid$row - 1) * b
  vertices <- cbind(x = x, y = y)

  idx <- function(r, c) (c - 1L) * nr + r
  creases <- list()
  for (r in seq_len(nr)) {
    for (c in seq_len(nc - 1L)) {
      creases[[length(creases) + 1L]] <- data.frame(
        from = idx(r, c), to = idx(r, c + 1L),
        assignment = if (r %% 2 == 1) "M" else "V",
        stringsAsFactors = FALSE
      )
    }
  }
  for (c in seq_len(nc)) {
    for (r in seq_len(nr - 1L)) {
      creases[[length(creases) + 1L]] <- data.frame(
        from = idx(r, c), to = idx(r + 1L, c),
        assignment = "V", stringsAsFactors = FALSE
      )
    }
  }

  structure(
    list(vertices = vertices, creases = do.call(rbind, creases), grid = grid,
         params = list(n_rows = n_rows, n_cols = n_cols, a = a, b = b),
         family = "yoshimura"),
    class = "crease_pattern"
  )
}


#' Construct a waterbomb crease pattern
#'
#' @inheritParams miura_ori
#' @return A `crease_pattern`, as described in [miura_ori()].
#' @family pattern constructors
#' @examples
#' p <- waterbomb(2, 2)
#' p$family
#' @export
waterbomb <- function(n_rows = 4, n_cols = 4, a = 1, b = 1) {
  p <- yoshimura(n_rows, n_cols, a, b)
  # Waterbomb differs from Yoshimura in crease assignment, not in the vertex
  # lattice: the horizontal creases alternate within a row rather than by row.
  p$creases$assignment <- ifelse(seq_len(nrow(p$creases)) %% 2 == 0, "M", "V")
  p$family <- "waterbomb"
  p$params <- list(n_rows = n_rows, n_cols = n_cols, a = a, b = b)
  p
}


#' @export
print.crease_pattern <- function(x, ...) {
  cat("<crease_pattern>", x$family, "\n")
  cat("  vertices:", nrow(x$vertices), "\n")
  cat("  creases: ", nrow(x$creases),
      sprintf("(M: %d, V: %d)\n", sum(x$creases$assignment == "M"),
              sum(x$creases$assignment == "V")))
  invisible(x)
}
