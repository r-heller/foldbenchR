#' Recover flat material coordinates
#'
#' Returns the coordinates before folding. Local surface lengths and angles
#' within each facet are preserved by the supported rigid map. Euclidean
#' distances between flat points are a benchmark target; they need not be
#' shortest paths constrained to a non-convex finite sheet.
#' @param pattern A `crease_pattern` or `fold_manifold`.
#' @return An `n x 2` matrix of flat material coordinates.
#' @family folding
#' @examples
#' unfold(miura_ori(2, 2))
#' @export
unfold <- function(pattern) {
  if (inherits(pattern, "fold_manifold")) return(pattern$truth)
  .check_pattern(pattern)
  pattern$vertices
}

#' Rigidly fold a Miura sheet
#'
#' Uses the closed-form Miura placement documented in `GEOMETRY.md`. All
#' pairwise distances within each facet are preserved. This is an ideal
#' zero-thickness kinematic model, without collision/contact mechanics.
#' @param pattern A Miura `crease_pattern`.
#' @param theta Normalised folding parameter in `[0, 1]`: zero is flat and one
#'   is fully folded, with overlapping facets. It is not a dihedral angle.
#'   For manifold benchmarks use `theta < 1`; the endpoint is not injective.
#' @return A `fold_manifold` with `coords`, `truth`, `pattern`, `theta` and
#'   `corrugation_angle` (radians, equal to `theta * alpha`).
#' @family folding
#' @examples
#' m <- fold(miura_ori(3, 3), theta = 0.5)
#' dim(m$coords)
#' @export
fold <- function(pattern, theta = 0) {
  .check_pattern(pattern); .check_theta(theta)
  if (pattern$family != "miura") {
    stop("No certified folding map for this family; only miura_ori() supports folding and sampling.", call. = FALSE)
  }
  p <- pattern$params
  phi <- theta * p$alpha
  i <- pattern$grid$col - 1L; j <- pattern$grid$row - 1L
  # The trigonometric product avoids cancellation at small alpha.
  Lx <- p$a * cos(phi); H <- p$a * sin(phi)
  dx <- p$b * cos(p$alpha) / cos(phi)
  Ly <- p$b * sqrt(max(0, sin(p$alpha + phi) * sin(p$alpha - phi))) / cos(phi)
  coords <- cbind(x = i * Lx + (j %% 2L) * dx, y = j * Ly, z = (i %% 2L) * H)
  if (theta == 0) coords <- cbind(pattern$vertices, z = 0)
  structure(list(coords = coords, truth = pattern$vertices, pattern = pattern,
                 theta = theta, corrugation_angle = phi), class = "fold_manifold")
}

# Map material points by their cell coordinates, preserving the affine facet map.
.map_flat <- function(pattern, truth, theta) {
  p <- pattern$params
  h <- p$b * sin(p$alpha); offset <- p$b * cos(p$alpha)
  ycell <- truth[, 2] / h
  row <- pmin(pmax(floor(ycell), 0), p$n_rows - 1L)
  v <- ycell - row
  shift <- ifelse(row %% 2L == 0L, v, 1 - v) * offset
  xcell <- (truth[, 1] - shift) / p$a
  col <- pmin(pmax(floor(xcell), 0), p$n_cols - 1L)
  u <- xcell - col
  idx <- col * (p$n_rows + 1L) + row + 1L
  vertices <- fold(pattern, theta)$coords
  vertices[idx, , drop = FALSE] +
    (vertices[idx + p$n_rows + 1L, , drop = FALSE] - vertices[idx, , drop = FALSE]) * u +
    (vertices[idx + 1L, , drop = FALSE] - vertices[idx, , drop = FALSE]) * v
}

.reflect_interval <- function(x, upper) {
  z <- x %% (2 * upper)
  ifelse(z <= upper, z, 2 * upper - z)
}

#' Sample a finite folded Miura sheet
#'
#' Selects equal-area facets uniformly and samples uniformly within each facet.
#' Coordinates are mapped through that same rigid facet map. A supplied seed
#' sets R's random-number state, as in `set.seed()`.
#' @param pattern A Miura `crease_pattern`.
#' @param theta Normalised folding parameter; see [fold()].
#' @param n Positive integer number of points.
#' @param noise List with `type` (`"none"`, `"ambient"`, `"manifold"`,
#'   `"outlier"`) and non-negative `sd`. Ambient Gaussian noise is in 3D.
#'   Manifold noise perturbs material coordinates and reflects them at the
#'   finite sheet boundaries, then folds; it is not unbounded Gaussian noise.
#'   Outlier noise displaces 2 percent of points (at least one) at `10 * sd`.
#'   For ambient/outlier noise, `truth` records the underlying clean material
#'   coordinates, not an exact inverse of the noisy observation.
#' @param seed Optional integer seed.
#' @return A list with `X` (`n x 3`) and `truth` (`n x 2`).
#' @family folding
#' @examples
#' sample_manifold(miura_ori(3, 3), theta = 0.5, n = 50, seed = 1)
#' @export
sample_manifold <- function(pattern, theta = 0, n = 500,
                            noise = list(type = "none", sd = 0), seed = NULL) {
  fold(pattern, theta)
  .check_count(n, "n")
  if (!is.list(noise)) stop("`noise` must be a list.", call. = FALSE)
  type <- match.arg(noise$type %||% "none", c("none", "ambient", "manifold", "outlier"))
  sd <- noise$sd %||% 0
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) {
    stop("noise `sd` must be a finite non-negative number.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  p <- pattern$params; flat <- pattern$vertices
  rows <- sample.int(p$n_rows, n, replace = TRUE)
  cols <- sample.int(p$n_cols, n, replace = TRUE)
  idx <- (cols - 1L) * (p$n_rows + 1L) + rows
  u <- flat[idx, , drop = FALSE] +
    (flat[idx + 1L, , drop = FALSE] - flat[idx, , drop = FALSE]) * stats::runif(n) +
    (flat[idx + p$n_rows + 1L, , drop = FALSE] - flat[idx, , drop = FALSE]) * stats::runif(n)
  if (type == "manifold" && sd > 0) {
    u <- u + matrix(stats::rnorm(2L * n, 0, sd), nrow = n)
    h <- p$b * sin(p$alpha)
    u[, 2] <- .reflect_interval(u[, 2], p$n_rows * h)
    ycell <- u[, 2] / h
    row <- pmin(floor(ycell), p$n_rows - 1L); v <- ycell - row
    shift <- ifelse(row %% 2L == 0L, v, 1 - v) * p$b * cos(p$alpha)
    u[, 1] <- shift + .reflect_interval(u[, 1] - shift, p$n_cols * p$a)
  }
  X <- .map_flat(pattern, u, theta)
  if (type == "ambient" && sd > 0) X <- X + matrix(stats::rnorm(length(X), 0, sd), nrow = n)
  if (type == "outlier" && sd > 0) {
    k <- max(1L, round(0.02 * n)); hit <- sample.int(n, k)
    X[hit, ] <- X[hit, ] + matrix(stats::rnorm(k * 3, 0, sd * 10), nrow = k)
  }
  list(X = X, truth = u)
}

#' Product of two folded manifolds
#'
#' Stacks the ambient coordinates of `k` folded sheets side by side. The result
#' has intrinsic dimension `2k` in ambient dimension `3k`, and its exact
#' unfolding is the product of the individual flat patterns — so ground truth
#' survives the construction. This is what lifts the benchmark off intrinsic
#' dimension 2.
#'
#' @param ... Two or more lists as returned by [sample_manifold()], all with
#'   the same number of rows.
#' @return A list with `X` and `truth`.
#' @family folding
#' @examples
#' a <- sample_manifold(miura_ori(2, 2), pi / 8, 40, seed = 1)
#' b <- sample_manifold(miura_ori(2, 2, alpha = pi / 4), pi / 8, 40, seed = 2)
#' dim(fold_product(a, b)$X)
#' @export
fold_product <- function(...) {
  parts <- lapply(list(...), .check_sample)
  if (length(parts) < 2L) stop("`fold_product()` needs at least two manifolds.",
                               call. = FALSE)
  n <- vapply(parts, function(p) nrow(p$X), integer(1))
  if (length(unique(n)) != 1L) {
    stop("All manifolds must have the same number of points.", call. = FALSE)
  }
  list(X = do.call(cbind, lapply(parts, `[[`, "X")),
       truth = do.call(cbind, lapply(parts, `[[`, "truth")))
}


#' Lift a manifold into a higher ambient dimension
#'
#' Applies a random orthonormal linear map, so ambient dimension can be set
#' arbitrarily while intrinsic geometry — and therefore the ground truth — is
#' untouched.
#'
#' @param manifold A list with `X` and `truth`, as returned by
#'   [sample_manifold()] or [fold_product()].
#' @param ambient_dim Target dimension; must exceed `ncol(manifold$X)`.
#' @param seed Optional integer seed.
#' @return A list with the lifted `X` and the unchanged `truth`.
#' @family folding
#' @examples
#' m <- sample_manifold(miura_ori(2, 2), pi / 8, 30, seed = 1)
#' dim(lift(m, 10, seed = 1)$X)
#' @export
lift <- function(manifold, ambient_dim, seed = NULL) {
  manifold <- .check_sample(manifold)
  .check_count(ambient_dim, "ambient_dim")
  d <- ncol(manifold$X)
  if (!is.numeric(ambient_dim) || length(ambient_dim) != 1L ||
      ambient_dim <= d) {
    stop("`ambient_dim` must be a single number greater than ncol(X).",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  A <- matrix(stats::rnorm(d * ambient_dim), nrow = d)
  A <- qr.Q(qr(t(A)))[, seq_len(d), drop = FALSE]   # orthonormal columns
  manifold$X <- manifold$X %*% t(A)
  manifold
}
