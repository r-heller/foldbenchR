#' Recover the exact flat embedding of a crease pattern
#'
#' Returns the pattern's flat coordinates. Because rigid folding is an
#' isometry, these coordinates are the true two-dimensional embedding of the
#' folded manifold at *every* fold angle — not an approximation of it. This is
#' the ground truth every metric in the package is measured against.
#'
#' @param pattern A `crease_pattern`, or a `fold_manifold` produced by [fold()].
#' @return An `n x 2` matrix of flat coordinates.
#' @family folding
#' @examples
#' p <- miura_ori(2, 2)
#' head(unfold(p))
#' @export
unfold <- function(pattern) {
  if (inherits(pattern, "fold_manifold")) return(pattern$truth)
  .check_pattern(pattern)
  pattern$vertices
}


#' Fold a crease pattern to a given fold angle
#'
#' Maps the flat pattern into `R^3` by rotating alternate rows out of plane by
#' the dihedral angle `theta`. Facets stay flat and lengths along the surface
#' are preserved, so the embedding is isometric and [unfold()] remains exact.
#'
#' @param pattern A `crease_pattern`.
#' @param theta Fold angle in radians, in `[0, pi/2)`. Zero returns the flat
#'   sheet lifted into `R^3`; larger values bring non-adjacent facets closer in
#'   ambient space while leaving geodesic distances unchanged.
#'
#' @return An object of class `fold_manifold`: a list with `coords` (`n x 3`),
#'   `truth` (`n x 2`, the exact unfolding), `pattern` and `theta`.
#' @family folding
#' @examples
#' p <- miura_ori(3, 3)
#' m <- fold(p, theta = pi / 6)
#' dim(m$coords)
#' @export
fold <- function(pattern, theta = 0) {
  .check_pattern(pattern)
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) ||
      theta < 0 || theta >= pi / 2) {
    stop("`theta` must be a single number in [0, pi/2).", call. = FALSE)
  }

  flat <- pattern$vertices
  row_id <- pattern$grid$row

  # Accordion fold along y. Each row band is rotated by +/- theta about the
  # crease line, alternating sign, so segment lengths in y are preserved:
  # a step of b in the flat sheet becomes b*cos(theta) in y and b*sin(theta)
  # in z. Arc length along the surface is therefore unchanged - the isometry
  # the whole package depends on.
  b <- if (!is.null(pattern$params$b)) pattern$params$b else 1
  steps <- row_id - 1L
  y <- steps * b * cos(theta)
  z <- ifelse(steps %% 2 == 0, 0, b * sin(theta))

  coords <- cbind(x = flat[, "x"], y = y, z = z)

  structure(
    list(coords = coords, truth = flat, pattern = pattern, theta = theta),
    class = "fold_manifold"
  )
}


#' Sample points on a folded manifold
#'
#' Samples uniformly on the facets of the folded surface and returns both the
#' ambient coordinates and the exact flat coordinates of the same points.
#'
#' @param pattern A `crease_pattern`.
#' @param theta Fold angle, see [fold()].
#' @param n Number of points.
#' @param noise A list with `type` (`"none"`, `"ambient"`, `"manifold"` or
#'   `"outlier"`) and `sd`. Ambient noise is added in `R^3` and moves points
#'   off the surface; on-manifold noise perturbs the flat coordinates before
#'   folding and therefore keeps the truth exact.
#' @param seed Optional integer seed.
#'
#' @return A list with `X` (`n x 3` ambient) and `truth` (`n x 2` exact).
#' @family folding
#' @examples
#' s <- sample_manifold(miura_ori(3, 3), theta = pi / 8, n = 50, seed = 1)
#' dim(s$X)
#' @export
sample_manifold <- function(pattern, theta = 0, n = 500,
                            noise = list(type = "none", sd = 0), seed = NULL) {
  .check_pattern(pattern)
  .check_count(n, "n")
  if (!is.null(seed)) set.seed(seed)
  type <- match.arg(noise$type %||% "none",
                    c("none", "ambient", "manifold", "outlier"))
  sd <- noise$sd %||% 0

  flat <- pattern$vertices
  rng_x <- range(flat[, "x"])
  rng_y <- range(flat[, "y"])
  u <- cbind(x = stats::runif(n, rng_x[1], rng_x[2]),
             y = stats::runif(n, rng_y[1], rng_y[2]))

  if (type == "manifold" && sd > 0) {
    u <- u + cbind(stats::rnorm(n, 0, sd), stats::rnorm(n, 0, sd))
  }

  b <- if (!is.null(pattern$params$b)) pattern$params$b else 1
  band <- floor(u[, "y"] / b)
  frac <- u[, "y"] - band * b
  y <- (band * b + frac) * cos(theta)
  z <- ifelse(band %% 2 == 0, frac, b - frac) * sin(theta)
  X <- cbind(x = u[, "x"], y = y, z = z)

  if (type == "ambient" && sd > 0) {
    X <- X + matrix(stats::rnorm(length(X), 0, sd), nrow = n)
  }
  if (type == "outlier" && sd > 0) {
    k <- max(1L, round(0.02 * n))
    hit <- sample.int(n, k)
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
#' b <- sample_manifold(yoshimura(2, 2), pi / 8, 40, seed = 2)
#' dim(fold_product(a, b)$X)
#' @export
fold_product <- function(...) {
  parts <- list(...)
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
