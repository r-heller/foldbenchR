`%||%` <- function(a, b) if (is.null(a)) b else a

# Set `seed` for the rest of the calling function and restore the caller's
# random-number state when that function exits.
.local_seed <- function(seed, envir = parent.frame()) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  restore <- function() {
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed)
  invisible(NULL)
}

.check_pattern <- function(x) {
  if (!inherits(x, "crease_pattern")) {
    stop("`pattern` must be a crease_pattern, as returned by miura_ori(), ",
         "yoshimura() or waterbomb().", call. = FALSE)
  }
  invisible(TRUE)
}

.check_count <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be a finite positive integer.", nm), call. = FALSE)
  }
  invisible(TRUE)
}

.check_positive <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive number.", nm), call. = FALSE)
  }
  invisible(TRUE)
}

.as_matrix <- function(x, nm) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || !all(is.finite(x)) ||
      nrow(x) < 2L || ncol(x) < 1L) {
    stop(sprintf("`%s` must be a finite numeric matrix with at least two rows and one column.", nm), call. = FALSE)
  }
  x
}

# Indices of the k nearest neighbours of each row, self excluded.
.knn <- function(X, k) {
  X <- .as_matrix(X, "X")
  .check_count(k, "k")
  if (k >= nrow(X)) stop("`k` must be smaller than the number of points.", call. = FALSE)
  d <- as.matrix(stats::dist(X))
  diag(d) <- Inf
  neighbours <- vapply(seq_len(nrow(d)), function(i) order(d[i, ])[seq_len(k)], integer(k))
  t(matrix(neighbours, nrow = k, ncol = nrow(X)))
}

# Shared engine for trustworthiness and continuity: both compare rank in one
# space against neighbourhood membership in the other, differing only in which
# space supplies which.
.rank_metric <- function(X, Y, k, direction = c("trust", "cont")) {
  direction <- match.arg(direction)
  X <- .as_matrix(X, "X"); Y <- .as_matrix(Y, "Y")
  n <- nrow(X)
  if (nrow(Y) != n) stop("`X` and `Y` must have the same number of rows.",
                         call. = FALSE)
  .check_count(k, "k")
  if (k >= n / 2) stop("`k` must be smaller than half the number of points for rank metrics.",
                   call. = FALSE)

  from <- if (direction == "trust") Y else X
  rank_in <- if (direction == "trust") X else Y

  dr <- as.matrix(stats::dist(rank_in)); diag(dr) <- Inf
  ranks <- t(apply(dr, 1, function(r) rank(r, ties.method = "first")))
  nb <- .knn(from, k)

  penalty <- 0
  for (i in seq_len(n)) {
    r <- ranks[i, nb[i, ]]
    penalty <- penalty + sum(pmax(0, r - k))
  }
  1 - (2 / (n * k * (2 * n - 3 * k - 1))) * penalty
}

.check_grid <- function(n_rows, n_cols, a, b) {
  .check_count(n_rows, "n_rows"); .check_count(n_cols, "n_cols")
  .check_positive(a, "a"); .check_positive(b, "b")
  if ((n_rows + 1) * (n_cols + 1) > 1e6) {
    stop("Pattern exceeds the supported limit of one million vertices.", call. = FALSE)
  }
}

.check_theta <- function(theta) {
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
      theta < 0 || theta > 1) {
    stop("`theta` is a normalised folding parameter in [0, 1], not an angle (see migration guide).", call. = FALSE)
  }
}

.check_sample <- function(x) {
  if (!is.list(x)) stop("A manifold must contain X and truth matrices.", call. = FALSE)
  x$X <- .as_matrix(x$X, "X"); x$truth <- .as_matrix(x$truth, "truth")
  if (nrow(x$X) != nrow(x$truth)) stop("X and truth must have matching rows.", call. = FALSE)
  x
}
