`%||%` <- function(a, b) if (is.null(a)) b else a

.check_pattern <- function(x) {
  if (!inherits(x, "crease_pattern")) {
    stop("`pattern` must be a crease_pattern, as returned by miura_ori(), ",
         "yoshimura() or waterbomb().", call. = FALSE)
  }
  invisible(TRUE)
}

.check_count <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1) {
    stop(sprintf("`%s` must be a single number >= 1.", nm), call. = FALSE)
  }
  invisible(TRUE)
}

.check_positive <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive number.", nm), call. = FALSE)
  }
  invisible(TRUE)
}

.as_matrix <- function(x, nm) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix.", nm), call. = FALSE)
  }
  x
}

# Indices of the k nearest neighbours of each row, self excluded.
.knn <- function(X, k) {
  d <- as.matrix(stats::dist(X))
  diag(d) <- Inf
  t(apply(d, 1, function(r) order(r)[seq_len(k)]))
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
  if (k >= n) stop("`k` must be smaller than the number of points.",
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
