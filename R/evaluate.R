#' Procrustes-aligned reconstruction error
#'
#' Root mean squared distance between an embedding and the known truth, after
#' removing the similarity transform (translation, rotation, reflection,
#' uniform scale) that no embedding method can be expected to recover. This is
#' the headline metric of the package: unlike trustworthiness or continuity it
#' is measured directly against the exact unfolding, not against a
#' neighbourhood proxy.
#'
#' @param embedding An `n x d` matrix produced by a dimension-reduction method.
#' @param truth An `n x d` matrix of true coordinates, from [unfold()].
#' @return A single non-negative number. Zero for an embedding that differs
#'   from the truth only by a similarity transform.
#' @family evaluation
#' @examples
#' p <- miura_ori(3, 3)
#' reconstruction_error(unfold(p), unfold(p))
#' @export
reconstruction_error <- function(embedding, truth) {
  embedding <- .as_matrix(embedding, "embedding")
  truth <- .as_matrix(truth, "truth")
  if (!identical(dim(embedding), dim(truth))) {
    stop("`embedding` and `truth` must have matching dimensions.",
         call. = FALSE)
  }
  if (sum(scale(embedding, scale = FALSE)^2) == 0 || sum(scale(truth, scale = FALSE)^2) == 0) {
    stop("Reconstruction requires non-constant embedding and truth.", call. = FALSE)
  }
  fit <- vegan::procrustes(truth, embedding, symmetric = FALSE)
  sqrt(mean(rowSums((fit$Yrot - fit$X)^2)))
}


#' Trustworthiness of an embedding
#'
#' Penalises points that are close in the embedding but far apart in the
#' original space. Field-standard, kNN-based, and known to be optimistic —
#' which is the point of measuring it here alongside
#' [reconstruction_error()].
#'
#' @param X Original coordinates, `n x p`.
#' @param Y Embedding, `n x d`.
#' @param k Positive integer neighbourhood size, less than half the number of
#'   points for trustworthiness and continuity, or less than the number of
#'   points for kNN preservation.
#' @return A number in `[0, 1]`; 1 is perfect.
#' @family evaluation
#' @examples
#' m <- sample_manifold(miura_ori(3, 3), pi / 8, 60, seed = 1)
#' trustworthiness(m$X, m$truth, k = 5)
#' @export
trustworthiness <- function(X, Y, k = 10) {
  .rank_metric(X, Y, k, direction = "trust")
}


#' Continuity of an embedding
#'
#' The dual of [trustworthiness()]: penalises points that are close in the
#' original space but pushed apart in the embedding.
#'
#' @inheritParams trustworthiness
#' @return A number in `[0, 1]`; 1 is perfect.
#' @family evaluation
#' @examples
#' m <- sample_manifold(miura_ori(3, 3), pi / 8, 60, seed = 1)
#' continuity(m$X, m$truth, k = 5)
#' @export
continuity <- function(X, Y, k = 10) {
  .rank_metric(X, Y, k, direction = "cont")
}


#' k-nearest-neighbour preservation
#'
#' Mean overlap between the `k` nearest neighbours in the original space and in
#' the embedding.
#'
#' @inheritParams trustworthiness
#' @return A number in `[0, 1]`.
#' @family evaluation
#' @examples
#' m <- sample_manifold(miura_ori(3, 3), pi / 8, 60, seed = 1)
#' knn_preservation(m$X, m$truth, k = 5)
#' @export
knn_preservation <- function(X, Y, k = 10) {
  X <- .as_matrix(X, "X"); Y <- .as_matrix(Y, "Y")
  if (nrow(X) != nrow(Y)) stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  nx <- .knn(X, k)
  ny <- .knn(Y, k)
  mean(vapply(seq_len(nrow(nx)),
              function(i) length(intersect(nx[i, ], ny[i, ])) / k,
              numeric(1)))
}


#' First folding parameter meeting a neighbourhood shortcut heuristic
#'
#' Uses one fixed material sample across the supplied parameter grid. A shortcut
#' is an ambient kNN edge whose flat Euclidean length exceeds the median of all
#' sampled pairwise flat distances. This is a sample- and k-dependent diagnostic,
#' not a theorem about failure of a dimension-reduction method. It may already
#' trigger at zero for sparse samples or large k.
#' @param pattern A Miura `crease_pattern`.
#' @param k Positive integer neighbourhood size less than `n`.
#' @param n Number of sampled material points.
#' @param grid Candidate normalised parameters in `[0, 1]`, sorted internally.
#' @param seed Optional integer seed.
#' @return First matching parameter or `NA_real_` when no grid value matches.
#' @family evaluation
#' @examples
#' short_circuit_theta(miura_ori(3, 3), k = 5, n = 60,
#'                     grid = seq(0, 0.95, length.out = 5), seed = 1)
#' @export
short_circuit_theta <- function(pattern, k = 10, n = 300,
                                grid = seq(0, 0.95, length.out = 25), seed = NULL) {
  .check_pattern(pattern); .check_count(n, "n"); .check_count(k, "k")
  if (k >= n) stop("`k` must be smaller than `n`.", call. = FALSE)
  if (!is.numeric(grid) || !length(grid) || any(!is.finite(grid)) || any(grid < 0 | grid > 1)) {
    stop("`grid` must contain normalised parameters in [0, 1].", call. = FALSE)
  }
  material <- sample_manifold(pattern, theta = 0, n = n, seed = seed)$truth
  d_true <- as.matrix(stats::dist(material))
  threshold <- stats::median(d_true[upper.tri(d_true)])
  for (theta in sort(unique(grid))) {
    nb <- .knn(.map_flat(pattern, material, theta), k)
    if (any(vapply(seq_len(n), function(i) any(d_true[i, nb[i, ]] > threshold), logical(1)))) return(theta)
  }
  NA_real_
}
