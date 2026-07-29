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
  if (nrow(embedding) != nrow(truth)) {
    stop("`embedding` and `truth` must have the same number of rows.",
         call. = FALSE)
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
#' @param k Neighbourhood size.
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
  nx <- .knn(X, k)
  ny <- .knn(Y, k)
  mean(vapply(seq_len(nrow(nx)),
              function(i) length(intersect(nx[i, ], ny[i, ])) / k,
              numeric(1)))
}


#' Fold angle at which the neighbourhood graph first short-circuits
#'
#' Finds the smallest fold angle on a grid at which the `k`-nearest-neighbour
#' graph of the folded manifold first connects two facets that are far apart
#' geodesically. This is the failure mode that breaks every kNN-graph method,
#' and knowing where it starts is what turns a curve over `theta` into a
#' diagnosis.
#'
#' @param pattern A `crease_pattern`.
#' @param k Neighbourhood size.
#' @param n Points sampled per candidate angle.
#' @param grid Vector of candidate fold angles.
#' @param seed Optional integer seed.
#' @return The first angle at which a short circuit appears, or `NA_real_` if
#'   none does on the supplied grid.
#' @family evaluation
#' @examples
#' short_circuit_theta(miura_ori(3, 3), k = 5, n = 60,
#'                     grid = seq(0, 1.2, length.out = 5), seed = 1)
#' @export
short_circuit_theta <- function(pattern, k = 10, n = 300,
                                grid = seq(0, pi / 2 - 0.05, length.out = 25),
                                seed = NULL) {
  .check_pattern(pattern)
  if (!is.null(seed)) set.seed(seed)
  for (th in grid) {
    m <- sample_manifold(pattern, theta = th, n = n)
    nb <- .knn(m$X, k)
    d_true <- as.matrix(stats::dist(m$truth))
    # A short circuit: ambient neighbours whose true separation exceeds the
    # median. The threshold is deliberately blunt — the transition is sharp,
    # so the exact cut-off does not move the answer much.
    thr <- stats::median(d_true[upper.tri(d_true)])
    hit <- any(vapply(seq_len(nrow(nb)),
                      function(i) any(d_true[i, nb[i, ]] > thr), logical(1)))
    if (hit) return(th)
  }
  NA_real_
}
