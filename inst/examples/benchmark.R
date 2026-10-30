library(foldbenchR)
pattern <- miura_ori(4, 6)
results <- do.call(rbind, lapply(c(0, 0.3, 0.6, 0.9), function(theta) {
  m <- sample_manifold(pattern, theta, n = 200, seed = 42)
  embedding <- prcomp(m$X)$x[, 1:2, drop = FALSE]
  data.frame(theta = theta,
             reconstruction_rmse = reconstruction_error(embedding, m$truth),
             trustworthiness = trustworthiness(m$truth, embedding, 10),
             continuity = continuity(m$truth, embedding, 10),
             neighbour_overlap = knn_preservation(m$truth, embedding, 10))
}))
stopifnot(all(is.finite(as.matrix(results))), results$reconstruction_rmse[1] < 1e-10)
print(results, row.names = FALSE)
# Every row uses the same 200 material points; only the folding parameter changes.
