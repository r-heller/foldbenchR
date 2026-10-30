test_that("reconstruction_error is zero for the identity embedding", {
  truth <- unfold(miura_ori(4, 4))
  expect_equal(reconstruction_error(truth, truth), 0, tolerance = 1e-8)
})

test_that("reconstruction_error ignores similarity transforms", {
  truth <- unfold(miura_ori(4, 4))
  ang <- pi / 5
  rot <- matrix(c(cos(ang), -sin(ang), sin(ang), cos(ang)), 2, 2)
  moved <- (truth %*% rot) * 3.7
  moved <- sweep(moved, 2, c(10, -4), "+")
  expect_equal(reconstruction_error(moved, truth), 0, tolerance = 1e-8)
})

test_that("rank metrics are 1 for a faithful embedding and bounded otherwise", {
  m <- sample_manifold(miura_ori(3, 3), theta = 0.3, n = 60, seed = 1)
  expect_equal(trustworthiness(m$truth, m$truth, k = 5), 1, tolerance = 1e-8)
  expect_equal(continuity(m$truth, m$truth, k = 5), 1, tolerance = 1e-8)
  expect_equal(knn_preservation(m$truth, m$truth, k = 5), 1)

  t_mixed <- trustworthiness(m$X, m$truth, k = 5)
  expect_gte(t_mixed, 0); expect_lte(t_mixed, 1)
})

test_that("shortcut sweep uses one material sample with an independently measured threshold", {
  p <- miura_ori(4, 4); grid <- seq(0, 0.95, length.out = 7); n <- 60; k <- 5
  truth <- sample_manifold(p, 0, n, seed = 23)$truth
  d <- as.matrix(dist(truth)); threshold <- median(d[upper.tri(d)])
  expected <- NA_real_
  for (th in grid) {
    X <- sample_manifold(p, th, n, seed = 23)$X
    ambient <- as.matrix(dist(X)); diag(ambient) <- Inf
    edges <- lapply(seq_len(n), function(i) order(ambient[i, ])[seq_len(k)])
    if (any(vapply(seq_len(n), function(i) any(d[i, edges[[i]]] > threshold), logical(1)))) {
      expected <- th; break
    }
  }
  expect_equal(short_circuit_theta(p, k, n, rev(grid), seed = 23), expected)
  expect_true(is.na(short_circuit_theta(p, 5, 60, grid = 0, seed = 1)))
  expect_error(reconstruction_error(matrix(1, 4, 2), matrix(2, 4, 2)), "non-constant")
  expect_error(reconstruction_error(matrix(1:8, 4), matrix(1:12, 4)), "dimensions")
})
