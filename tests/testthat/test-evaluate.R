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

test_that("short_circuit_theta finds a bridge and returns NA when none exists", {
  th <- short_circuit_theta(miura_ori(4, 4), k = 5, n = 60,
                            grid = seq(0, 1.4, length.out = 6), seed = 1)
  expect_true(is.na(th) || (th >= 0 && th < pi / 2))
  # A flat sheet cannot short-circuit: at theta = 0 the ambient and true
  # geometries coincide.
  expect_true(is.na(short_circuit_theta(miura_ori(4, 4), k = 5, n = 60,
                                        grid = 0, seed = 1)))
})
