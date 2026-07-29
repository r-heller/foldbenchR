test_that("unfold(fold(p, theta)) recovers unfold(p) for every theta", {
  p <- miura_ori(4, 4)
  for (th in seq(0, 1.5, length.out = 8)) {
    expect_equal(unfold(fold(p, th)), unfold(p), tolerance = 1e-12)
  }
})

test_that("folding is an isometry: geodesic distances are invariant", {
  # Arc length along the accordion is preserved by construction. Measured on
  # the vertex lattice: the summed step lengths along a column must not change
  # with theta.
  p <- miura_ori(6, 2)
  col1 <- which(p$grid$col == 1)
  arc <- function(m) sum(sqrt(rowSums(diff(m[col1, , drop = FALSE])^2)))
  flat <- arc(cbind(p$vertices, z = 0))
  for (th in c(0, 0.3, 0.7, 1.1)) {
    expect_equal(arc(fold(p, th)$coords), flat, tolerance = 1e-10)
  }
})

test_that("fold rejects theta outside [0, pi/2)", {
  p <- miura_ori(2, 2)
  expect_error(fold(p, -0.1), "theta")
  expect_error(fold(p, pi / 2), "theta")
})

test_that("sample_manifold returns matching ambient and true coordinates", {
  s <- sample_manifold(miura_ori(3, 3), theta = 0.4, n = 40, seed = 1)
  expect_equal(nrow(s$X), 40L)
  expect_equal(nrow(s$truth), 40L)
  expect_equal(ncol(s$X), 3L)
  expect_equal(ncol(s$truth), 2L)
})

test_that("product and lift raise dimension while truth stays exact", {
  a <- sample_manifold(miura_ori(2, 2), 0.3, 30, seed = 1)
  b <- sample_manifold(yoshimura(2, 2), 0.3, 30, seed = 2)
  pr <- fold_product(a, b)
  expect_equal(ncol(pr$X), 6L)       # ambient 3k
  expect_equal(ncol(pr$truth), 4L)   # intrinsic 2k
  lf <- lift(pr, 20, seed = 1)
  expect_equal(ncol(lf$X), 20L)
  expect_equal(lf$truth, pr$truth)   # lift must not touch the truth
})
