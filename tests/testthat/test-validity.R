test_that("one-neighbour metrics compare every point", {
  X <- matrix(c(0, 1, 4, 10), ncol = 1)
  Y <- matrix(c(0, 4, 1, 10), ncol = 1)
  expect_equal(knn_preservation(X, Y, 1), 0)
  expect_equal(trustworthiness(X, Y, 1), 0.5)
  expect_equal(continuity(X, Y, 1), 0.5)
})

test_that("metric and sampler invalid domains are rejected", {
  X <- matrix(seq_len(20), ncol = 2)
  for (k in list(0, 1.5, NA_real_, Inf, 5)) expect_error(trustworthiness(X, X, k), "k")
  expect_error(knn_preservation(X, X[-1, ], 1), "same number")
  X[1] <- NA_real_
  expect_error(knn_preservation(X, X, 1), "finite")
  expect_error(sample_manifold(miura_ori(), theta = Inf), "theta")
  expect_error(sample_manifold(miura_ori(), n = 2.5), "integer")
  expect_error(sample_manifold(miura_ori(), noise = list(sd = -1)), "sd")
})

test_that("sampled points lie within the staggered facets", {
  p <- miura_ori(3, 2, a = 1.2, b = 0.8)
  s <- sample_manifold(p, theta = 0, n = 1000, seed = 31)
  y <- s$truth[, 'y'] / (p$params$b * sin(p$params$alpha))
  row <- pmin(floor(y), p$params$n_rows - 1)
  frac <- y - row
  shift <- (ifelse(row %% 2 == 0, frac, 1 - frac)) * p$params$b * cos(p$params$alpha)
  relative_x <- s$truth[, 'x'] - shift
  expect_true(all(relative_x >= 0 & relative_x <= p$params$n_cols * p$params$a))
  expect_equal(unname(s$X[, 1:2]), unname(s$truth))
  expect_true(all(s$X[, 3] == 0))
})

test_that("Miura interior crease assignments satisfy the vertex rule", {
  p <- miura_ori(4, 4)
  inside <- which(p$grid$row > 1 & p$grid$row < 5 & p$grid$col > 1 & p$grid$col < 5)
  for (i in inside) {
    edges <- p$creases[p$creases$from == i | p$creases$to == i, ]
    expect_equal(nrow(edges), 4)
    expect_equal(abs(sum(edges$assignment == 'M') - sum(edges$assignment == 'V')), 2)
  }
})


test_that("manifold noise stays on the finite sheet and noise truth is explicit", {
  p <- miura_ori(2, 3)
  s <- sample_manifold(p, 0.7, 500, noise = list(type = "manifold", sd = 10), seed = 2)
  y <- s$truth[, 2] / (p$params$b * sin(p$params$alpha))
  row <- pmin(floor(y), p$params$n_rows - 1L); v <- y - row
  shift <- ifelse(row %% 2 == 0, v, 1-v) * p$params$b * cos(p$params$alpha)
  x <- (s$truth[, 1] - shift) / p$params$a
  expect_true(all(y >= 0 & y <= 2 & x >= 0 & x <= 3))
  expect_equal(s$X, foldbenchR:::.map_flat(p, s$truth, 0.7))
  clean <- sample_manifold(p, 0.5, 30, seed=4)
  noisy <- sample_manifold(p, 0.5, 30, list(type="ambient", sd=0.1), seed=4)
  expect_equal(noisy$truth, clean$truth)
  expect_gt(sum((noisy$X - clean$X)^2), 0)
})
