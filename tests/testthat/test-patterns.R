test_that("constructors expose actual facet meshes and unsupported mechanisms fail", {
  p <- miura_ori(3, 4)
  expect_equal(nrow(p$vertices), 20L)
  expect_length(p$facets, 12L)
  expect_true(all(lengths(p$facets) == 4L))
  expect_setequal(unique(p$creases$assignment), c("B", "M", "V"))
  y <- yoshimura(3, 4)
  expect_length(y$facets, 24L)
  expect_true(all(lengths(y$facets) == 3L))
  expect_setequal(unique(y$creases$assignment), c("B", "U"))
  expect_error(fold(y, 0.5), "No certified folding")
  expect_error(sample_manifold(y), "No certified folding")
  expect_error(waterbomb(), "not implemented")
})

test_that("Kawasaki and Maekawa are checked on actual interior incidence", {
  for (alpha in c(0.2, pi / 4, 1.4)) {
    p <- miura_ori(4, 3, alpha = alpha)
    inside <- which(p$grid$row > 1 & p$grid$row < 5 & p$grid$col > 1 & p$grid$col < 4)
    for (v in inside) {
      edges <- p$creases[p$creases$from == v | p$creases$to == v, ]
      neighbours <- ifelse(edges$from == v, edges$to, edges$from)
      rays <- sweep(p$vertices[neighbours, , drop = FALSE], 2, p$vertices[v, ], "-")
      bearings <- sort(atan2(rays[, 2], rays[, 1]))
      angles <- diff(c(bearings, bearings[1] + 2 * pi))
      expect_equal(sum(angles[c(1, 3)]), pi, tolerance = 1e-12)
      expect_equal(sum(angles[c(2, 4)]), pi, tolerance = 1e-12)
      expect_equal(abs(sum(edges$assignment == "M") - sum(edges$assignment == "V")), 2)
    }
  }
})

test_that("constructors reject invalid and excessive grids", {
  for (bad in list(0, NA_real_, Inf, 1.5)) expect_error(miura_ori(bad, 3), "n_rows")
  expect_error(miura_ori(3, 3, a = -1), "a")
  expect_error(miura_ori(3, 3, alpha = pi), "alpha")
  expect_error(miura_ori(1e6, 1e6), "limit")
})
