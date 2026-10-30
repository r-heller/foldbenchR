cross3 <- function(a, b) c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])

test_that("rigid folding preserves all facet distances and planarity across shapes", {
  for (nr in c(1, 3)) for (nc in c(1, 4)) for (alpha in c(0.05, pi/4, 1.5)) {
    p <- miura_ori(nr, nc, a = 1.7, b = 0.8, alpha = alpha)
    for (theta in c(0, 0.1, 0.5, 0.95, 1)) {
      m <- fold(p, theta)
      expect_equal(unfold(m), p$vertices)
      if (theta == 0) expect_equal(m$coords, cbind(p$vertices, z = 0))
      for (ids in p$facets) {
        q <- m$coords[ids, ]
        expect_equal(as.numeric(dist(q)), as.numeric(dist(p$vertices[ids, ])), tolerance = 1e-11)
        normal <- cross3(q[2, ] - q[1, ], q[3, ] - q[1, ])
        expect_equal(sum(normal * (q[4, ] - q[1, ])), 0, tolerance = 1e-11)
      }
    }
  }
})

test_that("every interior crease folds with the advertised sign", {
  for (alpha in c(0.2, pi/3, 1.4)) for (theta in c(0.1, 0.5, 0.9)) {
    p <- miura_ori(4, 5, alpha = alpha)
    V <- fold(p, theta)$coords
    for (e in seq_len(nrow(p$creases))) {
      i <- p$creases$from[e]; j <- p$creases$to[e]
      fs <- which(vapply(p$facets, function(f) all(c(i,j) %in% f), logical(1)))
      if (length(fs) == 1L) { expect_equal(p$creases$assignment[e], "B"); next }
      # A height relative to the oriented neighbouring facet determines side.
      # With CCW facets, the facet walking i->j has its interior to the left.
      walks <- vapply(p$facets[fs], function(f) any(f == i & c(f[-1],f[1]) == j), logical(1))
      expect_equal(sum(walks), 1L)
      f <- p$facets[[fs[walks]]]; g <- p$facets[[fs[!walks]]]
      off <- setdiff(f, c(i,j))[1]; opposite <- setdiff(g, c(i,j))[1]
      normal <- cross3(V[j, ] - V[i, ], V[off, ] - V[i, ])
      height <- sum(normal * (V[opposite, ] - V[i, ]))
      expect_gt(abs(height), 1e-9)
      expect_equal(if (height < 0) "M" else "V", p$creases$assignment[e])
    }
  }
})

test_that("normalised fold domain and rigid sample mapping are enforced", {
  p <- miura_ori(1, 1, a = 1.3, b = 0.7, alpha = 0.6)
  expect_error(fold(p, -0.1), "theta")
  expect_error(fold(p, 1.1), "theta")
  flat <- sample_manifold(p, 0, 100, seed = 1)
  for (theta in c(0.2, 0.7, 1)) {
    s <- sample_manifold(p, theta, 100, seed = 1)
    expect_equal(s$truth, flat$truth)
    expect_equal(as.numeric(dist(s$X)), as.numeric(dist(flat$truth)), tolerance = 1e-10)
  }
  p <- miura_ori(3, 4)
  expect_equal(foldbenchR:::.map_flat(p, p$vertices, 0.5), fold(p, 0.5)$coords, tolerance = 1e-12)
})

test_that("product and lift retain material rows and ambient distances", {
  a <- sample_manifold(miura_ori(2, 2), 0.3, 30, seed = 1)
  b <- sample_manifold(miura_ori(2, 2, alpha = 0.6), 0.7, 30, seed = 2)
  pr <- fold_product(a, b)
  expect_equal(dim(pr$X), c(30L, 6L)); expect_equal(dim(pr$truth), c(30L, 4L))
  lf <- lift(pr, 20, seed = 1)
  expect_equal(lf$truth, pr$truth)
  expect_equal(as.numeric(dist(lf$X)), as.numeric(dist(pr$X)), tolerance = 1e-10)
  expect_error(fold_product(a, list(X=b$X, truth=b$truth[-1, ])), "matching rows")
})
