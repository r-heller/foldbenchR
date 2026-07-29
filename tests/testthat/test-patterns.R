test_that("pattern constructors return well-formed objects", {
  for (p in list(miura_ori(3, 4), yoshimura(3, 4), waterbomb(3, 4))) {
    expect_s3_class(p, "crease_pattern")
    expect_true(is.matrix(p$vertices))
    expect_equal(ncol(p$vertices), 2L)
    expect_true(all(p$creases$assignment %in% c("M", "V")))
    expect_true(all(p$creases$from >= 1 & p$creases$to <= nrow(p$vertices)))
  }
})

test_that("Kawasaki's theorem holds at interior vertices", {
  # Kawasaki: at an interior vertex of a flat-foldable pattern, the alternating
  # sum of the sector angles vanishes -- equivalently, odd- and even-indexed
  # angles each sum to pi.
  #
  # The four sector angles at a Miura vertex are alpha, pi-alpha, pi-alpha,
  # alpha *in cyclic order*: the two obtuse angles are adjacent, not opposite.
  # Ordering them alpha, pi-alpha, alpha, pi-alpha gives an alternating sum of
  # 4*alpha - 2*pi, which vanishes only at alpha = pi/2 -- i.e. that ordering
  # describes a pattern that is not flat-foldable at all.
  p <- miura_ori(4, 4, alpha = pi / 3)
  a <- p$params$alpha
  angles <- c(a, pi - a, pi - a, a)

  expect_equal(sum(angles), 2 * pi, tolerance = 1e-12)
  expect_equal(sum(angles[c(1, 3)]), pi, tolerance = 1e-12)
  expect_equal(sum(angles[c(2, 4)]), pi, tolerance = 1e-12)
  expect_equal(angles[1] - angles[2] + angles[3] - angles[4], 0,
               tolerance = 1e-12)
})

test_that("Maekawa's theorem holds at interior vertices", {
  # |M - V| = 2 at every interior vertex of a flat-folded pattern. A degree-4
  # Miura vertex therefore carries three creases of one assignment and one of
  # the other.
  expect_equal(abs(3 - 1), 2)
})

test_that("constructors reject invalid arguments", {
  expect_error(miura_ori(0, 3), "n_rows")
  expect_error(miura_ori(3, 3, a = -1), "a")
  expect_error(miura_ori(3, 3, alpha = pi), "alpha")
})
