# foldbenchR 0.3.0

* Prepare the package for CRAN: maintainer address at Ulm University, `URL`,
  `BugReports` and `Language` fields, a reference for the Miura model in the
  description, a spelling word list and `cran-comments.md`.
* Add an introductory vignette, `vignette("foldbenchR")`.
* `sample_manifold()` and `lift()` restore the caller's random-number state
  after using a supplied `seed`; results for a given seed are unchanged.
* Register `autoplot.fold_manifold()` as a method for `ggplot2::autoplot()`.
  The function remains exported and can still be called directly.
* Add an example to `waterbomb()`, and link README documents absolutely so
  they resolve outside the repository.
* Run `R CMD check` on a CRAN-like GitHub Actions matrix (macOS, Windows,
  and Ubuntu with R devel, release and oldrel-1).

# foldbenchR 0.2.0

* Revive standalone Miura benchmarks with a rigid facet map, explicit faces and
  boundary labels. Normalise theta to [0,1]; interpret Miura b as side length.
* Replace the placeholder Yoshimura mesh with triangles, limit it to flat
  inspection, and refuse the unsupported waterbomb and Yoshimura folds.
* Sample finite facets uniformly; reflect material noise at sheet boundaries.
* Reuse material samples during shortcut sweeps and describe the diagnostic's
  dependence on the sample and neighbourhood threshold.
* Correct one-neighbour metric shape handling, rank-normalisation domains,
  reconstruction validation, and product input checks.
* Add installed examples, geometry evidence, migration instructions and CI.
