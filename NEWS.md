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
