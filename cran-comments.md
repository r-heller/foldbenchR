## Submission

This is a new submission.

foldbenchR is pure R code with CRAN dependencies only. It needs no network
access, external software or compiled code. Examples, tests and the vignette
use small synthetic samples and complete in a few seconds.

## Test environments

* local: macOS (aarch64), R 4.6.1, `R CMD check --as-cran --no-manual`
  (the PDF manual was built locally with TinyTeX (TeX Live 2026) via
  `R CMD Rd2pdf` without errors; the check run itself used `--no-manual`)
* GitHub Actions: macos-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1), with `--as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

Should the spell checker flag 'Miura', 'Yoshimura' or 'waterbomb' in the
description: these are names of origami crease patterns and are spelled
correctly.
