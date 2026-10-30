# Standalone revival and verification — 30 October 2026

Owner decision: revive the standalone package. The release scope and acceptance
criteria are in ROADMAP.md; first-use and intentional migration changes are in
README.md. Original source checkout and scientific source data were preserved.

## Findings and completed work

| Finding | Evidence and user impact | Resolution |
|---|---|---|
| High: incorrect mechanism claims | The old `R/fold.R` used the same accordion map for every family. Its Miura tests checked only a column arc; hard-coded vertex rules could not detect the wrong mechanism. | Implemented explicit Miura facets and the closed-form map, actual incidence/sector tests, full facet distances/planarity and independently derived crease signs. |
| High: unsupported families | Yoshimura had no triangular faces; waterbomb merely relabelled that grid. | Correct flat Yoshimura mesh, explicit unsupported folding errors, waterbomb refusal and documented research roadmap. |
| Medium: sample geometry | Bounding-box sampling could leave the staggered finite sheet; material noise could leave its boundary. | Uniform equal-area facet sampling and reflected material noise, verified against finite-sheet membership and affine rigid maps. |
| Medium: metric/sweep errors | k=1 collapsed a neighbour matrix; rank normalisation accepted invalid k; each sweep angle resampled truth. | Stable matrix shape, valid domains, matching dimensions/nonconstant reconstruction, and a fixed material sample across the sweep. Hand-computed k=1 expectations protect the metric repair. |
| Medium: distribution | No onboarding/CI; a linked worktree's `.git` file entered an initial source archive. | Installed example, migration/geometry docs, pinned CI, explicit `.git` exclusion; final archive check is clean. |

The version-0.1 results cannot be reinterpreted as Miura results. Version 0.2
requires recomputation and uses normalised theta and the stated facet side lengths.
No existing empirical result or book claim was rewritten to fit the implementation.

## Observed verification

Environment: macOS 27.0.1 ARM64, R 4.6.1. Dependencies installed in an isolated
maintenance R library; no user library was modified for this work.

- Full source suite and final packaged tests passed: 1,805 assertions, no failures,
  warnings or skips. These include full facet-distance and planarity sweeps,
  geometry-derived crease signs, actual Kawasaki/Maekawa incidence, finite
  sampling/noise, identical sweep samples, product/lift distances and metric
  failure domains.
- `R CMD build .` and `R CMD check --no-manual foldbenchR_0.2.0.tar.gz` passed with
  **Status: OK**. Examples and all packaged tests ran. PDF manual generation was
  not requested; no TeX qualification is claimed.
- Installed that archive into a fresh package library with declared dependencies
  available, changed out of the checkout, and ran the installed benchmark.
  Flat-sheet PCA reconstruction RMSE was 2.001275e-15; all four sweep rows were
  finite. This is a deterministic software demonstration, not external evidence
  about a dimension-reduction method's general performance.
- Rendered and visually checked the installed crease plot and folded x-z profile.
  Boundaries, signed crease legend and axis labels were legible.
- Checked the final diff and package contents. Two accidentally overlapping
  intermediate check attempts were discarded; only the separate final run above
  supports readiness. An initial install without the dependency-library path
  correctly failed for missing vegan; the stated installed check includes it.

Local logs, archive and plots are retained under
`/Users/rh/Documents/GitHub/.repository-audit/2026-10-30-rheller-refinement/foldbenchR/`:
`package-check-final/`, `installed-first-use-final.log`, `pattern.png`, `fold.png`.

## Readiness and remaining work

The defined standalone Miura scope is locally usable and source-package checked.
Windows/Linux CI is configured but has not run remotely. No CRAN submission,
GitHub publication or external scientific review was performed. Yoshimura and
waterbomb folding, finite-thickness/contact mechanics and large-n metric scaling
remain explicitly deferred. The package makes no universal shortcut or performance
claim. Next release action: review the 0.2 migration and candidate, run the configured
CI on an authorized branch, then prepare the chosen distribution submission.
