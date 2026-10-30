# Standalone package revival

Owner decision, 30 October 2026: revive the standalone scientific package.

The release target is an installable R package for reproducible dimension-
reduction experiments with known material coordinates. The supported mechanism
is a finite, zero-thickness Miura sheet with free boundaries.

1. Essential: replace the accordion map with Miura kinematics; verify full facet
   distances, planarity, signed crease assignments and flat endpoints across
   sizes and shape parameters. Expose explicit facets and boundary edges.
2. Everyday use: sample the actual finite sheet, use the same sample across a
   fold sweep, validate metrics and invalid inputs, and provide an installed
   first-use example, migration guide and reproducible CI.
3. Core enhancement: product and orthonormal lift workflows with distance
   checks, a compact reproducible benchmark, and documented provenance.
4. Deferred research: certify Yoshimura/waterbomb kinematics, collision/thickness
   mechanics, and scalable nearest-neighbour metrics. A flat pattern alone does
   not qualify a mechanism for the benchmark.

Completion means the geometry and workflow tests, source-package check, clean
installation and example pass locally. Remote CI, CRAN acceptance and external
scientific review are separate evidence. Version 0.2 changes the old prototype's
angle convention and side-length interpretation; the migration guide is required.
