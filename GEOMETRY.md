# Geometry and provenance

The revived Miura implementation follows the closed-form derivation in
`r-heller/folding-in-r`, commit `5f2fee3d81508498b17f86dcd58c07486f4c07ac`,
`R/patterns.R` and `R/folding.R`. That code is MIT, copyright 2026 R. Heller
(`LICENSE-CODE.md` in the book); this package retains the same copyright and
license. The package preserves its own row-major-in-column indexing and public
row/column argument names. No book narrative or empirical results are imported.

For material grid indices i (column) and j (row), the flat coordinates are
`(i*a + (j mod 2)*b*cos(alpha), j*b*sin(alpha))`. Put `phi=theta*alpha` and define:

```
Lx = a*cos(phi)                   H  = a*sin(phi)
dx = b*cos(alpha)/cos(phi)
Ly = b*sqrt(sin(alpha+phi)*sin(alpha-phi))/cos(phi)
V(i,j) = (i*Lx + (j mod 2)*dx, j*Ly, (i mod 2)*H)
```

The two facet side vectors obey `Lx^2+H^2=a^2`, `dx^2+Ly^2=b^2`, and
`Lx*dx=a*b*cos(alpha)`. Their lengths and included angle are therefore fixed.
Opposite sides coincide as vectors, giving planar rigid parallelograms.
The product-of-sines expression avoids subtracting nearly equal cosines.
Samples use the affine map of their own facet, including its staggered boundary.

The test suite checks every pairwise facet distance (including diagonals),
planarity, flat endpoints, and signed creases determined from the opposite
facet's displacement relative to an oriented facet normal. It measures cyclic
sectors and incidence at actual interior vertices for Kawasaki/Maekawa checks;
these necessary vertex conditions alone do not prove rigid foldability.

Context: Schenk and Guest (2013), [Geometry of Miura-folded metamaterials](https://doi.org/10.1073/pnas.1217998110), describes the ideal rigid-panel,
frictionless-hinge model and its Miura folding mode. Its angle notation differs
from this API. The equations above specify this implementation explicitly.

The book records no certified Yoshimura folding that activates every advertised
crease family. Its waterbomb implementation is also withheld. The package makes
no general impossibility claim about either origami family: those are unresolved
implementation and scientific-validation tasks here.
