# Team 3 - EECE 306 Toolbox

A Octave toolbox for EECE 306, Fundamentals of Electromagnetics Lab.

## Repository structure

- `+em/+const/` — physical constants (`eps0`, `mu0`, `c0`, `eta0`).
- `+em/+vec/` — vector magnitude, unit vectors, angles, and displacements.
- `+em/+coord/` — Cartesian, cylindrical, and spherical coordinate and vector-component conversions.
- `+em/+quad/` — one-, two-, and three-parameter numerical quadrature with geometric Jacobians.
- `+em/+src/` — discretized charge/current source constructors and source merging.
- `+em/+field/` — electric field, potential, electric flux density, flux, and circulation.
- `+em/+op/` — numerical gradient, divergence, curl, and Laplacian operators as implemented.
- `+em/+solve/` — Cartesian grid construction and a 2-D variable-permittivity Poisson solver with harmonic (default) or arithmetic face averaging.
- `+em/+viz/` — vector-field and scalar-potential plotting helpers.
- `+em/+test/` — assertion and numerical-convergence verification helpers.
- `tests/` — lab-specific regression checks (`test_lab01.m` through the latest included lab test).
- `runTests.m` — provided test runner; run it from the toolbox root.

## Running tests

Open Octave, change the current folder to `EECE306-Toolbox`, and run:

```octave
runTests
```

The runner executes the `tests/test_lab*.m` files it finds and summarizes
successes and failures. A passing run is meaningful only when each lab test
contains assertions for its required numerical results, dimensions, errors,
and warnings.

## Conventions

- Positions and vector fields use one point per row (`N`-by-`3`). Scalar fields use `N`-by-`1` outputs.
- Toolbox angles are in radians and quantities use SI units.
- Functions are namespaced under `+em`; do not add toolbox functions to the Octave path root.
- Each public function should document its signature, inputs, outputs, and important error conditions in its MATLAB help block.

## Current status

Lab08 tests are passing along with the previous lab tests. Lab06 tests are giving warnings reguarding the assumption that the electric potential is zero at infinity.
