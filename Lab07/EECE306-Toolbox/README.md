# Team 3 - EECE 306 Toolbox

A Octave toolbox for EECE 306, Fundamentals of Electromagnetics Lab.

## Repository structure

EECE306-Toolbox/
├── +em/                # Main MATLAB package
│   ├── +const/         # Physical constants
│   ├── +coord/         # Coordinate-system conversions
│   ├── +derive/        # Reserved for derived expressions
│   ├── +field/         # Electromagnetic field calculations
│   ├── +op/            # Gradient, divergence, and Laplacian
│   ├── +quad/          # Numerical quadrature and Jacobians
│   ├── +solve/         # Numerical PDE solvers
│   ├── +src/           # Charge/current source constructors
│   ├── +test/          # Numerical verification utilities
│   ├── +vec/           # Vector utilities
│   └── +viz/           # Field and potential visualization
├── tests/              # Lab-specific regression tests
├── Lab07_notebook.m    # Lab 7 examples and analysis
├── runTests.m          # Regression test runner
└── README.md           # Project documentation

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

The toolbox is still incomplete, and the toolbox should not be seen as functional until each test is validating as passed with Octave.
