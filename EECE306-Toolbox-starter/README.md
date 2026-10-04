# EECE 306 Toolbox

## Team 3 Members

- Carissa McWilliams
- Amelia Harris
- Mitch Nazareth

## Layout

```
EECE306-Toolbox/
├── +em/          the toolbox package, one subfolder per module
├── tests/        your test files, test_lab01.m onward
├── runTests.m    the test runner, provided, do not modify
└── README.md     this file, replace with your own documentation
```
## Modules

### em.const
contains electromagnetic physical constatns used throughout the toolbox.

### em.vec
Contains vector operations including magnitude, unit vectors, angles, distance, projection, and rejection.

### em.coord
Contains Cartesian, spherical, and cylindrical coordinate conversions for points and vector components. Functions accept nx3 arrays.

### em.src
Contains functions for creating and combining electromagnetic sources,
including point charges.

### em.field
Contains functions for calculating electromagnetic fields from source
structures at Nx3 observation points.

### em.viz
Contains functions for visualizing electromagnetic fields.

### em.test
Contains testing utilities used to verify toolbox functions.

% Create a point charge
s = em.src.pointCharge(1e-9, [0 0 0]);

% Define observation points
r = [1 0 0;
     2 0 0;
     3 0 0];

% Calculate the electric field
E = em.field.E(s, r);

