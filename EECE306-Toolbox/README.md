## EECE 306 Toolbox, team starter

Rename this folder to `EECE306-Toolbox` and keep it under version control or in
a dated backup. Everything you build this semester lives here.

## Layout

```
EECE306-Toolbox/
├── +em/          the toolbox package, one subfolder per module
├── tests/        your test files, test\\\\\\\_lab01.m onward
├── runTests.m    the test runner, provided, do not modify
└── README.md     this file, replace with your own documentation
```

## First steps

1. Put the folder containing `+em` on the Octave path with `addpath`. Do **not** add `+em` itself.
2. Delete each `CONTENTS.m` placeholder as you add real functions to that module.
3. Run `runTests` from the toolbox root. With no functions written yet it will
   report a failure, which is the correct starting state.

## What to put in this README by the end of the semester

Team name and members, a one line description of every module, and a short
example showing somebody else how to compute a field with your library.
The engineering quality portion of the grade is assessed partly by whether
another team could use your toolbox from this file alone.

## Team 3 - Amelia Harris, Carissa McWilliams and Mitch Nazareth

## Module Descriptions

Constant Functions
* c0 = em.const.c0() - returns the scalar value of the speed of light c. 
  Note: This constant is computed from eps0 and mu0.
* e0 = em.const.eps0() - returns the scalar value of the electric constant 
  epsilon_0 (permittivity of free space).
* eta = em.const.eta0() - returns the scalar value of the intrinsic 
  impedance of free space eta (eta_0). Note: This constant is
  computed from eps0 and mu0.
* m0 = em.const.mu0() - returns the scalar value of the magnetic constant 
  mu_0 (permeability of free space).

Coordinate Functions (Point Conversion)
  Note: 0 <= phi < 2pi and 0 <= theta <= pi
* p = em.coord.c2cyl(r) - returns the corresponding Nx3 cylindrical 
  coordinate (rho, phi, z) of an Nx3 Cartesian coordinate (x, y, z); takes 
  an Nx3 Cartesian coordinate input.
* s = em.coord.c2sph(r) - returns the corresponding Nx3 spherical 
  coordinate (r, theta, phi) of an Nx3 Cartesian coordinate (x, y, z); 
  takes an Nx3 Cartesian coordinate input.
* r = em.coord.cyl2c(p) - returns the corresponding Nx3 Cartesian 
  coordinate (x, y, z) of an Nx3 cylindrical coordinate (rho, phi, z); 
  takes an Nx3 cylindrical coordinate input.
* r = em.coord.sph2c(s) - returns the corresponding Nx3 Cartesian 
  coordinate (x, y, z) of an Nx3 spherical coordinate (r, theta, phi); 
  takes an Nx3 spherical coordinate input.

Coordinate Functions (Vector Component Conversion)
* Ac = em.coord.vecC2Cyl(A, r) - returns the corresponding Nx3 cylindrical 
  vector components (A_rho, A_phi, A_z) of an Nx3 Cartesian component set
  (A_x, A_y, A_z) along an Nx3 vector in the Cartesian plane; takes an 
  array of Nx3 Cartesian vector components and an Nx3 Cartesian position.
* As = em.coord.vecC2Sph(A, r) - returns the corresponding Nx3 spherical 
  vector components (A_r, A_theta, A_phi) of an Nx3 Cartesian component set
  (A_x, A_y, A_z) along an Nx3 vector in the Cartesian plane; takes an 
  array of Nx3 Cartesian vector components and an Nx3 Cartesian position.
* A = em.coord.vecCyl2C(Ac, r) - returns the corresponding Nx3 Cartesian
  vector components (A_x, A_y, A_z) of a set of Nx3 cylindrical vector
  component (A_rho, A_phi, A_z) along an Nx3 vector in the Cartesian plane;
  takes an array of Nx3 cylindrical vector components and an Nx3 Cartesian
  position.
* A = em.coord.vecSph2C(As, r) - returns the corresponding Nx3 Cartesian
  vector components (A_x, A_y, A_z) of a set of Nx3 spherical vector
  component (A_r, A_theta, A_phi) along an Nx3 vector in the Cartesian
  plane; takes an array of Nx3 cylindrical vector components and an Nx3 
  Cartesian position.

Field Functions
* C = em.field.circulation(F, curve, tspan, N, rule) returns the scalar 
  flux of a vector field around a curve; takes a vector field function 
  handle F, parametrized curve, span of t tspan, number of quadrature 
  points N, and quadrature rule.
* Dv = em.field.D(s, r, epsr) returns the Nx3 flux density produced by a 
  charge source or a set of charge sources; takes a discretized charge 
  source struct s, a position vector r, and the relative permittivity epsr.
* Ev = em.field.E(s, r) returns the Nx3 field produced by a constructed 
  point charge at a 1x3 position vector r; takes a discretized charge 
  source struct s and a position vector r.
* Phi = em.field.flux(F, surf, uspan, vspan, Nu, Nv, rule) returns the 
  scalar flux passing though a parametrized surface; takes a vector field 
  function handle F, parametrized surface, span of u uspan, span of v 
  vspan, the number of quadrature points of u (Nu) and v (Nv), and 
  quadrature rule.
* V = em.field.V(s, r) returns the Nx1 potentials of points inside an 
  electric field; takes a discretized charge source struct s and a position 
  vector r.

Differential Operators
* d = em.op.div(F, r, h) returns the Nx3 divergence of a vector field;
  takes a vector field function handle F, a position vector r, and a finite
  difference step size h.
* g = em.op.grad(f, r, h) returns the Nx3 gradient of a scalar field; takes 
  a vector field function handle F, a position vector r, and a finite
  difference step size h.
* L = em.op.laplacian(f, r, h) evaluates d2f/dx2+d2f/dy2 at rows of r; takes 
  a vector field function handle F, a position vector r, and a finite
  difference step size h.

Quadratures
* [pos, w] = em.quad.line(curve, tspan, N, rule) returns the Nx3 element 
  positions POS and Nx1 weights W with units of length; takes parametrized 
  curve, span of t tspan, number of quadrature points N, quadrature rule.
* [t, w] = em.quad.nodes(a, b, N, rule) returns a Nx1 column vector of 
  nodes (T) and an Nx1 column vector of weights (W); takes interval 
  boundaries a and b, number of quadrature points N, and quadrature rule.
* [pos, w] = em.quad.surf(surface, uspan, vspan, Nu, Nv, rule) returns Mx3
  element positions and Mx1 weights with units of area; takes parametrized 
  curve, span of u uspan, span of v vspan, the number of quadrature points 
  of u (Nu) and v (Nv), and quadrature rule.
* [pos, weight] = em.quad.vol(volume, uspan, vspan, wspan, Nu, Nv, Nw, rule) 
  returns Mx3 element positions and Mx1 weights with units of volume; takes 
  parametrized curve, spans of u (uspan), v (vspan), and w (wspan) the 
  number of quadrature points of u (Nu), v (Nv), and w (Nw), and quadrature
  rule.

Solvers
* G = em.solve.grid2d(xlim, ylim, Nx, Ny) returns a struct with G.x, G.y, 
  G.X, G.Y, G.hx and G.hy; takes x-axis limits, y-axis limits,  number of 
  quadrature points for x (Nx) and y (Ny).
* [V, info] = em.solve.poisson2d(G, rho, epsr, mask, Vfix, opts) solves 
  div(eps grad V) = -rho on a 2D grid; takes grid struct G from 
  em.solve.grid2d, relative permittivity epsr, boolean value mask, values 
  at the fixed nodes where mask is true (Vfix), and option struct opts.

Source Functions
* s = em.src.lineCharge(rhoL, curve, tspan, N, rule) returns a struct with 
  a discretized list of the differential elements of a line charge; takes 
  line charge density rhoL, parametrized curve, span of t tspan, number of
  quadrature points N, and the quadrature rule.
* s = merge(s1, s2) returns a struct that joins the element lists of two 
  point sources together; takes two point charge structs.
* s = pointCharge(Q, r0) eturns a struct with a discretized list of 
  differential elements of a point charge. It contains the type,
  position, quadrature weights, and charge density of a point source. Takes
  a charge density and position vectors.
* s = surfCharge(rhoS, surf, uspan, vspan, Nu, Nv, rule) returns a struct 
  with a discretized list of the differential elements of a surface charge;
  takes the surface charge density rhoS, parametrized surface, span of u
  uspan and v vspan, number of quadrature points for u (Nu) and v (Nv), and 
  quadrature rule.

Tests
* assertClose(actual, expected, tol, msg) compares with a relative 
  tolerance when expected is nonzero and an absolute tolerance when 
  expected is zero. It fails (throws an error) when the error is greater 
  than the tolerance and passes (is silent) otherwise.
* [p, h] = convergence(errFcn, Nlist, opts) calls errFcn(N) for each N in 
  Nlist, fits a straight line to log(error) against log(N), and returns the 
  negative of the slope as the observed order p. It also produces the 
  log-log plot. 
* rho = em.test mms2d(Vexact, epsr) implements the method of manufactured 
  solutions. 

Vector Functions
* th = em.vec.angle(A, B) - returns the Nx1 interior angle theta between 
  two Nx3 vectors A and B (in radians); takes two Nx3 vector inputs. 
  Note: cos(theta) is clamped to [-1, 1]
* d = em.vec.fromTo(P, Q) - returns the Nx3 displacement vector d (= Q - P)
  between two Nx3 vectors P and Q; takes two Nx3 vector inputs.
* m = em.vec.mag(A) - returns the Nx1 magnitude |A| of an Nx3 vector A; 
  takes an Nx3 vector input.
* u = em.vec.unit(A) - returns the normalized Nx3 vector of a non-zero Nx3
  vector (a zero vector returns NaN in all directions); takes an Nx3 vector
  input. Note: a zero vector input will return rows of NaN, due to divison 
  by 0.

Visualization Functions
* h = contourV(F, xlim, ylim, N, opts) returns a contour plot of a scalar 
  field; takes an x-range, y-range, unit scale, and customization options.
* h = em.viz.quiver2(F, xlim, ylim, N, opts) returns the vectorized plot of
  an electric field given a function handle to sample; takes a function 
  handle, x-range, y-range, unit scale, and customization options. Note: if
  false is instead set to true, it will normalize the plot and makes the 
  field directions easier to read. However, the relative magnitudes of the 
  field lines will be lost.
* h = streamlines(F, xlim, ylim, seeds, opts) sketches the field lines from 
  the seed points; takes a function handle, x-range, y-range, unit scale, 
  and customization options.


## Helpful Notes
* An NxM vector has N rows and M columns (Ex: an Nx1 vector has an 
  arbitrary number of rows N and one column)
* NaN means "Not a Number", which is a shorthand for indeterminate and 
  undefined values such as 0/0.

