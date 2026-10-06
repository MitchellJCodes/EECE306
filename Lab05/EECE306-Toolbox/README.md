## Team 3, Amelia, Mitch, and Carissa EECE 306 Toolbox

## Layout

```
EECE306-Toolbox/
├── +em/          the toolbox package, one subfolder per module
├── tests/        your test files, test\\\\\\\_lab01.m onward
├── runTests.m    the test runner, provided, do not modify
└── README.md     this file, replace with your own documentation
```

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

## Helpful Notes
* An NxM vector has N rows and M columns (Ex: an Nx1 vector has an 
  arbitrary number of rows N and one column)
* NaN means "Not a Number", which is a shorthand for indeterminate and 
  undefined values such as 0/0.

