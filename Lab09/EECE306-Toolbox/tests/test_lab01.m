function test_01 = test_lab01()
%TEST_LAB01 Lab 1 Test Function
%   Runs required tests for lab 1
clear; close all;
disp(version)
fprintf('eps0 = %.10e F/m\n', em.const.eps0());
fprintf('mu0  = %.10e H/m\n', em.const.mu0());
fprintf('c0   = %.10e m/s\n', em.const.c0());
fprintf('eta0 = %.6f ohm\n',  em.const.eta0());
em.test.assertClose(em.const.c0(), 1/sqrt(em.const.mu0()*em.const.eps0()), 1e-14, 'c0 derived');
em.test.assertClose(em.const.eta0(), sqrt(em.const.mu0()/em.const.eps0()), 1e-14, 'eta0 derived');
disp('Derived constant checks passed.')
fprintf('mag([3 4 0])  = %.15g   (expect 5 exactly)\n', em.vec.mag([3 4 0]));
disp('unit([3 4 0]) ='); disp(em.vec.unit([3 4 0]))
disp('unit([0 0 0]) (expect all NaN) ='); disp(em.vec.unit([0 0 0]))
fprintf('angle perpendicular = %.15g  (expect pi/2)\n', em.vec.angle([1 0 0],[0 1 0]));
ap = em.vec.angle([1 0 0],[-1 0 0]);
fprintf('angle antiparallel  = %.15g, real = %d  (expect pi, 1)\n', ap, isreal(ap));
disp('fromTo([1 1 1],[2 3 4]) ='); disp(em.vec.fromTo([1 1 1],[2 3 4]))
A = randn(100,3); B = randn(100,3);
fprintf('mag    input 100x3 output %s\n', mat2str(size(em.vec.mag(A))));
fprintf('unit   input 100x3 output %s\n', mat2str(size(em.vec.unit(A))));
fprintf('angle  input 100x3 output %s\n', mat2str(size(em.vec.angle(A,B))));
fprintf('fromTo input 100x3 output %s\n', mat2str(size(em.vec.fromTo(A,B))));