%% EECE 306 Lab 6 notebook. Electric Potential, the Gradient, and Path Independence
% *Team 03 Carissa McWilliams, Amelia Harris, and Mitch Nazareth.* 
%
% Fill in every TODO, run |publish('Lab06_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
e0 = em.const.eps0();
Q  = 2e-9;
s  = em.src.pointCharge(Q, [0 0 0]);
FE = @(r) em.field.E(s, r);
FV = @(r) em.field.V(s, r);

%% Potential against the closed form
radii = (1:10)';
Vnum = em.field.V(s, [radii zeros(10,2)]);
Vref = Q ./ (4*pi*e0*radii);
disp('   r (m)    relative error')
disp([radii abs(Vnum - Vref)./Vref])

%% E from the gradient, and the error map
% E computed directly and E recovered as minus the gradient of V must
% agree away from the source.
[X, Y] = meshgrid(linspace(0.3, 2, 40), linspace(-1, 1, 40));
P = [X(:) Y(:) zeros(numel(X),1)];
Ed = em.field.E(s, P);
Eg = -em.op.grad(FV, P);
errmap = reshape(em.vec.mag(Ed - Eg) ./ em.vec.mag(Ed), size(X));
figure; contourf(X, Y, log10(errmap), 20); colorbar
xlabel('x (m)'); ylabel('y (m)')
title('log10 relative error of E recovered from  -grad V')
% Error is largest near the edge of the region where the gradient
% is calculated because the gradient is approximated using finite
% difference. This creates relative error that gets bigger when
% the electric field is smaller.

%% Choosing the step size, the V shaped curve
% Truncation error falls as h squared, round off grows as 1 over h. The
% total has a minimum. Sweep h and find it for this machine.
p0 = [1 0.5 0];
Eexact = em.field.E(s, p0);
hs = logspace(-12, -1, 23)';
errh = zeros(size(hs));
for k = 1:numel(hs)
    Eh = -em.op.grad(FV, p0, hs(k));
    errh(k) = em.vec.mag(Eh - Eexact) / em.vec.mag(Eexact);
end
figure; loglog(hs, errh, 'o-'); grid on
xlabel('step h (m)'); ylabel('relative error of  -grad V')
title('The two error sources compete, the total has a minimum')
[emin, imin] = min(errh);
fprintf('best h on this machine about %.1e with error %.1e\n', hs(imin), emin);
fprintf('the API default h = 1e-5 sits on the flat bottom of this curve\n');

%% The field is conservative, with a control
loops = { @(t)[2*cos(t) 2*sin(t) zeros(size(t))], ...
          @(t)[1+0.5*cos(t) 0.5*sin(3*t) 0.3*sin(t)], ...
          @(t)[0.5*cos(t) 1+0.8*sin(t) 0.2*cos(2*t)] };
for k = 1:3
    C = em.field.circulation(FE, loops{k}, [0 2*pi], 400);
    fprintf('circulation of E on path %d = %.3e\n', k, C);
end
Fnc = @(r) [-r(:,2) r(:,1) zeros(size(r,1),1)];
Cnc = em.field.circulation(Fnc, loops{1}, [0 2*pi], 400);
fprintf('circulation of the control field = %.4f  (expect 8 pi, clearly nonzero)\n', Cnc);
% The control matters. A test that only ever confirms zero cannot tell
% correct code from code that always returns zero.

%% Potential difference two ways
pa = [1 0 0]; pb = [0.3 0.4 0];
dV = em.field.V(s, pb) - em.field.V(s, pa);
seg = @(t) pa + t.*(pb - pa);
work = -em.field.circulation(FE, seg, [0 1], 400);
fprintf('V(b) - V(a) from em.field.V        = %.6e V\n', dV);
fprintf('minus line integral of E, a to b   = %.6e V\n', work);

%% Seeing the potential
figure;
em.viz.contourV(@(r) em.field.V(s, r), [-2 2], [-2 2], 60, struct('title', 'Equipotentials of a point charge'));
hold on
em.viz.streamlines(FE, [-2 2], [-2 2], [cosd(0:45:315)' sind(0:45:315)' zeros(8,1)]*0.2, struct());
title('Equipotentials with field lines overlaid')

%% Interpretation
% The step size study implies that there is a size of h
% small enough to receive high accuracy without
% driving error higher due to round-off error.
% This error is caused by a loss in numerical precision
% when calculating the difference between two very close
% potential values. Driving h to be smaller is not always
% safer because eventually the two numerical values become
% close enough that the difference between them loses precision.

%% Problems encountered
% At first, I used midpoint quadrature and numerical derivative, but quickly found that my relative error was larger than my expected relative error. Switching to guass's equations resolved the issue.

%% Full test suite
runTests
