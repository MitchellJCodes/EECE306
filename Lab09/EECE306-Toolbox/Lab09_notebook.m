%% EECE 306 Lab 9 notebook. Current Sources and the Biot-Savart Law
% *Team 03 — Carissa McWilliams, Amelia Harris, and Mitch Nazareth.*

%% Setup
clear; close all;

%% How long em.field.H took
% em.field.H is the same loop as em.field.E with a different numerator.
minutes_for_H = 70;
% The biggest challenge was taking the structure from em.field.E and
% changing it to calculate the magnetic field using the Biot-Savart law.
% I had to make sure the cross product was handled correctly and that
% the current source elements and their tangent vectors were being used
% properly. I also had to make sure the function worked with multiple
% observation points and handled points where the field calculation
% could become singular. Reusing the source and integration code from
% the previous labs helped save a lot of time. I think this would have
% gone faster if the current source structure and tangent vector
% calculations had already been set up.



%% Circular loop against the closed form
a = 1; I = 1;
ring = @(t)[a*cos(t) a*sin(t) zeros(size(t))];
sr = em.src.lineCurrent(I, ring, [0 2*pi], 400);
z = linspace(0, 2, 10)';
Hnum = em.field.H(sr, [zeros(10,2) z]);
Href = I*a^2 ./ (2*(a^2 + z.^2).^1.5);
disp('    z (m)     relative error')
disp([z abs(Hnum(:,3) - Href)./Href])
Hc = em.field.H(sr, [0 0 0]);
fprintf('center H_z = %.6f   I over 2a = %.6f\n', Hc(3), I/(2*a));
fprintf('transverse component at center = %.3e  (consistent with zero)\n', hypot(Hc(1), Hc(2)));

%% Finite segment against the infinite wire
% The same study as the Lab 4 finite line, with the field engine swapped.
% The curve has the same shape, and the handout asks why.
rho = 1;
hor = logspace(0, 2, 15)';
ratio = zeros(size(hor));
for k = 1:numel(hor)
    h = hor(k) * rho;
    seg = @(t)[zeros(size(t)) zeros(size(t)) t];
    ss = em.src.lineCurrent(I, seg, [-h/2 h/2], 800);
    Hs = em.field.H(ss, [rho 0 0]);
    ratio(k) = abs(Hs(2)) / (I/(2*pi*rho));
end
figure; semilogx(hor, ratio, 'o-'); grid on
xlabel('h / rho'); ylabel('H finite over H infinite')
title('Finite segment approaching the infinite wire')
fprintf('ratio first exceeds 0.99 at h/rho about %.1f\n', hor(find(ratio > 0.99, 1)));

%% Solenoid, built by reusing the ring
aS = 0.4; hS = 8; Nturns = 25;
zk = linspace(-hS/2, hS/2, Nturns);
srcs = cell(1, Nturns);
for k = 1:Nturns
    zc = zk(k);
    srcs{k} = em.src.lineCurrent(I, @(t)[aS*cos(t) aS*sin(t) zc*ones(size(t))], [0 2*pi], 200);
end
sol = em.src.merge(srcs{:});
zax = linspace(-6, 6, 100)';
Hax = em.field.H(sol, [zeros(100,2) zax]);
figure; plot(zax, Hax(:,3), '-'); grid on
xlabel('z (m)'); ylabel('H_z (A/m)')
title('Axial profile of a 25 turn solenoid')
Hcen = em.field.H(sol, [0 0 0]);
fprintf('center H_z = %.4f   N I over h = %.4f   ratio = %.4f\n', ...
    Hcen(3), Nturns*I/hS, Hcen(3)/(Nturns*I/hS));
fprintf('aspect ratio h/a = %.0f\n', hS/aS);

%% Reversal, linearity, and the merge policy
sneg = em.src.lineCurrent(-I, ring, [0 2*pi], 400);
fprintf('max |H(-I) + H(I)|  = %.3e  (reversal)\n', ...
    max(em.vec.mag(em.field.H(sneg, [0.3 0.2 0.5]) + em.field.H(sr, [0.3 0.2 0.5]))));
s2I = em.src.lineCurrent(2*I, ring, [0 2*pi], 400);
fprintf('max |H(2I) - 2H(I)| = %.3e  (linearity)\n', ...
    max(em.vec.mag(em.field.H(s2I, [0.3 0.2 0.5]) - 2*em.field.H(sr, [0.3 0.2 0.5]))));
sm = em.src.merge(sr, s2I);
fprintf('merged current source element count = %d  (sum of parts = %d)\n', ...
    size(sm.pos,1), size(sr.pos,1) + size(s2I.pos,1));
try
    em.src.merge(sr, em.src.pointCharge(1e-9, [0 0 0]));
    disp('ERROR, this line must not be reached')
catch err
    fprintf('merging current with charge correctly threw. %s\n', err.message);
end

%% Timing at N = 5000
pts = randn(5000,3) + 5;
tic; em.field.H(sr, pts); t = toc;
fprintf('em.field.H on 5000 points took %.3f s  (requirement, under 2 s)\n', t);

%% Interpretation
% Lab 3 made this lab easier because the field engine was already set up to calculate fields at multiple observation points. Lab 4 also helped since I could reuse the source structures and numerical integration methods instead of starting over. Most of the work for em.field.H was changing the calculation to use the Biot-Savart law rather than having to build a completely new field engine. Going forward, I would try to design interfaces so they can be reused for different types of sources and field calculations. This would make adding new features easier and avoid having to rewrite code that already works.

%% Problems encountered
% None.

%% Full test suite
runTests
