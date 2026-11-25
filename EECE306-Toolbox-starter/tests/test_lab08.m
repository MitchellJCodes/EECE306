function test_lab08
%% Parallel plate capacitor

n = 81;
G = em.solve.grid2d([0 1], [0 1], n, n);

iy1 = round(0.25*(n-1)) + 1;
iy2 = round(0.75*(n-1)) + 1;

mask = false(n,n);
mask(iy1,:) = true;
mask(iy2,:) = true;

Vfix = zeros(n,n);
Vfix(iy2,:) = 1;

opts = struct('method', 'direct');
C = em.derive.capacitance(G, 1, mask, Vfix, opts);

d = 0.5;
A = 1;  % Width per unit depth
e0 = em.const.eps0();
Cref = e0*A/d;

relErr = abs(C-Cref)/Cref;
fprintf('Calculated C = %.12e F\n', C);
fprintf('Expected C   = %.12e F\n', Cref);
fprintf('Signed error = %.4f%%\n', 100*(C-Cref)/Cref);
fprintf('Test 1 error: %.4f%%\n', relErr*100);
assert(relErr <= 0.01, 'Test 1 failed');
%% Two dielectric layers in series
Ns = [42 82 162];
Ch = zeros(size(Ns));
Ca = zeros(size(Ns));
Crefs = zeros(size(Ns));

for k = 1:numel(Ns)

    n2 = Ns(k);
    G2 = em.solve.grid2d([0 1], [0 1], n2, n2);

    % Find electrode positions
    [~, iy1] = min(abs(G2.y - 0.25));
    [~, iy2] = min(abs(G2.y - 0.75));

    % Define electrode boundaries
    mask2 = false(n2,n2);
    mask2(iy1,:) = true;
    mask2(iy2,:) = true;

    Vfix2 = zeros(n2,n2);
    Vfix2(iy2,:) = 1;

    % Define dielectric layers
    epsr2 = ones(n2,n2);
    epsr2(G2.Y < 0.5) = 4;

    % Analytical capacitance
    d1 = 0.5 - G2.y(iy1);
    d2 = G2.y(iy2) - 0.5;

    Crefs(k) = e0 / (d1/4 + d2);

    % Harmonic and arithmetic averaging
    optsH = struct('method','direct','avgMethod','harmonic');
    optsA = struct('method','direct','avgMethod','arithmetic');

    Ch(k) = em.derive.capacitance(G2,epsr2,mask2,Vfix2,optsH);
    Ca(k) = em.derive.capacitance(G2,epsr2,mask2,Vfix2,optsA);

    % Print results for each grid
    relErr2 = abs(Ch(k)-Crefs(k))/Crefs(k);

    fprintf('Grid size: %d\n', n2);
    fprintf('Test 2 Numerical C = %.6e\n', Ch(k));
    fprintf('Test 2 Expected C  = %.6e\n', Crefs(k));
    fprintf('Test 2 Relative Error = %.4f%%\n\n', relErr2*100);

end

% Check harmonic capacitance accuracy
relErr2 = abs(Ch-Crefs)./Crefs;

assert(max(relErr2) <= 0.01, 'Test 2 failed');

% Calculate arithmetic convergence order
ea = abs(Ca-Crefs);

pa = log(ea(1:end-1)./ea(2:end)) ./ ...
     log((Ns(2:end)-1)./(Ns(1:end-1)-1));

pA = mean(pa);

fprintf('Arithmetic convergence order = %.3f\n', pA);
%% Coaxial capacitance

n = 201;
Gc = em.solve.grid2d([-1 1], [-1 1], n, n);

r = hypot(Gc.X,Gc.Y);
a = 0.3;
b = 0.9;

maskC = r <= a | r >= b;

VfixC = zeros(n,n);
VfixC(r <= a) = 1;

Cc = em.derive.capacitance(Gc, 1, maskC, VfixC, opts);
CcRef = 2*pi*e0/log(b/a);

err4 = abs(Cc-CcRef)/CcRef;

fprintf('Test 3 error: %.4f%%\n', 100*err4);
assert(err4 <= 0.02, 'Test 3 failed');
%% Voltage independence

VfixC100 = zeros(n,n);
VfixC100(r <= a) = 100;

Cc100 = em.derive.capacitance( ...
    Gc, 1, maskC, VfixC100, opts);

err5 = abs(Cc100-Cc)/abs(Cc);

fprintf('Test 4 error: %.4e\n', err5);
assert(err5 <= 1e-8, 'Test 4 failed');
%% Dielectric scaling

epsr = 3;

Cfill = em.derive.capacitance( ...
    Gc, epsr, maskC, VfixC, opts);

ratio = Cfill/Cc;

fprintf('Test 5 ratio: %.6f\n', ratio);
assert(abs(ratio-3)/3 <= 1e-8, 'Test 5 failed');
%% Resistance of a uniform bar

n = 81;
Gr = em.solve.grid2d([0 1], [0 1], n, n);

maskR = false(n,n);
maskR(:,1) = true;
maskR(:,end) = true;

VfixR = zeros(n,n);
VfixR(:,1) = 1;

sigma = 100;

Rbar = em.derive.resistance( ...
    Gr, sigma, maskR, VfixR, opts);

Rref = 1/(sigma*1);

err7 = abs(Rbar-Rref)/Rref;

fprintf('Test 6 error: %.4f%%\n', 100*err7);
assert(err7 <= 0.01, 'Test 6 failed');
end