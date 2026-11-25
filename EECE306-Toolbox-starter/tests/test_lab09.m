function ok = test_lab09()
% TEST_LAB09 Verify Lab 9 magnetic field calculations.

ok = false;

%% Test 1: Circular loop on-axis field

a = 1;
I = 1;

ring = @(t)[a*cos(t) a*sin(t) zeros(size(t))];

sr = em.src.lineCurrent(I, ring, [0 2*pi], 400);

z = linspace(0, 2, 10)';

Hnum = em.field.H(sr, [zeros(10,2) z]);

Href = I*a^2 ./ (2*(a^2 + z.^2).^1.5);

relErr = abs(Hnum(:,3) - Href) ./ abs(Href);

assert(max(relErr) <= 1e-6, ...
    'Test 1 failed: Circular loop field error too large');

fprintf('Test 1 passed: max relative error = %.3e\n', ...
    max(relErr));


%% Test 2: Magnetic field at the center of the loop

Hc = em.field.H(sr, [0 0 0]);

Hexpected = I/(2*a);

relErr = abs(Hc(3) - Hexpected)/abs(Hexpected);

assert(relErr <= 1e-6, ...
    'Test 2 failed: Incorrect center magnetic field');

fprintf('Test 2 passed: Center H_z = %.6f A/m\n', ...
    Hc(3));


%% Test 3: Transverse components vanish on the axis

% Hx and Hy should be zero due to circular symmetry.
% Use absolute tolerance because the expected value is zero.

absTol = 1e-10;

Htransverse = hypot(Hnum(:,1), Hnum(:,2));

assert(max(Htransverse) <= absTol, ...
    'Test 3 failed: Transverse components are not zero');

fprintf('Test 3 passed: Maximum transverse field = %.3e\n', ...
    max(Htransverse));


%% Test 4: Finite segment approaching infinite wire

rho = 1;

% h represents the half-length of the wire
hor = logspace(0, 2, 15)';

ratio = zeros(size(hor));

Hinf = I/(2*pi*rho);

seg = @(t)[zeros(size(t)) zeros(size(t)) t];

for k = 1:numel(hor)

    h = hor(k)*rho;

    ss = em.src.lineCurrent(I, seg, [-h h], 800);

    Hs = em.field.H(ss, [rho 0 0]);

    ratio(k) = Hs(2)/Hinf;

end

idx = find(ratio > 0.99, 1);

assert(~isempty(idx), ...
    'Test 4 failed: Finite wire did not reach 99 percent');

fprintf('Test 4 passed: h/rho = %.4f, ratio = %.6f\n', ...
    hor(idx), ratio(idx));


%% Test 5: Solenoid approaching the long limit

aS = 0.4;

Nturns = 25;

aspect = [2 4 6 8 10 12 15 20];

solRatio = zeros(size(aspect));

for j = 1:numel(aspect)

    hS = aspect(j)*aS;

    % Place turns at the centers of equally spaced intervals
    dz = hS/Nturns;

    zk = -hS/2 + ((1:Nturns)-0.5)*dz;

    srcs = cell(1,Nturns);

    for k = 1:Nturns

        zc = zk(k);

        loop = @(t)[aS*cos(t) ...
                    aS*sin(t) ...
                    zc*ones(size(t))];

        srcs{k} = em.src.lineCurrent( ...
            I, loop, [0 2*pi], 200);

    end

    sol = srcs{1};

    for k = 2:numel(srcs)
        sol = em.src.merge(sol,srcs{k});
    end
    
    Hcen = em.field.H(sol, [0 0 0]);

    Hlong = Nturns*I/hS;

    solRatio(j) = Hcen(3)/Hlong;

end

idx = find(solRatio >= 0.99, 1);

assert(~isempty(idx), ...
    'Test 5 failed: Solenoid did not reach 99 percent');

fprintf('Test 5 passed: h/a = %.2f, ratio = %.6f\n', ...
    aspect(idx), solRatio(idx));


%% Test 6: Reversing current reverses the field

pts = [0.3 0.2 0.5;
       0.5 0.4 1.0;
       0.2 0.6 0.8];

sneg = em.src.lineCurrent(-I, ring, [0 2*pi], 400);

Hpos = em.field.H(sr, pts);

Hneg = em.field.H(sneg, pts);

reverseErr = abs(Hneg + Hpos);

assert(max(reverseErr(:)) <= 1e-10, ...
    'Test 6 failed: Current reversal');

fprintf('Test 6 passed: Reversal error = %.3e\n', ...
    max(reverseErr(:)));


%% Test 7: Doubling current doubles the field

s2I = em.src.lineCurrent(2*I, ring, [0 2*pi], 400);

H2I = em.field.H(s2I, pts);

linearErr = abs(H2I - 2*Hpos);

assert(max(linearErr(:)) <= 1e-10, ...
    'Test 7 failed: Current linearity');

fprintf('Test 7 passed: Linearity error = %.3e\n', ...
    max(linearErr(:)));


%% Test 8: Source merging

% Part A: Merge two current sources

sm = em.src.merge(sr, s2I);

expectedCount = size(sr.pos,1) + size(s2I.pos,1);

assert(strcmp(sm.type,'current'), ...
    'Test 8 failed: Incorrect merged source type');

assert(size(sm.pos,1) == expectedCount, ...
    'Test 8 failed: Incorrect merged element count');

assert(size(sm.Idl,1) == expectedCount, ...
    'Test 8 failed: Incorrect current element count');

% Verify superposition

Hm = em.field.H(sm, pts);

mergeErr = abs(Hm - (Hpos + H2I));

assert(max(mergeErr(:)) <= 1e-10, ...
    'Test 8 failed: Merged field violates superposition');

% Part B: Reject merging current and charge sources

sq = em.src.pointCharge(1e-9, [0 0 0]);

errorThrown = false;

try

    em.src.merge(sr, sq);

catch

    errorThrown = true;

end

assert(errorThrown, ...
    'Test 8 failed: Current and charge sources were merged');

fprintf('Test 8 passed: Source merging and type rejection\n');


%% Test 9: Performance and output dimensions

rng(9);

Nobs = 5000;

pts = randn(Nobs,3) + 5;

tic;

Hlarge = em.field.H(sr, pts);

elapsed = toc;

assert(isequal(size(Hlarge), [Nobs 3]), ...
    'Test 9 failed: Incorrect output dimensions');

assert(all(isfinite(Hlarge(:))), ...
    'Test 9 failed: Non-finite magnetic field values');

assert(elapsed < 2, ...
    'Test 9 failed: Execution time exceeded 2 seconds');

fprintf('Test 9 passed: Output size = %d x %d\n', ...
    size(Hlarge,1), size(Hlarge,2));

fprintf('Execution time = %.4f seconds\n', elapsed);


%% All tests completed

ok = true;

fprintf('All nine Lab 9 tests passed!\n');

end
