function test_03 = test_lab03()
%TEST_LAB03 Lab 3 Test Function
%   Runs required tests for Lab 3
%% Setup
clear; close all;
rng(306);
e0 = em.const.eps0();

%% Closed form comparison at ten radii
Q = 2e-9;
s = em.src.pointCharge(Q, [0 0 0]);
radii = (1:10)';
Enum  = em.field.E(s, [radii zeros(10,2)]);
Eref  = Q ./ (4*pi*e0*radii.^2);
relerr = abs(Enum(:,1) - Eref) ./ Eref;
disp('   r (m)    relative error')
disp([radii relerr])

%% Radial field direction
Qpos = em.src.pointCharge(1e-9, [0 0 0]);
Qneg = em.src.pointCharge(-1e-9, [0 0 0]);
pts_dir = randn(50,3);
Epos = em.field.E(Qpos, pts_dir);
Eneg = em.field.E(Qneg, pts_dir);
rhat = em.vec.unit(pts_dir);
dot_pos = sum(Epos.* rhat, 2);
dot_neg = sum(Eneg.* rhat, 2);
fprintf(['Minimum charge of positive charges = %.3e ' ...
    '(should be greater than 0)\n'], min(dot_pos));
fprintf(['Maximum charge of negative charges = %.3e ' ...
    '(should be less than 0)\n'], max(dot_neg));

%% Symmetry on the bisector plane
q1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
q2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
s2 = em.src.merge(q1, q2);
pts = [zeros(200,1) randn(200,2)];
Ev  = em.field.E(s2, pts);
tol = 1e-6 * max(em.vec.mag(Ev)); 
em.test.assertClose(max(abs(Ev(:,1))), 0, tol, 'bisector symmetry');
fprintf('max |Ex| on bisector      = %.3e V/m\n', max(abs(Ev(:,1))));
fprintf('max |E| on the same plane = %.3e V/m\n', max(em.vec.mag(Ev)));

%% Merge Functionality
sum_count1 = em.src.merge(q1, q2);
fprintf('element count after first merge = %d (expect 2)\n', ...
    size(sum_count1.pos, 1))
crnt_struct = struct('type', 'current', 'pos', [0 0 0], 'I', 1);
try
    try_merge = em.src.merge(q1, crnt_struct); %#ok<*NASGU>
    fprintf(['Mismatch test failed: em.src.merge(q1, crnt_struct) ' ...
        'failed to raise an error.\n']);
catch
    fprintf(['Mismatch test passed: em.src.merge(q1, crnt_struct) ' ...
        'successfully raised an error.\n']);
end


%% Dipole far field
sd = em.src.merge(em.src.pointCharge(1e-9,[0 0 0.05]), ...
    em.src.pointCharge(-1e-9,[0 0 -0.05]));
E1 = em.vec.mag(em.field.E(sd, [0 0 10]));
E2 = em.vec.mag(em.field.E(sd, [0 0 20]));
fprintf('|E(r)| / |E(2r)| = %.4f   (expect near 8 for a dipole)\n', E1/E2);

%% The picture, raw and normalized
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', false, 'title', ['Two equal charges, ' ...
    'raw magnitudes']));
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', true, 'title', ['Two equal charges, ' ...
    'direction only']));

%% Evaluating the field at the point charge
try
    Einf = em.field.E(q, [0 0 0]);
    fprintf(['Singularity test failed: em.field.E(q, [0 0 0]) failed ' ...
        'to raise an error.\n']);
catch
    fprintf(['Singularity test passed: em.field.E(q, [0 0 0]) ' ...
        'successfully raised an error.\n']);
end

%% Timing at N = 5000
pts5k = randn(5000,3) + 5;
tic; em.field.E(s2, pts5k); t = toc;
fprintf(['em.field.E on 5000 points took %.3f s  (requirement, ' ...
    'under 2 s)\n'], t);

%% Field Engine Limits
s = em.src.pointCharge(1, [0 0 0]);
sec_max = 5.0;
for N = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8]
    r = randn(N, 3);
    tic
    Ev = em.field.E(s, r);
    t = toc;
    if t > sec_max
        fprintf(['Compution time surpasses time threshold at N=%d points: ' ...
            'took %.3f seconds\n'], N, t);
        break;
    end
    fprintf(['Observation points: N = %d, Source Charges: M = %d, Product = ' ...
        '%d, Time = %.3f s\n'], N, size(s.pos, 1), N*size(s.pos, 1), t);
end