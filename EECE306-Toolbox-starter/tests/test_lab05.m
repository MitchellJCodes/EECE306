function test_05 = test_lab05()
%TEST_LAB05 Tests for EECE 306 Lab 5
%
% Tests:
% 1. Enclosed point charge
% 2. Point charge outside sphere
% 3. Surface shape independence
% 4. Multiple charge superposition
% 5. Charge position independence
% 6. Divergence inside uniform volume charge
% 7. Divergence in charge-free region
% 8. Surface orientation sign

%% Constants

Q = 3e-9;

%% Test 1: Point charge inside sphere

Q = 3e-9;

s = em.src.pointCharge(Q,[0 0 0]);

FD = @(r) em.field.D(s,r);

sphere = @(th,ph) ...
    0.5*[sin(th).*cos(ph), ...
         sin(th).*sin(ph), ...
         cos(th)];

Phi = em.field.flux(FD,sphere,...
    [0 pi],[0 2*pi],60,120,'gauss');

em.test.assertClose(Phi,Q,1e-6,...
    'flux through enclosing sphere');
%% Test 2: Point charge outside sphere

s_out = em.src.pointCharge(Q,[2 0 0]);

FD_out = @(r) em.field.D(s_out,r);

Phi_out = em.field.flux(FD_out,sphere,...
    [0 pi],[0 2*pi],60,120);

% Absolute tolerance for nearly zero flux
tol_zero = 1e-12;

if abs(Phi_out) > tol_zero
    error('outside charge flux is not approximately zero');
end


%% Test 3: Sphere and cube give same enclosed flux

% Charge at origin
s = em.src.pointCharge(Q,[0 0 0]);
FD = @(r) em.field.D(s,r);

% Cube side coordinates from -0.5 to 0.5

% +x face
face_px = @(u,v) [0.5*ones(size(u)),u,v];

% -x face
face_nx = @(u,v) [-0.5*ones(size(u)),v,u];

% +y face
face_py = @(u,v) [v,0.5*ones(size(u)),u];

% -y face
face_ny = @(u,v) [u,-0.5*ones(size(u)),v];

% +z face
face_pz = @(u,v) [u,v,0.5*ones(size(u))];

% -z face
face_nz = @(u,v) [v,u,-0.5*ones(size(u))];

span = [-0.5 0.5];

Phi_cube = ...
    em.field.flux(FD,face_px,span,span,60,60) + ...
    em.field.flux(FD,face_nx,span,span,60,60) + ...
    em.field.flux(FD,face_py,span,span,60,60) + ...
    em.field.flux(FD,face_ny,span,span,60,60) + ...
    em.field.flux(FD,face_pz,span,span,60,60) + ...
    em.field.flux(FD,face_nz,span,span,60,60);

rel_cube = abs(Phi_cube-Phi)/abs(Phi);

if rel_cube > 1e-4
    error('cube and sphere flux do not agree');
end


%% Test 4: Three charges inside, two outside

s1 = em.src.pointCharge(1e-9,[0 0 0]);
s2 = em.src.pointCharge(2e-9,[0.1 0 0]);
s3 = em.src.pointCharge(3e-9,[0 0.1 0]);

s4 = em.src.pointCharge(4e-9,[2 0 0]);
s5 = em.src.pointCharge(5e-9,[0 2 0]);

smerge = em.src.merge(...
    em.src.merge(...
    em.src.merge(...
    em.src.merge(s1,s2),s3),s4),s5);

FDmerge = @(r) em.field.D(smerge,r);

Phi_merge = em.field.flux(FDmerge,sphere,...
    [0 pi],[0 2*pi],60,120,'gauss');

Q_inside = 1e-9 + 2e-9 + 3e-9;

em.test.assertClose(Phi_merge,Q_inside,1e-4,...
    'three charges inside sphere');


%% Test 5: Flux independent of charge position

positions = [ ...
     0.10  0.00  0.00;
     0.15  0.05  0.00;
    -0.10  0.10  0.05];

for k = 1:size(positions,1)

    sk = em.src.pointCharge(Q,positions(k,:));

    FDk = @(r) em.field.D(sk,r);

    Phi_k = em.field.flux(FDk,sphere,...
        [0 pi],[0 2*pi],60,120,'gauss');

    em.test.assertClose(Phi_k,Q,1e-4,...
        'flux independent of charge position');

end


%% Test 6: Divergence inside uniform volume charge

rhoV = 1e-9;

ball = @(r,th,ph) ...
    [r.*sin(th).*cos(ph), ...
     r.*sin(th).*sin(ph), ...
     r.*cos(th)];

sb = em.src.volCharge(rhoV,ball,...
    [0 0.5],[0 pi],[0 2*pi],...
    20,20,40,'midpoint');

FDb = @(r) em.field.D(sb,r);

% Point safely inside charged volume
r_inside = [0.1 0 0];

div_inside = em.op.div(FDb,r_inside,0.05);

rel_div = abs(div_inside-rhoV)/abs(rhoV);

if rel_div > 0.05
    error('divergence inside charge distribution differs by more than 5 percent');
end


%% Test 7: Divergence in charge-free region

r_free = [1 0 0];

div_free = em.op.div(FDb,r_free,0.05);

tol_div_zero = 1e-12;

if abs(div_free) > tol_div_zero
    error('divergence in charge-free region is not approximately zero');
end


%% Test 8: Reversing orientation flips flux sign

sphere_reverse = @(ph,th) ...
    0.5*[sin(th).*cos(ph), ...
         sin(th).*sin(ph), ...
         cos(th)];

Phi_forward = em.field.flux(FD,sphere,...
    [0 pi],[0 2*pi],60,120);

Phi_reverse = em.field.flux(FD,sphere_reverse,...
    [0 2*pi],[0 pi],120,60);

em.test.assertClose(Phi_reverse,-Phi_forward,1e-4,...
    'reversing surface parameters flips flux sign');


%% Test completed

test_05 = true;

end