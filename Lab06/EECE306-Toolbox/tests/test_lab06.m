function test_06 = test_lab06()
%TEST_LAB06 Lab 6 checks.
%   Tests electric potential, gradient, circulation, and contour plotting.
%

clear; close all;

%% Setup

tol = 1e-10;
e0 = em.const.eps0();

Q = 2e-9;
s = em.src.pointCharge(Q,[0 0 0]);

%% 1. Potential of a point charge

r = [1 0 0;
     0 2 0;
     0 0 3;
     1 2 3];

V = em.field.V(s,r);

Rmag = sqrt(sum(r.^2,2));
Vref = Q ./ (4*pi*e0*Rmag);

em.test.assertClose(V,Vref,tol, ...
    'Point charge potential');

%% 2. Potential obeys superposition

s1 = em.src.pointCharge(1e-9,[-0.1 0 0]);
s2 = em.src.pointCharge(2e-9,[ 0.1 0 0]);
sm = em.src.merge(s1,s2);

rtest = [0.3 0.2 0.1];

Vsum = em.field.V(sm,rtest);

V1 = em.field.V(s1,rtest);
V2 = em.field.V(s2,rtest);

em.test.assertClose(Vsum,V1+V2,tol, ...
    'Potential obeys superposition');

%% 3. Potential rejects observation point at source

caught = false;

try
    em.field.V(s,[0 0 0]);
catch
    caught = true;
end

assert(caught, ...
    'V must reject observation points at source locations');

%% 4. Gradient of a linear scalar field

% f(x,y,z) = 2x - 3y + 4z
% grad(f) = [2 -3 4]

flinear = @(r) ...
    2*r(:,1) - 3*r(:,2) + 4*r(:,3);

r0 = [0.37 -0.21 0.44];

g = em.op.grad(flinear,r0,1e-4);

gref = [2 -3 4];

em.test.assertClose(g,gref,1e-8, ...
    'Gradient of linear scalar field');

%% 5. Gradient of a constant scalar field

fconst = @(r) ...
    7*ones(size(r,1),1);

gzero = em.op.grad(fconst,[0.2 0.3 0.4],1e-4);

em.test.assertClose(gzero,[0 0 0],1e-10, ...
    'Gradient of constant scalar field');

%% 6. Electric field equals negative gradient of potential

p = [1 0.5 0];

E = em.field.E(s,p);
Eg = -em.op.grad(@(r) em.field.V(s,r),p,1e-5);

em.test.assertClose(E,Eg,1e-5, ...
    'E equals negative gradient of V');

%% 7. Circulation of a constant field along a closed path is zero

Fconst = @(r) ...
    repmat([2 -1 3],size(r,1),1);

circle = @(t) ...
    [2*cos(t) 2*sin(t) zeros(size(t))];

C = em.field.circulation( ...
    Fconst,circle,[0 2*pi],400);

em.test.assertClose(C,0,1e-10, ...
    'Circulation of constant field around closed path');

%% 8. Circulation of electrostatic field is zero

loops = { ...
    @(t)[2*cos(t) 2*sin(t) zeros(size(t))], ...
    @(t)[1+0.5*cos(t) 0.5*sin(3*t) 0.3*sin(t)], ...
    @(t)[0.5*cos(t) 1+0.8*sin(t) 0.2*cos(2*t)]};

for k = 1:3
    C = em.field.circulation( ...
        @(r) em.field.E(s,r), ...
        loops{k},[0 2*pi],400);

    em.test.assertClose(C,0,1e-8, ...
        'Circulation of electrostatic field');
end

%% 9. Non-conservative control field has nonzero circulation

% F = [-y, x, 0]
% Around x = 2 cos(t), y = 2 sin(t):
% integral F . dr = 8*pi

Fnc = @(r) ...
    [-r(:,2) r(:,1) zeros(size(r,1),1)];

Cnc = em.field.circulation( ...
    Fnc,circle,[0 2*pi],400);

em.test.assertClose(Cnc,8*pi,1e-8, ...
    'Circulation of nonconservative control field');

%% 10. Potential difference equals negative line integral of E

pa = [1 0 0];
pb = [0.3 0.4 0];

dV = em.field.V(s,pb) - em.field.V(s,pa);

seg = @(t) pa + t.*(pb-pa);

work = -em.field.circulation( ...
    @(r) em.field.E(s,r), ...
    seg,[0 1],400);

em.test.assertClose(work,dV,1e-8, ...
    'Potential difference equals negative line integral of E');

%% 11. contourV creates a figure and axes

figBefore = findall(0,'Type','figure');

em.viz.contourV( ...
    @(r) em.field.V(s,r), ...
    [-2 2],[-2 2],30, ...
    struct('title','Test potential'));

figAfter = findall(0,'Type','figure');

assert(numel(figAfter) > numel(figBefore), ...
    'contourV should create a figure');

ax = findall(gcf,'Type','axes');

assert(~isempty(ax), ...
    'contourV should create axes');

disp('  test_lab06 checks complete');

end
