
function test_06 = test_lab06()
%TEST_LAB06 Tests for electric potential, gradient,
% circulation, and potential visualization.

%% Constants
eps0 = em.const.eps0();
Q = 2e-9;
k = 1/(4*pi*eps0);

%% Test 1: Potential of a point charge

s = em.src.pointCharge(Q,[0 0 0]);

radii = linspace(0.1,2,10).';
r = [radii,zeros(10,2)];

V_calc = em.field.V(s,r);
V_expected = k*Q./radii;

em.test.assertClose(V_calc,V_expected,1e-6,...
    'Point charge potential');

fprintf('Test 1 passed: Point charge potential\n');


%% Additional Test: Potential obeys superposition

s1 = em.src.pointCharge(1e-9,[-0.1 0 0]);
s2 = em.src.pointCharge(2e-9,[0.1 0 0]);

sm = em.src.merge(s1,s2);

rtest = [0.3 0.2 0.1];

Vsum = em.field.V(sm,rtest);
V1 = em.field.V(s1,rtest);
V2 = em.field.V(s2,rtest);

em.test.assertClose(Vsum,V1+V2,1e-6,...
    'Potential obeys superposition');

fprintf('Additional test passed: Superposition\n');


%% Test 2: Electric field equals negative gradient

s1 = em.src.pointCharge(3e-9,[-0.2 0 0]);
s2 = em.src.pointCharge(-2e-9,[0.2 0 0]);

smerge = em.src.merge(s1,s2);

r = [0.5 0.5 0.2;
     1.0 0.5 0.3;
    -0.5 0.8 0.4];

E_direct = em.field.E(smerge,r);

Vfun = @(r) em.field.V(smerge,r);
E_grad = -em.op.grad(Vfun,r,1e-5);

em.test.assertClose(E_grad,E_direct,1e-4,...
    'Electric field and negative gradient');

fprintf('Test 2 passed: E = -grad(V)\n');


%% Test 3: Circulation on 3 paths and 3 sources

% Source 1: Point charges
s_points = smerge;

% Source 2: Finite line charge
line = @(t) [zeros(size(t)),zeros(size(t)),t];

s_line = em.src.lineCharge(1e-9,line,...
    [-0.3 0.3],100,'midpoint');

% Source 3: Uniform volume charge
ball = @(r,th,ph) ...
    [r.*sin(th).*cos(ph), ...
     r.*sin(th).*sin(ph), ...
     r.*cos(th)];

s_vol = em.src.volCharge(1e-6,ball,...
    [0 0.3],[0 pi],[0 2*pi],...
    10,20,10,'midpoint');

sources = {s_points,s_line,s_vol};

% Path 1: Circle
circle = @(t) ...
    [1.5*cos(t),1.5*sin(t),zeros(size(t))];

% Path 2: Square
square = @(t) squarePath06(t);

% Path 3: Nonplanar closed curve
custom = @(t) ...
    [0.4+1.4*cos(t), ...
     0.3+1.1*sin(t), ...
     0.5+0.2*cos(t)+0.3*sin(t)+0.1*sin(2*t)];

paths = {circle,square,custom};
spans = {[0 2*pi],[0 4],[0 2*pi]};

for i = 1:3

    F = @(r) em.field.E(sources{i},r);

    for j = 1:3

        C = em.field.circulation(F,paths{j},...
            spans{j},400,'midpoint');

        % Estimate characteristic field-line integral
        [pos,w] = em.quad.line(paths{j},...
            spans{j},400,'midpoint');

        pathLength = sum(w);
        scale = max(em.vec.mag(F(pos))) * pathLength;

        fprintf('Source %d, Path %d: C = %.6e\n',...
            i,j,C);

        if abs(C) > 1e-3*scale
            error('Circulation is not approximately zero');
        end

    end
end

fprintf('Test 3 passed: Nine closed-path circulations\n');


%% Test 4: Non-conservative field

Fbad = @(r) ...
    [-r(:,2),r(:,1),zeros(size(r,1),1)];

Cbad = em.field.circulation(Fbad,circle,...
    [0 2*pi],400,'midpoint');

C_expected = 2*pi*(1.5^2);

em.test.assertClose(Cbad,C_expected,1e-4,...
    'Non-conservative field circulation');

fprintf('Test 4 passed: Non-conservative field\n');


%% Test 5: Potential difference and line integral

s = em.src.pointCharge(Q,[0 0 0]);

F = @(r) em.field.E(s,r);

A = [0.5 0.2 0.1];
B = [1.0 0.4 0.2];

V_A = em.field.V(s,A);
V_B = em.field.V(s,B);

deltaV = V_B-V_A;

% Straight path from A to B
curve = @(t) ...
    [A(1)+(B(1)-A(1)).*t, ...
     A(2)+(B(2)-A(2)).*t, ...
     A(3)+(B(3)-A(3)).*t];

lineIntegral = em.field.circulation(F,curve,...
    [0 1],200,'gauss');

em.test.assertClose(deltaV,-lineIntegral,1e-5,...
    'Potential difference and line integral');

fprintf('Test 5 passed: Potential difference\n');



%% Test 6: Growing-source warning

% Create a finite line charge
lambda = 1e-9;
line = @(t) [zeros(size(t)), zeros(size(t)), t];

s_growing = em.src.lineCharge(lambda,line,...
    [-5 5],100,'midpoint');

% Indicate that the source extent is growing
s_growing.growingExtent = true;

% Reset the last warning
lastwarn('');

% Evaluate potential
V_growing = em.field.V(s_growing,[1 0 0]);

% Check whether the correct warning occurred
[msg,id] = lastwarn;

if strcmp(id,'em:field:InfiniteSource')
    fprintf('Test 6 passed: Growing-source warning\n');
else
    error('Test 6 failed: No growing-source warning detected');
end



%% Visualization: Dipole potential and electric field

s1 = em.src.pointCharge(1e-9,[-0.5 0 0]);
s2 = em.src.pointCharge(-1e-9,[0.5 0 0]);

sdipole = em.src.merge(s1,s2);

Vfun = @(r) em.field.V(sdipole,r);
Efun = @(r) em.field.E(sdipole,r);

% Use an even grid size to avoid evaluating exactly
% at the point-charge singularities.
em.viz.contourV(Vfun,[-2 2],[-2 2],60);

hold on;

opts.normalize = true;
em.viz.quiver2(Efun,[-2 2],[-2 2],16,opts);

title('Dipole Potential and Electric Field');

fprintf('Visualization completed\n');


%% Tests completed

% False until all six required tests are implemented.
test_06 = true;

fprintf('All six Lab 6 tests passed!\n');

end


%% Helper function: Square path

function pos = squarePath06(t)
% Four counterclockwise sides of a square
% centered at the origin, side length 3.

t = t(:);
pos = zeros(numel(t),3);

i1 = t >= 0 & t < 1;
i2 = t >= 1 & t < 2;
i3 = t >= 2 & t < 3;
i4 = t >= 3 & t <= 4;

pos(i1,:) = [-1.5+3*t(i1), ...
              -1.5*ones(sum(i1),1), ...
              zeros(sum(i1),1)];

pos(i2,:) = [1.5*ones(sum(i2),1), ...
              -1.5+3*(t(i2)-1), ...
              zeros(sum(i2),1)];

pos(i3,:) = [1.5-3*(t(i3)-2), ...
              1.5*ones(sum(i3),1), ...
              zeros(sum(i3),1)];

pos(i4,:) = [-1.5*ones(sum(i4),1), ...
              1.5-3*(t(i4)-3), ...
              zeros(sum(i4),1)];

end
