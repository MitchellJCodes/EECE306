function test_10 = test_lab10()
%TEST_LAB10 Lab 10: Ampere's law, magnetic Gauss law, curl, and field lines.
%
% The test deliberately reuses the Lab 5 flux integrator, Lab 6
% circulation integrator, and Lab 9 magnetic field engine unchanged.

clear; close all;

%% 1-2. Ampere's law for the four coaxial-cable regions
I = 1;
a = 0.5;
b = 1.0;
c = 1.3;
coax = makeCoax(I,a,b,c,0.1,40,400);
FH = @(r) em.field.H(coax,r);

rhoTest = [0.25; 0.75; 1.18; 1.5];
Iref = zeros(size(rhoTest));
Icalc = zeros(size(rhoTest));
for k = 1:numel(rhoTest)
    rho = rhoTest(k);
    loop = @(t) [rho*cos(t), rho*sin(t), zeros(size(t))];
    Icalc(k) = em.field.circulation(FH,loop,[0 2*pi],96);
    if rho <= a
        Iref(k) = I*rho^2/a^2;
    elseif rho < b
        Iref(k) = I;
    elseif rho <= c
        Iref(k) = I*(1-(rho^2-b^2)/(c^2-b^2));
    else
        Iref(k) = 0;
    end
end

% The filament bundle approximates a uniform volume current. The three
% nonzero-region values are required to agree within 5% of total current.
for k = 1:3
    assert(abs(Icalc(k)-Iref(k)) < 0.05*I, ...
        sprintf('Ampere circulation region %d: got %.6f A, expected %.6f A', ...
        k,Icalc(k),Iref(k)));
end

% Outside the return-current shield, the enclosed current is zero.
assert(abs(Icalc(4)) < 0.05, ...
    sprintf('Outside-coax circulation %.6g A is not consistent with zero',Icalc(4)));

fprintf('  Coax circulation [inside, dielectric, shield, outside] = ');
fprintf('%.4f ',Icalc);
fprintf('A\n');

%% 3. Circular, square, and irregular paths enclose the same current
rho = 0.75;
circle = em.field.circulation(FH, ...
    @(t)[rho*cos(t),rho*sin(t),zeros(size(t))],[0 2*pi],96);

half = 0.65;
sides = { ...
    @(t)[t,-half+zeros(size(t)),zeros(size(t))], ...
    @(t)[half+zeros(size(t)),t,zeros(size(t))], ...
    @(t)[-t,half+zeros(size(t)),zeros(size(t))], ...
    @(t)[-half+zeros(size(t)),-t,zeros(size(t))]};
square = 0;
for k = 1:4
    square = square + em.field.circulation(FH,sides{k},[-half half],48);
end

irregular = em.field.circulation(FH, ...
    @(t)[0.7*(1+0.1*sin(3*t)).*cos(t), ...
         0.68*(1+0.08*cos(2*t)).*sin(t),zeros(size(t))], ...
    [0 2*pi],128);

pathValues = [circle square irregular];
assert(max(abs(pathValues-I))/I < 0.01, ...
    sprintf('Ampere path-shape agreement exceeds 1%%: [%.6f %.6f %.6f] A', ...
    pathValues(1),pathValues(2),pathValues(3)));
fprintf('  Path circulations [circle square irregular] = %.5f %.5f %.5f A\n', ...
    pathValues(1),pathValues(2),pathValues(3));

%% 4. Magnetic Gauss law: two sources and three closed surfaces
% Both sources are wholly inside the sphere, upper hemisphere, and cube.
wire = em.src.lineCurrent(1, ...
    @(t)[t,0.15+zeros(size(t)),0.2+zeros(size(t))],[-0.3 0.3],100,'gauss');
loop = em.src.lineCurrent(1, ...
    @(t)[0.2*cos(t),0.2*sin(t),0.3+zeros(size(t))], ...
    [0 2*pi],128,'gauss');

sources = {wire,loop};
sourceNames = {'finite wire','circular loop'};
R = 0.8;
for isrc = 1:2
    FB = @(r) em.field.B(sources{isrc},r);
    scale = em.vec.mag(FB([0 0 0]))*pi*R^2;

    phiSphere = em.field.flux(FB, ...
        @(th,ph) R*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)], ...
        [0 pi],[0 2*pi],36,72);

    phiHemisphere = closedHemisphereFlux(FB,R);

    phiCube = cubeFlux(FB,R);

    fluxValues = [phiSphere phiHemisphere phiCube];
    assert(all(abs(fluxValues) < max(1e-12,0.05*scale)), ...
        sprintf('Magnetic Gauss law failed for %s: fluxes [%.3e %.3e %.3e], scale %.3e', ...
        sourceNames{isrc},fluxValues(1),fluxValues(2),fluxValues(3),scale));

    fprintf('  B flux, %s [sphere hemisphere+disk cube] = ',sourceNames{isrc});
    fprintf('%.3e ',fluxValues);
    fprintf('Wb (reference scale %.3e Wb)\n',scale);
end

%% 5. Magnetic Gauss law when the closed surface cuts the current loop
loopCut = em.src.lineCurrent(1, ...
    @(t)[0.2*cos(t),0.2*sin(t),0.3+zeros(size(t))], ...
    [0 2*pi],128,'gauss');
FBcut = @(r) em.field.B(loopCut,r);
Rcut = 0.12;
cutSphere = @(th,ph) [0.2+Rcut*sin(th).*cos(ph), ...
                      Rcut*sin(th).*sin(ph), ...
                      0.3+Rcut*cos(th)];
phiCut = em.field.flux(FBcut,cutSphere,[0 pi],[0 2*pi],48,96);
cutScale = em.vec.mag(FBcut([0.2 0.2 0.3]))*pi*Rcut^2;
assert(abs(phiCut) < max(1e-12,0.10*cutScale), ...
    sprintf('Magnetic flux through loop-cutting sphere %.3e Wb exceeds scale %.3e Wb', ...
    phiCut,cutScale));
fprintf('  B flux through sphere cutting current loop = %.3e Wb (scale %.3e Wb)\n', ...
    phiCut,cutScale);

%% 6-7. Point form of Ampere's law: curl H = J inside, zero outside
J = I/(pi*a^2);
curlInside = em.op.curl(FH,[0.25 0 0],0.15);
curlOutside = em.op.curl(FH,[0.75 0 0],0.1);

assert(abs(curlInside(3)-J)/J < 0.05, ...
    sprintf('curl H inside conductor %.6g differs from J %.6g by more than 5%%', ...
    curlInside(3),J));
assert(norm(curlInside(1:2)) < 0.05*J, ...
    'curl H inside the z-directed conductor should have negligible x/y components');
assert(norm(curlOutside) < 0.05*J, ...
    sprintf('curl H in the current-free dielectric is not near zero (%.6g)', ...
    norm(curlOutside)));

fprintf('  curl H inside = [%.4f %.4f %.4f], Jz = %.4f A/m^2\n', ...
    curlInside(1),curlInside(2),curlInside(3),J);
fprintf('  curl H in dielectric magnitude = %.3e A/m^2\n',norm(curlOutside));

%% 8. Magnetic divergence is zero inside and outside the conductor
FBcoax = @(r) em.field.B(coax,r);
divInside = em.op.div(FBcoax,[0.25 0 0],1e-4);
divOutside = em.op.div(FBcoax,[0.75 0 0],1e-4);
assert(abs(divInside) < 1e-8, ...
    sprintf('div B inside conductor is %.3e, not near zero',divInside));
assert(abs(divOutside) < 1e-8, ...
    sprintf('div B outside conductor is %.3e, not near zero',divOutside));

%% 9. Curl of an electrostatic field is zero
charge = em.src.pointCharge(2e-9,[0 0 0]);
FE = @(r) em.field.E(charge,r);
curlE = em.op.curl(FE,[0.8 0.4 0.6; -0.7 0.5 0.9],1e-4);
assert(max(abs(curlE(:))) < 1e-8, ...
    'Curl of an electrostatic point-charge field should vanish.');

% Also verify the curl operator against a field with a known nonzero curl.
Fknown = @(r)[3*r(:,2),4*r(:,3),5*r(:,1)];
knownCurl = em.op.curl(Fknown,[0.1 0.2 0.3; -0.4 0.5 0.6],1e-4);
em.test.assertClose(knownCurl,repmat([-4 -5 -3],2,1),1e-8, ...
    'Central-difference curl of a linear vector field');

%% Field-line visualization: one wire and an antiparallel pair
theta = (0:30:330)';
seedsOne = [0.3+0.45*cos(theta),0.45*sin(theta)];
fig = figure('visible','off');
linesOne = em.viz.streamlines(@(r) straightWireField(r,0.3,1), ...
    [-1 1],[-1 1],seedsOne, ...
    struct('stepSize',0.01,'title','One straight wire'));
assert(all(cellfun(@isClosedLine,linesOne)), ...
    'Every single-wire streamline should close on itself.');
close(fig);

seedsPair = [[0.3+0.12*cos(theta),0.12*sin(theta)]; ...
             [-0.3+0.12*cos(theta),0.12*sin(theta)]];
fig = figure('visible','off');
linesPair = em.viz.streamlines(@(r) ...
    straightWireField(r,0.3,1)+straightWireField(r,-0.3,-1), ...
    [-1.2 1.2],[-1.2 1.2],seedsPair, ...
    struct('stepSize',0.008,'title','Antiparallel wires'));
assert(all(cellfun(@isClosedLine,linesPair)), ...
    'Every antiparallel-wire streamline should close on itself.');
close(fig);

disp('  test_lab10 checks complete');
test_10 = true;
end

function s = makeCoax(I,a,b,c,spacing,L,Nline)
[x,y] = meshgrid(-1.3:spacing:1.3);
r = hypot(x,y);
inner = r <= a+1e-10;
shell = r >= b-1e-10 & r <= c+1e-10;
xin = x(inner); yin = y(inner);
xsh = x(shell); ysh = y(shell);
parts = cell(1,numel(xin)+numel(xsh));
for k = 1:numel(xin)
    x0 = xin(k); y0 = yin(k);
    parts{k} = em.src.lineCurrent(I/numel(xin), ...
        @(t)[x0+zeros(size(t)),y0+zeros(size(t)),t],[-L L],Nline);
end
for k = 1:numel(xsh)
    x0 = xsh(k); y0 = ysh(k);
    parts{numel(xin)+k} = em.src.lineCurrent(-I/numel(xsh), ...
        @(t)[x0+zeros(size(t)),y0+zeros(size(t)),t],[-L L],Nline);
end
s = em.src.merge(parts{:});
end

function phi = closedHemisphereFlux(F,R)
top = em.field.flux(F, ...
    @(th,ph) R*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)], ...
    [0 pi/2],[0 2*pi],28,56);
% The disk is parameterized with reversed azimuth so its normal is -z,
% closing the outward-oriented upper hemisphere.
disk = @(u,v)[u.*cos(v),-u.*sin(v),zeros(size(u))];
bottom = em.field.flux(F,disk,[0 R],[0 2*pi],24,56);
phi = top+bottom;
end

function phi = cubeFlux(F,R)
faces = { ...
    @(u,v)[u v  R*ones(size(u))], ...
    @(u,v)[v u -R*ones(size(u))], ...
    @(u,v)[v  R*ones(size(u)) u], ...
    @(u,v)[u -R*ones(size(u)) v], ...
    @(u,v)[ R*ones(size(u)) u v], ...
    @(u,v)[-R*ones(size(u)) v u]};
phi = 0;
for k = 1:6
    phi = phi + em.field.flux(F,faces{k},[-R R],[-R R],24,24);
end
end

function F = straightWireField(r,x0,I)
dx = r(:,1)-x0;
dy = r(:,2);
r2 = dx.^2+dy.^2;
F = I/(2*pi)*[-dy./r2,dx./r2,zeros(size(dx))];
end

function tf = isClosedLine(line)
tf = size(line,1) >= 20 && norm(line(1,:)-line(end,:)) < 1e-12;
end
