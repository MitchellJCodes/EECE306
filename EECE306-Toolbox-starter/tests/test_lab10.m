function passed = test_lab10()

%% Construct coaxial cable

I = 1;
a = 0.5;
b = 0.8;
c = 1.2;

spacing = 0.1;
L = 40;
Nwire = 400;

[X,Y] = meshgrid((-12:12)*spacing);
R = sqrt(X.^2 + Y.^2);

inner = R <= a + 1e-12;
shell = R >= b - 1e-12 & R <= c + 1e-12;

innerPos = [X(inner), Y(inner)];
shellPos = [X(shell), Y(shell)];

Ni = size(innerPos,1);
Ns = size(shellPos,1);

srcs = cell(Ni+Ns,1);

% Inner conductor: positive current
for k = 1:Ni

    x0 = innerPos(k,1);
    y0 = innerPos(k,2);

    wire = @(t) [x0+zeros(size(t)), ...
                 y0+zeros(size(t)), t];

    srcs{k} = em.src.lineCurrent( ...
        I/Ni,wire,[-L/2 L/2],Nwire,'gauss');

end

% Outer conductor: negative current
for k = 1:Ns

    x0 = shellPos(k,1);
    y0 = shellPos(k,2);

    wire = @(t) [x0+zeros(size(t)), ...
                 y0+zeros(size(t)), t];

    srcs{Ni+k} = em.src.lineCurrent( ...
        -I/Ns,wire,[-L/2 L/2],Nwire,'gauss');

end

coax = srcs{1};

for k = 2:numel(srcs)
    coax = em.src.merge(coax,srcs{k});
end

FH = @(r) fieldBatch(@(p) em.field.H(coax,p),r);
FB = @(r) fieldBatch(@(p) em.field.B(coax,p),r);

fprintf('Coaxial cable: %d inner and %d shell filaments\n',Ni,Ns);

%% Test 1: Circulation in all four regions

rho = [0.25 0.65 1.0 1.5];

expected = [I*rho(1)^2/a^2;
            I;
            I*(1-(rho(3)^2-b^2)/(c^2-b^2));
            0];

actual = zeros(4,1);

for k = 1:4

    radius = rho(k);

    path = @(t) [radius*cos(t), ...
                 radius*sin(t), zeros(size(t))];

    actual(k) = em.field.circulation( ...
        FH,path,[0 2*pi],400,'midpoint');

end

err = abs(actual-expected);

assert(all(err <= 0.05*abs(I)), ...
    'Test 1 failed: Coaxial circulation error exceeds 5%.');

fprintf('Test 1 passed: Circulation in all four regions\n');

for k = 1:4
    fprintf('rho = %.2f, numerical = %.6f, expected = %.6f A\n', ...
        rho(k),actual(k),expected(k));
end

%% Test 2: Circulation outside cable

tolOutside = 1e-3*abs(I);

assert(abs(actual(4)) <= tolOutside, ...
    'Test 2 failed: External circulation is not zero.');

fprintf('Test 2 passed: External circulation = %.3e A\n', ...
    actual(4));

%% Test 3: Independence of path shape

% Circular path
radius = 0.65;

circle = @(t) [radius*cos(t), ...
               radius*sin(t),zeros(size(t))];

Ccircle = em.field.circulation( ...
    FH,circle,[0 2*pi],400,'midpoint');

% Square path
d = 0.55;

squareVertices = [-d -d 0;
                   d -d 0;
                   d  d 0;
                  -d  d 0;
                  -d -d 0];

Csquare = polygonCirculation(FH,squareVertices);

% Irregular closed path
irregular = @(t) [ ...
    (0.65+0.04*cos(3*t)).*cos(t), ...
    (0.65+0.04*cos(3*t)).*sin(t), ...
    zeros(size(t))];

Cirregular = em.field.circulation( ...
    FH,irregular,[0 2*pi],400,'midpoint');

Cvalues = [Ccircle Csquare Cirregular];

pathError = (max(Cvalues)-min(Cvalues))/abs(I);

assert(pathError <= 0.01, ...
    'Test 3 failed: Path disagreement exceeds 1%.');

fprintf('Test 3 passed: Path disagreement = %.4f%%\n', ...
    100*pathError);

%% Construct magnetic sources for flux tests

% Finite straight current segment
wire = @(t) [zeros(size(t)),zeros(size(t)),t];

sWire = em.src.lineCurrent( ...
    I,wire,[-0.4 0.4],200,'gauss');

% Circular current loop
radius = 0.4;

ring = @(t) [radius*cos(t), ...
             radius*sin(t),zeros(size(t))];

sRing = em.src.lineCurrent( ...
    I,ring,[0 2*pi],400,'gauss');

sources = {sWire,sRing};
names = {'Straight wire','Current loop'};

%% Test 4: Magnetic flux through closed surfaces

Rsurf = 0.6;
center = [0.08 0.03 0.1];

for k = 1:2

    s = sources{k};

    Bfun = @(r) fieldBatch(@(p) em.field.B(s,p),r);

    % Closed sphere
    sphere = @(th,ph) center + Rsurf*[ ...
        sin(th).*cos(ph), ...
        sin(th).*sin(ph),cos(th)];

    PhiSphere = em.field.flux( ...
        Bfun,sphere,[0 pi],[0 2*pi],60,120);

    % Upper hemisphere
    hemisphere = sphere;

    PhiHemisphere = em.field.flux( ...
        Bfun,hemisphere,[0 pi/2],[0 2*pi],60,120);

    % Bottom disk with outward normal in -z direction
    disk = @(ph,rho) center + [ ...
        rho.*cos(ph),rho.*sin(ph),zeros(size(rho))];

    PhiDisk = em.field.flux( ...
        Bfun,disk,[0 2*pi],[0 Rsurf],120,60);

    PhiClosedHemisphere = PhiHemisphere + PhiDisk;

    % Closed cube
    PhiCube = cubeFlux(Bfun,center,Rsurf);

    % Reference magnetic flux scale
    Bref = norm(Bfun(center));
    scale = max(Bref*pi*Rsurf^2, ...
                em.const.mu0()*abs(I)*Rsurf);

    tolFlux = 0.01*scale;

    values = [PhiSphere PhiClosedHemisphere PhiCube];

    assert(all(abs(values) <= tolFlux), ...
        'Test 4 failed: Nonzero closed-surface magnetic flux.');

    fprintf('%s:\n',names{k});
    fprintf('  Sphere:     %.3e Wb\n',PhiSphere);
    fprintf('  Hemisphere: %.3e Wb\n',PhiClosedHemisphere);
    fprintf('  Cube:       %.3e Wb\n',PhiCube);
    fprintf('  Reference:  %.3e Wb\n',scale);

end

fprintf('Test 4 passed: All six closed-surface flux checks\n');

%% Test 5: Surface cutting through current loop

Bfun = @(r) fieldBatch(@(p) em.field.B(sRing,p),r);

centerCut = [0.3 0 0];
Rcut = 0.35;

sphereCut = @(th,ph) centerCut + Rcut*[ ...
    sin(th).*cos(ph), ...
    sin(th).*sin(ph),cos(th)];

PhiCut = em.field.flux( ...
    Bfun,sphereCut,[0 pi],[0 2*pi],100,200);

scaleCut = em.const.mu0()*abs(I)*Rcut;

assert(abs(PhiCut) <= 0.01*scaleCut, ...
    'Test 5 failed: Flux through surface cutting loop.');

fprintf('Test 5 passed: Cutting surface flux = %.3e Wb\n', ...
    PhiCut);

%% Test 6: Curl H inside conductor

h = 0.15;

rin = [0.25 0 0];

curlInside = em.op.curl(FH,rin,h);

J = I/(pi*a^2);
Jexpected = [0 0 J];

curlError = norm(curlInside-Jexpected)/abs(J);

assert(curlError <= 0.05, ...
    'Test 6 failed: Curl H differs from J by more than 5%.');

fprintf('Test 6 passed: Curl H inside conductor\n');
fprintf('  Calculated Jz = %.6f A/m^2\n',curlInside(3));
fprintf('  Expected Jz   = %.6f A/m^2\n',J);
fprintf('  Error         = %.4f%%\n',100*curlError);

%% Test 7: Curl H in dielectric

rout = [0.65 0 0];

curlOutside = em.op.curl(FH,rout,h);

tolCurl = 0.05*abs(J);

assert(norm(curlOutside) <= tolCurl, ...
    'Test 7 failed: Curl H is not zero in dielectric.');

fprintf('Test 7 passed: Curl H in dielectric = %.3e A/m^2\n', ...
    norm(curlOutside));

%% Test 8: Divergence B inside and outside

divInside = em.op.div(FB,rin,h);
divOutside = em.op.div(FB,rout,h);

tolDiv = 0.05*em.const.mu0()*abs(J);

assert(abs(divInside) <= tolDiv, ...
    'Test 8 failed: Divergence B inside conductor.');

assert(abs(divOutside) <= tolDiv, ...
    'Test 8 failed: Divergence B outside conductor.');

fprintf('Test 8 passed: Divergence B\n');
fprintf('  Inside  = %.3e T/m\n',divInside);
fprintf('  Outside = %.3e T/m\n',divOutside);

%% Test 9: Curl of electrostatic field

q = 1e-9;

sCharge = em.src.pointCharge(q,[0 0 0]);

FE = @(r) em.field.E(sCharge,r);

rTest = [1.0 0.5 0.3;
         0.4 1.0 0.6;
         0.7 0.8 1.0];

curlE = em.op.curl(FE,rTest,1e-5);

curlEmag = sqrt(sum(curlE.^2,2));

tolE = 1e-6;

assert(all(curlEmag <= tolE), ...
    'Test 9 failed: Electrostatic curl is not zero.');

fprintf('Test 9 passed: Maximum curl E = %.3e V/m^2\n', ...
    max(curlEmag));

%% All numerical tests completed

fprintf('All nine Lab 10 numerical tests passed!\n');

passed = true;

end


%% Helper function: Batched field evaluation

function Fval = fieldBatch(F,r)

N = size(r,1);
Fval = zeros(N,3);

batchSize = 50;

for k = 1:batchSize:N

    idx = k:min(k+batchSize-1,N);

    Fval(idx,:) = F(r(idx,:));

end

end


%% Helper function: Polygon circulation

function C = polygonCirculation(F,vertices)

C = 0;

for k = 1:size(vertices,1)-1

    p1 = vertices(k,:);
    p2 = vertices(k+1,:);

    path = @(t) p1 + t.*(p2-p1);

    C = C + em.field.circulation( ...
        F,path,[0 1],200,'midpoint');

end

end


%% Helper function: Closed cube flux

function Phi = cubeFlux(F,center,d)

% Positive x face
Sxp = @(u,v) center + [d+zeros(size(u)),u,v];

% Negative x face
Sxm = @(u,v) center + [-d+zeros(size(u)),v,u];

% Positive y face
Syp = @(u,v) center + [v,d+zeros(size(u)),u];

% Negative y face
Sym = @(u,v) center + [u,-d+zeros(size(u)),v];

% Positive z face
Szp = @(u,v) center + [u,v,d+zeros(size(u))];

% Negative z face
Szm = @(u,v) center + [v,u,-d+zeros(size(u))];

faces = {Sxp,Sxm,Syp,Sym,Szp,Szm};

Phi = 0;

for k = 1:6

    Phi = Phi + em.field.flux( ...
        F,faces{k},[-d d],[-d d],60,60);

end

end
