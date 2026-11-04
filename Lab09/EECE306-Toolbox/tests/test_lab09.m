function test_09 = test_lab09()
%TEST_LAB09 Magnetostatic source & Biot-Savart

mu0 = em.const.mu0();
I = 2.3;
a = 0.4;
circle = @(t) [a*cos(t), a*sin(t), zeros(size(t))];

%% Circular loop
loop = em.src.lineCurrent(I,circle,[0 2*pi],128,'gauss');
zvals = linspace(-1.2,1.2,10)';
rax = [zeros(size(zvals)),zeros(size(zvals)),zvals];
Hax = em.field.H(loop,rax);
Href = [zeros(size(zvals)),zeros(size(zvals)), ...
    I*a^2 ./ (2*(a^2+zvals.^2).^(3/2))];
em.test.assertClose(Hax,Href,1e-6,'Circular-loop axial field');

Hcenter = em.field.H(loop,[0 0 0]);
em.test.assertClose(Hcenter(3),I/(2*a),1e-6, ...
    'Circular-loop center field');
transverse = Hax(:,1:2);
assert(max(abs(transverse(:))) < 1e-11, ...
    'Transverse magnetic-field components should vanish on loop axis');

%% Magnetic flux density
Btest = em.field.B(loop,[0 0 0],2.5);
em.test.assertClose(Btest,2.5*mu0*Hcenter,1e-12, ...
    'B = mu0*mur*H');

%% Finite straight segment approaching the infinite-wire limit
rho = 0.2;
h = 8*rho;
wire = em.src.lineCurrent(I,@(t)[zeros(size(t)),zeros(size(t)),t], ...
    [-h h],128,'gauss');
Hwire = em.field.H(wire,[rho 0 0]);
Hinf = I/(2*pi*rho);
wireRatio = Hwire(2)/Hinf;
assert(wireRatio > 0.99, ...
    sprintf('Finite-wire ratio %.6f does not exceed 0.99 at h/rho=%.3f', ...
    wireRatio,h/rho));

%% Finite solenoid
Nturns = 100;
asol = 0.2;
hsol = 15*asol;
sol = [];
for k = 1:Nturns
    zk = -hsol/2 + (k-0.5)*hsol/Nturns;
    c = @(t) [asol*cos(t),asol*sin(t),zk+zeros(size(t))];
    sk = em.src.lineCurrent(I,c,[0 2*pi],64,'gauss');
    if isempty(sol), sol = sk; else, sol = em.src.merge(sol,sk); end
end
Hsol = em.field.H(sol,[0 0 0]);
Hlong = Nturns*I/hsol;
solRatio = Hsol(3)/Hlong;
assert(solRatio >= 0.99, ...
    sprintf('Solenoid center ratio %.6f is below 0.99 at h/a=%.3f', ...
    solRatio,hsol/asol));

%% Current reversal and linearity
Hforward = em.field.H(loop,[0.17 -0.11 0.31]);
loopReverse = em.src.lineCurrent(-I,circle,[0 2*pi],128,'gauss');
Hreverse = em.field.H(loopReverse,[0.17 -0.11 0.31]);
em.test.assertClose(Hreverse,-Hforward,1e-12, ...
    'Reversing current reverses each field component');
loopDouble = em.src.lineCurrent(2*I,circle,[0 2*pi],128,'gauss');
Hdouble = em.field.H(loopDouble,[0.17 -0.11 0.31]);
em.test.assertClose(Hdouble,2*Hforward,1e-12, ...
    'Doubling current doubles each field component');

%% Surface-current constructor
K = @(r) repmat([0 0 3],size(r,1),1);
plane = @(u,v)[u v zeros(size(u))];
ss = em.src.surfCurrent(K,plane,[0 1],[0 1],5,6,'gauss');
assert(isequal(size(ss.Idl),[30 3]), ...
    'surfCurrent must produce Nu*Nv by 3 current elements');
em.test.assertClose(sum(ss.Idl,1),[0 0 3],1e-12, ...
    'Uniform surface current integrates over unit square');

%% Source merge compatibility
loopA = em.src.lineCurrent(I,circle,[0 pi],32,'gauss');
loopB = em.src.lineCurrent(I,circle,[pi 2*pi],32,'gauss');
merged = em.src.merge(loopA,loopB);
assert(size(merged.pos,1)==64 && isequal(size(merged.Idl),[64 3]), ...
    'Merging current sources must concatenate all elements');
charge = em.src.pointCharge(1e-9,[0 0 0]);
mixedRejected = false;
try
    em.src.merge(loopA,charge);
catch
    mixedRejected = true;
end
assert(mixedRejected,'Merging current and charge sources must raise an error');

%% Quadrature tangent
[pos,w,tangent] = em.quad.line(circle,[0 2*pi],32,'gauss');
assert(isequal(size(pos),[32 3]) && isequal(size(w),[32 1]) && ...
    isequal(size(tangent),[32 3]), ...
    'em.quad.line must return positions, weights, and tangent');

%% Vectorization/performance
perfSource = em.src.lineCurrent(I,circle,[0 2*pi],96,'gauss');
obs = [linspace(-0.8,0.8,5000)',0.7*ones(5000,1),0.5*ones(5000,1)];
tic;
Hmany = em.field.H(perfSource,obs);
elapsed = toc;
assert(isequal(size(Hmany),[5000 3]), ...
    'em.field.H must return N-by-3 for 5000 observation points');
assert(elapsed < 2, ...
    sprintf('em.field.H took %.3f s for 5000 points (limit 2 s)',elapsed));

fprintf(['  Lab 9 loop center Hz %.6g A/m; wire ratio %.5f at h/rho=%.1f; ' ...
    'solenoid ratio %.5f at h/a=%.1f; 5000-point H %.3f s\n'], ...
    Hcenter(3),wireRatio,h/rho,solRatio,hsol/asol,elapsed);
disp('  test_lab09 checks complete');
end
