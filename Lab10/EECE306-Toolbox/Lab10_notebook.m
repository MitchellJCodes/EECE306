%% EECE 306 Lab 10 notebook. Ampere's Circuital Law, Curl, and the Magnetic Gauss Law
% *Team NN.* TODO replace with your team number and member names.
%
% Fill in every TODO, run |publish('Lab10_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages. The coax cell below
% is the slowest computation of the semester so far, expect about a
% minute.

%% Setup
clear; close all;
mu0 = em.const.mu0();

%% The coaxial cable as a bundle of filaments
% The API has no volume current constructor, and none is needed. A
% uniform volume current is a merge of many parallel line currents, each
% carrying its share of the total. The inner conductor carries +1 A, the
% shell returns it.
I = 1; a = 0.5; b = 1.0; c = 1.3; L = 40; nel = 400;
[gx, gy] = meshgrid(-1.3:0.1:1.3);
rin  = hypot(gx, gy) <= a;
rsh  = hypot(gx, gy) >= b & hypot(gx, gy) <= c;
xin = gx(rin);  yin = gy(rin);
xsh = gx(rsh);  ysh = gy(rsh);
srcs = cell(1, numel(xin) + numel(xsh));
for k = 1:numel(xin)
    x0 = xin(k); y0 = yin(k);
    srcs{k} = em.src.lineCurrent(I/numel(xin), @(t)[x0+0*t y0+0*t t], [-L L], nel);
end
for k = 1:numel(xsh)
    x0 = xsh(k); y0 = ysh(k);
    srcs{numel(xin)+k} = em.src.lineCurrent(-I/numel(xsh), @(t)[x0+0*t y0+0*t t], [-L L], nel);
end
coax = em.src.merge(srcs{:});
fprintf('%d filaments inside, %d in the shell, %d elements total\n', ...
    numel(xin), numel(xsh), size(coax.pos,1));
FH = @(r) em.field.H(coax, r);

%% Circulation against enclosed current, region by region
rho = [0.12 0.25 0.38 0.47 0.6 0.75 0.9 1.05 1.15 1.25 1.45 1.6]';
Ienc = zeros(size(rho));
Circ = zeros(size(rho));
IencRef = @(r) I*(r<=a).*(r.^2/a^2) + I*(r>a & r<b) + ...
    I*(r>=b & r<=c).*(1 - (r.^2 - b^2)/(c^2 - b^2)) + 0*(r>c);
for k = 1:numel(rho)
    rk = rho(k);
    loop = @(t)[rk*cos(t) rk*sin(t) zeros(size(t))];
    Circ(k) = em.field.circulation(FH, loop, [0 2*pi], 96);
    Ienc(k) = IencRef(rk);
end
figure; plot(rho, Circ, 'o', rho, Ienc, '-'); grid on
xlabel('loop radius (m)'); ylabel('A')
legend('circulation of H', 'analytic enclosed current')
title('Ampere''s law in all four regions of the coax')
fprintf('largest deviation = %.3f A\n', max(abs(Circ - Ienc)));
fprintf('outside the cable, circulation = %.3e A  (consistent with zero)\n', Circ(end));

%% One enclosed current, three path shapes
rk = 0.7;
circle = em.field.circulation(FH, @(t)[rk*cos(t) rk*sin(t) 0*t], [0 2*pi], 96);
sq = 0; s2 = 0.7;
sides = { @(t)[t -s2+0*t 0*t], @(t)[s2+0*t t 0*t], @(t)[-t s2+0*t 0*t], @(t)[-s2-0*t -t 0*t] };
for k = 1:4
    sq = sq + em.field.circulation(FH, sides{k}, [-s2 s2], 48);
end
irr = em.field.circulation(FH, ...
    @(t)[0.7*(1+0.1*sin(3*t)).*cos(t) 0.68*(1+0.08*cos(2*t)).*sin(t) 0*t], [0 2*pi], 128);
fprintf('circle    %.4f A\nsquare    %.4f A\nirregular %.4f A\n', circle, sq, irr);

%% The magnetic Gauss law
% Closed surface integral of B over three surfaces, for a current ring,
% including a surface the ring passes through.
ring = em.src.lineCurrent(1, @(t)[cos(t) sin(t) 0*t], [0 2*pi], 400);
FB = @(r) em.field.B(ring, r);
sphS = @(th,ph) 0.6*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
P1 = em.field.flux(FB, sphS, [0 pi], [0 2*pi], 60, 120);
cubeF = 0;
faces = { @(u,v)[u v  0.4*ones(size(u))], @(u,v)[v u -0.4*ones(size(u))], ...
          @(u,v)[v  0.4*ones(size(u)) u], @(u,v)[u -0.4*ones(size(u)) v], ...
          @(u,v)[ 0.4*ones(size(u)) u v], @(u,v)[-0.4*ones(size(u)) v u] };
for k = 1:6
    cubeF = cubeF + em.field.flux(FB, faces{k}, [-0.4 0.4], [-0.4 0.4], 50, 50);
end
sphCut = @(th,ph) [1 0 0] + 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
P3 = em.field.flux(FB, sphCut, [0 pi], [0 2*pi], 60, 120);
Bscale = em.vec.mag(em.field.B(ring, [0 0 0])) * pi * 0.6^2;
fprintf('sphere                    %.3e Wb\n', P1);
fprintf('cube                      %.3e Wb\n', cubeF);
fprintf('sphere cutting the ring   %.3e Wb\n', P3);
fprintf('scale for comparison, B(0) times a disk area = %.3e Wb\n', Bscale);

%% The point form
J = I / (pi*a^2);
cin  = em.op.curl(FH, [0.25 0 0], 0.15);
cout = em.op.curl(FH, [0.75 0 0], 0.1);
dB1  = em.op.div(FB, [0.3 0.2 0.4], 1e-4);
dB2  = em.op.div(FB, [2 1 1], 1e-4);
fprintf('curl H inside conductor  z component = %.4f   J = %.4f A/m^2\n', cin(3), J);
fprintf('curl H in the dielectric z component = %.4f  (consistent with zero)\n', cout(3));
fprintf('div B near the ring = %.3e, far away = %.3e\n', dB1, dB2);
% The derivative step for the curl must average over several filaments,
% which is why h = 0.15 m is passed rather than the default.

%% What field lines can and cannot show
wire1 = em.src.lineCurrent(1, @(t)[0.3+0*t 0*t t], [-40 40], 400);
wire2 = em.src.lineCurrent(-1, @(t)[-0.3+0*t 0*t t], [-40 40], 400);
pair = em.src.merge(wire1, wire2);
th = (0:30:330)';
figure;
em.viz.streamlines(@(r) em.field.H(wire1, r), [-1 1], [-1 1], ...
    [0.3+0.15*cosd(th) 0.15*sind(th)], struct('title', 'One wire, circles'));
figure;
em.viz.streamlines(@(r) em.field.H(pair, r), [-1.2 1.2], [-1.2 1.2], ...
    [[0.3+0.12*cosd(th) 0.12*sind(th)]; [-0.3+0.12*cosd(th) 0.12*sind(th)]], ...
    struct('title', 'Antiparallel pair'));
% TODO one paragraph. Every line in both figures closes on itself and
% none begins or ends anywhere. Connect that to div B = 0, and contrast
% it with the Lab 3 field, whose lines begin and end on charges because
% div D = rho. Also say what a field line plot fails to show, which the
% Lab 3 quiver comparison already taught.

%% Interpretation
% TODO three to six sentences. Ampere's law was checked in four regions
% with one geometry. State why the region outside the shield is the most
% informative of the four.

%% Problems encountered
% TODO honest account, or NONE.

%% Full test suite
runTests
