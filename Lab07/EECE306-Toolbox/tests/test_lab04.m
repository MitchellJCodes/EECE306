function test_04 = test_lab04()
% test_lab04.m  Lab 4 checks. Add to this file, do not delete existing checks.
% Run from the toolbox root with:  runTests

% 1-D nodes and weights

N = 16;
[t,w] = em.quad.nodes(0,1,N);

assert(isequal(size(t), [N 1]), ...
    'nodes must return N-by-1 nodes');
assert(isequal(size(w), [N 1]), ...
    'nodes must return N-by-1 weights');

em.test.assertClose(sum(w), 1, 1e-12, ...
    'midpoint weights sum to interval length');

[t2,w2] = em.quad.nodes(0,1,N,'midpoint');

em.test.assertClose(t,t2,1e-12, ...
    'default quadrature rule is midpoint');
em.test.assertClose(w,w2,1e-12, ...
    'default quadrature weights are midpoint');


% 2. Weight sums for all rules

rules = {'midpoint','trapz','simpson','gauss'};

for k = 1:numel(rules)
    [~,wk] = em.quad.nodes(-2,3,20,rules{k});

    em.test.assertClose(sum(wk),5,1e-11, ...
        [rules{k}, ' weights sum to interval length']);
end


% 3. Midpoint rule

[tm,wm] = em.quad.nodes(0,1,32,'midpoint');

Im = sum(exp(tm).*wm);
Iref = exp(1)-1;

assert(abs(Im-Iref) < 1e-3, ...
    'midpoint exponential integral is inaccurate');


% 4. Trapezoidal rule

[tt,wt] = em.quad.nodes(0,1,32,'trapz');

It = sum(exp(tt).*wt);

assert(abs(It-Iref) < 3e-4, ...
    'trapz exponential integral is inaccurate');


% 5. Simpson rule

[ts,ws] = em.quad.nodes(0,1,16,'simpson');

Is = sum(exp(ts).*ws);

assert(abs(Is-Iref) < 3e-7, ...
    'Simpson exponential integral is inaccurate');


% 6. Gauss rule

[tg,wg] = em.quad.nodes(0,1,8,'gauss');

Ig = sum(exp(tg).*wg);

assert(abs(Ig-Iref) < 1e-12, ...
    'Gauss exponential integral is inaccurate');


% 7. Gauss polynomial exactness

[tg,wg] = em.quad.nodes(-1,1,6,'gauss');

for p = 0:11
    Ip = sum((tg.^p).*wg);

    if mod(p,2) == 0
        Ipref = 2/(p+1);
    else
        Ipref = 0;
    end

    em.test.assertClose(Ip,Ipref,1e-10, ...
        sprintf('Gauss polynomial degree %d',p));
end


% 8. Line quadrature: straight line

lineC = @(t)[t zeros(size(t)) zeros(size(t))];

[pos,wline] = em.quad.line(lineC,[0 3],100);

assert(isequal(size(pos),[100 3]), ...
    'line positions must be N-by-3');
assert(isequal(size(wline),[100 1]), ...
    'line weights must be N-by-1');

em.test.assertClose(sum(wline),3,1e-10, ...
    'straight line length');


% 9. Line quadrature: unit circle

circ = @(t)[cos(t) sin(t) zeros(size(t))];

[pos,wcirc] = em.quad.line(circ,[0 2*pi],400);

em.test.assertClose(sum(wcirc),2*pi,1e-5, ...
    'unit circle circumference');


% 10. Line Jacobian: nonuniform parameterization

curve = @(t)[t.^2 zeros(size(t)) zeros(size(t))];

[pos,wcurve] = em.quad.line(curve,[0 1],400);

em.test.assertClose(sum(wcurve),1,1e-5, ...
    'line quadrature includes |dr/dt|');


% 11. Surface quadrature: unit square

sq = @(u,v)[u v zeros(size(u))];

[pos,wsq] = em.quad.surf(sq,[0 1],[0 1],20,20);

assert(isequal(size(pos),[400 3]), ...
    'surface positions must be Nu*Nv by 3');
assert(isequal(size(wsq),[400 1]), ...
    'surface weights must be Nu*Nv by 1');

em.test.assertClose(sum(wsq),1,1e-10, ...
    'unit square area');


% 12. Surface quadrature: unit sphere

sph = @(th,ph) ...
    [sin(th).*cos(ph) ...
     sin(th).*sin(ph) ...
     cos(th)];

[pos,wsph] = em.quad.surf( ...
    sph,[0 pi],[0 2*pi],60,120);

em.test.assertClose(sum(wsph),4*pi,3e-3, ...
    'unit sphere surface area');


% 13. Surface Jacobian must be positive

[pos,wsq2] = em.quad.surf( ...
    sq,[1 0],[0 1],20,20);

em.test.assertClose(sum(wsq2),1,1e-10, ...
    'surface uses absolute Jacobian');


% 14. Volume quadrature: unit cube

cube = @(u,v,q)[u v q];

[pos,wcube] = em.quad.vol( ...
    cube,[0 1],[0 1],[0 1],10,10,10);

assert(isequal(size(pos),[1000 3]), ...
    'volume positions must be Nu*Nv*Nw by 3');
assert(isequal(size(wcube),[1000 1]), ...
    'volume weights must be Nu*Nv*Nw by 1');

em.test.assertClose(sum(wcube),1,1e-10, ...
    'unit cube volume');


% 15. Volume Jacobian

box = @(u,v,q)[2*u 3*v 4*q];

[pos,wbox] = em.quad.vol( ...
    box,[0 1],[0 1],[0 1],8,8,8);

em.test.assertClose(sum(wbox),24,1e-10, ...
    'volume includes absolute Jacobian determinant');


% 16. Volume quadrature: unit ball

ball = @(u,v,q) ...
    [q.*sin(u).*cos(v) ...
     q.*sin(u).*sin(v) ...
     q.*cos(u)];

[pos,wball] = em.quad.vol( ...
    ball,[0 pi],[0 2*pi],[0 1],30,60,30);

em.test.assertClose(sum(wball),4*pi/3,6e-3, ...
    'unit ball volume');


% 17. Charge constructors and quadrature-rule forwarding
lineSource = em.src.lineCharge(2, @(t)[t zeros(size(t)) zeros(size(t))], ...
    [0 1], 12, 'gauss');
assert(numel(lineSource.q) == 12 && numel(lineSource.w) == 12, ...
    'lineCharge must preserve the requested number of elements');
em.test.assertClose(sum(lineSource.w),1,1e-10, ...
    'lineCharge forwards the requested quadrature rule');

surfaceSource = em.src.surfCharge(3, @(u,v)[u v zeros(size(u))], ...
    [0 1],[0 1],5,6,'gauss');
assert(numel(surfaceSource.q) == 30 && numel(surfaceSource.w) == 30, ...
    'surfCharge must preserve Nu*Nv elements');
em.test.assertClose(sum(surfaceSource.w),1,1e-10, ...
    'surfCharge forwards the requested quadrature rule');

volumeSource = em.src.volCharge(4, @(u,v,q)[u v q], ...
    [0 1],[0 1],[0 1],4,5,6,'gauss');
assert(numel(volumeSource.q) == 120 && numel(volumeSource.w) == 120, ...
    'volCharge must preserve Nu*Nv*Nw elements');
em.test.assertClose(sum(volumeSource.w),1,1e-10, ...
    'volCharge forwards the requested quadrature rule');

% 18. Scalar-only parameterizations using matrix multiplication (*).
% These fail on vector inputs and must be evaluated point by point.
scalarSurface = @(u,v) [u, v, u*v*0];
[posScalarSurface,wScalarSurface] = em.quad.surf( ...
    scalarSurface,[0 1],[0 1],12,10);
assert(isequal(size(posScalarSurface),[120 3]), ...
    'scalar-only surface parameterization must return Nu*Nv points');
em.test.assertClose(sum(wScalarSurface),1,1e-8, ...
    'scalar-only surface parameterization area');

scalarVolume = @(u,v,q) [u, v, u*q*0 + q];
[posScalarVolume,wScalarVolume] = em.quad.vol( ...
    scalarVolume,[0 1],[0 1],[0 1],6,5,4);
assert(isequal(size(posScalarVolume),[120 3]), ...
    'scalar-only volume parameterization must return Nu*Nv*Nq points');
em.test.assertClose(sum(wScalarVolume),1,1e-8, ...
    'scalar-only volume parameterization volume');

% Simpson must accept an odd requested node count as well.
[to,wo] = em.quad.nodes(0,1,9,'simpson');
em.test.assertClose(sum(wo),1,1e-12, ...
    'odd-node Simpson weights sum to interval length');
em.test.assertClose(sum(exp(to).*wo),Iref,3e-6, ...
    'odd-node Simpson exponential integral');

% 19. Public convergence helper returns the expected order for O(N^-2).
[pconv,hconv] = em.test.convergence(@(n) n.^(-2),[8 16 32 64], ...
    struct('title','Convergence helper regression test'));
em.test.assertClose(pconv,2,1e-12, ...
    'convergence helper estimates the observed order');
if isgraphics(hconv), close(hconv); end

% 19. Observed convergence order

Ns = [8 16 32 64 128]';

for k = 1:3

    rule = rules{k};
    err = zeros(size(Ns));

    for j = 1:numel(Ns)

        [tj,wj] = em.quad.nodes(0,1,Ns(j),rule);

        Ij = sum(exp(tj).*wj);

        err(j) = abs(Ij-Iref);
    end

    pfit = polyfit(log(Ns),log(err),1);
    p = -pfit(1);

    if strcmp(rule,'midpoint') || strcmp(rule,'trapz')
        assert(p > 1.5, ...
            [rule, ' should have approximately second-order convergence']);
    end

    if strcmp(rule,'simpson')
        assert(p > 3.0, ...
            'Simpson should have approximately fourth-order convergence');
    end
end


% Finished

disp('  test_lab04 checks complete');
