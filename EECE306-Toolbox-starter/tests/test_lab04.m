function test_04 = test_lab04()
% Unit Circle
circ = @(t) [cos(t), sin(t), zeros(size(t))];
[~, w] = em.quad.line(circ, [0 2*pi], 200,'midpoint');
em.test.assertClose(sum(w), 2*pi, 1e-6, 'unit circle arc length');
% Unit Sphere Surface
sphere = @(theta,phi)[sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
[~, w] = em.quad.surf(sphere,[0 pi],[0 2*pi],101,101,'gauss');
em.test.assertClose(sum(w),4*pi,1e-6,'sphere surface area');
% Unit Ball
ball = @(r,theta,phi ) [r.*sin(theta).*cos(phi), r.*sin(theta).*sin(phi), r.*cos(theta)];
[~, w] = em.quad.vol(ball, [0 1], [0 pi], [0 2*pi], 21,21,21,'gauss');
em.test.assertClose(sum(w), 4*pi/3, 1e-6, 'ball volume');
% Unit Square
square = @(u,v) [u, v, zeros(size(u))];
[~, w] = em.quad.surf(square, [0 1], [0 1], 50, 50,'midpoint');
em.test.assertClose(sum(w), 1, 1e-6, 'square area');
test_04 = true;
end 