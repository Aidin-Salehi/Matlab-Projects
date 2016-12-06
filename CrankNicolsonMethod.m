clear;clc;

theta = 1
C = 1/(pi^2);

T = 1;
nt = 20;
dt = T/nt;
t = 0:dt:T;

X = 1;
nx = 20;
dx = X/nx;
x = 0:dx:X;

r = C*dt/(dx)^2;
rl = r*theta;
rr = r*(1-theta);

e = ones(nx+1, 1);
Al = spdiags([-rl*e 1+2*rl*e -rl*e], -1:1, nx+1, nx+1);
Ar = spdiags([rr*e 1-2*rr*e rr*e], -1:1, nx+1, nx+1);

% u0 = @(x) 1/(pi)^2*sin(pi*x);
% u0 = @(x) 1-abs(1-2*x);
u0 = @(x) rand(size(x));

v = zeros(nx+1, nt+1);
v(:, 1) = u0(x);

B = (Al)^-1;
for i = 1:nt
    v(:, i+1) = B * Ar * v(:, i);
    v(1, i+1) = 0;
    v(nx+1, i+1) = 0;
end
v;

mesh(x,t,v')

title('Crank-Nicolson method')
xlabel('x')
ylabel('t')
zlabel('v')
zlim([0 max(max(v))])
% colormap autumn

% u = @(x,t) 1/(pi^2).*exp(-t).*sin(pi*x);
% u = u(x,T);
% d = u - v(:,nx+1)';
% e = d*d'
