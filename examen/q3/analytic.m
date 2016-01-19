r = sym('r');
e = sym('e');
t = sym('t');
n1 = sym('n1');
n2 = sym('n2');
n3 = sym('n3');
N = [n1; n2; n3];

theta = 0:2*pi/3:2*pi;
theta = theta(1:end-1);

x = (sin(r)*cos(theta))';
y = (sin(r)*sin(theta))';
z = sqrt(1 - x.^2 - y.^2);
L = [x, y, z];

xd = (sin(r)*cos(t))';
yd = (sin(r)*sin(t))';
zd = sqrt(1 - xd.^2 - yd.^2);
Ld = [xd, yd, zd];

b = L*N;


xe = (sin(r+e)*cos(theta))';
ye = (sin(r+e)*sin(theta))';
ze = sqrt(1 - xe.^2 - ye.^2);
Le = [xe, ye, ze];

r = Le\b;
%angle = acosd(normc((Le\b)')*N);

%%
clear x y z;
%syms x y z;
[x, y] = meshgrid(0:70, -30:30);

r1 = sind(x)./sind(y+x);
r3 = sqrt(1-sind(x).^2)./sqrt(1-sind(y+x).^2);
%ezplot(r1,[0 10], [0 10]);
surf(x, y, r1);
view(-155, 30);
xlabel('r'); ylabel('\epsilon'); zlabel('\Delta_\epsilon')
zlim([-2 5])
title('Components n_x and n_y')
export_fig('q3_analytic_nx_ny.pdf', '-transparent');

%%

surf(x, y, r3);
view(-55, 30);
xlabel('r'); ylabel('\epsilon'); zlabel('\Delta_\epsilon')
zlim([0, 5]);
title('Component n_z')
export_fig('q3_analytic_nz.pdf', '-transparent');