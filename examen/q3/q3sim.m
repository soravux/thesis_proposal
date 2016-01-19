function [angle, rc] = q3sim(r, N, e)
%% Perform render

theta = 0:2*pi/3:2*pi;
theta = theta(1:end-1);

x = (sin(r)*cos(theta))';
y = (sin(r)*sin(theta))';
z = sqrt(1 - x.^2 - y.^2);
L = [x, y, z];

b = max(0, L*N);


%% Perfom PS

xe = (sin(r+e)*cos(theta))';
ye = (sin(r+e)*sin(theta))';
ze = sqrt(1 - xe.^2 - ye.^2);
Le = [xe, ye, ze];

rc = rcond(Le);
if rc < eps
    angle = NaN;
else
    angle = acosd(min(1, normc(Le\b)'*N));
end
%angle = acosd((Le\b)'*N);
% This supposes the albedo is known to be unit
assert(isreal(angle));
end