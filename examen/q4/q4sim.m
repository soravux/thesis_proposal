t = 0:0.1:2*pi;
x = cos(t);
z = sin(t);
y = tand(10)*sin(t*2);
r = [1 0 0; 0 cosd(80) -sind(80); 0 sind(80) cosd(80)];

a = [x' y' z'];
a = normr(a);
b = r*a';

%plot3(a(:,1), a(:,2), a(:,3));
%hold on; 
plot3(b(1,:), b(2,:), b(3,:));
%hold off;
axis equal;