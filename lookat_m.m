function R = lookat_m(pos,target,up)

if nargin<3, up = [0 1 0]; end
up = up/norm(up);
d = target(1:3) - pos(1:3);
n = norm(d);
if n<1e-13, R = eye(3); return; end

z = d(:)/n;

% x = [z(3); 0; -z(1)];
x = [up(2)*z(3)-up(3)*z(2); up(3)*z(1)-up(1)*z(3); up(1)*z(2)-up(2)*z(1)];
nx = norm(x);
if nx<1e-13, x = [1; 0; 0];
else x = x/nx; end

y = [z(2)*x(3)-z(3)*x(2); z(3)*x(1)-z(1)*x(3); z(1)*x(2)-z(2)*x(1)];

R = [x y z]';