function h = arrow(xlim, ylim, zlim, varargin)
%
% Eric A. Mehiel
% Cal Poly, SLO
% Aerospace Engineering Department

s = [xlim(1); ylim(1); zlim(1)];
e = [xlim(2); ylim(2); zlim(2)];
l_vec = e-s;
l = norm(l_vec);
%if nargin ~= 4
    d = l/30;
%end

r = [d/2*ones(1,16), d:-d/3:0];

[x, y, z] = cylinder(r);

x = x + xlim(1);
y = y + ylim(1);
z = l*z + zlim(1);

if isempty(varargin)
    h = patch(surf2patch(x,y,z,'triangles'), 'edgeColor','none', 'faceColor', 'blue');
else
    h = patch(surf2patch(x,y,z,'triangles'), 'edgeColor','none', varargin{:});
end

lam = cross([0 0 1], l_vec)/l;
theta = 180/pi*acos(dot([0 0 1], l_vec/l));
if theta ~= 0
    rotate(h, lam, theta, s);
end

len = norm(diff([xlim; ylim; zlim],1,2));

set(h, 'userData', [xlim ylim zlim len]);