function h = cylin(rx, ry, l, n, varargin)
%
% h = cylin(rx, ry, l) forms a cylinder with radius rx by ry and
% length l.
%
% h = cylin(rx, ry, l, 'PropertyName',PropertyValue,...) sets any patch
% properties.
%
%
% Eric A. Mehiel
% Cal Poly, SLO
% Aerospace Engineering Department

r = ones(1,16);

[x, y, z] = cylinder(r, n);

x = rx*x;
y = ry*y;
z = l*z;

if (nargin == 4)
    hold on
    h = patch(surf2patch(x,y,z,'triangles'), 'edgeColor','none', 'faceColor', [.5 .5 .5]);
elseif (strcmp(varargin{1},'top'))
    th = 0:pi/8:2*pi;
    tx = rx*cos(th);
    ty = ry*sin(th);
    tz = l*ones(size(tx));
    bx = rx*cos(th);
    by = ry*sin(th);
    bz = zeros(size(tx));
    trit = delaunay(tx,ty);
    trib = delaunay(bx,by);
    hold on
    t = trimesh(trit,tx,ty,tz, 'edgeColor','none', 'faceColor', [.5 .5 .5]);
    hold on
    b = trimesh(trib,bx,by,bz, 'edgeColor','none', 'faceColor', [.5 .5 .5]);
    c = patch(surf2patch(x,y,z,'triangles'), 'edgeColor','none', 'faceColor', [.5 .5 .5]);
    h = CPGmerge(b,t, c);
    if (nargin > 4)
        set(h, varargin{2:end});
    end
else
    h = patch(surf2patch(x,y,z), 'edgeColor','none', varargin{:});
end