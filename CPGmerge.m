function h = CPGmerge(varargin)
%
% Eric A. Mehiel
% Cal Poly, SLO
% Aerospace Engineering Department

newfacecolor = get(varargin{1}, 'facecolor');
newedgecolor = get(varargin{1}, 'edgecolor');
newfacealpha = get(varargin{1}, 'facealpha');
faces = get(varargin{1}, 'faces');
vertices = get(varargin{1}, 'vertices');

delete(varargin{1});

for i = 2:length(varargin)
    temp_faces = get(varargin{i}, 'faces') + max(max(faces));
    temp_vertices = get(varargin{i}, 'vertices');
    
    faces = [faces; temp_faces];
    vertices = [vertices; temp_vertices];
    
    delete(varargin{i});
end

h = patch('faces', faces, 'vertices', vertices, 'edgeColor',newedgecolor, 'faceColor', newfacecolor, 'facealpha', newfacealpha);

reducepatch(h,1/length(varargin));

%faces = get(h, 'faces');
%verts = get(h, 'vertices');