function surf_from_scatter(x,y,z)

%%
% The problem is that the data is made up of individual (x,y,z)
% measurements. It isn't laid out on a rectilinear grid, which is what the
% SURF command expects. A simple plot command isn't very useful.

% plot3(x,y,z,'.-')

%% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.

tri = delaunay(x,y);
% plot(x,y,'.')

%%
% How many triangles are there?

% [r,c] = size(tri);
% disp(r)

%% Plot it with TRISURF

h = trisurf(tri, x, y, z);

% axis vis3d
% 
% %% Clean it up
% 
% axis off
% l = light('Position',[-50 -15 29])
% set(gca,'CameraPosition',[208 -50 7687])
% lighting phong
% shading interp
colorbar EastOutside
