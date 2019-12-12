%% For a given set of rays it changes ordering of columns in each cell according to a given order ord 
% it exists because Light Field Toolbox 4.0 uses different ordering of light field coordinates (s,t,u,v) from 
% ours (u,v,s,t)
% Input:	rays : matrix or cell array with at least one ray (row) and as many columns as length(ord)
% 		ord  : desired permutation of coordiantes
% Output	rays : new order
%
% (c) O. Johannsen, A. Sulc, B. Goldluecke, N. Marniok, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, N. Marniok, B. Goldluecke: Layered scene reconstruction from multiple light field camera views 
function [rays] = change_ray_coordinates(rays,ord)
is_cell = iscell(rays);

if ~is_cell
    rays = {rays};
end

for i=1:size(rays,1)
    for j=1:size(rays,2)
        rays{i,j} = rays{i,j}(:,ord);
    end
end

if ~is_cell
    rays = rays{1};
end
