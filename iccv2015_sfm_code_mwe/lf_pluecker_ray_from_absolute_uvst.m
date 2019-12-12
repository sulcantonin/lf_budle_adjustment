%% Transform absolute (u,v,s,t) coords into Pluecker ray coords
%% Requires a given focal length (distance between two planes)
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [q,m] = lf_pluecker_ray_from_absolute_uvst( ray, f )

    % project to (s,t) view
    q = [ ray(1) - ray(3); ray(2) - ray(4); f ];
%     q = q/norm(q);
    
    % m is the crossproduct c x q, where c is the center of projection
    m = cross( [ ray(3); ray(4); 0 ], q );

    % sanity check
    %[ok x1] = intersect_pluecker_ray_plane( q,m, [0; 0; 1; 0] )
    %[ok x2] = intersect_pluecker_ray_plane( q,m, [0; 0; 1; -f] )
    
end
