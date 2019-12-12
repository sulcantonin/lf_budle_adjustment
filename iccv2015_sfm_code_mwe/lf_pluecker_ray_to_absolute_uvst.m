%% Project a given ray in Pluecker world coordinates into a light field
%% Takes into account rotation and translation of the light field camera,
%% also requires focal length.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [ok ray_t] = pluecker_ray_to_uvst( q,m, R,t,f )

    % new version uses matrix transform to homogenous uvst
    [qt mt] = transform_pluecker_ray( q,m, R,t );
    ray = lf_ray_projection_matrix( f ) * [qt;mt];
    ok = ( abs( ray(5) ) ~= 0.0 );
    ray = ray / ray(5);
    ray_t = ray( 1:4 );

end
