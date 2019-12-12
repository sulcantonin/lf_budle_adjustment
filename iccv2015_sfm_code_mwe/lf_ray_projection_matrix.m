%% Compute ray projection matrix for a light field, equation (9)
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [P] = lf_ray_projection_matrix( f )
    P = [ f 0 0   0 -1 0; 0 f 0  1 0 0;  0 0 0  0 -1 0;   0 0 0  1 0 0; 0 0 1  0 0 0 ];
end
