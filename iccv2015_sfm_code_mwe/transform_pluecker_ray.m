%% Transform a given ray in Pluecker coords
%% according to a rotation and translation of the coordinate system
%% See equation (3) in the paper.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [qt mt] = transform_pluecker_ray( q,m, R,t )
    qt = R * q;
    mt = R*m + cross( t, qt );
end
