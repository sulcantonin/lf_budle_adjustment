%% create a random rotation matrix up to a maximum angle per axis
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [R] = random_rotation( max_angle )
    R = rotx( max_angle * rand() ) * roty( max_angle * rand() ) * rotz( max_angle * rand() );
end
