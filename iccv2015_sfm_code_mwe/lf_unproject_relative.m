%% Unproject 4D ray coordinates to light field camera coordinate frame
%% Requires disparity and focal length.
%% Note: this function uses relative two-plane parametrization
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [X] = lf_unproject_relative( ray, d, f )

    X = zeros(3,1);
    X(3) = f / d;
    X(1) = ray(1) * X(3) / f + ray(3);
    X(2) = ray(2) * X(3) / f + ray(4);

end
