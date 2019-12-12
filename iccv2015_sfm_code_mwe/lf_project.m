%% Project world point X into light field ray with given (s,t) coordinate
%% Light field camera defined by rotation, center, focal length
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [u,v] = lf_project_absolute( X, s,t, R,C,f )

    Xc = R * X + C;
    u = f*(Xc(1) - s) / Xc(3) + s;
    v = f*(Xc(2) - t) / Xc(3) + t;

end
