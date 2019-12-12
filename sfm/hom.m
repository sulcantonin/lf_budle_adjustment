%% Adds extra row to matrix V, to create homogenous coords
%
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras
function [v] = hom(v);
v(end+1,:) =1;
end
