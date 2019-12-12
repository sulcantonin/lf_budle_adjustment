%% transforms 3D points, poses and focal distance into a single vector for 
%% bundle adjustment
% Input:    X    : 3D points (X \in R^{3 x n}), 
%           pose : camrea poses (pose \in R^{3 x 4 x m}
%           f    : focal distance (f \in R) 
% Output    theta : [7*m; 1; 3*n] vector (ordered as cameras, f, 3d points)
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function [theta] = lf_cameras_points_2_vector(X,pose,f)
n_cam = size(pose,3);
% extrinsics, 4 params for rotation, 3 for translation for each camera
cam = zeros(7,n_cam);
for i=1:n_cam
    R = pose(1:3,1:3,i);
    t = pose(:,4,i);
    cam(:,i) = lf_camera_2_vector(R,t);
end
% intrinsics 1 for f, 
theta = cat(1,cam(:),f(:));
% 3d point cloud, 3 params per point
theta = cat(1,theta,X(:));
end