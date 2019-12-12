%% transforms a vector with parameters for BA into 3D points, poses and focal distance
% Input      theta : Vector with parameters (ordered as cameras, f, 3d
%                       points ~ [7*m; 1; 3*n])
% Output:    X    : 3D points (X \in R^{3 x n}), 
%            pose : camrea poses, pose(1:3,1:3) is rotation, pose(1:3,4) translation
%                       (pose \in R^{3 x 4 x m}_
%            f    : focal distance (f \in R) 
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function [X, pose,f] = lf_vector_2_cameras_points(theta,n_cam,n_point)
    pose = zeros(3,4,n_cam);

    % 7 parameters for camera (4 rotation, 3 translation),
    % 3 parameters per point, 1 focal length for all cameras
    
    % extrinsics -  7 params per camera
    cam = reshape(theta(1:7*n_cam),7,[]);
    for i=1:n_cam
        [pose(1:3,1:3,i),pose(:,4,i)] = lf_vector_2_camera(cam(:,i));
    end
    offset = 7*n_cam;
    % focal distance - 1 param
    f = repmat( theta(offset+1),n_cam,1);
    offset = offset + 1;

    % 3D point cloud - 3 params per point 
    X(1:3,:) = reshape(	theta(offset+1 :end),3,[]);

    assert(size(X,1) == 3); 
    assert(size(X,2) == n_point);
end