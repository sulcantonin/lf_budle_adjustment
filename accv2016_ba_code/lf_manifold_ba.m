%% Bundle adjustment energy, the fuction is minimized by LSQNONLIN 
%% Matlab builtin, see "Layered scene reconstruction from multiple light field camera views" eq. (1)
% Input:    theta:  Vector with parameters to be optimized (see lf_cameras_points_2_vector)
%           rays :  Cell array (n x m) where each non-empty cell (i,j)
%                       represents a set of rays of a i-th 3D point visible in j-th camera 
%                       in light field coordinates [u,v,s,t]
% Output    e :     Energy (1) in "Layered scene reconstruction from multiple light 
%                       field camera views" before ^2
%
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function e = lf_manifold_ba(theta,rays)
    n_cam = size(rays,2);
    n_points = size(rays,1);
    n_projections = arrayfun(@(i) size(rays{i},1),1:numel(rays));

    e = zeros(sum(n_projections),1);
    X = ones(4,n_points);

    [X(1:3,:), pose,f] = lf_vector_2_cameras_points(theta,n_cam,n_points);
    
    % transforming entire point cloud into coordinate frame of
    % corresponding camera
    X_t = ones(3,n_points,n_cam);
    for i=1:n_cam
        X_t(:,:,i) = pose(:,:,i) * X;
    end

    
    % transforming the 3D points into 2D subspace, checkout eq. (10) of "On
    % Linear Structure from Motion for Light Field Cameras" or (1) of
    % "Layered scene reconstruction from multiple light field camera views"
    % for more details
    M = lf_M_from_X(X_t,f);

    
    k = 1;
    for j=1:n_cam
        for i=1:n_points
            if ~isempty(rays{i,j})
                e_ = M(:,:,i,j) * (hom(rays{i,j}'));
                e(k:k+numel(e_)-1)=e_(:);
                k = k+numel(e_);
            end
        end
    end

    e(~isfinite(e(:))) = 0;

