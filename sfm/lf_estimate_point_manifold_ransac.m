% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras
function [M,best_inliers] = lf_estimate_point_manifold_ransac(rays,inlier_threshold,p)
e = 0.5; 
max_iter = uint32(log(1-p)/log(1-(1-e)^2));
K = size(rays,1);

% inliers = false(K);
best_inliers = false(1,K);
iter = 1;
while iter < max_iter 
    ind = randperm(K,2); % minimal set, 2 rays
    M = lf_estimate_point_manifold(rays(ind,:)); % model
    err = sum(( M * hom(rays') ).^2); % error 
    inliers = err < inlier_threshold; % inliers
    if sum(inliers) > sum(best_inliers) % 
        best_inliers = inliers;
        % updating outlier percentage and adapt number of iterations
        e = 1-sum(inliers) / length(inliers);
        max_iter = uint32(log(1-p)/log(1-(1-e)^2));
    end    
    
    iter = iter + 1;
    
end

M = lf_estimate_point_manifold(rays(best_inliers,:));
