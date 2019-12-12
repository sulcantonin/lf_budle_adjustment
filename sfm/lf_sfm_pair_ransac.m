%% Light Field Structure from Motion of camera pair
% Input: 		matches : 	A pair of matches
%			f 	: 	focal length
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function [R,t,best_inliers] = lf_sfm_pair_ransac(matches,f,varargin)

% parsing input parameters
p = inputParser;
addRequired (p,'matches');
addRequired(p,'f',@isnumeric);
addParameter(p,'inlier_threshold',0.0001,@isnumeric);
addParameter(p,'p',0.999,@isnumeric);

parse(p,matches,f,varargin{:});
% inlier threshold for SFM RANSAC
inlier_threshold = p.Results.inlier_threshold;
% p for SFM RANSAC
p = p.Results.p;


% inlier probability
e = 0.5; 
% minimum number of points
s = 3; 
% maxium number of iterations (check Multiple View Geometry in Computer
% Vision, formula 4.18)
N = log(1-p)/log(1-(1-e)^2);

npoints = matches.npoints;
iter = 1;

best_inliers = false(npoints,1);
err = zeros(npoints,1);
while iter < N
    subindices = false(npoints,1);
    subindices(randperm(npoints,s)) = true;
    submatches = lf_select_points(matches,subindices);    
    [R,t] = lf_sfm(submatches,f,f,20);
    
    % calculating error for individual paris to find inliers
    for i=1:npoints
        indices = false(npoints,1);
        indices(i) = true;
        pair = lf_select_points(matches,indices);
        err(i) = lf_forward_warp_error(pair,R,t,f,f);
    end
    
    inliers=err<inlier_threshold;
    if sum(inliers) > sum(best_inliers)
        best_inliers = inliers;
        % updating outlier percentage and adapt number of iterations
        e = 1-sum(inliers) / length(inliers);
        N = log(1-p)/log(1-(1-e)^s);
        
        fprintf('iter %i maxiter %i inliers %i of %i err %f \n',uint16(iter), uint16(N), uint16(sum(inliers)), uint16(length(inliers)),sum(err));
    end
    iter = iter + 1;
    
end

matches = lf_select_points(matches,best_inliers);
[R,t] = lf_sfm(matches ,f,f,50);
