%% Light Field Structure from Motion of camera sequence
% Input: 		matches : 	A set of matches a feature detections
%			f 	: 	
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras
function [pose,matches] =  lf_sfm_sequence_ransac(matches,f,H,varargin)

    % parsing input parameters 
    p = inputParser;
    addRequired (p,'matches');
    addRequired(p,'f',@isnumeric);
    addRequired(p,'H',@isnumeric);
    addParameter(p,'inlier_threshold',0.0001,@isnumeric);
    addParameter(p,'p',0.999,@isnumeric);

    parse(p,matches,f,H,varargin{:});
    f = p.Results.f;
    % calibration matrix
    H = p.Results.H;
    % inlier threshold for SFM RANSAC
    inlier_threshold = p.Results.inlier_threshold;
    % p for SFM RANSAC
    p = p.Results.p;

    
    N = length(matches);
    pose_pairwise = zeros(3,4,N);
    for comb=1:N
        fprintf('\ncombination (%i,%i)\n',matches(comb).lf1,matches(comb).lf2);        
        [pose_pairwise(1:3,1:3,comb),pose_pairwise(1:3,4,comb),inliers] = lf_sfm_pair_ransac(matches(comb),f,'inlier_threshold',inlier_threshold,'p',p);
        matches(comb) = lf_select_points(matches(comb),inliers);
    end

    
    pose = lf_pairwise_pose_to_sequence(pose_pairwise,H);
  
end
