%% function which finds features and their matches betweeen corresponding center views
% some outlining features which doesn't fit to 2D manifold M are filtered
% out by RANSAC
% Input :  LF                   : cell array with light fields, coordinates should
%                                   be ordered (n,m,y,x) (the same ordering as
%                                   Light Field Toolbox 4.0)
%          H                    : intrinsic calibration matrix (pixel coordiantes
%                                   to light field coordinates)
%          (inlier_threshold)   : threshold for error of a projection of the same
%                                   3D point for RANSANC
%          (p)                  : probablity p for RANSAC
%          (support)            : minimum number of feature rays per light
%                                   field feature             
%          (cross)              : true if features are detected only in
%                                   crosshair, otherwise in all subaperture views
% 
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function [matches,R] = lf_extract_matches(LF,H,varargin)
vec = @(i) i(:); % vectorize function 

% parsing input parameters 
argparser = inputParser;
addRequired (argparser,'LF');
addRequired (argparser,'H');
addParameter(argparser,'inlier_threshold',0.0001);
addParameter(argparser,'p',0.999);
addParameter(argparser,'support',2);
addParameter(argparser,'cross',false);

parse(argparser,LF,H,varargin{:});
inlier_threshold = argparser.Results.inlier_threshold;
p = argparser.Results.p;
support = argparser.Results.support;
cross = argparser.Results.cross;

N = length(LF); 
R = cell(N,1); % features [m,n,x,y] (in pixel coordinates)
D = cell(N,1); % center view descriptors


% feature detection and descriptor extraction
fprintf('Feature detection [')
for j=1:length(LF)
    [R{j},D{j}] = lf_extract_features(LF{j},'support',support,'cross',cross);
    fprintf('.')
end
fprintf(']')
%%
combinations = reshape([1;vec(repmat(2:length(LF)-1,1,2)); length(LF)]',[],2);
ncomb = size(combinations,1);
matches = [];
for comb=1:ncomb
    m1 = combinations(comb,1);
    m2 = combinations(comb,2);
    matches(comb).matches = vl_ubcmatch(D{m1},D{m2})';
    matches(comb).lf1 = m1;
    matches(comb).lf2 = m2;
end

for comb=1:ncomb

    m1 = matches(comb).matches(:,1);
    m2 = matches(comb).matches(:,2);
    lf1 = matches(comb).lf1;
    lf2 = matches(comb).lf2;
    
    n = length(m1);
    matches(comb).npoints = n;
    
    matches(comb).rays1 = R{lf1}(m1);
    matches(comb).rays2 = R{lf2}(m2);
    matches(comb).disparity1 = nan(n,1);
    matches(comb).disparity2 = nan(n,1);
    matches(comb).M1s = cell(n,1);
    matches(comb).M2s = cell(n,1);
    
    keep_rays = true(n,1);
    
    for i=1:n
        r1 = matches(comb).rays1{i};
        r2 = matches(comb).rays2{i};
        % pixel coordinates to light field coordinateas
        r1 = (H(:,:,lf1) * hom(r1'))';
        r1 = r1(:,1:4);    
        r2 = (H(:,:,lf2) * hom(r2'))';
        r2 = r2(:,1:4);
        
        % [s,t,u,v] to [u,v,s,t]
        r1 = change_ray_coordinates(r1,[3 4 1 2]);
        r2 = change_ray_coordinates(r2,[3 4 1 2]);        
        
        [M1,inliers1] = lf_estimate_point_manifold_ransac(r1,inlier_threshold,p);
        [M2,inliers2] = lf_estimate_point_manifold_ransac(r2,inlier_threshold,p);
        if sum(inliers1) < support || sum(inliers2) < support
            keep_rays(i) = false;
        end
        
        % storing inlier rays
        matches(comb).rays1{i} = r1(inliers1,:);
        matches(comb).rays2{i} = r2(inliers2,:);
        % storing inlier manifold
        matches(comb).M1s{i} = M1;
        matches(comb).M2s{i} = M2;
        % extracting and storing disparity 
        M1 = M1 ./ M1(1,1);
        M2 = M2 ./ M2(1,1);
        matches(comb).disparity1(i) = M1(1,3) + 1.0;
        matches(comb).disparity2(i) = M2(1,3) + 1.0;
    end
    
    % removing those who doesn't have enough rays (support parameter)
    matches(comb) = lf_select_points(matches(comb),keep_rays); 
    
    
end
