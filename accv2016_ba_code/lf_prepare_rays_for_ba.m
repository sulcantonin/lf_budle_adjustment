%% Function builds structures (3D point cloud, rays) for bundle adjustment
% Inputs :      matches :   a set of feature matches, preferably only inliers
%                               after SfM (RANSAC), bundle adjustment is very
%                               sensitive to outliers!, to see structure
%                               checkout function "lf_extract_matches"
%               pose :      camera poses where the first camera is assumed to
%                               be origin. Inverted (pose \in R^{3 x 4 x m})
%               f :         the focal length, distance between two planes,
%                               needs to be defined for unprojection
% Outputs:      X :         a 3D point cloud of all cameras
%               rays :  Cell array (n x m) where each non-empty cell (i,j)
%                           represents a set of rays of a i-th 3D point visible in 
%                           j-th camera in light field coordinates [u,v,s,t]
% 
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function [X,rays,disparity] = lf_prepare_rays_for_ba(matches,pose,f)

    n = length(matches);
    n_cams = length(matches)+1;
    n_matches = arrayfun(@(i) matches(i).npoints, 1:n);

    % true if point was projected with another feature
    point_projected = cell(n_cams,1);
    % for each feature in each camera we have row with indices that point to 
    % feature (row) endoed by column as camera 
    feature_projection = cell(n_cams,1);

    for c1=1:n_cams-1
        point_projected{c1} = false(n_matches(c1),1);
        feature_projection{c1} = nan(n_matches(c1),n_cams);
    end
    point_projected{end} = false(n_matches(end),1);
    feature_projection{end} = nan(n_matches(c1),n_cams);


    for c1=1:n % for each set of matches
        for i=1:n_matches(c1) % for each match 
            pair = matches(c1).matches(i,:);
            for c2 = c1 + 1 : n
                ind = find(matches(c2).matches(:,1) == pair(2));
                if ~isempty(ind) && ~point_projected{c2}(ind)
                    pair = matches(c2).matches(ind,:);
                    point_projected{c2}(ind) = true;
                    feature_projection{c1}(i,c2+1) = ind;
                else
                    break
                end
            end
        end
    end

    rays = [];
    disparity = [];
    for j=1:n
        r = cell(n_matches(j),n_cams);
        d = nan(n_matches(j),n_cams);
        r(:,j) =   matches(j).rays1;
        r(:,j+1) = matches(j).rays2;
        d(:,j) =   matches(j).disparity1;
        d(:,j+1) = matches(j).disparity2;

        [feature, camera] = find(~isnan(feature_projection{j}));
        for i=1:length(feature)
            ind = feature_projection{j}(feature(i),camera(i));
            r(feature(i),camera(i)) =  matches(camera(i)-1).rays2(ind);
            d(feature(i),camera(i)) =  matches(camera(i)-1).disparity2(ind);
        end

        r = r(~point_projected{j},:);
        d = d(~point_projected{j},:);

        rays = cat(1,rays,r);
        disparity = cat(1,disparity,d);
    end
    
    % unprojecting rays into 3d points
    X = zeros(3,size(rays,1));
    for i=1:size(rays,1)
        r = rays(i,:);
        d = disparity(i,:);
        X_ = [];
        for j=1:length(r)
            if ~isempty( r{j} )
                Xray = lf_get_3DP_from_raysnd(r(j),d(j),f);
                [Ri,ti] = invert_Rt(pose(1:3,1:3,j),pose(:,4,j));
                Xray = [Ri,ti] * hom(Xray{1});
                X_ = cat(2,X_,Xray);
            end
        end
        X(:,i) = mean(X_')';
    end