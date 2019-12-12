function [ f,error_inlier ] = lf_estimate_f( matches,rastering )
    error_inlier = [];
    for i=1:length(rastering)
        fprintf('trying f=%f of (%i of %i)\n',rastering(i),i,length(rastering));
        [ R, t,inliers] = lf_sfm_pair_ransac( matches, rastering(i));
        submatches = lf_select_points(matches,inliers);
        error_inlier(i) = lf_forward_warp_error( matches, R, t, rastering(i),rastering(i));
    end
    [~,ind] = min(error_inlier);
    f = rastering(ind);
    plot(rastering, error_inlier); xlabel('f'); ylabel('error')
    
end

