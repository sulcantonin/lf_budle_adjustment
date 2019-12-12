function X = lf_unproject_to_point_cloud(rays,disparity,pose)    
    N = size(rays_for_ba,1);

    % unprojecting 3d points 
    X = zeros(3,N);
    for i=1:N
        r = rays(i,:);
        d = disparity(i,:);
        X_ = [];
        for j=1:length(r)
            if ~isempty( r{j} )
                X_ray = lf_get_3DP_from_raysnd(r(j),d(j),f);
                [Ri,ti] = invert_Rt(pose{j}.R,pose{j}.t);
                X_ray = [Ri,ti] * hom(X_ray{1});
                X_ = cat(2,X_,X_ray);
            end
        end
        X(:,i) = mean(X_')';
    end
end