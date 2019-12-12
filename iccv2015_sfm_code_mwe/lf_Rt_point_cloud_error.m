%% Total squared distance error for 3D point cloud estimates
%% computed from first and second light field,
%% given feature matches with disparity estimates,
%% rotation, translation and both focal lengths.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [err] = lf_Rt_point_cloud_error( matches, R, t, f1, f2 )

    npoints = length( matches.rays1 );
    assert( length( matches.rays2 ) == npoints );

    [Ri ti] = invert_Rt( R,t );
    err = 0.0;
    
    % for each feature match:
    for i=1:npoints
        % compute location of 3D points from all feature lists,
        % add up all pairwise errors
        % TODO: better, do a real triangulation and check backprojection
        % errors.
        nprojections1 = size( matches.rays1{i}, 1 );
        nprojections2 = size( matches.rays2{i}, 1 );
        assert( size( matches.rays1{i}, 2 ) == 4 );
        assert( size( matches.rays2{i}, 2 ) == 4 );
        
        for j=1:nprojections1

            X2 = lf_unproject_absolute( matches.rays2{i} (j,:), matches.disparity2( i ), f2 );
            X2t = Ri*X2 + ti;
            for k=1:nprojections2
                X1 = lf_unproject_absolute( matches.rays1{i} (k,:), matches.disparity1( i ), f1 );
                err = err + norm( X2t - X1 )^2;
            end
        end
    end

end
